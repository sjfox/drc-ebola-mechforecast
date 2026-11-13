"""
Prepare the inputs the R/Stan joint nowcast model needs from the Python side.

    uv run python -m fcast.bdbv_inputs --hub ../BDBV-Modeling-Hub

Writes:
- `data/processed/bdbv_zones.csv`: one row per hub health zone (`CD-IT-10`, ...) with
  its GRID3 node ID, province, and population. Zones are matched by DHIS2 code first
  (`health-zone-code-mapping.csv`), then by normalized name within the province.
- `data/processed/bdbv_kernel.csv.gz`: the dense gravity movement kernel over all 519
  zones (row = origin node ID, column = destination node ID, rows sum to 1), from
  `fcast.movement.gravity_kernel` with the same defaults as the epymorph `Gravity`
  model. It is too diffuse to sparsify (keeping 95% of each zone's travel takes ~260
  destinations on average), so the model restricts it to the zones it fits and
  renormalizes.
- `data/raw/bdbv/`: a snapshot of the hub's auxiliary files the model reads (reporting
  delays, including archived versions, daily counts, and the code mapping), so a run
  doesn't depend on the state of the hub checkout. Delay releases that the hub never
  archived (10-09) are recovered from its git history.
"""

import argparse
import io
import re
import shutil
import subprocess
import unicodedata
from pathlib import Path

import numpy as np
import pandas as pd

from fcast.movement import gravity_kernel

ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "data" / "processed"
RAW_BDBV = ROOT / "data" / "raw" / "bdbv"

# Hub zone names that don't normalize to the GRID3 name. Keys are hub codes.
NAME_OVERRIDES = {
    "CD-BC-20": "Massa",
    "CD-BC-23": "Muanda",
    "CD-EQ-11": "Mampoko",
    "CD-HK-25": "Ruashi",
    "CD-HL-13": "Malemba Nkulu",
    "CD-IT-14": "Gety",
    "CD-KC-02": "Bena Tshadi",
    "CD-KE-07": "Kabeya Kamwanga",
    "CD-KE-19": "Citenge",
    "CD-KG-07": "Kimbao",
    "CD-KG-14": "Wamba Luadi",
    "CD-MN-01": "Banzow Moke",
    "CD-MN-13": "Pendjwa",
    "CD-SA-03": "Djalo Ndjeka",
    "CD-SK-04": "Haut Plateau",
    "CD-SU-02": "Bogosenubea",
    "CD-TA-05": "Kiambi",
    "CD-TU-04": "Busanga",
    "CD-TU-12": "Yalifafu",
}

AUX_FILES = [
    "weekly-reporting-delays-by-province.csv",
    "daily-confirmed-case.csv",
    "daily-missing-symptom-by-report-date.csv",
    "health-zone-code-mapping.csv",
    "locations.csv",
]

DELAY_FILE = "auxiliary-data/weekly-reporting-delays-by-province.csv"
DELAY_COLS = ["week_ending", "location_name", "days_sym_to_report", "as_of"]


def _norm(s: str) -> str:
    s = unicodedata.normalize("NFKD", str(s)).encode("ascii", "ignore").decode().lower()
    s = re.sub(r"\bii\b", "2", s)
    s = re.sub(r"\bi\b", "1", s)
    return re.sub(r"[^a-z0-9]", "", s)


def snapshot_aux(hub: Path) -> None:
    """Copy the hub auxiliary files the model reads into data/raw/bdbv/."""
    aux = hub / "auxiliary-data"
    (RAW_BDBV / "archived-delay-data").mkdir(parents=True, exist_ok=True)
    for name in AUX_FILES:
        shutil.copy2(aux / name, RAW_BDBV / name)
    for f in sorted((aux / "archived-delay-data").glob("*.csv")):
        shutil.copy2(f, RAW_BDBV / "archived-delay-data" / f.name)
    recover_delay_vintages(hub)


def recover_delay_vintages(hub: Path) -> None:
    """Write delay-file releases that only survive in the hub's git history.

    The hub archives each delay release when the next one replaces it, but archiving
    started with 10-17: the 10-09 release is only in the commits before it. For every
    `as_of` in the delay file's history with no archived copy, write the last committed
    version of it in the archive's format.
    """
    def git(*args: str) -> str:
        return subprocess.run(["git", "-C", str(hub), *args], capture_output=True,
                              text=True, check=True).stdout

    archive = RAW_BDBV / "archived-delay-data"
    have = {f.stem.rsplit("-province-", 1)[1] for f in archive.glob("*.csv")}
    have |= set(pd.read_csv(RAW_BDBV / "weekly-reporting-delays-by-province.csv")["as_of"])
    try:
        commits = git("log", "--all", "--format=%H", "--", DELAY_FILE).split()
    except (subprocess.CalledProcessError, FileNotFoundError):
        print(f"{hub} is not a git checkout; unarchived delay releases not recovered.")
        return
    # Newest commit first, so the first version seen of each release is its last.
    for commit in commits:
        try:
            text = git("show", f"{commit}:{DELAY_FILE}")
        except subprocess.CalledProcessError:
            continue
        df = pd.read_csv(io.StringIO(text))
        # Early versions called the province name column `location`.
        if "location_name" not in df:
            df = df.rename(columns={"location": "location_name"})
        for as_of, rel in df.groupby("as_of"):
            if as_of in have:
                continue
            rel[DELAY_COLS].to_csv(
                archive / f"weekly-reporting-delays-by-province-{as_of}.csv", index=False
            )
            have.add(as_of)
            print(f"Recovered delay release {as_of} ({len(rel)} cases) from {commit[:7]}")


def build_zones(mapping: pd.DataFrame, geo: pd.DataFrame) -> pd.DataFrame:
    """Match each hub health zone to a GRID3 zone (DHIS2 code, then name)."""
    hz = mapping[mapping["location_type"] == "health zone"].copy()
    hz["prov_code"] = hz["location"].str[:5]

    # Province names differ only in accents and hyphens.
    geo = geo.assign(prov_key=geo["province"].map(_norm), name_key=geo["name"].map(_norm))
    hz["prov_key"] = hz["parent"].map(_norm)
    unknown = set(hz["prov_key"]) - set(geo["prov_key"])
    if unknown:
        raise ValueError(f"Hub provinces not in GRID3: {sorted(unknown)}")

    # A DHIS2 code listed for two hub zones (Bandalungwa and Kinshasa) is wrong for at
    # least one of them, so those fall through to name matching.
    by_dhis2 = dict(zip(geo["node_id"], geo["node_id"]))
    node = hz["Code_DHIS2"].map(by_dhis2).where(~hz["Code_DHIS2"].duplicated(keep=False))
    matched_by = np.where(node.notna(), "dhis2", None)

    name = hz["location"].map(NAME_OVERRIDES).fillna(hz["location_name"]).map(_norm)
    by_name = geo.set_index(["prov_key", "name_key"])["node_id"]
    for i in np.flatnonzero(node.isna().to_numpy()):
        key = (hz["prov_key"].iloc[i], name.iloc[i])
        if key in by_name.index:
            node.iloc[i] = by_name.loc[key]
            matched_by[i] = "name"
    hz["node_id"] = node
    hz["matched_by"] = matched_by

    dup = hz["node_id"].dropna().duplicated(keep=False)
    if dup.any():
        raise ValueError(f"GRID3 zones matched twice: {hz.loc[dup[dup].index, 'location']}")

    out = hz.merge(
        geo[["node_id", "province", "province_id", "population"]], on="node_id", how="left"
    )
    return out[
        ["location", "location_name", "prov_code", "parent", "node_id", "matched_by",
         "province", "province_id", "population"]
    ].rename(columns={"parent": "prov_name"}).sort_values("location")


def build_kernel(geo: pd.DataFrame, travel_path: Path) -> pd.DataFrame:
    """Dense gravity kernel over all GRID3 zones, indexed by node ID."""
    with np.load(travel_path, allow_pickle=False) as data:
        ids = data["node_ids"]
        tt = data["travel_time_h"]
    pop = geo.set_index("node_id").loc[ids, "population"].to_numpy(dtype=float)
    return pd.DataFrame(gravity_kernel(pop, tt), index=pd.Index(ids, name="from"), columns=ids)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--hub", type=Path, default=ROOT.parent / "BDBV-Modeling-Hub")
    args = parser.parse_args()

    snapshot_aux(args.hub)
    geo = pd.read_csv(PROCESSED / "drc_health_zones.csv", dtype={"node_id": str})
    mapping = pd.read_csv(RAW_BDBV / "health-zone-code-mapping.csv")

    zones = build_zones(mapping, geo)
    print(zones["matched_by"].value_counts(dropna=False).to_string())
    unmatched = zones[zones["node_id"].isna()]
    if len(unmatched):
        print("Unmatched hub zones (no population or movement):")
        print(unmatched[["location", "location_name", "prov_name"]].to_string(index=False))
    zones.to_csv(PROCESSED / "bdbv_zones.csv", index=False)

    kernel = build_kernel(geo, PROCESSED / "drc_travel_times.npz")
    kernel.to_csv(PROCESSED / "bdbv_kernel.csv.gz", float_format="%.6g")
    print(f"Wrote {len(zones)} zones and a {kernel.shape[0]}x{kernel.shape[1]} kernel.")


if __name__ == "__main__":
    main()
