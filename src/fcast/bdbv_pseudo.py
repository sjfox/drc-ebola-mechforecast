"""
Rebuild weekly pseudo-releases of the hub's province counts from its delay files.

    uv run python -m fcast.bdbv_pseudo

The hub has published only five releases of the case counts (08-18, then weekly from
10-09), but each delay-file release lists every case with an observed onset: province,
onset week, and days from onset to report. So for any date t before that release, we
can count the cases already reported by t. Those counts would be pseudo-releases at
any t, going back to May.

Taken at face value they run ahead of the real releases. A case's report date is not
when it appears in the published data: cases are often entered days or weeks later and
back-dated. Checked against the real releases, counts reconstructed from report dates
alone overstate the newest week by 20-110% and weeks 1-4 weeks old by 5-25%. So we
also estimate a publication lag G (days from report to publication) from the real
releases, and a case reported on day r counts as published by t with probability
G(t - r) / G(a - r), where a is the delay release it came from (every case in it was
published by a).

The delay files cover cases with observed onsets only (they match the hub's
`daily-confirmed-case.csv` counts with `imputed = FALSE` exactly). The pseudo-releases
are of those cases, and are checked against those counts. Imputed-onset cases (~15%)
are left out.

Writes:
- `data/processed/bdbv_pseudo_vintages.csv`: one row per delay release `source_as_of`,
  pseudo-release date `as_of` (every Saturday before the source release), province
  `location` and onset week `target_end_date`. `observation` is the expected count
  published by `as_of` (not an integer), `naive` the count reported by `as_of`, and
  `source_count` the count in the source release.
- `data/processed/bdbv_pseudo_validation.csv`: the leave-one-release-out check against
  the real releases, by week age.
- `data/processed/bdbv_pseudo_lag.json`: the fitted publication lag.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats

ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "data" / "processed"
RAW_BDBV = ROOT / "data" / "raw" / "bdbv"

AGE_BINS = [-1, 6, 13, 27, 55, 10_000]
AGE_LABELS = ["0-6", "7-13", "14-27", "28-55", "56+"]


def _saturday(d: pd.Series) -> pd.Series:
    """Epi-week end (Saturday) of each date."""
    return d + pd.to_timedelta((5 - d.dt.dayofweek) % 7, unit="D")


def load_delays() -> pd.DataFrame:
    """All delay-file releases: one row per case with an observed onset."""
    files = sorted((RAW_BDBV / "archived-delay-data").glob("*.csv"))
    files.append(RAW_BDBV / "weekly-reporting-delays-by-province.csv")
    # Rows are cases, so identical rows are distinct cases. A release in two files (the
    # latest may also be archived) is read from the first only.
    parts = {}
    for f in files:
        df = pd.read_csv(f)
        for as_of, rel in df.groupby("as_of"):
            parts.setdefault(as_of, rel)
    dl = pd.concat(parts.values(), ignore_index=True)
    codes = pd.read_csv(RAW_BDBV / "locations.csv")
    codes = codes[codes["location_type"] == "province"]
    dl["location"] = dl["location_name"].map(dict(zip(codes["location_name"], codes["location"])))
    if dl["location"].isna().any():
        raise ValueError(f"Unknown provinces: {dl.loc[dl['location'].isna(), 'location_name'].unique()}")
    # Negative delays (onset after report) are data errors; treat as same-day reports.
    dl["delay"] = dl["days_sym_to_report"].clip(lower=0)
    return dl.assign(week=pd.to_datetime(dl["week_ending"]), as_of=pd.to_datetime(dl["as_of"]))[
        ["as_of", "location", "week", "delay"]
    ]


def load_real_observed() -> pd.DataFrame:
    """Real releases: observed-onset counts by province and onset week."""
    dc = pd.read_csv(RAW_BDBV / "daily-confirmed-case.csv")
    dc = dc[(dc["location_type"] == "province") & ~dc["imputed"]]
    dc = dc.assign(as_of=pd.to_datetime(dc["as_of"]),
                   week=_saturday(pd.to_datetime(dc["target_end_date"])))
    return dc.groupby(["as_of", "location", "week"], as_index=False)["observation"].sum()


def lag_cdf(x: np.ndarray, theta: np.ndarray) -> np.ndarray:
    """Publication lag G(x) = P(published within x days of report).

    A share `pi` is published by the end of the report day; the rest after a
    discretized lognormal lag (median exp(m) days, log-spread s).
    """
    pi = 1 / (1 + np.exp(-theta[0]))
    m, s = theta[1], np.exp(theta[2])
    x = np.asarray(x, dtype=float)
    tail = stats.norm.cdf((np.log(np.maximum(x, 0) + 1) - m) / s)
    return np.where(x < 0, 0.0, pi + (1 - pi) * tail)


def _onset_report(delays: pd.DataFrame) -> pd.DataFrame:
    """One row per case and possible onset day (7 per case, weight 1/7 each).

    The delay file only gives the onset week; onsets are spread evenly over its days.
    """
    rep = delays.loc[delays.index.repeat(7)].reset_index(drop=True)
    rep["report"] = rep["week"] - pd.to_timedelta(np.tile(np.arange(6, -1, -1), len(delays)), unit="D") \
        + pd.to_timedelta(rep["delay"], unit="D")
    return rep


def published_by(rep: pd.DataFrame, t: pd.Timestamp, theta: np.ndarray | None) -> pd.Series:
    """Expected count published by t from one source release, by province and week.

    rep: `_onset_report` rows of one release, with `as_of` = the release date a. With
    theta=None, a case counts once it is reported (the naive reconstruction).
    """
    to_t = (t - rep["report"]).dt.days.to_numpy()
    if theta is None:
        p = (to_t >= 0).astype(float)
    else:
        to_a = (rep["as_of"] - rep["report"]).dt.days.to_numpy()
        p = lag_cdf(to_t, theta) / np.maximum(lag_cdf(to_a, theta), 1e-12)
    return pd.Series(p / 7, index=rep.index).groupby([rep["location"], rep["week"]]).sum()


def _pairs(delays: pd.DataFrame, real: pd.DataFrame, releases: list) -> list:
    """For each real release v, the next delay release after it as the source."""
    sources = sorted(delays["as_of"].unique())
    out = []
    for v in releases:
        later = [a for a in sources if a > v]
        if not later:
            continue
        a = later[0]
        rep = _onset_report(delays[delays["as_of"] == a].reset_index(drop=True))
        rv = real[real["as_of"] == v].set_index(["location", "week"])["observation"]
        n = delays[delays["as_of"] == a].groupby(["location", "week"]).size()
        out.append((v, a, rep, rv, n))
    return out


def fit_lag(pairs: list) -> np.ndarray:
    """Fit G by binomial likelihood: real count at v out of the source release count."""
    def nll(theta):
        total = 0.0
        for v, _, rep, rv, n in pairs:
            e = published_by(rep, v, theta)
            # Only weeks that had started by v; real counts above the source count (onset
            # corrections) are capped.
            e = e[e.index.get_level_values("week") - pd.Timedelta(days=6) <= v]
            k = rv.reindex(e.index, fill_value=0).clip(upper=n.reindex(e.index)).to_numpy()
            nn = n.reindex(e.index).to_numpy()
            p = np.clip(e.to_numpy() / nn, 1e-9, 1 - 1e-9)
            total -= np.sum(k * np.log(p) + (nn - k) * np.log1p(-p))
        return total

    res = optimize.minimize(nll, x0=np.array([0.0, np.log(5), 0.0]), method="Nelder-Mead",
                            options={"xatol": 1e-4, "fatol": 1e-4, "maxiter": 2000})
    return res.x


def describe_lag(theta: np.ndarray) -> dict:
    g = lag_cdf(np.array([0, 3, 7, 14, 28, 56]), theta)
    return {
        "pi_same_day": float(1 / (1 + np.exp(-theta[0]))),
        "tail_median_days": float(np.exp(theta[1])),
        "tail_log_spread": float(np.exp(theta[2])),
        "G": {f"{d}d": round(float(x), 3) for d, x in zip([0, 3, 7, 14, 28, 56], g)},
    }


def validate(delays: pd.DataFrame, real: pd.DataFrame) -> pd.DataFrame:
    """Leave one real release out: fit G on the others, reconstruct the held-out one."""
    releases = sorted(real["as_of"].unique())
    pairs = _pairs(delays, real, releases)
    rows = []
    for i, (v, a, rep, rv, n) in enumerate(pairs):
        theta = fit_lag(pairs[:i] + pairs[i + 1:])
        for method, th in [("naive", None), ("lag-adjusted", theta)]:
            e = published_by(rep, v, th).rename("reconstructed")
            j = pd.concat([e, rv.rename("real"), n.rename("source_count")], axis=1).fillna(0)
            j = j.reset_index().assign(release=v, source_as_of=a, method=method)
            j["age"] = (v - j["week"]).dt.days
            rows.append(j[j["age"] >= 0])
    out = pd.concat(rows, ignore_index=True)
    out["age_bin"] = pd.cut(out["age"], AGE_BINS, labels=AGE_LABELS)
    out = out.groupby(["release", "source_as_of", "method", "age_bin"], observed=True)[
        ["real", "reconstructed", "source_count"]
    ].sum().reset_index()
    out["ratio"] = out["reconstructed"] / out["real"]
    return out


def build(delays: pd.DataFrame, theta: np.ndarray) -> pd.DataFrame:
    """Pseudo-releases on every Saturday before each delay release."""
    rows = []
    for a, d in delays.groupby("as_of"):
        rep = _onset_report(d.reset_index(drop=True))
        n = d.groupby(["location", "week"]).size().rename("source_count")
        first = d["week"].min()
        for t in pd.date_range(first, a - pd.Timedelta(days=1), freq="W-SAT"):
            e = published_by(rep, t, theta).rename("observation")
            naive = published_by(rep, t, None).rename("naive")
            j = pd.concat([e, naive, n], axis=1).reset_index()
            j = j[j["week"] - pd.Timedelta(days=6) <= t]
            rows.append(j.assign(source_as_of=a, as_of=t))
    out = pd.concat(rows, ignore_index=True).rename(columns={"week": "target_end_date"})
    out = out[["source_as_of", "as_of", "location", "target_end_date", "observation", "naive",
               "source_count"]]
    for c in ["source_as_of", "as_of", "target_end_date"]:
        out[c] = out[c].dt.date
    return out.round({"observation": 3, "naive": 3})


def main() -> None:
    delays = load_delays()
    real = load_real_observed()

    check = validate(delays, real)
    check.to_csv(PROCESSED / "bdbv_pseudo_validation.csv", index=False, float_format="%.3f")
    wide = check.pivot_table(index=["release", "age_bin"], columns="method", values="ratio",
                             observed=True)
    print("Reconstructed / real (observed-onset cases), leave-one-release-out:")
    print(wide.round(3).to_string())

    pairs = _pairs(delays, real, sorted(real["as_of"].unique()))
    theta = fit_lag(pairs)
    lag = describe_lag(theta)
    print("Publication lag (all releases):", json.dumps(lag))
    (PROCESSED / "bdbv_pseudo_lag.json").write_text(
        json.dumps({"theta": theta.tolist(), **lag}, indent=2) + "\n")

    out = build(delays, theta)
    out.to_csv(PROCESSED / "bdbv_pseudo_vintages.csv", index=False)
    print(f"Wrote {len(out)} rows: {out['as_of'].nunique()} pseudo-release dates from "
          f"{out['source_as_of'].nunique()} delay releases.")


if __name__ == "__main__":
    main()
