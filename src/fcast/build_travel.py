"""
Build a zone-to-zone road travel-time matrix for the DRC health zones.

Driving times come from the OSRM table service (OpenStreetMap road network, car
profile), queried between the population-weighted centroids written by
`fcast.build_geo`. Travel time for each pair is

    road time between the snapped points + time to reach the road at each end

OSRM snaps each centroid to the nearest road, which can be tens of km away in remote
zones; that gap is covered at `ACCESS_SPEED_KMH`. Some pairs have no road route (a
centroid snaps onto a disconnected road fragment); those fall back to great-circle
distance times `DETOUR_FACTOR` at `FALLBACK_SPEED_KMH`, and are flagged in `routed`.

OSRM's car speeds assume roads in typical condition, so absolute times are
optimistic for much of the DRC. Relative times are what the gravity model uses.

Run with:  uv run python -m fcast.build_travel
"""

import argparse
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from epymorph.util import pairwise_haversine

from fcast.build_geo import PROCESSED, RAW

OSRM_URL = "https://router.project-osrm.org/table/v1/driving/"
# The public OSRM server allows at most 100 coordinates per table request.
BLOCK = 50
ACCESS_SPEED_KMH = 10.0
FALLBACK_SPEED_KMH = 20.0
DETOUR_FACTOR = 1.4

ZONES_CSV = PROCESSED / "drc_health_zones.csv"
RAW_DURATIONS = RAW / "osrm_durations.npz"
OUT = PROCESSED / "drc_travel_times.npz"


def _osrm_block(coords: list[str], n_src: int) -> dict:
    """One table request: the first `n_src` coordinates are sources, the rest destinations."""
    sources = ";".join(str(i) for i in range(n_src))
    destinations = ";".join(str(i) for i in range(n_src, len(coords)))
    for attempt in range(5):
        resp = requests.get(
            OSRM_URL + ";".join(coords),
            params={
                "sources": sources,
                "destinations": destinations,
                "annotations": "duration",
            },
            timeout=120,
        )
        if resp.status_code == 429 or resp.status_code >= 500:
            time.sleep(2**attempt)
            continue
        resp.raise_for_status()
        body = resp.json()
        if body["code"] != "Ok":
            raise RuntimeError(f"OSRM error: {body['code']}: {body.get('message')}")
        return body
    resp.raise_for_status()
    raise RuntimeError("OSRM request failed after retries.")


def fetch_osrm_durations(
    zones: pd.DataFrame, dest: Path, refresh: bool = False
) -> tuple[np.ndarray, np.ndarray]:
    """
    Road durations (seconds, NaN where unroutable) between all zone centroids, and each
    centroid's snapping distance to the road network (meters).
    """
    if dest.exists() and not refresh:
        cached = np.load(dest, allow_pickle=False)
        if list(cached["node_ids"]) == zones["node_id"].to_list():
            return cached["durations_s"], cached["snap_m"]
    coords = [f"{x:.6f},{y:.6f}" for x, y in zip(zones["centroid_lon"], zones["centroid_lat"])]
    n = len(coords)
    durations = np.full((n, n), np.nan)
    snap = np.full(n, np.nan)
    blocks = range(0, n, BLOCK)
    for i0 in blocks:
        src = coords[i0 : i0 + BLOCK]
        for j0 in blocks:
            dst = coords[j0 : j0 + BLOCK]
            body = _osrm_block(src + dst, len(src))
            block = np.array(body["durations"], dtype=float)  # None -> nan
            durations[i0 : i0 + len(src), j0 : j0 + len(dst)] = block
            snap[i0 : i0 + len(src)] = [s["distance"] for s in body["sources"]]
            time.sleep(1.0)  # stay polite on the public server
        print(f"  rows {i0}-{i0 + len(src) - 1} of {n}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        dest, node_ids=zones["node_id"].to_numpy(dtype=np.str_), durations_s=durations, snap_m=snap
    )
    return durations, snap


def build(refresh: bool = False) -> dict:
    zones = pd.read_csv(ZONES_CSV, dtype={"node_id": str})
    durations_s, snap_m = fetch_osrm_durations(zones, RAW_DURATIONS, refresh)

    access_h = snap_m / 1000 / ACCESS_SPEED_KMH
    road_h = durations_s / 3600
    centroid = np.array(
        list(zip(zones["centroid_lon"], zones["centroid_lat"])),
        dtype=[("longitude", np.float64), ("latitude", np.float64)],
    )
    # pairwise_haversine returns miles.
    gc_km = pairwise_haversine(centroid) * 1.609344
    fallback_h = gc_km * DETOUR_FACTOR / FALLBACK_SPEED_KMH

    routed = np.isfinite(road_h)
    travel_h = np.where(routed, road_h + access_h[:, None] + access_h[None, :], fallback_h)
    np.fill_diagonal(travel_h, 0.0)
    np.fill_diagonal(routed, True)

    np.savez_compressed(
        OUT,
        node_ids=zones["node_id"].to_numpy(dtype=np.str_),
        travel_time_h=travel_h,
        routed=routed,
    )
    n = len(zones)
    unrouted_zones = zones.loc[(~routed).sum(axis=1) > n / 2, "name"].tolist()
    meta = {
        "built_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n_zones": n,
        "units": "hours",
        "layout": "travel_time_h[i, j] is from zone i to zone j; rows/cols follow node_ids",
        "source": {"name": "OSRM table service, car profile", "url": OSRM_URL},
        "access_speed_kmh": ACCESS_SPEED_KMH,
        "fallback": {
            "speed_kmh": FALLBACK_SPEED_KMH,
            "detour_factor": DETOUR_FACTOR,
            "pairs": int((~routed).sum()),
            "zones_mostly_unrouted": unrouted_zones,
        },
        "snap_km": {
            "median": round(float(np.median(snap_m)) / 1000, 2),
            "max": round(float(np.max(snap_m)) / 1000, 2),
        },
    }
    OUT.with_suffix(".meta.json").write_text(json.dumps(meta, indent=2))
    return meta


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--refresh", action="store_true", help="Re-query OSRM.")
    args = parser.parse_args()
    meta = build(refresh=args.refresh)
    print(f"Wrote {meta['n_zones']}x{meta['n_zones']} travel times to {OUT}")
    print(
        f"Unrouted pairs (great-circle fallback): {meta['fallback']['pairs']:,}; "
        f"zones mostly unrouted: {meta['fallback']['zones_mostly_unrouted']}"
    )


if __name__ == "__main__":
    main()
