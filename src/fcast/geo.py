"""
Load the processed health-zone geography as epymorph inputs.

    from fcast.geo import load_geo
    geo = load_geo(provinces=["Ituri", "Nord-Kivu"])
    rume = SingleStrataRUME.build(..., scope=geo.scope, params={**geo.params, ...})

`geo.params` provides the per-node attributes epymorph expects:
- "population": total population (int, shape N)
- "centroid": population-weighted (longitude, latitude) (CentroidType, shape N);
  used by the built-in `Centroids` movement model
- "label": node ID (str, shape N); used by initializers like `LabeledLocations`
- "travel_time": road travel time in hours, row = origin, column = destination
  (float, shape NxN); used by `fcast.movement.Gravity`. Present when
  `fcast.build_travel` has been run.

Node IDs are DHIS2 health zone UIDs (GRID3 `zs_uid`). `CustomScope` sorts its node IDs,
so every array here, and `geo.table`, is in that same sorted order.
"""

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from epymorph.data_type import CentroidDType
from epymorph.geography.custom import CustomScope

DEFAULT_PATH = (
    Path(__file__).resolve().parents[2] / "data" / "processed" / "drc_health_zones.csv"
)
TRAVEL_PATH = DEFAULT_PATH.with_name("drc_travel_times.npz")


@dataclass(frozen=True)
class DRCGeo:
    scope: CustomScope
    params: dict[str, np.ndarray]
    table: pd.DataFrame
    """One row per node, in scope order."""

    def ids_for(self, names: Sequence[str], province: str | None = None) -> list[str]:
        """Look up node IDs by health zone name (case-insensitive)."""
        t = self.table
        if province is not None:
            t = t[t["province"].str.casefold() == province.casefold()]
        ids = []
        for name in names:
            match = t[t["name"].str.casefold() == name.casefold()]
            if len(match) != 1:
                where = f" in {province}" if province else ""
                raise KeyError(f"{len(match)} health zones named {name!r}{where}.")
            ids.append(match["node_id"].iloc[0])
        return ids


def load_geo(
    path: Path | str = DEFAULT_PATH,
    provinces: Sequence[str] | None = None,
    node_ids: Sequence[str] | None = None,
    travel_path: Path | str | None = TRAVEL_PATH,
) -> DRCGeo:
    """
    Build an epymorph scope and node parameters from the processed health zones.

    Parameters
    ----------
    path :
        The CSV written by `fcast.build_geo`.
    provinces :
        Keep only health zones in these provinces (case-insensitive).
    node_ids :
        Keep only these health zones.
    travel_path :
        The travel-time matrix written by `fcast.build_travel`. Skipped if None or
        the file doesn't exist.
    """
    table = pd.read_csv(path, dtype={"node_id": str})
    if provinces is not None:
        wanted = {p.casefold() for p in provinces}
        unknown = wanted - set(table["province"].str.casefold())
        if unknown:
            raise KeyError(f"Unknown provinces: {sorted(unknown)}")
        table = table[table["province"].str.casefold().isin(wanted)]
    if node_ids is not None:
        missing = set(node_ids) - set(table["node_id"])
        if missing:
            raise KeyError(f"Unknown node IDs: {sorted(missing)}")
        table = table[table["node_id"].isin(node_ids)]

    scope = CustomScope(table["node_id"].to_list())
    table = table.set_index("node_id").loc[scope.node_ids].reset_index()

    centroid = np.array(
        list(zip(table["centroid_lon"], table["centroid_lat"])), dtype=CentroidDType
    )
    params = {
        "population": table["population"].to_numpy(dtype=np.int64),
        "centroid": centroid,
        "label": table["node_id"].to_numpy(dtype=np.str_),
    }
    if travel_path is not None and Path(travel_path).exists():
        params["travel_time"] = _load_travel_time(travel_path, scope.node_ids)
    return DRCGeo(scope=scope, params=params, table=table)


def _load_travel_time(path: Path | str, node_ids: Sequence[str]) -> np.ndarray:
    """The travel-time matrix, subset and reordered to `node_ids`."""
    with np.load(path, allow_pickle=False) as data:
        index = {nid: i for i, nid in enumerate(data["node_ids"])}
        missing = [nid for nid in node_ids if nid not in index]
        if missing:
            raise KeyError(f"{len(missing)} node IDs missing from {path}; rebuild it.")
        idx = np.array([index[nid] for nid in node_ids])
        return data["travel_time_h"][np.ix_(idx, idx)]
