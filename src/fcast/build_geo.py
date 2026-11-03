"""
Build the health-zone geography for the DRC metapopulation model.

Everything comes from GRID3:
- Health zone boundaries: GRID3 COD - Health Zones v9.0 (CIESIN / Ministère de la
  Santé Publique), pulled from GRID3's ArcGIS feature service.
- Population: GRID3 COD - Gridded Population Estimates v4.4 (~100 m grid), pulled from
  WorldPop, which hosts the GRID3 release.

For each health zone we sum the gridded population inside its boundary and compute a
population-weighted centroid (where people live, rather than the polygon's geometric
center, which matters for large rural zones).

Run with:  uv run python -m fcast.build_geo
"""

import argparse
import json
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import rasterio.mask
import requests

HZ_SERVICE = (
    "https://services3.arcgis.com/BU6Aadhn6tbBEdyk/arcgis/rest/services/"
    "GRID3_COD_health_zones_v9_0/FeatureServer/0"
)
POP_URL = (
    "https://data.worldpop.org/repo/wopr/COD/population/v4.4/"
    "COD_population_v4_4_gridded.zip"
)
EQUAL_AREA_CRS = "ESRI:102022"  # Africa Albers Equal Area Conic

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"


def fetch_health_zones(dest: Path, refresh: bool = False) -> gpd.GeoDataFrame:
    """Download all health zone polygons from the GRID3 feature service."""
    if dest.exists() and not refresh:
        return gpd.read_file(dest)
    features, offset, page = [], 0, 100
    while True:
        resp = requests.get(
            f"{HZ_SERVICE}/query",
            params={
                "where": "1=1",
                "outFields": "*",
                "orderByFields": "OBJECTID",
                "resultOffset": offset,
                "resultRecordCount": page,
                "outSR": 4326,
                "f": "geojson",
            },
            timeout=120,
        )
        resp.raise_for_status()
        batch = resp.json()["features"]
        features.extend(batch)
        if len(batch) < page:
            break
        offset += page
    gdf = gpd.GeoDataFrame.from_features(features, crs="EPSG:4326")
    dest.parent.mkdir(parents=True, exist_ok=True)
    gdf.to_file(dest, driver="GPKG")
    return gdf


def _find_raster(dest_dir: Path, kind: str) -> Path | None:
    """Find the release's `gridded`, `lower` or `upper` raster (names vary in case)."""
    for path in sorted(dest_dir.glob("**/*.tif")):
        if path.name.lower() == f"cod_population_v4_4_{kind}.tif":
            return path
    return None


def fetch_population_raster(dest_dir: Path, refresh: bool = False) -> dict[str, Path]:
    """Download and unzip the GRID3 v4.4 population rasters (estimate and bounds)."""
    kinds = ("gridded", "lower", "upper")
    existing = {k: _find_raster(dest_dir, k) for k in kinds}
    if all(existing.values()) and not refresh:
        return existing
    zip_path = dest_dir.parent / "COD_population_v4_4_gridded.zip"
    if refresh or not zip_path.exists():
        with requests.get(POP_URL, stream=True, timeout=600) as resp:
            resp.raise_for_status()
            with open(zip_path, "wb") as f:
                for chunk in resp.iter_content(chunk_size=1 << 20):
                    f.write(chunk)
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(dest_dir)
    found = {k: _find_raster(dest_dir, k) for k in kinds}
    missing = [k for k, v in found.items() if v is None]
    if missing:
        raise FileNotFoundError(f"Population rasters not found in zip: {missing}")
    return found


def _masked_sum(src: rasterio.DatasetReader, geom) -> float:
    data, _ = rasterio.mask.mask(src, [geom], crop=True, filled=True, nodata=np.nan)
    return float(np.nansum(np.clip(data[0], 0, None)))


def zonal_population(
    zones: gpd.GeoDataFrame, rasters: dict[str, Path]
) -> pd.DataFrame:
    """Population sum, bounds, and population-weighted centroid for each zone."""
    rows = []
    with (
        rasterio.open(rasters["gridded"]) as src,
        rasterio.open(rasters["lower"]) as lower,
        rasterio.open(rasters["upper"]) as upper,
    ):
        zones = zones.to_crs(src.crs)
        for geom in zones.geometry:
            data, transform = rasterio.mask.mask(
                src, [geom], crop=True, filled=True, nodata=np.nan
            )
            pop = np.where(np.isfinite(data[0]), data[0], 0.0)
            pop[pop < 0] = 0.0
            total = float(pop.sum())
            if total > 0:
                r, c = np.nonzero(pop)
                xs, ys = rasterio.transform.xy(transform, r, c)
                w = pop[r, c]
                lon = float(np.average(xs, weights=w))
                lat = float(np.average(ys, weights=w))
            else:
                lon = lat = np.nan
            rows.append(
                (total, _masked_sum(lower, geom), _masked_sum(upper, geom), lon, lat)
            )
    return pd.DataFrame(
        rows,
        columns=["population", "pop_lower", "pop_upper", "centroid_lon", "centroid_lat"],
        index=zones.index,
    )


def raster_total(raster_path: Path) -> float:
    """National population total in the raster, for a coverage check."""
    total = 0.0
    with rasterio.open(raster_path) as src:
        for _, window in src.block_windows(1):
            block = src.read(1, window=window, masked=True)
            total += float(block.sum() or 0.0)
    return total


def build(refresh: bool = False) -> pd.DataFrame:
    zones = fetch_health_zones(RAW / "grid3_cod_health_zones_v9_0.gpkg", refresh)
    rasters = fetch_population_raster(RAW / "COD_population_v4_4", refresh)

    stats = zonal_population(zones, rasters)
    geom_centroid = zones.to_crs(EQUAL_AREA_CRS).geometry.centroid.to_crs(4326)
    area_km2 = zones.to_crs(EQUAL_AREA_CRS).geometry.area / 1e6

    out = gpd.GeoDataFrame(
        {
            "node_id": zones["zs_uid"],
            "name": zones["zonesante"],
            "province": zones["province"],
            "province_id": zones["prov_uid"],
            "antenne": zones["antenne"],
            "grid3id": zones["grid3id"],
            "population": stats["population"].round().astype("int64"),
            # Sums of per-cell credible bounds: wider than a true zone-level interval.
            "population_lower": stats["pop_lower"].round().astype("int64"),
            "population_upper": stats["pop_upper"].round().astype("int64"),
            # Fall back to the geometric centroid for any zone with no gridded pop.
            "centroid_lon": stats["centroid_lon"].fillna(geom_centroid.x),
            "centroid_lat": stats["centroid_lat"].fillna(geom_centroid.y),
            "geom_centroid_lon": geom_centroid.x,
            "geom_centroid_lat": geom_centroid.y,
            "area_km2": area_km2,
        },
        geometry=zones.geometry,
        crs=zones.crs,
    )
    out["density_per_km2"] = out["population"] / out["area_km2"]
    out = out.sort_values("node_id").reset_index(drop=True)

    if out["node_id"].duplicated().any() or out["node_id"].isna().any():
        raise ValueError("Health zone IDs (zs_uid) are missing or not unique.")

    PROCESSED.mkdir(parents=True, exist_ok=True)
    out.to_file(PROCESSED / "drc_health_zones.gpkg", driver="GPKG")
    out.drop(columns="geometry").to_csv(PROCESSED / "drc_health_zones.csv", index=False)

    national = raster_total(rasters["gridded"])
    meta = {
        "built_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "node_id": "zs_uid (DHIS2 health zone UID)",
        "n_zones": int(len(out)),
        "population_total_zones": int(out["population"].sum()),
        "population_total_raster": round(national),
        "zones_with_zero_population": out.loc[out["population"] == 0, "name"].tolist(),
        "centroid": "population-weighted; geometric centroid where population is 0",
        "sources": {
            "health_zones": {
                "name": "GRID3 COD - Health Zones v9.0",
                "url": HZ_SERVICE,
                "hdx": "https://data.humdata.org/dataset/grid3-cod-health-zones-v9-0",
            },
            "population": {
                "name": "GRID3 COD - Gridded Population Estimates v4.4",
                "url": POP_URL,
                "hdx": "https://data.humdata.org/dataset/"
                "grid3-cod-gridded-population-estimates-v4-4",
            },
        },
    }
    (PROCESSED / "drc_health_zones.meta.json").write_text(json.dumps(meta, indent=2))
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--refresh", action="store_true", help="Re-download source data."
    )
    args = parser.parse_args()
    out = build(refresh=args.refresh)
    meta = json.loads((PROCESSED / "drc_health_zones.meta.json").read_text())
    print(f"Wrote {len(out)} health zones to {PROCESSED}")
    print(
        f"Population in zones: {meta['population_total_zones']:,} "
        f"(raster total {meta['population_total_raster']:,})"
    )


if __name__ == "__main__":
    main()
