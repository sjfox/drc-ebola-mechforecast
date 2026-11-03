# DRC outbreak forecasting

Metapopulation forecasting for the ongoing DRC outbreak, built on
[epymorph](https://github.com/NAU-CCL/Epymorph). Each node is a health zone.

## Geography and population

All health-zone geography and population come from a single source, GRID3:

| Input | Source |
|---|---|
| Health zone boundaries (519 zones, 26 provinces) | [GRID3 COD - Health Zones v9.0](https://data.humdata.org/dataset/grid3-cod-health-zones-v9-0), pulled from GRID3's ArcGIS feature service |
| Population | [GRID3 COD - Gridded Population Estimates v4.4](https://data.humdata.org/dataset/grid3-cod-gridded-population-estimates-v4-4) (~100 m grid), pulled from WorldPop, which hosts the release |

Node IDs are DHIS2 health zone UIDs (GRID3 `zs_uid`), so case data from DHIS2 or line
lists that carry the UID join directly.

### Build

```bash
uv sync
uv run python -m fcast.build_geo
```

This downloads the source data to `data/raw/` (gitignored, ~250 MB) and writes:

- `data/processed/drc_health_zones.csv`: one row per health zone: `node_id`, `name`,
  `province`, `province_id`, `antenne`, `grid3id`, `population` (with
  `population_lower`/`population_upper`, the sums of GRID3's per-cell credible bounds,
  which are wider than a true zone-level interval), population-weighted
  `centroid_lon`/`centroid_lat`, geometric centroid, `area_km2`, `density_per_km2`
- `data/processed/drc_health_zones.gpkg`: the same, with boundaries
- `data/processed/drc_health_zones.meta.json`: sources, build time, and a population
  coverage check against the raster total

Pass `--refresh` to re-download the source data.

### Use in epymorph

```python
from epymorph.kit import *
from fcast.geo import load_geo

geo = load_geo(provinces=["Ituri", "Nord-Kivu"])  # or load_geo() for all 519 zones
seeds = geo.ids_for(["Bunia", "Mongbwalu", "Rwampara"], province="Ituri")

rume = SingleStrataRUME.build(
    ipm=ipm.SEIRS(),
    mm=mm.Centroids(),
    init=init.LabeledLocations(labels=seeds, seed_size=30),
    scope=geo.scope,
    time_frame=TimeFrame.of("2026-05-01", 90),
    params={**geo.params, "beta": 0.25, "sigma": 1 / 9, "gamma": 1 / 10, "xi": 0.0},
)
out = BasicSimulator(rume).run()
```

`geo.params` supplies `population`, `centroid`, and `label`. `geo.table` is the zone
table in the same order as `geo.scope.node_ids`. `CustomScope` sorts IDs, so always
line up other per-node data against `geo.table`, not the CSV's order.

`scripts/smoke_test_epymorph.py` runs this end to end with placeholder parameters.

### Notes

- **Connectivity.** epymorph's `Centroids` movement model needs only `population` and
  `centroid`. It uses an exponential distance kernel (`phi`) and a commuting fraction
  (`commuter_proportion`), both of which need fitting. Mobile-phone flows from
  Flowminder are a different source and are not used here.
- **Provincial updates.** GRID3 has published newer emergency boundary packages for
  Ituri and Nord-Kivu. v9.0 is used everywhere for consistency.
- **Population year.** GRID3 v4.4 estimates cover 2022–2024 and are not projected to
  2026.
- **Population differs by source.** The national total is 109.7M. GRID3 puts Bunia at
  663k, while Flowminder's phone-based estimate is 228k. Kinshasa comes out at 7.7M,
  well below most other estimates. For zones that drive the forecast, check these
  numbers against local denominators.
- **Centroids.** Population-weighted centroids can fall outside non-convex zones (13 of
  519). They are only used for distances.
