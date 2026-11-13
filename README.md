# DRC outbreak forecasting

Metapopulation forecasting for the ongoing DRC outbreak, built on
[epymorph](https://github.com/NAU-CCL/Epymorph). Each node is a health zone.

## Setup on a new machine

The steps below get the joint R/Stan model (the forecasting model) running from a fresh
clone. Later sections are the reference for each piece.

**Prerequisites**

- [uv](https://docs.astral.sh/uv/) (it installs Python 3.12, which `pyproject.toml`
  requires).
- R ≥ 4.4.
- A C++ toolchain for CmdStan: Xcode command-line tools on macOS
  (`xcode-select --install`), `build-essential` on Linux, RTools on Windows.
- git, and access to the hub repo
  ([sjfox/BDBV-Modeling-Hub](https://github.com/sjfox/BDBV-Modeling-Hub)).

**Steps**

1. Clone this repo and the hub side by side. The hub must be a full git clone, not a
   downloaded copy: `fcast.bdbv_inputs` recovers the 10-09 delay release from the hub's
   git history.

   ```bash
   git clone https://github.com/sjfox/drc-ebola-mechforecast.git
   git clone https://github.com/sjfox/BDBV-Modeling-Hub.git
   cd drc-ebola-mechforecast
   ```

2. Install the Python environment:

   ```bash
   uv sync
   ```

3. Install the R packages and CmdStan (once per machine; CmdStan takes several minutes
   to build):

   ```r
   install.packages(c("posterior", "dplyr", "tidyr", "readr", "ggplot2", "patchwork",
                      "ggrepel", "scales"))
   install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
   cmdstanr::install_cmdstan()
   ```

4. Copy the hub data in and prepare the model inputs. Repeat this after each new hub
   release (`git -C ../BDBV-Modeling-Hub pull` first).

   ```bash
   mkdir -p data/raw
   cp ../BDBV-Modeling-Hub/target-data/time-series.csv data/raw/time-series.csv
   uv run python -m fcast.bdbv_inputs --hub ../BDBV-Modeling-Hub
   ```

   This writes `data/raw/bdbv/` (the hub's auxiliary files) and refreshes
   `data/processed/bdbv_zones.csv` and `bdbv_kernel.csv.gz`.

5. Check that everything runs with a short fit (2 chains, 200 iterations; the posterior
   isn't reliable, but it exercises the whole pipeline and compiles the Stan model):

   ```bash
   Rscript scripts/run_bdbv.R 2026-11-01 quick
   ```

   It writes a hub sample file, diagnostics and plots to `output/bdbv/2026-11-01/`.

6. Run the full model on the latest release with `Rscript scripts/run_bdbv.R`, or the
   backtests in "Finalized-data backtests" below.

**Optional**

- The geography build (`fcast.build_geo`, about 250 MB of downloads) and the OSRM
  travel times (`fcast.build_travel`) are only needed to change the zone geography or
  movement inputs. Their outputs are committed in `data/processed/`.
- `fcast.bdbv_pseudo` is a diagnostic and isn't used by the model.

**Run time and cores.** Fits use 4 chains (4 cores) by default. Set `BDBV_CHAINS` to
use fewer, at the cost of fewer posterior samples. Backtest fits take 2–25 minutes per
forecast date, growing with the number of weeks fitted, so a full 13-date backtest
takes a few hours.

**Not in git.** `output/` (fits, samples, scores) is gitignored. Results from another
machine have to be copied across (`output/final/` for the backtests) or regenerated.

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

### Road travel times

```bash
uv run python -m fcast.build_travel
```

Queries the public [OSRM](https://project-osrm.org/) table service (OpenStreetMap
roads, car profile) between the population-weighted centroids and writes
`data/processed/drc_travel_times.npz`: a 519×519 matrix of travel times in hours
(`travel_time_h[i, j]` is from zone i to zone j), plus a `routed` mask. Each time is the
road time plus the time to reach the road from each centroid at 10 km/h. Pairs with no
road route fall back to great-circle distance × 1.4 at 20 km/h. Raw OSRM responses are
cached in `data/raw/osrm_durations.npz`; pass `--refresh` to re-query. Takes a few
minutes (121 requests, since the public server allows 100 coordinates per request).

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

`geo.params` supplies `population`, `centroid`, `label`, and (once
`fcast.build_travel` has run) `travel_time`. `geo.table` is the zone
table in the same order as `geo.scope.node_ids`. `CustomScope` sorts IDs, so always
line up other per-node data against `geo.table`, not the CSV's order.

`scripts/smoke_test_epymorph.py` runs this end to end with placeholder parameters.

### Notes

- **Connectivity.** Two movement models, both daily commute-and-return with a
  commuting fraction (`commuter_proportion`):
  - epymorph's `Centroids`: exponential kernel over straight-line distance (`phi`).
  - `fcast.movement.Gravity`: a gravity kernel over road travel time,
    `K[i, j] ∝ population[j]^gravity_pop_exp / (travel_time[i, j] + gravity_time_offset)^gravity_time_exp`,
    with no self-travel, so `commuter_proportion` is the share leaving each zone daily.
    Defaults are 1, 2, and 1 hour; all need fitting. `fcast.movement.gravity_kernel`
    computes the same matrix outside a simulation.

  Swap in `mm=Gravity()` (from `fcast.movement import Gravity`) for `mm.Centroids()`.
  Mobile-phone flows from Flowminder are a different source and are not used here.
- **Travel times are optimistic.** OSRM assumes typical road speeds for each road
  class; many DRC roads are much slower, and seasonal. Relative times matter more than
  absolute ones here, but `gravity_time_offset` is in the same hours.
- **Comparison with INRB-UMIE's OSRM matrix.** [BDBV2026-Data](https://github.com/INRB-UMIE/BDBV2026-Data/tree/main/data/osrm)
  publishes a 519×519 car travel-time matrix from the same public OSRM server (regenerated
  2026-07-10 on GRID3 v8.0 boundaries; its `ZSCode` matches our `node_id`). Compared on
  2026-11-03, the two agree closely for longer trips (Spearman 0.99 for pairs ≥ 12 h,
  median ratio 1.00), but less for short ones (pairs < 3 h: Spearman 0.89, 16% differ
  by more than 1.5×), which is where the gravity kernel puts most weight: with default
  parameters about 7% of each zone's travel goes to different destinations (median). The
  difference comes from the routing points, not our road-access time. INRB routes from
  `st_point_on_surface`, an arbitrary interior point; we route from the
  population-weighted centroid, a median 9 km away (> 28 km for 10% of zones). For
  Mongbwalu their point is about 48 km west of the town, and for Kowe (Lubumbashi) about
  226 km off. Their matrix is also missing every route to or from Idjwi and Bokoro. We
  keep ours; theirs is a drop-in sensitivity check for comparing with the INRB
  invasion-risk model. Both use OSRM's optimistic car speeds (theirs implies a median
  55 km/h on the road).
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

## Joint nowcast/forecast model (R + Stan)

A joint model of transmission (deterministic metapopulation SEIR with gravity mixing)
and reporting delays, fit with Stan to the BDBV Modeling Hub data. Design decisions are
logged in [docs/bdbv_model_decisions.md](docs/bdbv_model_decisions.md).

Requires R with `cmdstanr` (CmdStan ≥ 2.36), `posterior`, `dplyr`, `tidyr`, `readr`,
`ggplot2`, `patchwork`, `ggrepel`, `scales`, plus a checkout of the hub next to this
repo. Developed with R 4.4.1, CmdStan 2.36.0, cmdstanr 0.9.0, posterior 1.6.1,
dplyr 1.2.0, tidyr 1.3.1, readr 2.2.0, ggplot2 4.0.2, patchwork 1.3.2, ggrepel 0.9.6 and
scales 1.4.0. The compiled model (`stan/bdbv_joint`) is gitignored and builds on the
first run.

```r
install.packages(c("posterior", "dplyr", "tidyr", "readr", "ggplot2", "patchwork",
                   "ggrepel", "scales"))
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()
```

```bash
cp ../BDBV-Modeling-Hub/target-data/time-series.csv data/raw/time-series.csv
uv run python -m fcast.bdbv_inputs --hub ../BDBV-Modeling-Hub
uv run python -m fcast.bdbv_pseudo               # diagnostic only, not used by the model (P1-P4)
Rscript scripts/run_bdbv.R                       # latest vintage
Rscript scripts/run_bdbv.R 2026-10-24            # an earlier vintage
Rscript scripts/backtest_bdbv.R 2026-10-17 2026-10-24
Rscript scripts/synthetic_check_bdbv.R output/bdbv/2026-11-01
```

Each run writes the hub sample file, diagnostics and plots to `output/bdbv/<as_of>/`.

Forecasts damp each zone's log R toward 0 (R toward 1) by a factor of 0.5 per forecast
week (`forecast_R = "damped"`, `damp_R = 0.5` in `bdbv_config()`), chosen by the
finalized-data backtests below; see "Forecast R damping" in the decisions doc.

### Finalized-data backtests

Weekly forecast dates from 2026-05-09 to 2026-08-01, each fit to complete counts up to
that date and scored with WIS against the hub flatline baseline, using only weeks at
least 60 days old in the latest release as truth.

```bash
Rscript scripts/backtest_final.R E0                # fits (a few hours with 4 chains)
Rscript scripts/backtest_final.R E0_damp050        # damped forecasts from E0's fits (minutes)
Rscript scripts/score_final.R E0_damp050           # output/final/scores.csv, summary.csv
Rscript scripts/plot_forecasts_final.R E0_damp050  # figures in output/final/review/
uv run python scripts/build_zone_report.py forecast_review   # reports/forecast_review.html
```

`scripts/backtest_final.R` lists the variants (E0, the damping grid E0_dampNNN, E3_ar
and F1 with in-fit AR(1) persistence). Variants marked "generated quantities only"
reuse another variant's fits, so run that one first. Set `BDBV_CHAINS` to use fewer
cores. Results are not tracked (`output/` is gitignored), so a new machine needs the E0
fits rerun, or `output/final/` copied across.
