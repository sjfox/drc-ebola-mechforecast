# Build the Stan data for the joint BDBV nowcast/forecast model from:
# - data/raw/time-series.csv: hub target data, all vintages (as_of)
# - data/raw/bdbv/: hub auxiliary snapshot (code mapping), written by
#   `uv run python -m fcast.bdbv_inputs`
# - data/processed/bdbv_zones.csv, bdbv_kernel.csv.gz: zone populations and the gravity
#   kernel, also written by fcast.bdbv_inputs
#
# All dates are epi-week end dates (Saturdays) unless noted.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

bdbv_config <- function(as_of = NULL, ...) {
  cfg <- list(
    root = ".",
    as_of = as_of,                 # data vintage to fit; NULL = latest in time-series
    min_province_cases = 10,       # provinces modeled = those with at least this many
    lead_in_weeks = 3,             # model starts this many weeks before the outbreak
    seed_weeks = 3,                # zones with cases in the outbreak's first weeks are seeded
    knot_weeks = 4,                # spline knot spacing, counted back from the last week
    delay = "lognormal",           # reporting delay: "lognormal" (stan/bdbv_joint.stan) or
                                   # "hazard" (weekly hazard, stan/bdbv_joint_hazard.stan)
    d_max = 120,                   # hazard: eventual count = cases published within
                                   # d_max days of onset (decision R6)
    truth_as_of = NULL,            # backtests: a later release to predict counts for
    report_lag_days = 0,           # days between the data cutoff and the as_of date
    forecast_weeks_after_ref = 3,  # hub: forecasts through reference date + 3 weeks
    output_weeks_before_ref = 5,   # hub submissions cover reference date - 5 weeks
    mean_incubation = 7.5,         # days; Bundibugyo, Uganda 2007 (decision T2)
    mean_infectious = 8,           # days after onset; gives a 15.5-day generation
                                   # interval, the pooled Ebola serial interval (T2)
    eps = 0.01,                    # background expected cases per zone-week
    # Finalized-data backtests (scripts/backtest_final.R): fit the latest release's
    # counts for onset weeks up to `origin` as if fully reported (no reporting model),
    # and forecast `horizon_weeks` weeks past it.
    finalized = FALSE,
    origin = NULL,                 # last data week (a Saturday) when finalized
    horizon_weeks = 4,
    provinces = NULL,              # fixed province codes to model; NULL = by case count
    # Province / zone R settings (defaults = the current model).
    rho_zone_prior = c(4, 2),      # beta prior on the zone AR(1) persistence
    sd_ar_zone_scale = 0.3,        # half-normal scale of the zone AR(1) step sd
    sd_rw_R_scale = 0.2,           # half-normal scale of the province weekly R step sd
    forecast_R = "damped",         # province log R in forecast weeks: "rw" continues the
                                   # random walk, "flat" holds the last week's value,
                                   # "damped" decays each zone's log R by damp_R per week;
                                   # with ar_R, "rw" continues the AR(1) and "mean" uses
                                   # its mean decay phi_R^h
    damp_R = 0.5,                  # chosen by the finalized-data backtests (2026-11-06)
    ar_R = FALSE,                  # province log R: AR(1) toward 0 with fitted phi_R
    phi_R_prior = c(2, 2),         # beta prior on phi_R (weekly persistence)
    zone_weekly = FALSE            # zone deviation: weekly AR(1) instead of the spline
                                   # (rho_zone_prior then applies per week)
  )
  modifyList(cfg, list(...))
}

read_vintages <- function(cfg) {
  read_csv(file.path(cfg$root, "data/raw/time-series.csv"), show_col_types = FALSE) |>
    mutate(as_of = as.Date(as_of), target_end_date = as.Date(target_end_date))
}

# Epi-week end (Saturday) of the week containing `date`.
week_end_of <- function(date) {
  date + (6 - as.integer(format(date, "%w")))
}

build_bdbv_data <- function(cfg) {
  vint <- read_vintages(cfg)
  as_of <- if (is.null(cfg$as_of)) max(vint$as_of) else as.Date(cfg$as_of)
  if (!as_of %in% vint$as_of) stop("as_of ", as_of, " is not a vintage in time-series.csv")
  snap <- filter(vint, as_of == !!as_of)
  if (cfg$finalized) {
    origin <- as.Date(cfg$origin)
    snap <- filter(snap, target_end_date <= origin)
  }

  zones_all <- read_csv(file.path(cfg$root, "data/processed/bdbv_zones.csv"),
                        show_col_types = FALSE)

  # Outbreak start: first week of a run of two weeks with cases nationally. Isolated
  # earlier cases (Feb, Mar) are before the model starts and are left out.
  country <- snap |> filter(location_type == "country") |> arrange(target_end_date)
  pos <- country$observation > 0
  first_run <- which(pos & dplyr::lead(pos, default = FALSE))[1]
  outbreak_start <- country$target_end_date[first_run]
  model_start <- outbreak_start - 7 * cfg$lead_in_weeks
  last_week <- max(snap$target_end_date)
  weeks <- seq(model_start, last_week, by = 7)

  hz <- snap |> filter(location_type == "health zone")
  prov_totals <- hz |>
    mutate(prov_code = substr(location, 1, 5)) |>
    group_by(prov_code) |>
    summarise(cases = sum(observation), .groups = "drop")
  provinces <- prov_totals |>
    filter(cases >= cfg$min_province_cases) |>
    arrange(desc(cases)) |>
    pull(prov_code)
  if (!is.null(cfg$provinces)) provinces <- cfg$provinces

  zones <- zones_all |>
    filter(prov_code %in% provinces, !is.na(node_id)) |>
    mutate(prov = match(prov_code, provinces)) |>
    arrange(prov, location)

  # Gravity kernel restricted to the modeled zones; travel elsewhere is dropped.
  kern <- read_csv(file.path(cfg$root, "data/processed/bdbv_kernel.csv.gz"),
                   show_col_types = FALSE)
  K <- as.matrix(kern[match(zones$node_id, kern$from), zones$node_id])
  kept_movement <- rowSums(K)
  K <- K / kept_movement
  dimnames(K) <- NULL

  cases <- hz |>
    filter(location %in% zones$location, target_end_date %in% weeks) |>
    select(location, target_end_date, observation) |>
    complete(location = zones$location, target_end_date = weeks,
             fill = list(observation = 0)) |>
    pivot_wider(names_from = target_end_date, values_from = observation) |>
    arrange(match(location, zones$location))
  case_mat <- as.matrix(cases[, as.character(weeks)])
  storage.mode(case_mat) <- "integer"

  # Seeds: zones with cases in the first `seed_weeks` weeks of the outbreak.
  seed_weeks <- outbreak_start + 7 * (0:(cfg$seed_weeks - 1))
  seed_idx <- which(rowSums(case_mat[, as.character(seed_weeks), drop = FALSE]) > 0)

  # Reporting triangle: counts for the same province-week (summed over modeled zones) in
  # every release up to as_of. Each earlier release's count is modeled given the next
  # release's, so the triangle identifies when cases appear in the published data.
  # Province totals, not zone counts: cases are reassigned between zones of a province
  # from one release to the next. At zone level that shows up as extra noise (triangle
  # concentration ~7 vs ~36) and, once drops are capped, as spurious late additions;
  # together they gave v3 a heavy delay tail (decision R5). The few province counts that
  # drop between releases are capped at the later count.
  vint_dates <- sort(unique(vint$as_of[vint$as_of <= as_of]))
  V <- length(vint_dates)
  age_v <- t(sapply(seq_len(V), function(k) {
    as.integer(vint_dates[k] - weeks) - cfg$report_lag_days
  }))
  tri <- vint |>
    filter(as_of %in% vint_dates, location %in% zones$location,
           target_end_date %in% weeks) |>
    transmute(p = zones$prov[match(location, zones$location)],
              w = match(target_end_date, weeks), v = match(as_of, vint_dates),
              count = as.integer(observation)) |>
    group_by(p, w, v) |>
    summarise(count = sum(count), .groups = "drop") |>
    complete(p = seq_along(provinces), w = seq_along(weeks), v = seq_len(V)) |>
    filter(age_v[cbind(v, w)] >= -6) |>
    mutate(count = coalesce(count, 0L)) |>
    arrange(p, w, desc(v)) |>
    group_by(p, w) |>
    mutate(count = cummin(count), later = lag(count)) |>
    ungroup() |>
    filter(!is.na(later), later > 0)
  # With a d_max cap, weeks already older than d_max at the earlier release can't grow
  # under the model; their (small) changes carry no information about the delay.
  if (cfg$delay == "hazard") tri <- filter(tri, age_v[cbind(v, w)] < cfg$d_max)
  # Finalized data: one "release" in which every week is complete, and no triangle.
  if (cfg$finalized) {
    vint_dates <- as_of
    V <- 1L
    age_v <- matrix(100000L, 1, length(weeks))
    tri <- tri[0, ]
  }

  # Forecast horizon: the hub's reference date is the Saturday ending the epi week in
  # which the forecast is submitted (Sunday to Tuesday), i.e. the week of the first
  # Sunday on or after as_of. Weeks after the last data week are forecast weeks.
  first_sunday <- as_of + (7 - as.integer(format(as_of, "%w"))) %% 7
  ref_date <- week_end_of(first_sunday)
  H <- as.integer((ref_date + 7 * cfg$forecast_weeks_after_ref - last_week) / 7)
  if (cfg$finalized) H <- as.integer(cfg$horizon_weeks)
  all_weeks <- c(weeks, last_week + 7 * seq_len(H))

  # Cubic B-spline basis over all weeks (observed + forecast), with knots every
  # knot_weeks weeks counted back from the last observed week, so the most recent
  # segment has the same length in every fit. The knot at the last observed week and
  # the forecast segment after it keep the end of the fit from being pinned by the
  # least complete weeks.
  n_all <- length(all_weeks)
  W <- length(weeks)
  knots <- rev(seq(W, 3, by = -cfg$knot_weeks))
  Bs <- splines::bs(seq_len(n_all), knots = knots, degree = 3, intercept = TRUE,
                    Boundary.knots = c(1, n_all))
  Bs <- matrix(as.numeric(Bs), nrow(Bs))

  # Backtests: each week's age at a later release, to predict the count published by
  # then. Without one, the ages at as_of (y_truth = the reported counts).
  truth_date <- if (is.null(cfg$truth_as_of)) as_of else as.Date(cfg$truth_as_of)
  truth_age <- as.integer(truth_date - all_weeks) - cfg$report_lag_days

  stan_data <- list(
    Z = nrow(zones), P = length(provinces), W = length(weeks), H = H,
    prov = zones$prov, pop = as.numeric(zones$population), K = K,
    cases = case_mat,
    S_n = length(seed_idx), seed_idx = as.array(seed_idx),
    NB = ncol(Bs), Bs = Bs,
    V = V, age_v = age_v,
    N_tri = nrow(tri), tri_prov = tri$p, tri_week = tri$w, tri_v = tri$v,
    tri_early = tri$count, tri_late = tri$later,
    mean_incubation = cfg$mean_incubation, mean_infectious = cfg$mean_infectious,
    eps = cfg$eps, do_gq = 0L, truth_age = truth_age, d_max = as.integer(cfg$d_max),
    perfect_reporting = as.integer(cfg$finalized),
    rho_zone_prior = cfg$rho_zone_prior, sd_ar_zone_scale = cfg$sd_ar_zone_scale,
    sd_rw_R_scale = cfg$sd_rw_R_scale, forecast_R_flat = as.integer(cfg$forecast_R == "flat"),
    forecast_R_damped = as.integer(cfg$forecast_R == "damped"), damp_R = cfg$damp_R,
    ar_R = as.integer(cfg$ar_R), phi_R_prior = cfg$phi_R_prior,
    forecast_R_mean = as.integer(cfg$forecast_R == "mean"), knot_weeks = as.integer(cfg$knot_weeks),
    zone_weekly = as.integer(cfg$zone_weekly)
  )

  list(
    stan_data = stan_data, cfg = cfg, as_of = as_of, ref_date = ref_date,
    weeks = weeks, forecast_weeks = last_week + 7 * seq_len(H),
    zones = zones, provinces = provinces, seed_zones = zones$location[seed_idx],
    kept_movement = kept_movement, truth_date = truth_date,
    snapshot = snap, vintages = vint, vint_dates = vint_dates,
    outbreak_start = outbreak_start
  )
}
