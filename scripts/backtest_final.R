# Finalized-data backtests of the transmission model: fit the latest release's counts for
# onset weeks up to each forecast origin, treated as complete (no reporting model), and
# forecast 4 weeks ahead. scripts/score_final.R scores the forecasts.
#
#   Rscript scripts/backtest_final.R E0                       # every origin
#   Rscript scripts/backtest_final.R E1_zone_tight 2026-05-16 2026-06-13
#
# Origins default to every Saturday from 2026-05-09 to 2026-08-01: the last one whose
# 4-week horizon ends on a week that is at least 60 days old in the latest release.
# Runs one fit at a time with 4 chains on 4 cores (BDBV_CHAINS, BDBV_ITER,
# BDBV_TREEDEPTH override). Existing origins are skipped. BDBV_FINAL_OUT changes the
# output root (default output/final).
#
# Variants (province / zone R only; everything else is the current model):
#   E0             current model: zone AR(1) spline, province weekly random walk
#   E1_zone_tight  tighter zone spline: rho_zone ~ beta(8, 2), sd_ar_zone ~ N+(0, 0.15)
#   E2_prov_tight  tighter province walk: sd_rw_R ~ N+(0, 0.1)
#   E0_flat        E0's fits, with province log R held at its last value in forecast
#                  weeks (generated quantities only; needs E0 run first)
#   E3_ar          province log R a mean-reverting AR(1) toward R = 1, persistence phi_R
#                  ~ beta(2, 2) fitted at each origin; forecasts continue the AR(1), zone
#                  deviations decay at the fitted zone persistence
#   E3_ar_mean     E3_ar's fits, forecast log R following the AR(1) mean phi_R^h (no
#                  innovations; generated quantities only)
#   E0_dampNNN     E0's fits, each zone's log R multiplied by NNN/100 per forecast week
#                  (NNN = 100, 090, 080, 070, 060, 050, 030, 000; generated quantities only)
#   F1             E3_ar with the zone deviation a weekly AR(1) toward the province
#                  (rho_zone ~ beta(2, 2) per week) instead of the spline; forecasts
#                  continue both AR(1)s, so no separate forecast rule
#   F1_mean        F1's fits, forecasts following both AR(1) means (generated quantities only)
#
# Writes <output root>/<variant>/<origin>/: forecast_samples.csv.gz (zone samples for the
# forecast weeks), diagnostics.txt, fit.rds (not for E0_flat).

args <- commandArgs(trailingOnly = TRUE)
variant <- args[1]
origins <- if (length(args) > 1) as.Date(args[-1]) else
  seq(as.Date("2026-05-09"), as.Date("2026-08-01"), by = 7)

variants <- list(
  E0 = list(),
  E1_zone_tight = list(rho_zone_prior = c(8, 2), sd_ar_zone_scale = 0.15),
  E2_prov_tight = list(sd_rw_R_scale = 0.1),
  E0_flat = list(forecast_R = "flat", from = "E0")
)
variants$E3_ar <- list(ar_R = TRUE)
variants$E3_ar_mean <- list(ar_R = TRUE, forecast_R = "mean", from = "E3_ar")
variants$F1 <- list(ar_R = TRUE, zone_weekly = TRUE, rho_zone_prior = c(2, 2))
variants$F1_mean <- list(ar_R = TRUE, zone_weekly = TRUE, rho_zone_prior = c(2, 2),
                         forecast_R = "mean", from = "F1")
# Damped forecast R from E0's fits: E0_damp080 multiplies each zone's log R by 0.8 per
# forecast week (generated quantities only).
for (phi in c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.3, 0)) {
  variants[[sprintf("E0_damp%03d", round(100 * phi))]] <-
    list(forecast_R = "damped", damp_R = phi, from = "E0")
}
if (!variant %in% names(variants)) stop("Unknown variant ", variant)
# The variants were run when the default forecast was the random walk; keep it unless
# a variant says otherwise.
v <- modifyList(list(forecast_R = "rw"), variants[[variant]])
provinces <- c("CD-IT", "CD-NK", "CD-HU", "CD-TO", "CD-BU")
chains <- as.integer(Sys.getenv("BDBV_CHAINS", "4"))
iter <- as.integer(Sys.getenv("BDBV_ITER", "300"))
treedepth <- as.integer(Sys.getenv("BDBV_TREEDEPTH", "9"))
n_samples <- 400
out_root <- Sys.getenv("BDBV_FINAL_OUT", "output/final")

suppressWarnings(suppressPackageStartupMessages({
  source("R/bdbv_data.R")
  source("R/bdbv_fit.R")
}))

for (o in as.list(origins)) {
  out_dir <- file.path(out_root, variant, o)
  out_file <- file.path(out_dir, "forecast_samples.csv.gz")
  if (file.exists(out_file)) next
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cfg <- do.call(bdbv_config, c(list(as_of = NULL, finalized = TRUE, origin = o,
                                     provinces = provinces), v[setdiff(names(v), "from")]))
  d <- build_bdbv_data(cfg)
  model <- bdbv_model(delay = "lognormal")
  t0 <- Sys.time()
  if (!is.null(v$from)) {
    fit <- readRDS(file.path(out_root, v$from, o, "fit.rds"))
  } else {
    message(sprintf("%s %s: %d zones, %d weeks", variant, o, d$stan_data$Z, d$stan_data$W))
    fit <- fit_bdbv(d, model, chains = chains, iter_warmup = iter, iter_sampling = iter,
                    max_treedepth = treedepth, output_dir = out_dir)
    fit$save_object(file.path(out_dir, "fit.rds"))
  }
  minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  samples <- generate_samples(fit, d, model,
                              n_samples = min(n_samples, posterior::ndraws(fit$draws("lp__"))))
  readr::write_csv(filter(samples, target_end_date %in% d$forecast_weeks), out_file)

  dg <- diagnostics_table(fit)
  sink(file.path(out_dir, "diagnostics.txt"))
  cat("variant:", variant, " origin:", as.character(o), " fit minutes:", round(minutes, 1), "\n")
  cat("divergences:", dg$divergences, " max treedepth hits:", dg$max_treedepth_hits,
      " E-BFMI:", round(dg$ebfmi, 2), "\n\n")
  print(dg$summary, n = 50, width = 120)
  sink()
  # The draws are in fit.rds; drop cmdstan's CSVs.
  unlink(list.files(out_dir, pattern = "^bdbv_joint-.*\\.csv$", full.names = TRUE))
  message(sprintf("%s %s done (%.1f min fit, %d divergences)", variant, o, minutes,
                  dg$divergences))
}
