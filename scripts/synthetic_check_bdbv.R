# Synthetic recovery check: simulate data from the model at the posterior median of a
# real fit, refit, and report how often the true values fall in the posterior intervals.
#
#   Rscript scripts/synthetic_check_bdbv.R output/bdbv/2026-11-01
#
# Writes recovery.txt and synthetic_fit.png to <run dir>/synthetic/.

args <- commandArgs(trailingOnly = TRUE)
run_dir <- if (length(args)) args[1] else stop("Pass a run directory from run_bdbv.R")

suppressWarnings(suppressPackageStartupMessages({
  source("R/bdbv_data.R")
  source("R/bdbv_fit.R")
  source("R/bdbv_plots.R")
  source("R/bdbv_synthetic.R")
}))

real_fit <- readRDS(file.path(run_dir, "fit.rds"))
as_of <- basename(normalizePath(run_dir))
d <- build_bdbv_data(bdbv_config(as_of = as_of))
model <- bdbv_model()

# Posterior median of every parameter as the truth. Medians of the random-walk
# innovations are consistent with each other, so this is a plausible single trajectory.
draws <- real_fit$draws(format = "draws_matrix")
truth <- as_draws_matrix(t(apply(draws, 2, median)))
sim <- simulate_bdbv(d, truth, model)

out_dir <- file.path(run_dir, "synthetic")
dir.create(out_dir, showWarnings = FALSE)
fit <- fit_bdbv(sim, model, output_dir = out_dir)

vars <- c("log_R0", "sd_prov_R", "sd_zone_R", "sd_rw_R", "phi_move", "inv_sqrt_nb",
          "delay_mu", "delay_sigma", "pi_late", "late_m", "late_s", "log_R_prov", "delay_m", "F_week")
rec <- recovery_table(fit, truth, vars)
samples <- generate_samples(fit, sim, model)
now <- nowcast_recovery(samples, sim, d)
dg <- diagnostics_table(fit)

sink(file.path(out_dir, "recovery.txt"))
cat("Synthetic recovery check from", run_dir, "\n")
cat("divergences:", dg$divergences, " max treedepth hits:", dg$max_treedepth_hits, "\n\n")
cat("Coverage of true values by posterior intervals (expect ~0.5 and ~0.9):\n")
print(rec, n = 50)
cat("\nProvince nowcasts for the last 6 weeks (true = simulated final count):\n")
print(now, n = 60, width = 120)
cat("\n90% coverage of province final counts:", mean(now$in90), "\n")
sink()

# Plot the synthetic final counts as if they were a single "true" vintage.
truth_long <- as.data.frame(sim$truth$final) |>
  setNames(as.character(d$weeks)) |>
  mutate(location = substr(d$zones$location, 1, 5)) |>
  tidyr::pivot_longer(-location, names_to = "target_end_date", values_to = "observation") |>
  mutate(target_end_date = as.Date(target_end_date)) |>
  group_by(location, target_end_date) |>
  summarise(observation = sum(observation), .groups = "drop")
rep_long <- as.data.frame(sim$stan_data$cases) |>
  setNames(as.character(d$weeks)) |>
  mutate(location = substr(d$zones$location, 1, 5)) |>
  tidyr::pivot_longer(-location, names_to = "target_end_date", values_to = "observation") |>
  mutate(target_end_date = as.Date(target_end_date)) |>
  group_by(location, target_end_date) |>
  summarise(observation = sum(observation), .groups = "drop")
sim_d <- d
sim_d$vintages <- bind_rows(
  mutate(truth_long, as_of = d$as_of + 365),
  mutate(rep_long, as_of = d$as_of)
)
sim_d$as_of <- d$as_of + 365
p <- plot_fit(aggregate_samples(samples, d), sim_d, d$provinces) +
  labs(title = "Synthetic check: samples vs simulated truth (dark) and simulated reported data (light)") +
  coord_cartesian(xlim = range(c(d$weeks, d$forecast_weeks)))
ggsave(file.path(out_dir, "synthetic_fit.png"), p, width = 10, height = 6, dpi = 120)
message("Wrote ", out_dir)
