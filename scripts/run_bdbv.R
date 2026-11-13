# Fit the joint BDBV nowcast/forecast model to one data vintage and write hub-format
# samples, diagnostics, and plots.
#
#   Rscript scripts/run_bdbv.R                 # latest vintage in data/raw/time-series.csv
#   Rscript scripts/run_bdbv.R 2026-10-24      # an earlier vintage (backtest)
#   Rscript scripts/run_bdbv.R 2026-11-01 quick
#
# "quick" runs 2 chains x 200/200 iterations to check the pipeline end to end; its
# posterior is not reliable. Run `uv run python -m fcast.bdbv_inputs` first.
#
# Environment variables:
#   BDBV_TEAM_MODEL   hub team-model name (default NAU-jointseir)
#   BDBV_OUT          output root (default output/bdbv), for comparing model variants
#   BDBV_TRUTH        backtests: a later release; also writes samples_truth.csv.gz, the
#                     predicted counts as published in that release
#   BDBV_RESUME       1 = reuse <output root>/<as_of>/fit.rds instead of refitting
#   BDBV_CHAINS, BDBV_ITER  chains and warmup/sampling iterations each (default 4, 500)
#   BDBV_TREEDEPTH    maximum tree depth (default 10)
#   BDBV_DELAY        reporting delay: lognormal (default) or hazard
#
# Writes to <output root>/<as_of>/:
#   <reference_date>-<team_model>.csv   hub sample file (all 546 locations, 9 weeks)
#   samples_all_weeks.csv.gz            modeled-zone samples for every model week
#   diagnostics.txt, fit.rds, *.png

args <- commandArgs(trailingOnly = TRUE)
as_of <- if (length(args) >= 1 && args[1] != "latest") args[1] else NULL
quick <- "quick" %in% args
team_model <- Sys.getenv("BDBV_TEAM_MODEL", "NAU-jointseir")
out_root <- Sys.getenv("BDBV_OUT", "output/bdbv")
truth_as_of <- if (nzchar(Sys.getenv("BDBV_TRUTH"))) Sys.getenv("BDBV_TRUTH") else NULL
chains <- as.integer(Sys.getenv("BDBV_CHAINS", "4"))
iter <- as.integer(Sys.getenv("BDBV_ITER", "500"))
resume <- Sys.getenv("BDBV_RESUME", "0") == "1"

suppressWarnings(suppressPackageStartupMessages({
  source("R/bdbv_data.R")
  source("R/bdbv_fit.R")
  source("R/bdbv_plots.R")
  library(patchwork)
}))

d <- build_bdbv_data(bdbv_config(as_of = as_of, truth_as_of = truth_as_of,
                                 delay = Sys.getenv("BDBV_DELAY", "lognormal")))
out_dir <- file.path(out_root, as.character(d$as_of))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message(sprintf(
  "as_of %s, reference date %s: %d zones in %s, %d weeks + %d forecast weeks, %d releases, %d triangle pairs",
  d$as_of, d$ref_date, d$stan_data$Z, paste(d$provinces, collapse = "/"),
  d$stan_data$W, d$stan_data$H, d$stan_data$V, d$stan_data$N_tri
))

model <- bdbv_model(delay = d$cfg$delay)
fit <- if (resume) {
  readRDS(file.path(out_dir, "fit.rds"))
} else if (quick) {
  fit_bdbv(d, model, chains = 2, iter_warmup = 200, iter_sampling = 200,
           output_dir = out_dir)
} else {
  fit_bdbv(d, model, chains = chains, iter_warmup = iter, iter_sampling = iter,
           max_treedepth = as.integer(Sys.getenv("BDBV_TREEDEPTH", "10")),
           output_dir = out_dir)
}
if (!resume) fit$save_object(file.path(out_dir, "fit.rds"))

dg <- diagnostics_table(fit, d$cfg$delay)
sink(file.path(out_dir, "diagnostics.txt"))
cat("as_of:", as.character(d$as_of), " reference_date:", as.character(d$ref_date), "\n")
cat("chains:", fit$num_chains(), " draws/chain:", fit$metadata()$iter_sampling, "\n")
cat("divergences:", dg$divergences, " max treedepth hits:", dg$max_treedepth_hits, "\n")
cat("E-BFMI:", round(dg$ebfmi, 2), "\n\n")
print(dg$summary, n = 50, width = 120)
all_summary <- fit$summary()
cat("\nWorst R-hat / bulk ESS over all parameters:\n")
print(all_summary |> arrange(desc(rhat)) |> head(5), width = 120)
print(all_summary |> arrange(ess_bulk) |> head(5), width = 120)
sink()

samples <- generate_samples(fit, d, model)
readr::write_csv(samples, file.path(out_dir, "samples_all_weeks.csv.gz"))
if (!is.null(truth_as_of)) {
  readr::write_csv(attr(samples, "truth"), file.path(out_dir, "samples_truth.csv.gz"))
}
hub <- hub_samples(samples, d)
hub_file <- file.path(out_dir, sprintf("%s-%s.csv", d$ref_date, team_model))
readr::write_csv(hub, hub_file)

agg <- aggregate_samples(samples, d)
top_zones <- d$zones$location[order(-rowSums(d$stan_data$cases))][1:12]
ggsave(file.path(out_dir, "fit_country_provinces.png"),
       plot_fit(agg, d, c("CD", d$provinces)), width = 10, height = 6, dpi = 120)
ggsave(file.path(out_dir, "fit_top_zones.png"),
       plot_fit(samples, d, top_zones), width = 12, height = 8, dpi = 120)
ggsave(file.path(out_dir, "reporting.png"), plot_reporting(fit, d),
       width = 8, height = 6, dpi = 120)
ggsave(file.path(out_dir, "R_provinces.png"), plot_R(fit, d), width = 8, height = 4, dpi = 120)

message("Wrote ", hub_file, " (", nrow(hub), " rows) and plots to ", out_dir)
