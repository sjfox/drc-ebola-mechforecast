# Score earlier-vintage runs against the latest data vintage.
#
#   Rscript scripts/backtest_bdbv.R 2026-10-17 2026-10-24
#   BDBV_OUT=output/bdbv_zrw Rscript scripts/backtest_bdbv.R 2026-10-17 2026-10-24
#
# For each as_of given (each must already have <root>/<as_of>/ from run_bdbv.R; root is
# BDBV_OUT, default output/bdbv), compares its final-count samples with the latest
# vintage in time-series.csv, for the country, provinces, and health zones in modeled
# provinces, for weeks that are at least `min_truth_age` days old in the latest vintage
# (so the "truth" is mostly complete). Reports CRPS, 50%/90% coverage, and the absolute
# error of the sample median against the naive nowcast (reported counts at face value).
#
# Writes <root>/backtest.csv, backtest_provinces.png and backtest_zones.png (plus
# backtest_published.csv/.png for runs fit with BDBV_TRUTH), and prints
# a summary. scripts/compare_backtests.R compares roots.
#
# Final-data scenarios: the latest release is itself incomplete for recent weeks, which
# favours models that expect little late reporting. backtest_scenarios.csv rescores
# against truth = latest count x (1 + s x g), s = 0, 0.25, ..., 1, where g is the growth
# after the latest release that a reference run expects for each province and week
# (mean eventual / mean count as published in the latest release, from its
# samples_truth.csv.gz). The reference is BDBV_GROWTH_REF (an output root; its run for
# the last as_of given is used), default this root when it has samples_truth. Use the
# same reference for every root you compare.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Pass one or more as_of dates with existing runs")
min_truth_age <- 14
root <- Sys.getenv("BDBV_OUT", "output/bdbv")

suppressWarnings(suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  source("R/bdbv_plots.R")
}))

week_end_of <- function(date) date + (6 - as.integer(format(date, "%w")))

# Same rule as R/bdbv_data.R: the week of the first Sunday on or after as_of.
reference_date_of <- function(as_of) {
  week_end_of(as_of + (7 - as.integer(format(as_of, "%w"))) %% 7)
}

crps_sample <- function(x, y) {
  mean(abs(x - y)) - 0.5 * mean(abs(outer(x, x, "-")))
}

vint <- read_csv("data/raw/time-series.csv", show_col_types = FALSE) |>
  mutate(as_of = as.Date(as_of), target_end_date = as.Date(target_end_date))
latest <- max(vint$as_of)
truth <- vint |>
  filter(as_of == latest, target_end_date <= latest - min_truth_age) |>
  select(location, target_end_date, truth = observation)

score_run <- function(a, truth = get("truth", envir = globalenv())) {
  a <- as.Date(a)
  dir <- file.path(root, a)
  hub <- list.files(dir, pattern = "^\\d{4}-\\d{2}-\\d{2}-.*\\.csv$", full.names = TRUE)
  if (!length(hub)) stop("No hub file in ", dir)
  naive <- vint |>
    filter(as_of == a) |>
    select(location, target_end_date, reported = observation)
  modeled_provs <- read_csv(file.path(dir, "samples_all_weeks.csv.gz"),
                            show_col_types = FALSE, col_select = "location") |>
    pull(location) |> substr(1, 5) |> unique()
  read_csv(hub[1], show_col_types = FALSE) |>
    filter(substr(location, 1, 5) %in% c("CD", modeled_provs)) |>
    inner_join(truth, by = c("location", "target_end_date")) |>
    group_by(location, target_end_date, truth) |>
    summarise(
      crps = crps_sample(value, truth[1]),
      median = median(value),
      in50 = truth[1] >= quantile(value, 0.25) & truth[1] <= quantile(value, 0.75),
      in90 = truth[1] >= quantile(value, 0.05) & truth[1] <= quantile(value, 0.95),
      .groups = "drop"
    ) |>
    left_join(naive, by = c("location", "target_end_date")) |>
    mutate(
      as_of = a,
      level = case_when(location == "CD" ~ "country", nchar(location) == 5 ~ "province",
                        TRUE ~ "zone"),
      weeks_from_ref = as.integer(target_end_date - reference_date_of(a)) / 7,
      ae_model = abs(median - truth),
      ae_naive = abs(coalesce(reported, 0) - truth)
    )
}

scores <- bind_rows(lapply(args, score_run))
write_csv(scores, file.path(root, "backtest.csv"))

growth_ref <- Sys.getenv("BDBV_GROWTH_REF", root)
ref_dir <- file.path(growth_ref, max(as.Date(args)))
if (file.exists(file.path(ref_dir, "samples_truth.csv.gz"))) {
  prov_mean <- function(f) {
    read_csv(f, show_col_types = FALSE) |>
      mutate(prov = substr(location, 1, 5)) |>
      group_by(prov, target_end_date, sample) |> summarise(value = sum(value), .groups = "drop") |>
      group_by(prov, target_end_date) |> summarise(value = mean(value), .groups = "drop")
  }
  growth <- inner_join(prov_mean(file.path(ref_dir, "samples_all_weeks.csv.gz")),
                       prov_mean(file.path(ref_dir, "samples_truth.csv.gz")),
                       by = c("prov", "target_end_date"), suffix = c("_final", "_pub")) |>
    transmute(prov, target_end_date, g = pmax(value_final / value_pub - 1, 0))
  zone_truth <- truth |> filter(nchar(location) > 5) |>
    mutate(prov = substr(location, 1, 5)) |>
    left_join(growth, by = c("prov", "target_end_date")) |>
    mutate(g = coalesce(g, 0))
  scenarios <- bind_rows(lapply(c(0, 0.25, 0.5, 0.75, 1), function(sc) {
    z <- mutate(zone_truth, truth = truth * (1 + sc * g))
    truth_s <- bind_rows(
      z |> group_by(target_end_date) |> summarise(truth = sum(truth), .groups = "drop") |>
        mutate(location = "CD"),
      z |> group_by(location = prov, target_end_date) |> summarise(truth = sum(truth), .groups = "drop"),
      select(z, location, target_end_date, truth)
    )
    bind_rows(lapply(args, score_run, truth = truth_s)) |>
      select(location, target_end_date, as_of, level, truth, crps, in90) |>
      mutate(s = sc)
  })) |>
    mutate(growth_ref = growth_ref)
  write_csv(scenarios, file.path(root, "backtest_scenarios.csv"))
  cat("Final-data scenarios: truth = latest x (1 + s x growth expected by", growth_ref, ")\n")
  scenarios |>
    group_by(as_of, level, s) |>
    summarise(crps = if (level[1] == "zone") sum(crps) else mean(crps), cover90 = mean(in90),
              .groups = "drop") |>
    print(n = 40)
} else {
  message("No samples_truth.csv.gz in ", ref_dir, "; skipping final-data scenarios ",
          "(set BDBV_GROWTH_REF to a root with BDBV_TRUTH runs).")
}

# Runs fit with BDBV_TRUTH set also predict each count as published in the latest
# release (samples_truth.csv.gz). Scoring those against it is fair for every week,
# including recent weeks that are still incomplete in the latest release. Zones outside
# the modeled provinces keep their counts at as_of, as in the hub file.
score_published <- function(a) {
  a <- as.Date(a)
  f <- file.path(root, a, "samples_truth.csv.gz")
  if (!file.exists(f)) return(NULL)
  zone <- read_csv(f, show_col_types = FALSE) |>
    filter(target_end_date >= reference_date_of(a) - 7 * 5, target_end_date <= latest)
  other <- vint |>
    filter(as_of == a, location_type == "health zone",
           !location %in% unique(zone$location), target_end_date %in% zone$target_end_date) |>
    group_by(target_end_date) |> summarise(other = sum(observation), .groups = "drop")
  samples <- bind_rows(
    zone |> group_by(target_end_date, sample) |> summarise(value = sum(value), .groups = "drop") |>
      left_join(other, by = "target_end_date") |>
      mutate(value = value + coalesce(other, 0), location = "CD") |> select(-other),
    zone |> mutate(location = substr(location, 1, 5)) |>
      group_by(location, target_end_date, sample) |> summarise(value = sum(value), .groups = "drop"),
    zone
  )
  truth_all <- vint |> filter(as_of == latest) |>
    select(location, target_end_date, truth = observation)
  samples |>
    inner_join(truth_all, by = c("location", "target_end_date")) |>
    group_by(location, target_end_date, truth) |>
    summarise(
      crps = crps_sample(value, truth[1]), median = median(value),
      in50 = truth[1] >= quantile(value, 0.25) & truth[1] <= quantile(value, 0.75),
      in90 = truth[1] >= quantile(value, 0.05) & truth[1] <= quantile(value, 0.95),
      .groups = "drop"
    ) |>
    mutate(as_of = a,
           level = case_when(location == "CD" ~ "country", nchar(location) == 5 ~ "province",
                             TRUE ~ "zone"),
           truth_age = as.integer(latest - target_end_date))
}
published <- bind_rows(lapply(args, score_published))
if (nrow(published)) {
  write_csv(published, file.path(root, "backtest_published.csv"))
  cat("Counts as published in", as.character(latest), "(all weeks):\n")
  published |>
    group_by(as_of, level) |>
    summarise(n = n(), crps = if (level[1] == "zone") sum(crps) else mean(crps),
              cover50 = mean(in50), cover90 = mean(in90), .groups = "drop") |>
    print()
}

cat("Truth: vintage", as.character(latest), "weeks at least", min_truth_age,
    "days old.\n\n")
scores |>
  filter(level != "zone") |>
  group_by(as_of, weeks_from_ref) |>
  summarise(n = n(), crps = mean(crps), cover50 = mean(in50), cover90 = mean(in90),
            mae_model = mean(ae_model), mae_naive = mean(ae_naive), .groups = "drop") |>
  print(n = 50)
cat("\nCountry only:\n")
scores |>
  filter(location == "CD") |>
  select(as_of, target_end_date, truth, reported, median, in90, crps) |>
  print(n = 50)

runs <- setNames(lapply(args, function(a) {
  read_csv(file.path(root, a, "samples_all_weeks.csv.gz"), show_col_types = FALSE)
}), args)
provs <- c("CD", "CD-IT", "CD-NK", "CD-HU", "CD-TO")
ggsave(file.path(root, "backtest_provinces.png"), plot_backtest(runs, vint, provs),
       width = 4 + 3.5 * length(args), height = 10, dpi = 120)
top_zones <- vint |>
  filter(as_of == latest, location_type == "health zone",
         target_end_date > latest - 7 * 8) |>
  group_by(location) |> summarise(n = sum(observation)) |>
  slice_max(n, n = 6) |> pull(location)
ggsave(file.path(root, "backtest_zones.png"), plot_backtest(runs, vint, top_zones),
       width = 4 + 3.5 * length(args), height = 11, dpi = 120)

# Figure: predicted counts as published in the latest release (median, 50%/90%) against
# that release, for the country and the four largest provinces.
if (nrow(published)) {
  locs <- c("CD", "CD-IT", "CD-NK", "CD-HU", "CD-TO")
  bands <- bind_rows(lapply(args, function(a) {
    zone <- read_csv(file.path(root, a, "samples_truth.csv.gz"), show_col_types = FALSE) |>
      filter(target_end_date >= as.Date(a) - 7 * 8, target_end_date <= latest)
    bind_rows(
      zone |> group_by(target_end_date, sample) |> summarise(value = sum(value), .groups = "drop") |>
        mutate(location = "CD"),
      zone |> mutate(location = substr(location, 1, 5)) |>
        group_by(location, target_end_date, sample) |> summarise(value = sum(value), .groups = "drop")
    ) |>
      filter(location %in% locs) |>
      sample_bands() |>
      mutate(run = a)
  }))
  obs <- vint |>
    filter(location %in% locs, as_of %in% c(as.Date(args), latest)) |>
    cross_join(tibble(run = args)) |>
    filter(as_of == as.Date(run) | as_of == latest,
           target_end_date >= as.Date(run) - 7 * 8, target_end_date <= latest) |>
    mutate(data = ifelse(as_of == latest, paste("Release", latest), "Release fit to"))
  p <- ggplot(bands, aes(target_end_date)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = "#2b6cb0", alpha = 0.2) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = "#2b6cb0", alpha = 0.35) +
    geom_line(aes(y = median), colour = "#2b6cb0") +
    geom_point(data = obs, aes(y = observation, shape = data), size = 1.2) +
    geom_vline(aes(xintercept = as.Date(run)), linetype = "dashed", colour = "grey40") +
    scale_shape_manual(values = c(16, 1)) +
    facet_grid(location ~ run, scales = "free_y",
               labeller = labeller(run = function(x) paste("Fit as of", x))) +
    labs(x = "Week of symptom onset (ending)", y = "Confirmed cases", shape = NULL,
         title = sprintf("Backtest: predicted counts as published in release %s (50%%/90%%)", latest),
         subtitle = "Includes cases from the cut-off week onward that are not yet published.") +
    theme_minimal(base_size = 10) + theme(legend.position = "bottom")
  ggsave(file.path(root, "backtest_published.png"), p, width = 9, height = 11, dpi = 120)
}
