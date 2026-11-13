# Compare backtest scores of model variants, each scored by backtest_bdbv.R.
#
#   Rscript scripts/compare_backtests.R output/bdbv_v1 output/bdbv output/bdbv_zrw
#
# Scores only the location-weeks every variant has. Zone CRPS is summed over zones
# (comparable across variants; a mean would be dominated by zones with no cases).
# Writes output/backtest_comparison.csv and output/backtest_comparison.png.

roots <- commandArgs(trailingOnly = TRUE)
if (length(roots) < 2) stop("Pass two or more output roots")

suppressWarnings(suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
}))

scores <- bind_rows(lapply(roots, function(r) {
  s <- read_csv(file.path(r, "backtest.csv"), show_col_types = FALSE)
  if (!"level" %in% names(s)) {
    s$level <- ifelse(s$location == "CD", "country",
                      ifelse(nchar(s$location) == 5, "province", "zone"))
  }
  mutate(s, model = basename(r))
})) |>
  group_by(as_of, location, target_end_date) |>
  filter(n_distinct(model) == length(roots)) |>
  ungroup()

summary <- scores |>
  group_by(model, level, as_of) |>
  summarise(n = n(), crps = if (first(level) == "zone") sum(crps) else mean(crps),
            cover50 = mean(in50), cover90 = mean(in90),
            mae_model = mean(ae_model), mae_naive = mean(ae_naive), .groups = "drop") |>
  mutate(model = factor(model, basename(roots))) |>
  arrange(level, as_of, model)
write_csv(summary, "output/backtest_comparison.csv")
cat("CRPS: mean over locations (country, province) or summed over zones; lower is better.\n\n")
print(summary, n = 100, width = 120)

# Final-data scenarios (backtest_bdbv.R), when every root has them with the same
# growth reference.
scen_files <- file.path(roots, "backtest_scenarios.csv")
if (all(file.exists(scen_files))) {
  scen <- bind_rows(lapply(seq_along(roots), function(i) {
    read_csv(scen_files[i], show_col_types = FALSE) |> mutate(model = basename(roots[i]))
  }))
  if (n_distinct(scen$growth_ref) == 1) {
    scen_summary <- scen |>
      semi_join(scores, by = c("as_of", "location", "target_end_date")) |>
      group_by(model, level, as_of, s) |>
      summarise(crps = if (first(level) == "zone") sum(crps) else mean(crps), .groups = "drop") |>
      tidyr::pivot_wider(names_from = model, values_from = crps) |>
      arrange(level, as_of, s)
    write_csv(scen_summary, "output/backtest_comparison_scenarios.csv")
    cat("\nCRPS against truth = latest x (1 + s x growth expected by", scen$growth_ref[1], "):\n")
    print(scen_summary, n = 100, width = 120)
  } else {
    message("Scenario files use different growth references; skipping scenario comparison.")
  }
}

# Figure: each variant's 90% interval and median for the country and the four largest
# provinces, against the data each run was fit to and the latest release.
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
vint <- read_csv("data/raw/time-series.csv", show_col_types = FALSE) |>
  mutate(as_of = as.Date(as_of))
latest <- max(vint$as_of)
locs <- c("CD", "CD-IT", "CD-NK", "CD-HU", "CD-TO")
runs <- sort(unique(scores$as_of))
bands <- bind_rows(lapply(roots, function(r) bind_rows(lapply(runs, function(a) {
  s <- read_csv(file.path(r, a, "samples_all_weeks.csv.gz"), show_col_types = FALSE) |>
    filter(target_end_date >= as.Date(a) - 7 * 8)
  bind_rows(
    s |> group_by(target_end_date, sample) |> summarise(value = sum(value), .groups = "drop") |>
      mutate(location = "CD"),
    s |> mutate(location = substr(location, 1, 5)) |>
      group_by(location, target_end_date, sample) |> summarise(value = sum(value), .groups = "drop")
  ) |>
    filter(location %in% locs) |>
    group_by(location, target_end_date) |>
    summarise(lo = quantile(value, 0.05), mid = median(value), hi = quantile(value, 0.95),
              .groups = "drop") |>
    mutate(run = as.character(a), model = basename(r))
}))))
obs <- vint |>
  filter(location %in% locs, as_of %in% c(runs, latest)) |>
  cross_join(tibble(run = as.character(runs))) |>
  filter(as_of == as.Date(run) | as_of == latest,
         target_end_date >= as.Date(run) - 7 * 8, target_end_date <= as.Date(run) + 7 * 4) |>
  mutate(data = ifelse(as_of == latest, paste("Release", latest), "Release fit to"))
# Cap bands at 3x the largest observed count per location so wide forecast tails don't
# flatten the panels.
caps <- obs |> group_by(location) |> summarise(cap = 3 * max(observation))
bands <- bands |> left_join(caps, by = "location") |> mutate(hi = pmin(hi, cap))
p <- ggplot(bands, aes(target_end_date)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = model), alpha = 0.15) +
  geom_line(aes(y = mid, colour = model)) +
  geom_point(data = obs, aes(y = observation, shape = data), size = 1) +
  geom_vline(aes(xintercept = as.Date(run)), linetype = "dashed", colour = "grey40") +
  scale_shape_manual(values = c(16, 1)) +
  facet_grid(location ~ run, scales = "free_y",
             labeller = labeller(run = function(x) paste("Fit as of", x))) +
  labs(x = "Week of symptom onset (ending)", y = "Confirmed cases", colour = "Model",
       fill = "Model", shape = NULL,
       title = "Backtest: eventual-count forecasts (median, 90%) by model",
       subtitle = sprintf("Points: data the run was fit to (open) and release %s (filled), itself incomplete in recent weeks.\nBands are cut at 3x the largest observed count.", latest)) +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave("output/backtest_comparison.png", p, width = 11, height = 12, dpi = 120)
cat("\nWrote output/backtest_comparison.csv and output/backtest_comparison.png\n")
