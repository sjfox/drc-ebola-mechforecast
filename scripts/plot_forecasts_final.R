# Retrospective forecasts from the finalized-data backtests (scripts/backtest_final.R)
# against the finalized counts, for the country, provinces and every health zone, and
# their WIS against the hub flatline baseline (from scripts/score_final.R's scores.csv).
#
#   Rscript scripts/plot_forecasts_final.R              # E0_damp050 (the default model)
#   Rscript scripts/plot_forecasts_final.R E3_ar_mean
#
# Writes to <root>/review/<variant>/ (root: BDBV_FINAL_OUT, default output/final):
#   forecasts_country_province.png  each origin's 1-4 week forecast (median, 50%, 90%)
#   forecasts_zones.pdf             the same for every zone, 20 per page, by province and
#                                   cases (also forecasts_zones_<page>.png)
#   rel_wis_horizon.png, rel_wis_origin.png, rel_wis_zones.png   WIS vs the flatline
#   location_scores.csv, level_scores.csv   WIS, relative WIS and coverage

args <- commandArgs(trailingOnly = TRUE)
variant <- if (length(args)) args[1] else "E0_damp050"
root <- Sys.getenv("BDBV_FINAL_OUT", "output/final")
out <- file.path(root, "review", variant)
provinces <- c("CD-IT", "CD-NK", "CD-HU", "CD-TO", "CD-BU")
prov_names <- c("CD-IT" = "Ituri", "CD-NK" = "Nord-Kivu", "CD-HU" = "Haut-Uele",
                "CD-TO" = "Tshopo", "CD-BU" = "Bas-Uele")
truth_end <- as.Date("2026-08-29")   # weeks at least 60 days old in the 11-01 release
start <- as.Date("2026-04-04")
per_page <- 20

suppressWarnings(suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2)
}))
dir.create(out, recursive = TRUE, showWarnings = FALSE)

zones <- read_csv("data/processed/bdbv_zones.csv", show_col_types = FALSE) |>
  filter(substr(location, 1, 5) %in% provinces)
vint <- read_csv("data/raw/time-series.csv", show_col_types = FALSE) |>
  mutate(as_of = as.Date(as_of), target_end_date = as.Date(target_end_date))
truth <- vint |>
  filter(as_of == max(as_of), target_end_date >= start, target_end_date <= truth_end,
         location %in% c("CD", provinces, zones$location)) |>
  select(location, target_end_date, observation)
# Zones missing from the release have no cases.
truth <- truth |>
  complete(location = c("CD", provinces, zones$location),
           target_end_date = seq(start, truth_end, by = 7), fill = list(observation = 0))

# Forecast quantiles for every origin, location and horizon (zone samples summed up).
origins <- sort(as.Date(basename(list.dirs(file.path(root, variant), recursive = FALSE))))
origins <- origins[file.exists(file.path(root, variant, origins, "forecast_samples.csv.gz"))]
fc <- bind_rows(lapply(origins, function(o) {
  z <- read_csv(file.path(root, variant, o, "forecast_samples.csv.gz"), show_col_types = FALSE) |>
    mutate(target_end_date = as.Date(target_end_date))
  bind_rows(
    z |> group_by(target_end_date, sample) |> summarise(value = sum(value), .groups = "drop") |>
      mutate(location = "CD"),
    z |> mutate(location = substr(location, 1, 5)) |>
      group_by(location, target_end_date, sample) |> summarise(value = sum(value), .groups = "drop"),
    z
  ) |>
    group_by(location, target_end_date) |>
    summarise(lo90 = quantile(value, 0.05), lo50 = quantile(value, 0.25), med = median(value),
              hi50 = quantile(value, 0.75), hi90 = quantile(value, 0.95), .groups = "drop") |>
    mutate(origin = o, horizon = as.integer(target_end_date - o) / 7)
}))

# Per-location scores against the flatline baseline.
sc <- read_csv(file.path(root, "scores.csv"), show_col_types = FALSE) |>
  filter(model %in% c(variant, "flatline"), origin %in% origins) |>
  group_by(location, model) |> summarise(wis = mean(wis), .groups = "drop") |>
  pivot_wider(names_from = model, values_from = wis) |>
  mutate(rel_wis = .data[[variant]] / flatline)
write_csv(sc, file.path(out, "location_scores.csv"))

cases <- truth |> group_by(location) |> summarise(cases = sum(observation), .groups = "drop")
name_of <- function(loc) {
  nm <- zones$location_name[match(loc, zones$location)]
  ifelse(loc == "CD", "DRC (5 provinces)", ifelse(is.na(nm), prov_names[loc], nm))
}
labels <- cases |>
  left_join(sc, by = "location") |>
  mutate(label = sprintf("%s %s\n%s cases%s", location, name_of(location),
                         format(cases, big.mark = ","),
                         ifelse(is.na(rel_wis) | !is.finite(rel_wis), "",
                                sprintf(", rel. WIS %.2f", rel_wis))))
lab <- setNames(labels$label, labels$location)

log1p_y <- scale_y_continuous(trans = "log1p", breaks = c(0, 1, 3, 10, 30, 100, 300, 1000, 3000),
                              labels = scales::comma)
base_theme <- theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 7.5, hjust = 0),
        legend.position = "bottom")
truth_layers <- function(d) list(
  geom_line(data = d, aes(y = observation), colour = "grey20", linewidth = 0.3),
  geom_point(data = d, aes(y = observation), colour = "grey20", size = 0.8)
)

# Each origin's 1-4 week forecast, starting at its last observed count.
ramp <- colorRampPalette(c("#a9c9f0", "#2a78d6", "#0b2f6b"))(length(origins))
trajectory_plot <- function(locs, title, ncol = 5) {
  anchor <- truth |> filter(location %in% locs) |>
    inner_join(tibble(origin = origins), by = c(target_end_date = "origin")) |>
    transmute(location, target_end_date, origin = target_end_date,
              lo90 = observation, lo50 = observation, med = observation,
              hi50 = observation, hi90 = observation)
  d <- bind_rows(fc |> filter(location %in% locs), anchor) |>
    mutate(location = factor(location, locs), origin = factor(origin, origins))
  tr <- truth |> filter(location %in% locs) |> mutate(location = factor(location, locs))
  ggplot(d, aes(target_end_date)) +
    geom_ribbon(aes(ymin = lo90, ymax = hi90, fill = origin, group = origin), alpha = 0.08) +
    geom_ribbon(aes(ymin = lo50, ymax = hi50, fill = origin, group = origin), alpha = 0.16) +
    geom_line(aes(y = med, colour = origin, group = origin), linewidth = 0.55) +
    truth_layers(tr) +
    scale_colour_manual(values = ramp, aesthetics = c("colour", "fill"), drop = FALSE) +
    log1p_y +
    facet_wrap(~location, ncol = ncol, scales = "free_y", labeller = as_labeller(lab),
               drop = FALSE) +
    labs(x = "Week of symptom onset (ending)", y = "Confirmed cases (log(1 + count) scale)",
         colour = "Forecast date", fill = "Forecast date", title = title,
         subtitle = paste("Each forecast date's 1-4 week forecast: median line from the last",
                          "observed count, 50% and 90% intervals shaded. Points: finalized counts",
                          "(11-01 release).\nRelative WIS vs the hub flatline baseline over all",
                          "forecast dates and horizons (below 1 = better).")) +
    base_theme + theme(legend.position = "right")
}
top <- c("CD", provinces)
ggsave(file.path(out, "forecasts_country_province.png"),
       trajectory_plot(top, sprintf("%s: country and provinces", variant), 3),
       width = 13, height = 7.5, dpi = 130)

ordered <- cases |> filter(location %in% zones$location) |>
  mutate(p = match(substr(location, 1, 5), provinces)) |>
  arrange(p, desc(cases)) |> pull(location)
pages <- split(ordered, ceiling(seq_along(ordered) / per_page))
zone_title <- function(pg) {
  provs <- unique(prov_names[substr(pages[[pg]], 1, 5)])
  sprintf("%s: health zones in %s, by cases (page %d of %d)",
          variant, paste(provs, collapse = ", "), pg, length(pages))
}
pdf(file.path(out, "forecasts_zones.pdf"), width = 14, height = 10)
for (pg in seq_along(pages)) print(trajectory_plot(pages[[pg]], zone_title(pg)))
dev.off()
for (pg in seq_along(pages))
  ggsave(file.path(out, sprintf("forecasts_zones_%d.png", pg)),
         trajectory_plot(pages[[pg]], zone_title(pg)), width = 14, height = 10, dpi = 100)

# --- WIS against the flatline ---------------------------------------------------------
lev_cols <- c(country = "#2a78d6", province = "#eb6834", zone = "#1baf7a")
s_all <- read_csv(file.path(root, "scores.csv"), show_col_types = FALSE) |>
  filter(model %in% c(variant, "flatline"), origin %in% origins)
rel_by <- function(...) s_all |>
  group_by(level, model, ...) |> summarise(wis = mean(wis), .groups = "drop") |>
  pivot_wider(names_from = model, values_from = wis) |>
  mutate(rel_wis = .data[[variant]] / flatline)
level_scores <- s_all |> group_by(level, model) |>
  summarise(n = n(), wis = mean(wis), cover50 = mean(in50), cover90 = mean(in90),
            .groups = "drop")
write_csv(level_scores, file.path(out, "level_scores.csv"))
rel_plot <- function(d, x, xlab, title) {
  ggplot(d, aes({{ x }}, rel_wis, colour = level)) +
    geom_hline(yintercept = 1, colour = "grey45", linewidth = 0.4) +
    geom_line(linewidth = 0.7) + geom_point(size = 2.2) +
    scale_colour_manual(values = lev_cols) +
    scale_y_log10(breaks = c(0.25, 0.5, 0.75, 1, 1.5, 2, 3)) +
    labs(x = xlab, y = "Relative WIS (log scale)", colour = NULL, title = title,
         subtitle = "Mean WIS of the model / mean WIS of the flatline baseline; below 1 = better than the flatline.") +
    base_theme + theme(legend.position = "top")
}
by_h <- rel_by(horizon)
write_csv(by_h, file.path(out, "rel_wis_horizon.csv"))
ggsave(file.path(out, "rel_wis_horizon.png"),
       rel_plot(by_h, horizon, "Horizon (weeks ahead)", sprintf("%s: relative WIS by horizon", variant)),
       width = 7, height = 4.2, dpi = 140)
by_o <- rel_by(origin)
write_csv(by_o, file.path(out, "rel_wis_origin.csv"))
ggsave(file.path(out, "rel_wis_origin.png"),
       rel_plot(by_o, origin, "Forecast date", sprintf("%s: relative WIS by forecast date (all horizons)", variant)) +
         scale_x_date(breaks = origins, date_labels = "%m-%d"),
       width = 10, height = 4.2, dpi = 140)

zs <- sc |> filter(location %in% zones$location, is.finite(rel_wis)) |>
  left_join(cases, by = "location") |>
  mutate(province = factor(prov_names[substr(location, 1, 5)], prov_names),
         name = name_of(location))
pcols <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"), prov_names)
p <- ggplot(zs, aes(cases, rel_wis, colour = province)) +
  geom_hline(yintercept = 1, colour = "grey45", linewidth = 0.4) +
  geom_point(size = 2.4, alpha = 0.85) +
  ggrepel::geom_text_repel(data = filter(zs, rel_wis > 1.1 | cases > 500), aes(label = name),
                           size = 2.8, colour = "grey25", seed = 1, min.segment.length = 0) +
  scale_colour_manual(values = pcols, drop = TRUE) +
  scale_x_log10(labels = scales::comma) + scale_y_log10(breaks = c(0.3, 0.5, 0.75, 1, 1.5)) +
  labs(x = "Finalized cases, April-August (log scale)", y = "Relative WIS (log scale)",
       colour = NULL, title = sprintf("%s: relative WIS by health zone", variant),
       subtitle = paste0("The ", nrow(zs), " zones with cases; all forecast dates and horizons. ",
                        "Below 1 = better than the flatline.")) +
  base_theme + theme(legend.position = "top")
ggsave(file.path(out, "rel_wis_zones.png"), p, width = 9, height = 5.5, dpi = 140)
message("Wrote ", out)
