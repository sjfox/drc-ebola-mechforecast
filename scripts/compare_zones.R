# Compare two model variants' backtests zone by zone, for every modeled health zone.
#
#   Rscript scripts/compare_zones.R output/bdbv_v1 output/bdbv_v5 2026-10-17 2026-10-24
#
# Each root needs <root>/<as_of>/samples_all_weeks.csv.gz (eventual counts) from
# run_bdbv.R. Truth is the latest release in time-series.csv. Writes, to output/:
#   zone_comparison_<as_of>.pdf  country, provinces, then every zone (20 per page):
#                                median and 90% band of each model, the release each
#                                run was fit to, and the latest release
#   zone_crps.csv, zone_crps.png per-zone CRPS of each model on weeks at least 14 days
#                                old in the latest release, summed over weeks

args <- commandArgs(trailingOnly = TRUE)
roots <- args[1:2]
as_ofs <- as.Date(args[-(1:2)])
if (length(as_ofs) == 0) stop("Pass two output roots and one or more as_of dates")
min_truth_age <- 14
weeks_before <- 8
per_page <- 20

suppressWarnings(suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
}))
source("R/bdbv_plots.R")

model_names <- basename(roots)
# Categorical slots 1 and 2 of the reference palette (validated: CVD dE 24.7).
model_cols <- setNames(c("#2a78d6", "#eb6834"), model_names)

crps_sample <- function(x, y) mean(abs(x - y)) - 0.5 * mean(abs(outer(x, x, "-")))

vint <- read_csv("data/raw/time-series.csv", show_col_types = FALSE) |>
  mutate(as_of = as.Date(as_of), target_end_date = as.Date(target_end_date))
latest <- max(vint$as_of)
zones <- read_csv("data/processed/bdbv_zones.csv", show_col_types = FALSE) |>
  select(location, location_name)
prov_order <- c("CD-IT", "CD-NK", "CD-HU", "CD-TO", "CD-BU")

read_run <- function(root, a) {
  read_csv(file.path(root, a, "samples_all_weeks.csv.gz"), show_col_types = FALSE) |>
    filter(target_end_date >= a - 7 * weeks_before)
}

# Zone samples plus province and country sums (within each sample).
with_aggregates <- function(s) {
  bind_rows(
    s |> group_by(target_end_date, sample) |> summarise(value = sum(value), .groups = "drop") |>
      mutate(location = "CD"),
    s |> mutate(location = substr(location, 1, 5)) |>
      group_by(location, target_end_date, sample) |> summarise(value = sum(value), .groups = "drop"),
    s
  )
}

crps_rows <- list()
for (i in seq_along(as_ofs)) {
  a <- as_ofs[i]
  runs <- lapply(roots, function(r) with_aggregates(read_run(r, a)))
  names(runs) <- model_names
  locs <- Reduce(intersect, lapply(runs, function(s) unique(s$location)))

  obs <- vint |>
    filter(location %in% locs, as_of %in% c(a, latest), target_end_date >= a - 7 * weeks_before) |>
    mutate(data = ifelse(as_of == latest, paste("Release", latest), paste("Release", a, "(fit to)")))
  # Cap bands at 3x the largest count seen in the window so one explosive forecast
  # doesn't flatten the panel; capped bands are marked in the subtitle.
  cap <- obs |> group_by(location) |> summarise(cap = 3 * max(observation, 5), .groups = "drop")
  bands <- bind_rows(lapply(model_names, function(m) {
    sample_bands(filter(runs[[m]], location %in% locs)) |> mutate(model = m)
  })) |>
    left_join(cap, by = "location") |>
    mutate(capped = q95 > cap, q95 = pmin(q95, cap), median = pmin(median, cap))

  # Order: country, provinces, then zones by province and by cases in the window.
  recent <- obs |> filter(as_of == latest) |> group_by(location) |>
    summarise(n = sum(observation), .groups = "drop")
  zone_locs <- tibble(location = setdiff(locs, c("CD", prov_order))) |>
    left_join(recent, by = "location") |>
    mutate(prov = factor(substr(location, 1, 5), prov_order), n = coalesce(n, 0)) |>
    arrange(prov, desc(n)) |>
    pull(location)
  ordered <- c("CD", intersect(prov_order, locs), zone_locs)
  labels <- tibble(location = ordered) |>
    left_join(zones, by = "location") |>
    mutate(label = ifelse(is.na(location_name), location, paste(location, location_name)))
  lab <- setNames(labels$label, labels$location)

  pages <- split(ordered, ceiling(seq_along(ordered) / per_page))
  pdf(file.path("output", sprintf("zone_comparison_%s.pdf", a)), width = 14, height = 10)
  for (pg in seq_along(pages)) {
    page_locs <- pages[[pg]]
    b <- filter(bands, location %in% page_locs) |>
      mutate(location = factor(location, page_locs))
    o <- filter(obs, location %in% page_locs) |> mutate(location = factor(location, page_locs))
    n_capped <- b |> filter(capped) |> distinct(location, model) |> nrow()
    p <- ggplot(b, aes(target_end_date)) +
      geom_ribbon(aes(ymin = q05, ymax = q95, fill = model), alpha = 0.18) +
      geom_line(aes(y = median, colour = model), linewidth = 0.6) +
      geom_vline(xintercept = a, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
      geom_point(data = o, aes(y = observation, shape = data), size = 1.3, colour = "grey15") +
      scale_colour_manual(values = model_cols, aesthetics = c("colour", "fill")) +
      scale_shape_manual(values = setNames(c(1, 16), c(paste("Release", a, "(fit to)"),
                                                        paste("Release", latest)))) +
      facet_wrap(~location, scales = "free_y", ncol = 5, labeller = as_labeller(lab)) +
      labs(x = "Week of symptom onset (ending)", y = "Confirmed cases (eventual)",
           colour = "Model (median, 90%)", fill = "Model (median, 90%)", shape = NULL,
           title = sprintf("Backtest as of %s: eventual counts, %s vs %s (page %d of %d)",
                           a, model_names[1], model_names[2], pg, length(pages)),
           subtitle = paste0(
             "Dashed line: as_of. Points: release the runs were fit to (open) and release ",
             latest, " (filled), itself incomplete for recent weeks. Bands and medians cut at ",
             "3x the largest count in the window",
             if (n_capped) sprintf(" (%d model-panels cut on this page).", n_capped) else ".")) +
      theme_minimal(base_size = 9) +
      theme(legend.position = "bottom", panel.grid.minor = element_blank(),
            strip.text = element_text(size = 7.5, hjust = 0))
    print(p)
  }
  dev.off()

  truth <- vint |>
    filter(as_of == latest, location %in% locs, target_end_date >= a - 7 * weeks_before,
           target_end_date <= latest - min_truth_age) |>
    select(location, target_end_date, truth = observation)
  crps_rows[[i]] <- bind_rows(lapply(model_names, function(m) {
    runs[[m]] |>
      inner_join(truth, by = c("location", "target_end_date")) |>
      group_by(location, target_end_date) |>
      summarise(crps = crps_sample(value, truth[1]), .groups = "drop") |>
      mutate(model = m)
  })) |>
    mutate(as_of = a)
}

crps <- bind_rows(crps_rows) |>
  group_by(as_of, location, model) |>
  summarise(crps = sum(crps), .groups = "drop") |>
  tidyr::pivot_wider(names_from = model, values_from = crps) |>
  mutate(level = case_when(location == "CD" ~ "country", nchar(location) == 5 ~ "province",
                           TRUE ~ "zone"))
write_csv(crps, "output/zone_crps.csv")

z <- filter(crps, level == "zone") |>
  rename(a = !!model_names[1], b = !!model_names[2]) |>
  mutate(better = case_when(b < a * 0.9 ~ paste(model_names[2], "better"),
                            a < b * 0.9 ~ paste(model_names[1], "better"),
                            TRUE ~ "within 10%"))
totals <- z |> group_by(as_of) |>
  summarise(a = sum(a), b = sum(b),
            n_b = sum(better == paste(model_names[2], "better")),
            n_a = sum(better == paste(model_names[1], "better")), n = n(), .groups = "drop") |>
  mutate(text = sprintf("Summed CRPS: %s %.0f, %s %.0f\n%s better in %d zones, %s in %d (of %d)",
                        model_names[1], a, model_names[2], b, model_names[2], n_b,
                        model_names[1], n_a, n))
lim <- c(0.05, max(c(z$a, z$b)) * 1.2)
p <- ggplot(z, aes(pmax(a, 0.05), pmax(b, 0.05))) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(aes(colour = better), size = 2, alpha = 0.85) +
  geom_text(data = totals, aes(x = lim[1], y = lim[2], label = text), inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 3, colour = "grey20") +
  scale_x_log10(limits = lim) + scale_y_log10(limits = lim) +
  scale_colour_manual(values = setNames(c(model_cols[[2]], model_cols[[1]], "grey65"),
                                        c(paste(model_names[2], "better"),
                                          paste(model_names[1], "better"), "within 10%"))) +
  facet_wrap(~as_of, labeller = labeller(as_of = function(x) paste("Fit as of", x))) +
  labs(x = paste(model_names[1], "CRPS (summed over weeks)"),
       y = paste(model_names[2], "CRPS (summed over weeks)"), colour = NULL,
       title = sprintf("Per-zone backtest score, %s vs %s (lower is better)", model_names[1], model_names[2]),
       subtitle = sprintf("Eventual-count samples scored against release %s, weeks at least %d days old. Below the diagonal: %s better.",
                          latest, min_truth_age, model_names[2])) +
  coord_equal() +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave("output/zone_crps.png", p, width = 10, height = 5.8, dpi = 130)

cat("Per-level CRPS (country/province: summed over weeks; zone: summed over zones and weeks):\n")
crps |> group_by(as_of, level) |>
  summarise(across(all_of(model_names), sum), n = n(), .groups = "drop") |> print()
