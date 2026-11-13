# Finalized weekly case series for every modeled health zone, by province, and how
# similar the zones' dynamics are within and between provinces.
#
#   Rscript scripts/zone_series.R
#
# Data: the latest release in time-series.csv. Weeks at least `final_age` days old in it
# are treated as final (as in scripts/score_final.R); later weeks are drawn but not
# analysed. Zones: the modeled zones of the five provinces (bdbv_zones.csv), zero-filled.
#
# Correlations use zones with at least `min_cases` final cases, over onset weeks from
# `start` to the last final week, on two scales:
#   level   log(1 + weekly cases): mostly the shared rise and plateau
#   growth  week-on-week change in log(1 + 3-week mean cases): the dynamics once the
#           shared level is removed
# Within- vs between-province differences are tested by permuting province labels.
#
# Writes to output/zone_series/: zone_timeseries_overview.png, zone_series_<province>.png
# (one per province), zone_corr_matrix.png, zone_corr_pairs.png, zone_vs_province.png, and
# zone_corr_pairs.csv, zone_corr_summary.csv, zone_vs_province.csv.

final_age <- 60
min_cases <- 20
start <- as.Date("2026-04-25")
display_start <- as.Date("2026-04-04")
provinces <- c("CD-IT", "CD-NK", "CD-HU", "CD-TO", "CD-BU")
prov_names <- c("CD-IT" = "Ituri", "CD-NK" = "Nord-Kivu", "CD-HU" = "Haut-Uele",
                "CD-TO" = "Tshopo", "CD-BU" = "Bas-Uele")
# Categorical slots 1-5 of the reference palette (validated in the dataviz skill).
prov_cols <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"), prov_names)
out <- "output/zone_series"
n_perm <- 2000

suppressWarnings(suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
}))
set.seed(20261105)
dir.create(out, showWarnings = FALSE)

vint <- read_csv("data/raw/time-series.csv", show_col_types = FALSE) |>
  mutate(as_of = as.Date(as_of), target_end_date = as.Date(target_end_date))
latest <- max(vint$as_of)
last_final <- max(vint$target_end_date[vint$target_end_date <= latest - final_age])
zones <- read_csv("data/processed/bdbv_zones.csv", show_col_types = FALSE) |>
  filter(prov_code %in% provinces, !is.na(node_id)) |>
  select(location, location_name, prov_code)
weeks <- seq(display_start, max(vint$target_end_date[vint$as_of == latest]), by = 7)

series <- vint |>
  filter(as_of == latest, location %in% zones$location, target_end_date %in% weeks) |>
  select(location, target_end_date, cases = observation) |>
  complete(location = zones$location, target_end_date = weeks, fill = list(cases = 0)) |>
  left_join(zones, by = "location") |>
  mutate(province = factor(prov_names[prov_code], prov_names),
         final = target_end_date <= last_final)
totals <- series |> filter(final) |> group_by(location) |>
  summarise(total = sum(cases), peak = target_end_date[which.max(cases)],
            first = suppressWarnings(min(target_end_date[cases > 0])), .groups = "drop")
series <- left_join(series, totals, by = "location")
write_csv(select(series, location, location_name, province, target_end_date, cases, final),
          file.path(out, "zone_series.csv"))

theme_zs <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.title.position = "plot")

# --- Overview: every zone's series, one panel per province ----------------------------

# Final weeks solid; weeks still being reported dashed, joined at the last final week.
recent_part <- function(d) filter(d, target_end_date >= last_final)
ov <- series |> filter(total > 0 | !final) |>
  group_by(location) |> filter(sum(cases) > 0) |> ungroup()
top_labels <- ov |> filter(final) |> distinct(location, location_name, province, total) |>
  group_by(province) |> slice_max(total, n = 4, with_ties = FALSE) |> ungroup() |>
  left_join(ov |> filter(final) |> group_by(location) |>
              slice_max(cases, n = 1, with_ties = FALSE) |>
              select(location, target_end_date, cases), by = "location")
n_with <- ov |> distinct(province, location) |> count(province, name = "n_with")
n_all <- zones |> count(province = factor(prov_names[prov_code], prov_names), name = "n_all")
strip <- left_join(n_all, n_with, by = "province") |>
  mutate(n_with = coalesce(n_with, 0L),
         label = sprintf("%s: %d of %d zones with cases", province, n_with, n_all))
strip_lab <- setNames(strip$label, as.character(strip$province))
p <- ggplot(ov, aes(target_end_date, cases, group = location, colour = province)) +
  annotate("rect", xmin = last_final, xmax = max(weeks), ymin = -Inf, ymax = Inf,
           fill = "grey94") +
  geom_line(data = filter(ov, final), linewidth = 0.45, alpha = 0.75) +
  geom_line(data = recent_part(ov), linewidth = 0.4, alpha = 0.5, linetype = "22") +
  ggrepel::geom_text_repel(data = top_labels, aes(label = location_name), size = 2.8,
                           colour = "grey20", min.segment.length = 0, seed = 1,
                           box.padding = 0.3, max.overlaps = 20) +
  scale_colour_manual(values = prov_cols, guide = "none") +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 3, 10, 30, 100)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  facet_wrap(~province, ncol = 1, labeller = as_labeller(strip_lab)) +
  labs(x = "Week of symptom onset (ending)", y = "Weekly confirmed cases (log(1 + x) scale)",
       title = sprintf("Weekly confirmed cases by health zone, release %s", latest),
       subtitle = sprintf(paste0("One line per zone with any cases; the 4 zones with most final cases ",
                                 "in each province are labelled.\nShaded and dashed: weeks less than %d ",
                                 "days old, still being reported."), final_age)) +
  theme_zs + theme(strip.text = element_text(hjust = 0, face = "bold"))
ggsave(file.path(out, "zone_timeseries_overview.png"), p, width = 10, height = 12, dpi = 130)

# --- Small multiples: one figure per province, zones by total cases -------------------

for (pr in provinces) {
  # Zones with no cases at all are left out (the overview counts them).
  d <- filter(series, prov_code == pr) |>
    group_by(location) |> filter(sum(cases) > 0) |> ungroup() |>
    mutate(label = factor(sprintf("%s %s (%d)", location, location_name, total)),
           label = reorder(label, -total))
  prov <- d |> group_by(target_end_date, final) |> summarise(cases = sum(cases), .groups = "drop") |>
    mutate(label = sprintf("%s total (%d)", prov_names[pr], sum(cases[final])))
  d <- bind_rows(prov, d) |> mutate(label = factor(label, unique(c(prov$label[1], levels(d$label)))))
  col <- prov_cols[[prov_names[pr]]]
  n_panels <- nlevels(d$label)
  p <- ggplot(d, aes(target_end_date, cases)) +
    annotate("rect", xmin = last_final, xmax = max(weeks), ymin = -Inf, ymax = Inf,
             fill = "grey94") +
    geom_line(data = filter(d, final), colour = col, linewidth = 0.6) +
    geom_point(data = filter(d, final), colour = col, size = 0.9) +
    geom_line(data = recent_part(d), colour = col, linewidth = 0.5, linetype = "22", alpha = 0.6) +
    facet_wrap(~label, scales = "free_y", ncol = 6) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b") +
    labs(x = "Week of symptom onset (ending)", y = "Weekly confirmed cases",
         title = sprintf("%s: weekly confirmed cases by health zone, release %s", prov_names[pr], latest),
         subtitle = sprintf("Zones with any cases, by final cases (n). Shaded and dashed: weeks less than %d days old, still being reported. Free y scales.",
                            final_age)) +
    theme_zs + theme(strip.text = element_text(size = 7, hjust = 0))
  ggsave(file.path(out, sprintf("zone_series_%s.png", pr)), p, width = 13,
         height = 1.2 + 1.75 * ceiling(n_panels / 6), dpi = 120)
}

# --- Correlations ---------------------------------------------------------------------

active <- totals |> filter(total >= min_cases) |> left_join(zones, by = "location") |>
  mutate(province = factor(prov_names[prov_code], prov_names)) |>
  arrange(province, first, desc(total))
cs <- series |> filter(location %in% active$location, target_end_date >= start, final) |>
  arrange(location, target_end_date) |>
  group_by(location) |>
  mutate(level = log1p(cases),
         ma3 = (cases + lag(cases, default = 0) + lag(cases, 2, default = 0)) / 3,
         growth = log1p(ma3) - lag(log1p(ma3))) |>
  ungroup()
mat <- function(var) {
  m <- cs |> select(location, target_end_date, value = all_of(var)) |>
    pivot_wider(names_from = location, values_from = value) |>
    arrange(target_end_date) |> select(-target_end_date) |> as.matrix()
  m[, active$location]
}
M <- list(level = mat("level"), growth = mat("growth"))
C <- lapply(M, function(m) cor(m, use = "pairwise.complete.obs"))

pairs <- bind_rows(lapply(names(C), function(s) {
  as_tibble(as.data.frame.table(C[[s]], responseName = "r")) |>
    rename(a = Var1, b = Var2) |>
    mutate(across(c(a, b), as.character), scale = s)
})) |>
  filter(a < b) |>
  mutate(pa = active$province[match(a, active$location)],
         pb = active$province[match(b, active$location)],
         type = ifelse(pa == pb, "Within province", "Between provinces"),
         pair = ifelse(pa == pb, as.character(pa),
                       paste(pmin(as.character(pa), as.character(pb)),
                             pmax(as.character(pa), as.character(pb)), sep = " - ")))
write_csv(pairs, file.path(out, "zone_corr_pairs.csv"))

# Permutation test: shuffle province labels among the zones, keeping group sizes.
perm_test <- function(cm) {
  ut <- upper.tri(cm)
  stat <- function(lab) {
    same <- outer(lab, lab, "==")[ut]
    mean(cm[ut][same]) - mean(cm[ut][!same])
  }
  obs <- stat(as.character(active$province))
  null <- replicate(n_perm, stat(sample(as.character(active$province))))
  c(diff = obs, p = mean(c(null, obs) >= obs))
}
tests <- lapply(C, perm_test)
summary_tab <- pairs |>
  group_by(scale, type, pair) |>
  summarise(n_pairs = n(), mean_r = mean(r), median_r = median(r), .groups = "drop") |>
  bind_rows(pairs |> group_by(scale, type) |>
              summarise(pair = "All", n_pairs = n(), mean_r = mean(r), median_r = median(r),
                        .groups = "drop")) |>
  left_join(tibble(scale = names(tests), within_minus_between = sapply(tests, `[`, "diff"),
                   perm_p = sapply(tests, `[`, "p")), by = "scale") |>
  arrange(scale, desc(type), pair)
write_csv(summary_tab, file.path(out, "zone_corr_summary.csv"))

# Each zone against the rest of its own province and against the other provinces
# (summed counts, growth scale), so a zone's own cases don't inflate the comparison.
agg_growth <- function(locs) {
  x <- series |> filter(location %in% locs, target_end_date >= start, final) |>
    group_by(target_end_date) |> summarise(cases = sum(cases), .groups = "drop") |>
    arrange(target_end_date)
  ma3 <- (x$cases + lag(x$cases, default = 0) + lag(x$cases, 2, default = 0)) / 3
  log1p(ma3) - lag(log1p(ma3))
}
all_zone_locs <- zones$location
vs_prov <- bind_rows(lapply(seq_len(nrow(active)), function(i) {
  z <- active$location[i]; pr <- active$prov_code[i]
  own <- agg_growth(setdiff(zones$location[zones$prov_code == pr], z))
  other <- agg_growth(zones$location[zones$prov_code != pr])
  g <- M$growth[, z]
  tibble(location = z, location_name = active$location_name[i], province = active$province[i],
         total = active$total[i],
         `Rest of own province` = cor(g, own, use = "complete.obs"),
         `Other provinces` = cor(g, other, use = "complete.obs"))
}))
write_csv(vs_prov, file.path(out, "zone_vs_province.csv"))

# --- Correlation figures --------------------------------------------------------------

lab_zone <- setNames(sprintf("%s (%s)", active$location_name, substr(active$prov_code, 4, 5)),
                     active$location)
cm_long <- bind_rows(lapply(names(C), function(s) {
  as_tibble(as.data.frame.table(C[[s]], responseName = "r")) |>
    mutate(scale = s)
})) |>
  mutate(Var1 = factor(lab_zone[as.character(Var1)], lab_zone),
         Var2 = factor(lab_zone[as.character(Var2)], rev(lab_zone)),
         scale = factor(ifelse(scale == "level", "Level: log(1 + cases)",
                               "Growth: weekly change in log(1 + 3-week mean)"),
                        c("Level: log(1 + cases)", "Growth: weekly change in log(1 + 3-week mean)")))
breaks_at <- cumsum(table(active$province))
breaks_at <- breaks_at[breaks_at < nrow(active)] + 0.5
p <- ggplot(cm_long, aes(Var1, Var2, fill = r)) +
  geom_tile() +
  geom_vline(xintercept = breaks_at, colour = "grey20", linewidth = 0.4) +
  geom_hline(yintercept = nrow(active) + 1 - breaks_at, colour = "grey20", linewidth = 0.4) +
  scale_fill_gradient2(low = "#eb6834", mid = "#f2f2f0", high = "#2a78d6", midpoint = 0,
                       limits = c(-1, 1), name = "Correlation") +
  facet_wrap(~scale) +
  coord_equal() +
  labs(x = NULL, y = NULL,
       title = "Correlation of weekly case series between health zones",
       subtitle = sprintf(paste0("%d zones with at least %d final cases, weeks %s to %s, ordered by province ",
                                 "(lines) and first case. IT Ituri, NK Nord-Kivu, HU Haut-Uele, TO Tshopo."),
                          nrow(active), min_cases, start, last_final)) +
  theme_zs + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
                   axis.text.y = element_text(size = 6), panel.grid = element_blank())
ggsave(file.path(out, "zone_corr_matrix.png"), p, width = 15, height = 8.5, dpi = 130)

pair_order <- pairs |> distinct(type, pair) |> arrange(desc(type), pair) |> pull(pair)
pp <- pairs |>
  mutate(pair = factor(pair, rev(pair_order)),
         scale = factor(ifelse(scale == "level", "Level: log(1 + cases)",
                               "Growth: weekly change in log(1 + 3-week mean)"),
                        levels(cm_long$scale)))
ptxt <- tibble(scale = factor(c("Level: log(1 + cases)",
                                "Growth: weekly change in log(1 + 3-week mean)"), levels(cm_long$scale)),
               text = sprintf("Within minus between: %+.2f (permutation p = %.3f)",
                              sapply(tests, `[`, "diff"), sapply(tests, `[`, "p")))
pn <- pp |> group_by(scale, pair) |> summarise(n = n(), .groups = "drop")
p <- ggplot(pp, aes(r, pair)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_boxplot(aes(fill = type), outlier.shape = NA, width = 0.6, alpha = 0.35, colour = "grey35") +
  geom_jitter(aes(colour = type), height = 0.15, size = 1, alpha = 0.6) +
  geom_text(data = pn, aes(x = 1.02, label = paste0("n=", n)), hjust = 0, size = 2.8, colour = "grey30") +
  geom_text(data = ptxt, aes(x = -1, y = Inf, label = text), hjust = 0, vjust = 1.3, size = 3,
            colour = "grey20", inherit.aes = FALSE) +
  scale_colour_manual(values = c("Within province" = "#2a78d6", "Between provinces" = "#eb6834"),
                      aesthetics = c("colour", "fill")) +
  scale_x_continuous(limits = c(-1, 1.15), breaks = seq(-1, 1, 0.5)) +
  facet_wrap(~scale) +
  labs(x = "Pairwise correlation between zones", y = NULL, colour = NULL, fill = NULL,
       title = "Zone-pair correlations within and between provinces",
       subtitle = "One point per pair of zones; boxes: median and quartiles. n = number of zone pairs.") +
  theme_zs + theme(legend.position = "bottom")
ggsave(file.path(out, "zone_corr_pairs.png"), p, width = 12, height = 6, dpi = 130)

vp <- vs_prov |> pivot_longer(c(`Rest of own province`, `Other provinces`), names_to = "versus",
                              values_to = "r") |>
  mutate(label = reorder(sprintf("%s (%d)", location_name, total), total))
p <- ggplot(vp, aes(r, label)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_line(aes(group = label), colour = "grey75", linewidth = 0.4) +
  geom_point(aes(colour = versus, shape = versus), size = 2.2) +
  scale_colour_manual(values = c("Rest of own province" = "#2a78d6", "Other provinces" = "#eb6834")) +
  scale_shape_manual(values = c("Rest of own province" = 16, "Other provinces" = 17)) +
  facet_grid(province ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Correlation of zone growth with the aggregate's growth", y = NULL, colour = NULL, shape = NULL,
       title = "Does a zone's growth track its own province or the other provinces?",
       subtitle = "Growth: weekly change in log(1 + 3-week mean cases). Rest of own province leaves the zone out. (n) = final cases.") +
  theme_zs + theme(legend.position = "bottom", strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"))
ggsave(file.path(out, "zone_vs_province.png"), p, width = 9, height = 9, dpi = 130)

cat(sprintf("Release %s; final weeks to %s; %d zones, %d with >= %d cases\n\n", latest, last_final,
            nrow(zones), nrow(active), min_cases))
print(summary_tab, n = 50, width = 120)
cat("\nZone vs aggregates (growth), mean over zones by province:\n")
vs_prov |> group_by(province) |>
  summarise(n = n(), own = mean(`Rest of own province`), other = mean(`Other provinces`)) |> print()
