# Score finalized-data backtests (scripts/backtest_final.R) with the weighted interval
# score (WIS), against baselines fit to the same data.
#
#   Rscript scripts/score_final.R                  # every variant in output/final
#   Rscript scripts/score_final.R E0 E1_zone_tight
#
# Truth: the latest release's counts for onset weeks at least `final_age` days old in it.
# Forecasts are scored for horizons 1-4 weeks past each origin, for the country, the five
# modeled provinces and their 117 health zones, on the hub's 23 quantile levels. Model
# provinces and country are sums of zone samples (the country leaves out the few cases
# outside the modeled provinces).
#
# Baselines, each fit separately to every location's series up to the origin:
#   flatline   the hub baseline: median = last count; spread from the symmetrized
#              week-to-week changes over the whole series, summed over the horizon and
#              truncated at 0 (as in the FluSight / hubverse quantile baseline)
#   loglinear  quasi-Poisson log-linear trend over the last 4 weeks, extrapolated from
#              the fitted last-week count with trend uncertainty (weekly growth capped at
#              x2 / x0.5) and negative binomial noise
#   renewal    renewal equation: R from the last 2 weeks (Cori et al. gamma posterior,
#              prior mean 5, sd 5) with the model's generation interval (incubation 7.5
#              + infectious 8 days, exponential), projected with Poisson noise
#
# Writes to <root> (BDBV_FINAL_OUT, default output/final): scores.csv (one row per
# model, origin, location, week), summary.csv, relative_wis.png, coverage.png, wis_by_origin.png,
# forecasts.png.

args <- commandArgs(trailingOnly = TRUE)
root <- Sys.getenv("BDBV_FINAL_OUT", "output/final")
final_age <- 60
horizon_weeks <- 4
quantile_levels <- round(c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99), 3)
qi <- function(p) which(quantile_levels == p)
provinces <- c("CD-IT", "CD-NK", "CD-HU", "CD-TO", "CD-BU")
n_baseline <- 4000

suppressWarnings(suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
}))
set.seed(20261105)

variants <- if (length(args)) args else list.dirs(root, recursive = FALSE, full.names = FALSE)
variants <- variants[file.exists(file.path(root, variants))]

vint <- read_csv("data/raw/time-series.csv", show_col_types = FALSE) |>
  mutate(as_of = as.Date(as_of), target_end_date = as.Date(target_end_date))
latest <- max(vint$as_of)
release <- vint |>
  filter(as_of == latest,
         location == "CD" | location %in% provinces | substr(location, 1, 5) %in% provinces) |>
  select(location, target_end_date, observation) |>
  mutate(level = case_when(location == "CD" ~ "country", nchar(location) == 5 ~ "province",
                           TRUE ~ "zone"))
truth <- release |> filter(target_end_date <= latest - final_age)

# --- Scoring --------------------------------------------------------------------------

# WIS from quantiles: the mean over quantile levels of twice the quantile (pinball)
# loss, which equals (1 / (K + 1/2)) (|y - median| / 2 + sum_k alpha_k / 2 * IS_alpha_k).
wis <- function(q, y) mean(2 * ((y < q) - quantile_levels) * (q - y))

score_quantiles <- function(qtab) {
  # qtab: model, origin, location, target_end_date, q (list column of quantiles)
  qtab |>
    inner_join(truth, by = c("location", "target_end_date")) |>
    rowwise() |>
    mutate(wis = wis(q, observation),
           median = q[qi(0.5)],
           in50 = observation >= q[qi(0.25)] & observation <= q[qi(0.75)],
           in90 = observation >= q[qi(0.05)] & observation <= q[qi(0.95)]) |>
    ungroup() |>
    mutate(horizon = as.integer(target_end_date - origin) / 7) |>
    select(model, origin, location, level, target_end_date, horizon, truth = observation,
           median, wis, in50, in90, q)
}

to_quantiles <- function(samples) {
  # samples: location, target_end_date, value (many rows per location-week)
  samples |>
    group_by(location, target_end_date) |>
    summarise(q = list(unname(quantile(value, quantile_levels, type = 7))), .groups = "drop")
}

# --- Model forecasts -----------------------------------------------------------------

model_quantiles <- function(variant) {
  dirs <- list.dirs(file.path(root, variant), recursive = FALSE)
  bind_rows(lapply(dirs, function(dir) {
    f <- file.path(dir, "forecast_samples.csv.gz")
    if (!file.exists(f)) return(NULL)
    z <- read_csv(f, show_col_types = FALSE)
    s <- bind_rows(
      z |> group_by(target_end_date, sample) |> summarise(value = sum(value), .groups = "drop") |>
        mutate(location = "CD"),
      z |> mutate(location = substr(location, 1, 5)) |>
        group_by(location, target_end_date, sample) |> summarise(value = sum(value), .groups = "drop"),
      z
    )
    to_quantiles(s) |> mutate(model = variant, origin = as.Date(basename(dir)))
  }))
}

# --- Baselines -----------------------------------------------------------------------

# Weekly generation interval: onset-to-onset time is an infectious-period draw plus an
# incubation draw (exponential, means 8 and 7.5 days); week k collects days
# [7k - 3.5, 7k + 3.5), with week 0 folded into week 1.
gi_weekly <- local({
  a <- 1 / 7.5; b <- 1 / 8
  cdf <- function(t) ifelse(t <= 0, 0, 1 - (b * exp(-a * t) - a * exp(-b * t)) / (b - a))
  k <- 1:12
  w <- cdf(7 * k + 3.5) - cdf(7 * k - 3.5)
  w[1] <- w[1] + cdf(3.5)
  w / sum(w)
})

baseline_samples <- function(x, n) {
  # x: weekly counts up to the origin (oldest first). Returns a list of n x horizon
  # sample matrices, one per baseline.
  last <- x[length(x)]
  H <- horizon_weeks

  d <- diff(x)
  sym <- c(d, -d)
  flat <- if (length(sym) == 0 || all(sym == 0)) {
    matrix(last, n, H)
  } else {
    steps <- matrix(sample(sym, n * H, replace = TRUE), n, H)
    pmax(last + t(apply(steps, 1, cumsum)), 0)
  }

  recent <- tail(x, 4)
  loglin <- if (sum(recent) == 0) {
    matrix(0, n, H)
  } else {
    tt <- seq_along(recent)
    fit <- suppressWarnings(glm(recent ~ tt, family = quasipoisson()))
    phi <- max(summary(fit)$dispersion, 1)
    # Start from the fitted count for the last week and draw only the weekly trend,
    # capped at doubling / halving: sparse series (0, 0, 1, 3) have huge coefficient
    # standard errors and otherwise extrapolate to absurd counts.
    slope <- rnorm(n, coef(fit)[2], sqrt(vcov(fit)[2, 2]))
    slope <- pmin(pmax(slope, -log(2)), log(2))
    mu <- fitted(fit)[length(recent)] * exp(slope %o% (1:H))
    if (phi > 1.001) {
      matrix(rnbinom(n * H, mu = mu, size = mu / (phi - 1)), n, H)
    } else {
      matrix(rpois(n * H, mu), n, H)
    }
  }

  lam <- function(series) {
    k <- seq_along(gi_weekly)
    sum(gi_weekly[k <= length(series)] * rev(series)[k[k <= length(series)]])
  }
  renew <- {
    t_idx <- (length(x) - 1):length(x)
    Lambda <- sapply(t_idx, function(t) lam(x[seq_len(t - 1)]))
    if (sum(Lambda) == 0) {
      matrix(0, n, H)
    } else {
      R <- rgamma(n, shape = 1 + sum(x[t_idx]), rate = 1 / 5 + sum(Lambda))
      out <- matrix(0, n, H)
      for (i in seq_len(n)) {
        series <- x
        for (h in 1:H) {
          series <- c(series, rpois(1, R[i] * lam(series)))
        }
        out[i, ] <- tail(series, H)
      }
      out
    }
  }
  list(flatline = flat, loglinear = loglin, renewal = renew)
}

baseline_quantiles <- function(origins) {
  locs <- unique(release$location)
  bind_rows(lapply(origins, function(o) {
    hist <- release |> filter(target_end_date <= o, target_end_date >= as.Date("2026-04-04"))
    bind_rows(lapply(locs, function(l) {
      x <- hist |> filter(location == l) |> arrange(target_end_date) |> pull(observation)
      if (length(x) < 2) return(NULL)
      bs <- baseline_samples(x, n_baseline)
      bind_rows(lapply(names(bs), function(m) {
        tibble(model = m, origin = o, location = l,
               target_end_date = o + 7 * seq_len(horizon_weeks),
               q = lapply(seq_len(horizon_weeks),
                          function(h) unname(quantile(bs[[m]][, h], quantile_levels))))
      }))
    }))
  }))
}

# --- Run -----------------------------------------------------------------------------

model_q <- bind_rows(lapply(variants, model_quantiles))
origins <- sort(unique(model_q$origin))
message(sprintf("Scoring %d variants over %d origins (%s to %s)", length(variants),
                length(origins), min(origins), max(origins)))
base_q <- baseline_quantiles(origins)
scores <- score_quantiles(bind_rows(model_q, base_q))
write_csv(select(scores, -q), file.path(root, "scores.csv"))

# Only location-week-origins every model has.
n_models <- n_distinct(scores$model)
common <- scores |> group_by(origin, location, target_end_date) |>
  filter(n_distinct(model) == n_models) |> ungroup()

summary <- common |>
  group_by(model, level, horizon) |>
  summarise(n = n(), wis = mean(wis), cover50 = mean(in50), cover90 = mean(in90),
            .groups = "drop") |>
  group_by(level, horizon) |>
  mutate(rel_wis = wis / wis[model == "flatline"]) |>
  ungroup()
overall <- common |>
  group_by(model, level) |>
  summarise(n = n(), wis = mean(wis), cover50 = mean(in50), cover90 = mean(in90),
            .groups = "drop") |>
  group_by(level) |>
  mutate(rel_wis = wis / wis[model == "flatline"]) |>
  ungroup() |>
  mutate(horizon = NA_integer_)
write_csv(bind_rows(summary, overall), file.path(root, "summary.csv"))
cat("WIS relative to the flatline baseline (all horizons; lower is better):\n")
print(overall |> select(level, model, n, wis, rel_wis, cover50, cover90) |>
        arrange(factor(level, c("country", "province", "zone")), rel_wis), n = 60)

# --- Figures -------------------------------------------------------------------------

baselines <- c("flatline", "loglinear", "renewal")
model_levels <- c(setdiff(sort(unique(scores$model)), baselines), baselines)
# Variants take categorical slots in order (validated palette); baselines are grey with
# distinct line types.
variant_cols <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300")
v_names <- setdiff(model_levels, baselines)
cols <- c(setNames(variant_cols[seq_along(v_names)], v_names),
          flatline = "grey20", loglinear = "grey45", renewal = "grey65")
ltys <- c(setNames(rep("solid", length(v_names)), v_names),
          flatline = "solid", loglinear = "dashed", renewal = "dotted")
level_order <- c("country", "province", "zone")
theme_bt <- theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

p <- summary |>
  mutate(model = factor(model, model_levels), level = factor(level, level_order)) |>
  ggplot(aes(horizon, rel_wis, colour = model, linetype = model)) +
  geom_hline(yintercept = 1, colour = "grey70", linewidth = 0.3) +
  geom_line(linewidth = 0.7) + geom_point(size = 2) +
  scale_colour_manual(values = cols) + scale_linetype_manual(values = ltys) +
  scale_y_log10() +
  facet_wrap(~level) +
  labs(x = "Horizon (weeks past the last data week)", y = "Mean WIS / flatline mean WIS (log scale)",
       colour = NULL, linetype = NULL,
       title = "Finalized-data backtest: WIS relative to the hub flatline baseline",
       subtitle = sprintf("%d origins, %s to %s; truth = release %s, weeks at least %d days old. Below 1: better than flatline.",
                          length(origins), min(origins), max(origins), latest, final_age)) +
  theme_bt
ggsave(file.path(root, "relative_wis.png"), p, width = 10, height = 4.8, dpi = 130)

p <- summary |>
  select(model, level, horizon, cover50, cover90) |>
  pivot_longer(c(cover50, cover90), names_to = "interval", values_to = "coverage") |>
  mutate(nominal = ifelse(interval == "cover50", 0.5, 0.9),
         interval = ifelse(interval == "cover50", "50% interval", "90% interval"),
         model = factor(model, model_levels), level = factor(level, level_order)) |>
  ggplot(aes(horizon, coverage, colour = model, linetype = model)) +
  geom_hline(aes(yintercept = nominal), colour = "grey70", linewidth = 0.3) +
  geom_line(linewidth = 0.7) + geom_point(size = 2) +
  scale_colour_manual(values = cols) + scale_linetype_manual(values = ltys) +
  facet_grid(interval ~ level) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Horizon (weeks)", y = "Share of truths inside the interval", colour = NULL,
       linetype = NULL, title = "Interval coverage by horizon (grey line: nominal)") +
  theme_bt
ggsave(file.path(root, "coverage.png"), p, width = 10, height = 6, dpi = 130)

p <- common |>
  filter(level != "zone") |>
  group_by(model, origin, level) |>
  summarise(wis = mean(wis), .groups = "drop") |>
  group_by(origin, level) |>
  mutate(rel_wis = wis / wis[model == "flatline"]) |>
  ungroup() |>
  mutate(model = factor(model, model_levels), level = factor(level, level_order)) |>
  ggplot(aes(origin, rel_wis, colour = model, linetype = model)) +
  geom_hline(yintercept = 1, colour = "grey70", linewidth = 0.3) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.6) +
  scale_colour_manual(values = cols) + scale_linetype_manual(values = ltys) +
  scale_y_log10() +
  facet_wrap(~level, ncol = 1) +
  labs(x = "Forecast origin (last data week)", y = "WIS / flatline WIS (log scale)",
       colour = NULL, linetype = NULL,
       title = "Relative WIS by origin (all horizons)") +
  theme_bt
ggsave(file.path(root, "wis_by_origin.png"), p, width = 10, height = 6.5, dpi = 130)

# Forecast fans for the country and the four largest provinces: every origin's median and
# 90% interval against the truth, one column per model.
fan <- scores |>
  filter(location %in% c("CD", "CD-IT", "CD-NK", "CD-HU", "CD-TO"),
         model %in% c(v_names, "flatline")) |>
  rowwise() |>
  mutate(lo = q[qi(0.05)], hi = q[qi(0.95)]) |>
  ungroup() |>
  mutate(model = factor(model, model_levels),
         location = factor(location, c("CD", "CD-IT", "CD-NK", "CD-HU", "CD-TO")))
obs <- truth |> filter(location %in% levels(fan$location), target_end_date >= as.Date("2026-04-25")) |>
  mutate(location = factor(location, levels(fan$location)))
cap <- obs |> group_by(location) |> summarise(cap = 3 * max(observation), .groups = "drop")
fan <- fan |> left_join(cap, by = "location") |> mutate(hi = pmin(hi, cap), median = pmin(median, cap))
p <- ggplot(fan, aes(target_end_date, group = origin)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = model), alpha = 0.12) +
  geom_line(aes(y = median, colour = model), linewidth = 0.4) +
  geom_line(data = obs, aes(target_end_date, observation), inherit.aes = FALSE,
            colour = "grey15", linewidth = 0.4) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  facet_grid(location ~ model, scales = "free_y") +
  labs(x = "Week of symptom onset (ending)", y = "Confirmed cases", colour = NULL, fill = NULL,
       title = "Forecasts from every origin (median and 90% interval, 1-4 weeks) against finalized counts",
       subtitle = "Black line: finalized counts. Bands cut at 3x the largest count.") +
  theme_bt + theme(legend.position = "none")
ggsave(file.path(root, "forecasts.png"), p, width = 3 + 2.6 * nlevels(droplevels(fan$model)),
       height = 9, dpi = 130)
message("Wrote scores.csv, summary.csv and figures to ", root)
