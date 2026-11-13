# Diagnostic plots for a fitted BDBV joint model.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(posterior)
})

sample_bands <- function(samples) {
  samples |>
    group_by(location, target_end_date) |>
    summarise(
      q05 = quantile(value, 0.05), q25 = quantile(value, 0.25),
      median = median(value),
      q75 = quantile(value, 0.75), q95 = quantile(value, 0.95), .groups = "drop"
    )
}

# Observed counts from every vintage up to as_of (lighter = older) against the model's
# final-count samples. A gap between the latest vintage and the bands in recent weeks is
# the nowcast: cases the model expects still to be reported.
plot_fit <- function(samples, d, locations) {
  obs <- d$vintages |>
    filter(location %in% locations, as_of <= d$as_of,
           target_end_date >= min(d$weeks)) |>
    mutate(vintage = factor(as_of))
  bands <- sample_bands(filter(samples, location %in% locations))
  ggplot(bands, aes(target_end_date)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = "#2b6cb0", alpha = 0.2) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = "#2b6cb0", alpha = 0.35) +
    geom_line(aes(y = median), colour = "#2b6cb0") +
    geom_line(data = obs, aes(y = observation, group = vintage, alpha = vintage)) +
    geom_point(data = filter(obs, as_of == d$as_of), aes(y = observation), size = 0.8) +
    geom_vline(xintercept = d$as_of, linetype = "dashed", colour = "grey40") +
    scale_alpha_discrete(range = c(0.15, 1)) +
    facet_wrap(~location, scales = "free_y") +
    labs(x = "Week of symptom onset (ending)", y = "Confirmed cases",
         alpha = "Data vintage",
         title = sprintf("Final-count samples (50%%/90%%) vs reported data, as of %s",
                         d$as_of)) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")
}

# Province-level aggregation of zone samples (modeled zones only) plus the country.
aggregate_samples <- function(samples, d) {
  prov <- samples |>
    mutate(location = substr(location, 1, 5)) |>
    group_by(location, target_end_date, sample) |>
    summarise(value = sum(value), .groups = "drop")
  country <- samples |>
    group_by(target_end_date, sample) |>
    summarise(value = sum(value), .groups = "drop") |>
    mutate(location = "CD")
  bind_rows(country, prov)
}

# Fitted week completeness F, delay median and spread per province, from the posterior.
plot_reporting <- function(fit, d) {
  P <- d$stan_data$P
  W <- d$stan_data$W
  tidy <- function(var, label) {
    fit$draws(var, format = "draws_matrix") |>
      apply(2, quantile, c(0.05, 0.5, 0.95)) |>
      t() |>
      as_tibble(.name_repair = ~ c("lo", "mid", "hi")) |>
      mutate(prov = rep(seq_len(P), W), w = rep(seq_len(W), each = P), what = label)
  }
  df <- tidy("F_week", "Fraction of week reported by as-of date")
  if (d$cfg$delay == "lognormal") {
    df <- bind_rows(
      df,
      tidy("delay_m", "Delay median (days)") |> mutate(across(lo:hi, exp)),
      tidy("delay_s", "Delay spread (log scale)")
    )
  }
  df <- df |>
    mutate(province = d$provinces[prov], week = d$weeks[w]) |>
    filter(week >= max(d$weeks) - 7 * 12)
  ggplot(df, aes(week, mid, colour = province, fill = province)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line() +
    facet_wrap(~what, scales = "free_y", ncol = 1) +
    labs(x = "Week of symptom onset (ending)", y = NULL, colour = NULL, fill = NULL,
         title = "Reporting completeness and delay (posterior median, 90%)") +
    theme_minimal(base_size = 10)
}

plot_R <- function(fit, d) {
  P <- d$stan_data$P
  W <- d$stan_data$W
  df <- fit$draws("log_R_prov", format = "draws_matrix") |>
    exp() |>
    apply(2, quantile, c(0.05, 0.5, 0.95)) |>
    t() |>
    as_tibble(.name_repair = ~ c("lo", "mid", "hi")) |>
    mutate(province = d$provinces[rep(seq_len(P), W)], week = d$weeks[rep(seq_len(W), each = P)])
  ggplot(df, aes(week, mid, colour = province, fill = province)) +
    geom_hline(yintercept = 1, colour = "grey50") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line() +
    labs(x = "Week (ending)", y = "R (province baseline, before zone offsets)",
         colour = NULL, fill = NULL,
         title = "Province transmission baseline (posterior median, 90%)") +
    theme_minimal(base_size = 10)
}

# Backtest: each earlier run's final-count samples (50%/90%) against the data it was fit
# to (grey) and the latest vintage (black). Recent weeks of the latest vintage are still
# incomplete, so it is a lower bound there. `runs` is a named list of sample tibbles
# (names = as_of dates), each from samples_all_weeks.csv.gz.
plot_backtest <- function(runs, vintages, locations, weeks_before = 10) {
  latest <- max(vintages$as_of)
  bands <- bind_rows(lapply(names(runs), function(a) {
    s <- runs[[a]] |> filter(target_end_date >= as.Date(a) - 7 * weeks_before)
    s_agg <- bind_rows(
      s |> group_by(target_end_date, sample) |>
        summarise(value = sum(value), .groups = "drop") |> mutate(location = "CD"),
      s |> mutate(location = substr(location, 1, 5)) |>
        group_by(location, target_end_date, sample) |>
        summarise(value = sum(value), .groups = "drop"),
      s
    )
    sample_bands(filter(s_agg, location %in% locations)) |> mutate(run = a)
  }))
  window <- bands |> group_by(run) |>
    summarise(from = min(target_end_date), to = max(target_end_date))
  obs <- vintages |>
    filter(location %in% locations) |>
    cross_join(window) |>
    filter(target_end_date >= from, target_end_date <= to,
           as_of == as.Date(run) | as_of == latest) |>
    mutate(data = ifelse(as_of == latest, paste("Latest vintage", latest),
                         "Data the run was fit to"))
  ggplot(bands, aes(target_end_date)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = "#2b6cb0", alpha = 0.2) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = "#2b6cb0", alpha = 0.35) +
    geom_line(aes(y = median), colour = "#2b6cb0") +
    geom_line(data = obs, aes(y = observation, colour = data)) +
    geom_point(data = obs, aes(y = observation, colour = data), size = 0.8) +
    geom_vline(data = distinct(bands, run), aes(xintercept = as.Date(run)),
               linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = c("grey60", "black")) +
    facet_grid(location ~ run, scales = "free_y",
               labeller = labeller(run = function(x) paste("Run as of", x))) +
    labs(x = "Week of symptom onset (ending)", y = "Confirmed cases", colour = NULL,
         title = "Backtest: final-count samples (50%/90%) vs later data",
         subtitle = "Dashed line = run's as-of date. Latest vintage is itself incomplete for its last ~4 weeks.") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")
}
