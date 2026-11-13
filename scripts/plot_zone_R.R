# Zone reproduction numbers from the finalized-data backtest fits (scripts/backtest_final.R),
# with posterior uncertainty, for every origin.
#
#   Rscript scripts/plot_zone_R.R            # E0 fits, province R held flat (as E0_flat)
#   Rscript scripts/plot_zone_R.R E1_zone_tight
#
# Zone R in week w is exp(province log R + zone spline deviation), the R the zone's
# transmission rate would give in a fully susceptible population (no depletion, no
# movement). In forecast weeks the province log R stays at its last value (E0_flat); the
# zone deviation follows the spline, as in the forecasts.
#
# Writes to <root>/zone_R/ (root: BDBV_FINAL_OUT, default output/final):
#   zone_R_<origin>.pdf    province and zone R (median, 50%, 90%), 20 panels per page
#   zone_R_vintages.png    median zone R from every origin, for the zones with most cases
#   zone_R_summary.csv     per origin and zone: R at the last data week and 4 weeks on
#   zone_R_quantiles.csv.gz  every origin, location and week's R quantiles

args <- commandArgs(trailingOnly = TRUE)
variant <- if (length(args)) args[1] else "E0"
root <- Sys.getenv("BDBV_FINAL_OUT", "output/final")
out <- file.path(root, "zone_R")
provinces <- c("CD-IT", "CD-NK", "CD-HU", "CD-TO", "CD-BU")
n_draws <- 400
per_page <- 20

suppressWarnings(suppressPackageStartupMessages({
  source("R/bdbv_data.R")
  library(ggplot2)
  library(cmdstanr)
}))
dir.create(out, showWarnings = FALSE)
zone_names <- read_csv("data/processed/bdbv_zones.csv", show_col_types = FALSE) |>
  select(location, location_name)
prov_names <- c("CD-IT" = "Ituri", "CD-NK" = "Nord-Kivu", "CD-HU" = "Haut-Uele",
                "CD-TO" = "Tshopo", "CD-BU" = "Bas-Uele")

r_quantiles <- function(m) {
  # m: draws x cells. Returns median, 50% and 90% bounds per cell.
  q <- apply(m, 2, quantile, c(0.05, 0.25, 0.5, 0.75, 0.95))
  tibble(lo90 = q[1, ], lo50 = q[2, ], med = q[3, ], hi50 = q[4, ], hi90 = q[5, ])
}

origins <- sort(as.Date(basename(list.dirs(file.path(root, variant), recursive = FALSE))))
origins <- origins[file.exists(file.path(root, variant, origins, "fit.rds"))]
all_q <- list()
for (o in as.list(origins)) {
  d <- build_bdbv_data(bdbv_config(finalized = TRUE, origin = o, provinces = provinces))
  sd <- d$stan_data
  Z <- sd$Z; P <- sd$P; W <- sd$W; H <- sd$H; NB <- sd$NB
  fit <- readRDS(file.path(root, variant, o, "fit.rds"))
  dr <- posterior::as_draws_matrix(fit$draws(c("log_R_prov", "a_zone")))
  dr <- dr[round(seq(1, nrow(dr), length.out = min(n_draws, nrow(dr)))), ]
  lrp <- dr[, grep("^log_R_prov", colnames(dr))]  # column-major: [p, w]
  az <- dr[, grep("^a_zone", colnames(dr))]       # column-major: [z, k]
  all_weeks <- c(d$weeks, d$forecast_weeks)
  zone_R <- array(NA_real_, c(nrow(dr), Z, W + H))
  prov_R <- array(NA_real_, c(nrow(dr), P, W + H))
  for (i in seq_len(nrow(dr))) {
    lr <- matrix(lrp[i, ], P, W)
    lr <- cbind(lr, matrix(lr[, W], P, H))
    dev <- matrix(az[i, ], Z, NB) %*% t(sd$Bs)
    zone_R[i, , ] <- exp(lr[sd$prov, ] + dev)
    prov_R[i, , ] <- exp(lr)
  }
  cells <- function(locs, n) expand.grid(i = seq_along(locs), w = seq_len(W + H)) |>
    as_tibble() |> mutate(location = locs[i], target_end_date = all_weeks[w])
  q <- bind_rows(
    cells(d$provinces, P) |> bind_cols(r_quantiles(matrix(prov_R, nrow(dr)))) |>
      mutate(level = "province"),
    cells(d$zones$location, Z) |> bind_cols(r_quantiles(matrix(zone_R, nrow(dr)))) |>
      mutate(level = "zone")
  ) |>
    select(-i, -w) |>
    mutate(origin = o, forecast = target_end_date > o)
  cases <- tibble(location = d$zones$location,
                  cases = rowSums(sd$cases),
                  recent = rowSums(sd$cases[, (W - 3):W, drop = FALSE]))
  cases <- bind_rows(cases, cases |> mutate(location = substr(location, 1, 5)) |>
                       group_by(location) |> summarise(across(c(cases, recent), sum)))
  q <- left_join(q, cases, by = "location")
  all_q[[as.character(o)]] <- q

  # One PDF per origin: provinces, then zones by province and cases to date.
  ordered <- c(d$provinces,
               cases |> filter(nchar(location) > 5) |>
                 mutate(p = match(substr(location, 1, 5), provinces)) |>
                 arrange(p, desc(cases)) |> pull(location))
  lab <- tibble(location = ordered) |>
    left_join(zone_names, by = "location") |>
    left_join(cases, by = "location") |>
    mutate(name = coalesce(location_name, prov_names[location]),
           label = sprintf("%s %s\n%d cases, %d in last 4 wk", location, name, cases, recent))
  lab <- setNames(lab$label, lab$location)
  pages <- split(ordered, ceiling(seq_along(ordered) / per_page))
  pdf(file.path(out, sprintf("zone_R_%s.pdf", o)), width = 14, height = 10)
  for (pg in seq_along(pages)) {
    b <- filter(q, location %in% pages[[pg]]) |>
      mutate(location = factor(location, pages[[pg]]),
             across(c(lo90, lo50, med, hi50, hi90), ~ pmin(pmax(.x, 0.1), 10)))
    p <- ggplot(b, aes(target_end_date)) +
      annotate("rect", xmin = o, xmax = max(all_weeks) + 3, ymin = -Inf, ymax = Inf,
               fill = "grey92") +
      geom_hline(yintercept = 1, colour = "grey50", linewidth = 0.3) +
      geom_ribbon(aes(ymin = lo90, ymax = hi90), fill = "#2a78d6", alpha = 0.18) +
      geom_ribbon(aes(ymin = lo50, ymax = hi50), fill = "#2a78d6", alpha = 0.35) +
      geom_line(aes(y = med), colour = "#1d5aa3", linewidth = 0.6) +
      scale_y_log10(breaks = c(0.1, 0.3, 1, 3, 10), labels = c("0.1", "0.3", "1", "3", "10")) +
      coord_cartesian(ylim = c(0.1, 10)) +
      facet_wrap(~location, ncol = 5, labeller = as_labeller(lab)) +
      labs(x = "Week of symptom onset (ending)", y = "R (log scale)",
           title = sprintf("%s, origin %s: province and zone R, median with 50%% and 90%% intervals (page %d of %d)",
                           variant, o, pg, length(pages)),
           subtitle = paste("Shaded: forecast weeks (province R held at its last value; zone deviation follows the spline).",
                            "Province panels: province random-walk R before zone deviations. Cut at 0.1 and 10.")) +
      theme_minimal(base_size = 9) +
      theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 7, hjust = 0))
    print(p)
  }
  dev.off()
  message(sprintf("%s %s: wrote zone_R_%s.pdf", variant, o, o))
}

q <- bind_rows(all_q)
write_csv(q, file.path(out, "zone_R_quantiles.csv.gz"))

# Summary: each zone's R at the last data week and 4 weeks later.
summ <- q |>
  group_by(origin, location) |>
  filter(target_end_date %in% c(origin[1], max(target_end_date))) |>
  arrange(target_end_date, .by_group = TRUE) |>
  summarise(level = first(level), cases = first(cases), recent = first(recent),
            R_last = first(med), R_last_lo90 = first(lo90), R_last_hi90 = first(hi90),
            R_h4 = last(med), R_h4_lo90 = last(lo90), R_h4_hi90 = last(hi90), .groups = "drop")
write_csv(summ, file.path(out, "zone_R_summary.csv"))

# Vintages: median R from every origin for the zones with most cases, plus provinces.
top <- q |> filter(level == "zone", origin == max(origin)) |> distinct(location, cases) |>
  arrange(desc(cases)) |> slice(1:15) |> pull(location)
v <- q |> filter(location %in% c(provinces, top)) |>
  mutate(location = factor(location, c(provinces, top)), origin = factor(origin))
vlab <- tibble(location = levels(v$location)) |> left_join(zone_names, by = "location") |>
  mutate(label = paste(location, coalesce(location_name, prov_names[location])))
ramp <- colorRampPalette(c("#a9c9f0", "#2a78d6", "#0b2f6b"))(nlevels(v$origin))
p <- ggplot(v, aes(target_end_date, pmin(pmax(med, 0.1), 10), colour = origin, group = origin)) +
  geom_hline(yintercept = 1, colour = "grey50", linewidth = 0.3) +
  geom_line(data = filter(v, !forecast), linewidth = 0.5) +
  geom_line(data = filter(v, target_end_date >= as.Date(as.character(origin))),
            linewidth = 0.5, linetype = "22") +
  scale_colour_manual(values = ramp) +
  scale_y_log10(breaks = c(0.3, 1, 3, 10)) +
  facet_wrap(~location, ncol = 5, labeller = as_labeller(setNames(vlab$label, vlab$location))) +
  labs(x = "Week of symptom onset (ending)", y = "Median R (log scale)", colour = "Origin",
       title = sprintf("%s: median R estimated at each origin (solid: fitted weeks; dashed: 4 forecast weeks)", variant),
       subtitle = "Provinces (province random-walk R) and the 15 zones with most cases. Cut at 0.1 and 10.") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "right", panel.grid.minor = element_blank(),
        strip.text = element_text(size = 7.5, hjust = 0))
ggsave(file.path(out, "zone_R_vintages.png"), p, width = 14, height = 9, dpi = 120)
message("Wrote ", out)
