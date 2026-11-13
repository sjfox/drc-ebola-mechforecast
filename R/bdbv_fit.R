# Fit the joint BDBV model with cmdstanr and turn posterior draws into hub-format
# sample trajectories. Use with R/bdbv_data.R; scripts/run_bdbv.R runs both end to end.

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(tidyr)
})

bdbv_model <- function(root = ".", delay = "lognormal") {
  file <- if (delay == "hazard") "stan/bdbv_joint_hazard.stan" else "stan/bdbv_joint.stan"
  cmdstan_model(file.path(root, file))
}

delay_params <- function(delay) {
  if (delay == "hazard") c("log_haz", "sd_prov_haz", "sd_step_haz")
  else c("delay_mu", "sd_prov_m", "sd_step_m", "log_s_mu", "sd_prov_s", "sd_step_s")
}

# Start every chain near the prior median, with small random effects, so the SEIR
# starts in a sensible regime (random inits can give explosive or extinct epidemics).
bdbv_inits <- function(sd, delay = "lognormal") {
  ncol_zone <- if (sd$zone_weekly) sd$W else sd$NB
  function() {
    init <- list(
      log_R0 = log(1.3) + rnorm(1, 0, 0.05),
      sd_prov_R = 0.1, z_prov_R = rnorm(sd$P, 0, 0.1),
      log_R_later = matrix(log(1.3) + rnorm(sd$P * (sd$W - 1), 0, 0.02), sd$P),
      sd_rw_R = 0.1, phi_R = 0.7,
      rho_zone = 0.8, sd_ar_zone = 0.2,
      eps_zone = matrix(rnorm(sd$Z * ncol_zone, 0, 0.1), sd$Z),
      phi_move = 0.03, log_seed = rnorm(sd$S_n, log(2), 0.1), inv_sqrt_nb = 0.5,
      inv_sqrt_kappa = 0.3
    )
    delay_init <- if (delay == "hazard") {
      list(log_haz = log(c(0.1, 0.08, 0.05, 0.03, 0.02, 0.01)), sd_prov_haz = 0.1,
           z_prov_haz = matrix(rnorm(sd$P * 6, 0, 0.1), sd$P), sd_step_haz = 0.1,
           eps_haz = matrix(rnorm(sd$P * (sd$NB - 1), 0, 0.1), sd$P))
    } else {
      list(delay_mu = log(10), sd_prov_m = 0.1, z_prov_m = rnorm(sd$P, 0, 0.1),
           sd_step_m = 0.1, eps_m = matrix(rnorm(sd$P * (sd$NB - 1), 0, 0.1), sd$P),
           log_s_mu = 0, sd_prov_s = 0.05, z_prov_s = rnorm(sd$P, 0, 0.1),
           sd_step_s = 0.05, eps_s = matrix(rnorm(sd$P * (sd$NB - 1), 0, 0.1), sd$P))
    }
    c(init, delay_init)
  }
}

fit_bdbv <- function(d, model = bdbv_model(d$cfg$root, d$cfg$delay), chains = 4, iter_warmup = 500,
                     iter_sampling = 500, adapt_delta = 0.9, max_treedepth = 10,
                     seed = 20261101, use_pathfinder = TRUE, ...) {
  sd <- d$stan_data
  init <- bdbv_inits(sd, d$cfg$delay)
  if (use_pathfinder) {
    # Pathfinder gives approximate posterior draws to start the chains from. If it
    # fails, fall back to the prior-median inits.
    pf <- tryCatch(
      model$pathfinder(data = sd, init = replicate(chains, init(), simplify = FALSE),
                       num_paths = chains, psis_resample = FALSE, seed = seed,
                       refresh = 0, show_messages = FALSE),
      error = function(e) {
        message("Pathfinder failed (", conditionMessage(e), "); using default inits.")
        NULL
      }
    )
    # Paths that collapse to one point give too few distinct draws to start from.
    if (!is.null(pf) && nrow(unique(pf$draws(format = "draws_matrix"))) >= chains) {
      init <- pf
    }
  }
  model$sample(
    data = sd, init = init, chains = chains, parallel_chains = chains,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling, adapt_delta = adapt_delta,
    max_treedepth = max_treedepth, seed = seed, refresh = 100, ...
  )
}

# Draw-level nowcast/forecast samples for `n_samples` posterior draws.
generate_samples <- function(fit, d, model = bdbv_model(d$cfg$root, d$cfg$delay), n_samples = 100,
                             seed = 1) {
  sd <- d$stan_data
  sd$do_gq <- 1L
  draws <- fit$draws(format = "draws_matrix")
  set.seed(seed)
  keep <- sort(sample(ndraws(draws), n_samples))
  gq <- model$generate_quantities(fitted_params = draws[keep, ], data = sd, seed = seed)

  all_weeks <- c(d$weeks, d$forecast_weeks)
  tidy <- function(var) {
    y <- gq$draws(var, format = "draws_matrix")
    idx <- do.call(rbind, regmatches(colnames(y), gregexpr("[0-9]+", colnames(y))))
    tibble(
      sample = rep(seq_len(nrow(y)), times = ncol(y)),
      location = rep(d$zones$location[as.integer(idx[, 1])], each = nrow(y)),
      target_end_date = rep(all_weeks[as.integer(idx[, 2])], each = nrow(y)),
      value = as.vector(y)
    )
  }
  # Eventual counts; with a truth release set (backtests), also the counts as they will
  # appear in it.
  out <- tidy("y_final")
  attr(out, "truth") <- tidy("y_truth")
  out
}

# Hub sample file covering every location. Modeled zones come from the model; zones
# outside the modeled provinces keep their reported counts for past weeks and get 0 for
# future weeks (they have had at most a handful of cases). Provinces and the country
# are sums over zones within each sample, so all levels agree.
hub_samples <- function(samples, d, forecast_date = Sys.Date()) {
  out_weeks <- d$ref_date + 7 * seq(-d$cfg$output_weeks_before_ref,
                                    d$cfg$forecast_weeks_after_ref)
  n <- max(samples$sample)
  modeled <- filter(samples, target_end_date %in% out_weeks)

  other_zones <- d$snapshot |>
    filter(location_type == "health zone", !location %in% d$zones$location)
  fallback <- expand_grid(location = unique(other_zones$location),
                          target_end_date = out_weeks, sample = seq_len(n)) |>
    left_join(select(other_zones, location, target_end_date, value = observation),
              by = c("location", "target_end_date")) |>
    mutate(value = coalesce(value, 0))

  zone <- bind_rows(modeled, fallback)
  province <- zone |>
    mutate(location = substr(location, 1, 5)) |>
    group_by(location, target_end_date, sample) |>
    summarise(value = sum(value), .groups = "drop")
  country <- zone |>
    group_by(target_end_date, sample) |>
    summarise(value = sum(value), .groups = "drop") |>
    mutate(location = "CD")

  bind_rows(country, province, zone) |>
    transmute(
      reference_date = d$ref_date, target = "wk inc confirmed cases", target_end_date,
      location, output_type = "sample", output_type_id = sample, value = as.numeric(value),
      forecast_date = as.Date(forecast_date)
    ) |>
    arrange(location, output_type_id, target_end_date)
}

diagnostics_table <- function(fit, delay = "lognormal") {
  s <- fit$summary(
    variables = intersect(
      c("log_R0", "sd_prov_R", "sd_rw_R", "phi_R", "rho_zone", "sd_ar_zone", "phi_move",
        "inv_sqrt_nb", delay_params(delay), "inv_sqrt_kappa", "lp__"),
      c(fit$metadata()$stan_variables, "lp__"))
  )
  d <- fit$diagnostic_summary(quiet = TRUE)
  list(summary = s, divergences = sum(d$num_divergent),
       max_treedepth_hits = sum(d$num_max_treedepth), ebfmi = d$ebfmi)
}
