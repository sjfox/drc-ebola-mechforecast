# Synthetic-data recovery check for the joint BDBV model.
#
# Given "true" parameter values, run the model's own SEIR (via generate_quantities), then
# simulate the reporting process case by case: final counts ~ NegBin, each case gets an
# onset day within its week, a lognormal delay D and, with probability pi (by province),
# an extra late-entry lag L. A case is counted in each vintage it was published by.
# Fitting the result should recover the true values.

suppressPackageStartupMessages({
  library(dplyr)
  library(posterior)
})

# True values: a fitted draw (e.g. the posterior median) as a one-row draws_matrix of
# the model's parameters.
simulate_bdbv <- function(d, truth, model, seed = 42) {
  set.seed(seed)
  sd <- d$stan_data
  sd$do_gq <- 1L
  gq <- model$generate_quantities(fitted_params = truth, data = sd, seed = seed)
  mu <- matrix(gq$draws("expected_onsets", format = "draws_matrix")[1, ], sd$Z)
  mu <- mu[, seq_len(sd$W)]
  get <- function(v) as.numeric(truth[1, v])
  phi <- 1 / get("inv_sqrt_nb")^2
  delay_m <- matrix(get(grep("^delay_m\\[", colnames(truth), value = TRUE)), sd$P)
  delay_sigma <- get("delay_sigma")
  pi_late <- get(grep("^pi_late\\[", colnames(truth), value = TRUE))
  late_m <- get("late_m")
  late_s <- get("late_s")

  final <- matrix(rnbinom(length(mu), mu = mu + sd$eps, size = phi), sd$Z)
  counts <- array(0L, c(sd$Z, sd$W, sd$V))
  delays <- list()
  for (z in seq_len(sd$Z)) {
    p <- sd$prov[z]
    for (w in seq_len(sd$W)) {
      n <- final[z, w]
      if (n == 0) next
      days_before_end <- sample(0:6, n, replace = TRUE)
      # The final count is cases published within final_age days of the week's end
      # (as in the model), so publication times are drawn conditional on that.
      draw <- function(k) {
        list(delay = floor(rlnorm(k, delay_m[p, w], delay_sigma)),
             lag = ifelse(runif(k) < pi_late[p], floor(rlnorm(k, late_m, late_s)), 0))
      }
      x <- draw(n)
      for (try in 1:100) {
        redo <- which(x$delay + x$lag > sd$final_age + days_before_end)
        if (!length(redo)) break
        y <- draw(length(redo))
        x$delay[redo] <- y$delay
        x$lag[redo] <- y$lag
      }
      delay <- x$delay
      lag <- x$lag
      for (v in seq_len(sd$V))
        counts[z, w, v] <- sum(delay + lag <= sd$age_v[v, w] + days_before_end)
      seen <- delay + lag <= sd$age_v[sd$V, w] + days_before_end
      if (any(seen)) delays[[length(delays) + 1]] <- tibble(prov = p, w = w, d = delay[seen])
    }
  }
  max_delay <- d$cfg$max_delay
  del <- bind_rows(delays) |>
    filter(d <= max_delay) |>
    count(prov, w, d, name = "n")

  # Triangle pairs, as in build_bdbv_data (vintages that include the week).
  tri <- expand.grid(z = seq_len(sd$Z), w = seq_len(sd$W), v = seq_len(sd$V - 1)) |>
    filter(sd$age_v[cbind(v, w)] >= -6) |>
    mutate(early = counts[cbind(z, w, v)], later = counts[cbind(z, w, v + 1)]) |>
    filter(later > 0, sd$age_v[cbind(v, w)] < sd$final_age)

  sim <- d
  sim$stan_data$cases <- counts[, , sd$V]
  sim$stan_data$N_tri <- nrow(tri)
  sim$stan_data$tri_zone <- tri$z
  sim$stan_data$tri_week <- tri$w
  sim$stan_data$tri_v <- tri$v
  sim$stan_data$tri_early <- tri$early
  sim$stan_data$tri_late <- tri$later
  sim$stan_data$N_del <- nrow(del)
  sim$stan_data$del_prov <- del$prov
  sim$stan_data$del_week <- del$w
  sim$stan_data$del_d <- del$d
  sim$stan_data$del_n <- del$n
  sim$truth <- list(final = final, expected = mu, draw = truth)
  sim
}

# Fraction of true values inside the posterior 50% and 90% intervals, per variable.
recovery_table <- function(fit, truth, variables) {
  s <- fit$summary(variables, ~ quantile(.x, c(0.05, 0.25, 0.75, 0.95)))
  names(s) <- c("variable", "q05", "q25", "q75", "q95")
  s$truth <- as.numeric(truth[1, s$variable])
  s |>
    mutate(group = sub("\\[.*", "", variable),
           in50 = truth >= q25 & truth <= q75, in90 = truth >= q05 & truth <= q95) |>
    group_by(group) |>
    summarise(n = n(), cover50 = mean(in50), cover90 = mean(in90), .groups = "drop")
}

# Coverage of the true final counts by the final-count samples, for the last k weeks,
# summed to province level.
nowcast_recovery <- function(samples, sim, d, k = 6) {
  weeks <- tail(d$weeks, k)
  truth <- as.data.frame(sim$truth$final[, seq_along(d$weeks)]) |>
    setNames(as.character(d$weeks)) |>
    mutate(location = substr(d$zones$location, 1, 5)) |>
    tidyr::pivot_longer(-location, names_to = "target_end_date", values_to = "true") |>
    mutate(target_end_date = as.Date(target_end_date)) |>
    filter(target_end_date %in% weeks) |>
    group_by(location, target_end_date) |>
    summarise(true = sum(true), .groups = "drop")
  reported <- as.data.frame(sim$stan_data$cases) |>
    setNames(as.character(d$weeks)) |>
    mutate(location = substr(d$zones$location, 1, 5)) |>
    tidyr::pivot_longer(-location, names_to = "target_end_date", values_to = "reported") |>
    mutate(target_end_date = as.Date(target_end_date)) |>
    group_by(location, target_end_date) |>
    summarise(reported = sum(reported), .groups = "drop")
  samples |>
    filter(target_end_date %in% weeks) |>
    mutate(location = substr(location, 1, 5)) |>
    group_by(location, target_end_date, sample) |>
    summarise(value = sum(value), .groups = "drop") |>
    group_by(location, target_end_date) |>
    summarise(q05 = quantile(value, 0.05), median = median(value),
              q95 = quantile(value, 0.95), .groups = "drop") |>
    left_join(truth, by = c("location", "target_end_date")) |>
    left_join(reported, by = c("location", "target_end_date")) |>
    mutate(in90 = true >= q05 & true <= q95)
}
