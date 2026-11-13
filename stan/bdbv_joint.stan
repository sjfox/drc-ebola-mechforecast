// Joint nowcast/forecast model for BDBV confirmed cases by week of symptom onset.
// (The previous version, with the delay file and a late-entry mixture, is
// stan/bdbv_joint_triangle.stan.)
//
// Transmission: deterministic daily SEIR in each health zone, coupled by a gravity
// mixing kernel (people spend a fraction phi_move of their time in other zones).
// Log R is a weekly random walk per province (starting from a national mean) plus a
// zone deviation that varies smoothly in time: a cubic B-spline whose coefficients
// follow an AR(1) towards 0, i.e. towards the province average. The spline covers the
// forecast weeks too, so the deviation reverts towards the province in forecasts.
// With zone_weekly, the zone deviation is instead a weekly AR(1) towards 0 (no spline),
// and forecasts continue that same AR(1).
//
// Reporting: a case appears in the published data T days after onset, T a discretized
// lognormal with no upper limit. Its log-median and log-spread vary smoothly by
// province and week (cubic B-splines on the same knots, random-walk coefficients).
// F(age) = P(T <= age) is the fraction of a week's eventual count published by a
// release `age` days after the week's end. Two likelihood terms:
//   1. Counts in the fitted (latest) release ~ NegBin(mean = onsets * F).
//   2. Reporting triangle on province totals: each earlier release's count given the
//      next release's ~ BetaBinomial(later, F(earlier) / F(later)), overdispersed for
//      reporting batches.

functions {
  // P(T <= x) for integer x, T = floor(continuous lognormal) in days.
  real delay_cumprob(int x, real m, real s) {
    if (x < 0) return 0;
    return Phi((log(x + 1.0) - m) / s);
  }

  // Fraction of a week's onsets published within `age` days of the week's last day,
  // with onsets spread evenly over the week's 7 days.
  real week_completeness(int age, real m, real s) {
    real total = 0;
    for (j in 0:6) total += delay_cumprob(age + j, m, s);
    return total / 7;
  }

  // Coefficients of a non-centered random walk [P, K] with a partially pooled start.
  matrix rw_coefs(real mu, real sd_start, vector z_start, real sd_step, matrix eps) {
    int P = rows(eps);
    int K = cols(eps) + 1;
    matrix[P, K] c;
    c[, 1] = mu + sd_start * z_start;
    for (k in 2:K) c[, k] = c[, k - 1] + sd_step * eps[, k - 1];
    return c;
  }

  // Zone deviations at the spline knots (or weeks, with zone_weekly) [Z, K]: stationary
  // AR(1) with mean 0.
  matrix ar1_coefs(real rho, real sd_step, matrix eps) {
    int Z = rows(eps);
    int K = cols(eps);
    matrix[Z, K] a;
    a[, 1] = sd_step / sqrt(1 - square(rho)) * eps[, 1];
    for (k in 2:K) a[, k] = rho * a[, k - 1] + sd_step * eps[, k];
    return a;
  }

  // Run the SEIR for n_weeks from state (S, E, I) and return weekly onsets [Z, n_weeks].
  matrix seir_weekly(int n_weeks, matrix beta_zw, vector S0, vector E0, vector I0,
                     vector N_star, matrix K, matrix Kt, real phi, real p_onset,
                     real p_rec) {
    int Z = rows(S0);
    matrix[Z, n_weeks] onsets = rep_matrix(0, Z, n_weeks);
    vector[Z] S = S0;
    vector[Z] E = E0;
    vector[Z] I = I0;
    for (w in 1:n_weeks) {
      vector[Z] b = beta_zw[, w];
      for (day in 1:7) {
        // Infectious people present in each zone, then the force of infection on
        // residents of each zone across the zones they spend time in.
        vector[Z] prev = ((1 - phi) * I + phi * (Kt * I)) ./ N_star;
        vector[Z] bp = b .* prev;
        vector[Z] lambda = (1 - phi) * bp + phi * (K * bp);
        vector[Z] new_inf = S .* -expm1(-lambda);
        vector[Z] ons = E * p_onset;
        vector[Z] rec = I * p_rec;
        S -= new_inf;
        E += new_inf - ons;
        I += ons - rec;
        onsets[, w] += ons;
      }
    }
    return onsets;
  }
}

data {
  int<lower=1> Z;                       // health zones
  int<lower=1> P;                       // provinces
  int<lower=1> W;                       // observed weeks (reference weeks)
  int<lower=0> H;                       // forecast weeks after the last observed week
  array[Z] int<lower=1, upper=P> prov;
  vector<lower=1>[Z] pop;
  matrix<lower=0>[Z, Z] K;              // row = home zone, column = zone visited

  array[Z, W] int<lower=0> cases;       // counts in the fitted (latest) release

  int<lower=4> NB;                      // spline basis functions
  matrix[W + H, NB] Bs;                 // cubic B-spline basis over all weeks

  int<lower=1> V;                       // releases up to and including the fitted one
  array[V, W] int age_v;                // release date - week end date, in days
  int<lower=0> N_tri;                   // triangle pairs (province, week, earlier release)
  array[N_tri] int<lower=1, upper=P> tri_prov;
  array[N_tri] int<lower=1, upper=W> tri_week;
  array[N_tri] int<lower=1, upper=V - 1> tri_v;
  array[N_tri] int<lower=0> tri_early;  // count in release tri_v
  array[N_tri] int<lower=1> tri_late;   // count in release tri_v + 1

  int<lower=0> S_n;                     // zones seeded at model start
  array[S_n] int<lower=1, upper=Z> seed_idx;

  real<lower=0> mean_incubation;        // days, onset at end of E
  real<lower=0> mean_infectious;        // days, onset to end of infectiousness
  real<lower=0> eps;                    // background expected cases per zone-week

  int<lower=0, upper=1> do_gq;          // 1 to generate nowcast/forecast samples
  array[W + H] int truth_age;           // ages at a later release, for y_truth (gq)

  // Finalized-data backtests: every observed week is complete (F = 1); the delay
  // parameters then sample from their priors and don't affect the fit.
  int<lower=0, upper=1> perfect_reporting;
  // Province / zone R settings.
  vector<lower=0>[2] rho_zone_prior;    // beta(a, b) prior on rho_zone
  real<lower=0> sd_ar_zone_scale;       // half-normal scale for sd_ar_zone
  real<lower=0> sd_rw_R_scale;          // half-normal scale for sd_rw_R
  int<lower=0, upper=1> forecast_R_flat; // 1: province log R stays at its last value
  // 1: each zone's log R in forecast week h is damp_R^h times its value in the last data
  // week (R decays geometrically toward 1; damp_R = 1 holds it, 0 sets R = 1)
  int<lower=0, upper=1> forecast_R_damped;
  real<lower=0, upper=1> damp_R;
  // 1: province log R is a mean-reverting AR(1) toward 0 (R = 1) with persistence phi_R
  // estimated in the fit, instead of a random walk. Forecast weeks then continue the
  // AR(1), and zone deviations decay at the fitted zone persistence (rho_zone per
  // knot_weeks weeks).
  int<lower=0, upper=1> ar_R;
  vector<lower=0>[2] phi_R_prior;       // beta(a, b) prior on phi_R
  int<lower=0, upper=1> forecast_R_mean; // ar_R: log R decays as phi_R^h, no innovations
  int<lower=1> knot_weeks;
  // 1: zone deviations are a weekly AR(1) (rho_zone per week) instead of the spline;
  // forecast weeks continue it (its mean rho_zone^h with forecast_R_mean).
  int<lower=0, upper=1> zone_weekly;
}

transformed data {
  matrix[Z, Z] Kt = K';
  vector[Z] Ktpop = Kt * pop;
  real rec_rate = 1 / mean_infectious;
  // Daily exit probabilities 1/mean give geometric stage durations with exactly these
  // means.
  real p_onset = 1 / mean_incubation;
  real p_rec = 1 / mean_infectious;
  matrix[NB, W] Bs_obs_t = Bs[1:W]';
  // Each province's zones, for the soft sum-to-zero on zone deviations.
  matrix[P, Z] prov_mean = rep_matrix(0, P, Z);
  for (z in 1:Z) prov_mean[prov[z], z] = 1;
  for (p in 1:P) prov_mean[p] /= sum(prov_mean[p]);
}

parameters {
  // Transmission. Province log R is a centered weekly random walk (or, with ar_R, an
  // AR(1) toward 0): the counts pin down each week's value directly. Week 1 (before the
  // outbreak, little data) is non-centered around the national mean.
  real log_R0;
  real<lower=0> sd_prov_R;
  vector[P] z_prov_R;
  matrix[P, W - 1] log_R_later;
  real<lower=0> sd_rw_R;
  real<lower=0, upper=1> phi_R;         // AR(1) persistence (used only with ar_R)
  // Zone deviation spline: AR(1) coefficients (non-centered); with zone_weekly, the
  // weekly deviations themselves.
  real<lower=0, upper=0.99> rho_zone;
  real<lower=0> sd_ar_zone;
  matrix[Z, zone_weekly ? W : NB] eps_zone;
  real<lower=0, upper=1> phi_move;
  vector[S_n] log_seed;
  real<lower=0> inv_sqrt_nb;

  // Reporting delay: spline coefficients for the log-median and log-spread, random walks
  // across knots per province (non-centered), each starting around a national mean.
  real delay_mu;
  real<lower=0> sd_prov_m;
  vector[P] z_prov_m;
  real<lower=0> sd_step_m;
  matrix[P, NB - 1] eps_m;
  real log_s_mu;
  real<lower=0> sd_prov_s;
  vector[P] z_prov_s;
  real<lower=0> sd_step_s;
  matrix[P, NB - 1] eps_s;
  real<lower=0> inv_sqrt_kappa;         // triangle overdispersion (0 = binomial)
}

transformed parameters {
  matrix[P, W] log_R_prov = append_col(log_R0 + sd_prov_R * z_prov_R, log_R_later);
  matrix[Z, zone_weekly ? W : NB] a_zone = ar1_coefs(rho_zone, sd_ar_zone, eps_zone);
  matrix[P, NB] c_m = rw_coefs(delay_mu, sd_prov_m, z_prov_m, sd_step_m, eps_m);
  matrix[P, NB] c_s = rw_coefs(log_s_mu, sd_prov_s, z_prov_s, sd_step_s, eps_s);
  matrix[P, W] delay_m = c_m * Bs_obs_t;
  matrix[P, W] delay_s = exp(c_s * Bs_obs_t);
  // F_v[v][p, w]: fraction of week w's eventual count (province p) published by
  // release v.
  array[V] matrix[P, W] F_v;
  for (v in 1:V)
    for (p in 1:P)
      for (w in 1:W)
        F_v[v][p, w] = perfect_reporting
                         ? 1 : week_completeness(age_v[v, w], delay_m[p, w], delay_s[p, w]);
  matrix[P, W] F_week = F_v[V];
}

model {
  matrix[Z, W] zone_dev;
  if (zone_weekly) zone_dev = a_zone;
  else zone_dev = a_zone * Bs_obs_t;
  matrix[Z, W] beta_zw;
  for (z in 1:Z)
    beta_zw[z] = exp(log_R_prov[prov[z]] + zone_dev[z]) * rec_rate;

  vector[Z] E0 = rep_vector(0, Z);
  E0[seed_idx] = exp(log_seed);
  vector[Z] I0 = 0.5 * E0;
  vector[Z] N_star = (1 - phi_move) * pop + phi_move * Ktpop;
  matrix[Z, W] onsets = seir_weekly(W, beta_zw, pop - E0 - I0, E0, I0, N_star, K, Kt,
                                    phi_move, p_onset, p_rec);

  real nb_phi = inv_square(inv_sqrt_nb);
  for (z in 1:Z)
    cases[z] ~ neg_binomial_2((onsets[z] + eps) .* F_week[prov[z]], nb_phi);

  real kappa = inv_square(inv_sqrt_kappa);
  for (i in 1:N_tri) {
    int p = tri_prov[i];
    int w = tri_week[i];
    int v = tri_v[i];
    real ratio = fmin(fmax(F_v[v][p, w] / F_v[v + 1][p, w], 1e-9), 1 - 1e-9);
    tri_early[i] ~ beta_binomial(tri_late[i], ratio * kappa, (1 - ratio) * kappa);
  }

  // Priors (defaults in bdbv_config(): rho_zone ~ beta(4, 2), sd_ar_zone ~ N+(0, 0.3),
  // sd_rw_R ~ N+(0, 0.2))
  log_R0 ~ normal(log(1.5), 0.5);
  sd_prov_R ~ normal(0, 0.3);
  z_prov_R ~ std_normal();
  {
    vector[P * (W - 1)] prev = to_vector(log_R_prov[, 1:(W - 1)]);
    to_vector(log_R_prov[, 2:W]) ~ normal(ar_R ? phi_R * prev : prev, sd_rw_R);
  }
  sd_rw_R ~ normal(0, sd_rw_R_scale);
  phi_R ~ beta(phi_R_prior[1], phi_R_prior[2]);
  // Zone deviations: stationary sd sd_ar_zone / sqrt(1 - rho^2), about 0.5 when
  // rho = 0.8 and sd_ar_zone = 0.3. A loose sum-to-zero within each province at each
  // knot keeps the province random walk near the zones' average (a tight one leaves a
  // narrow ridge that saturates the sampler's tree depth).
  rho_zone ~ beta(rho_zone_prior[1], rho_zone_prior[2]);
  sd_ar_zone ~ normal(0, sd_ar_zone_scale);
  to_vector(eps_zone) ~ std_normal();
  to_vector(prov_mean * a_zone) ~ normal(0, 0.2);
  phi_move ~ beta(2, 58);
  log_seed ~ normal(0, 1.5);
  inv_sqrt_nb ~ normal(0, 1);

  // Delay: median around 10 days; spread around 1 on the log scale. The spread changes
  // little between knots unless the data require it.
  delay_mu ~ normal(log(10), 0.7);
  sd_prov_m ~ normal(0, 0.3);
  z_prov_m ~ std_normal();
  sd_step_m ~ normal(0, 0.3);
  to_vector(eps_m) ~ std_normal();
  log_s_mu ~ normal(0, 0.4);
  sd_prov_s ~ normal(0, 0.2);
  z_prov_s ~ std_normal();
  sd_step_s ~ normal(0, 0.1);
  to_vector(eps_s) ~ std_normal();
  inv_sqrt_kappa ~ normal(0, 0.3);
}

generated quantities {
  // Eventual counts for every zone and week: observed weeks get the reported count plus
  // a draw of the cases still to be published; forecast weeks get a full draw. Province
  // log R continues its random walk into the forecast weeks (or stays at its last value
  // with forecast_R_flat); zone deviations follow the spline, whose last coefficients
  // come from the AR(1). With forecast_R_damped, each zone's log R instead decays
  // geometrically from its last-data-week value. With ar_R, province log R continues
  // its fitted AR(1) (or its mean, phi_R^h, with forecast_R_mean) and zone deviations
  // decay at the fitted zone persistence.
  // y_truth: the count as it will appear in a later release (ages truth_age), for
  // backtests.
  array[do_gq ? Z : 0, do_gq ? W + H : 0] int y_final;
  array[do_gq ? Z : 0, do_gq ? W + H : 0] int y_truth;
  matrix[do_gq ? Z : 0, do_gq ? W + H : 0] expected_onsets;
  matrix[do_gq ? P : 0, do_gq ? W + H : 0] R_prov;
  if (do_gq) {
    real nb_phi = inv_square(inv_sqrt_nb);
    matrix[P, W + H] log_R_all;
    log_R_all[, 1:W] = log_R_prov;
    for (p in 1:P) {
      real lr = log_R_prov[p, W];
      for (h in 1:H) {
        if (ar_R) {
          lr = forecast_R_mean ? phi_R^h * log_R_prov[p, W]
                               : phi_R * lr + normal_rng(0, sd_rw_R);
        } else if (forecast_R_damped) {
          lr = damp_R^h * log_R_prov[p, W];
        } else if (!forecast_R_flat) {
          lr += normal_rng(0, sd_rw_R);
        }
        log_R_all[p, W + h] = lr;
      }
    }
    matrix[Z, W + H] zone_dev_all;
    if (zone_weekly) {
      // Weekly AR(1): forecasts continue it (its mean with forecast_R_mean).
      zone_dev_all[, 1:W] = a_zone;
      vector[Z] d = a_zone[, W];
      for (h in 1:H) {
        d = forecast_R_mean ? rho_zone^h * a_zone[, W]
                            : rho_zone * d + to_vector(normal_rng(rep_vector(0, Z), sd_ar_zone));
        zone_dev_all[, W + h] = d;
      }
    } else {
      zone_dev_all = a_zone * Bs';
      // Damped: zone deviations decay from their last-data-week value too, instead of
      // following the spline (whose last coefficient has no data behind it).
      if (ar_R) {
        for (h in 1:H)
          zone_dev_all[, W + h] = rho_zone^(h * 1.0 / knot_weeks) * zone_dev_all[, W];
      } else if (forecast_R_damped) {
        for (h in 1:H) zone_dev_all[, W + h] = damp_R^h * zone_dev_all[, W];
      }
    }
    matrix[Z, W + H] beta_zw;
    for (z in 1:Z)
      beta_zw[z] = exp(log_R_all[prov[z]] + zone_dev_all[z]) * rec_rate;

    vector[Z] E0 = rep_vector(0, Z);
    E0[seed_idx] = exp(log_seed);
    vector[Z] I0 = 0.5 * E0;
    vector[Z] N_star = (1 - phi_move) * pop + phi_move * Ktpop;
    expected_onsets = seir_weekly(W + H, beta_zw, pop - E0 - I0, E0, I0, N_star, K, Kt,
                                  phi_move, p_onset, p_rec);
    R_prov = exp(log_R_all);

    // Delay for every week (forecast weeks use the spline's continuation).
    matrix[P, W + H] m_all = c_m * Bs';
    matrix[P, W + H] s_all = exp(c_s * Bs');
    for (z in 1:Z) {
      int p = prov[z];
      for (w in 1:(W + H)) {
        real mu = expected_onsets[z, w] + eps;
        real F_now = w <= W ? F_week[p, w] : 0;
        real F_truth = perfect_reporting
                         ? 1
                         : fmax(week_completeness(truth_age[w], m_all[p, w], s_all[p, w]),
                                F_now);
        int obs = w <= W ? cases[z, w] : 0;
        // Cases published by the truth release, then the rest, as one split of the
        // remaining eventual count.
        int add_truth = mu * (F_truth - F_now) > 1e-8
                          ? neg_binomial_2_rng(mu * (F_truth - F_now), nb_phi) : 0;
        int add_rest = mu * (1 - F_truth) > 1e-8
                         ? neg_binomial_2_rng(mu * (1 - F_truth), nb_phi) : 0;
        y_truth[z, w] = obs + add_truth;
        y_final[z, w] = obs + add_truth + add_rest;
      }
    }
  }
}
