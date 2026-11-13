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
//
// Reporting (this file; stan/bdbv_joint.stan has the lognormal delay): a case appears
// in the published data T days after onset. T has a constant daily publication rate
// within each week since onset (0, 1, 2, 3, 4, 5+): log rate = national rate for that
// week + a province x week-since-onset deviation + a province trend over onset weeks
// (cubic B-spline, random-walk coefficients starting at 0) that speeds up or slows
// down all of a province's rates together. A week's eventual count is the cases
// published within d_max days of onset, so F(age) = P(T <= age | T <= d_max) is the
// fraction of it published by a release `age` days after the week's end. Two
// likelihood terms:
//   1. Counts in the fitted (latest) release ~ NegBin(mean = onsets * F).
//   2. Reporting triangle on province totals: each earlier release's count given the
//      next release's ~ BetaBinomial(later, F(earlier) / F(later)), overdispersed for
//      reporting batches.

functions {
  // Cumulative publication hazard through day x after onset; r = daily rates for
  // weeks 0-4 since onset and 5+.
  real cum_hazard(int x, vector r) {
    int a = min(x %/% 7, 5);
    real H = 0;
    for (k in 1:a) H += 7 * r[k];
    return H + (x - 7 * a + 1) * r[a + 1];
  }

  // P(T <= x | T <= d_max) for integer x.
  real delay_cumprob(int x, vector r, int d_max) {
    if (x < 0) return 0;
    return expm1(-cum_hazard(min(x, d_max), r)) / expm1(-cum_hazard(d_max, r));
  }

  // Fraction of a week's onsets published within `age` days of the week's last day,
  // with onsets spread evenly over the week's 7 days.
  real week_completeness(int age, vector r, int d_max) {
    real total = 0;
    for (j in 0:6) total += delay_cumprob(age + j, r, d_max);
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

  // Zone deviations at the spline knots [Z, K]: stationary AR(1) with mean 0.
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
  int<lower=1> d_max;                   // days after onset that define the eventual count
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
  // Transmission. Province log R is a centered weekly random walk: the counts pin down
  // each week's value directly. Week 1 (before the outbreak, little data) is
  // non-centered around the national mean.
  real log_R0;
  real<lower=0> sd_prov_R;
  vector[P] z_prov_R;
  matrix[P, W - 1] log_R_later;
  real<lower=0> sd_rw_R;
  // Zone deviation spline: AR(1) coefficients (non-centered).
  real<lower=0, upper=0.99> rho_zone;
  real<lower=0> sd_ar_zone;
  matrix[Z, NB] eps_zone;
  real<lower=0, upper=1> phi_move;
  vector[S_n] log_seed;
  real<lower=0> inv_sqrt_nb;

  // Reporting delay: log daily publication rate by week since onset (national), province
  // deviations (non-centered), and a province trend whose spline coefficients are a
  // random walk across knots starting at 0 (non-centered).
  vector[6] log_haz;
  real<lower=0> sd_prov_haz;
  matrix[P, 6] z_prov_haz;
  real<lower=0> sd_step_haz;
  matrix[P, NB - 1] eps_haz;
  real<lower=0> inv_sqrt_kappa;         // triangle overdispersion (0 = binomial)
}

transformed parameters {
  matrix[P, W] log_R_prov = append_col(log_R0 + sd_prov_R * z_prov_R, log_R_later);
  matrix[Z, NB] a_zone = ar1_coefs(rho_zone, sd_ar_zone, eps_zone);
  matrix[P, 6] log_haz_prov = rep_matrix(log_haz', P) + sd_prov_haz * z_prov_haz;
  matrix[P, NB] c_haz = rw_coefs(0, 0, rep_vector(0, P), sd_step_haz, eps_haz);
  matrix[P, W] haz_trend = c_haz * Bs_obs_t;
  // F_v[v][p, w]: fraction of week w's eventual count (province p) published by
  // release v.
  array[V] matrix[P, W] F_v;
  for (v in 1:V)
    for (p in 1:P)
      for (w in 1:W)
        F_v[v][p, w] = week_completeness(age_v[v, w],
                                         exp(log_haz_prov[p]' + haz_trend[p, w]), d_max);
  matrix[P, W] F_week = F_v[V];
}

model {
  matrix[Z, W] zone_dev = a_zone * Bs_obs_t;
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

  // Priors
  log_R0 ~ normal(log(1.5), 0.5);
  sd_prov_R ~ normal(0, 0.3);
  z_prov_R ~ std_normal();
  to_vector(log_R_prov[, 2:W]) ~ normal(to_vector(log_R_prov[, 1:(W - 1)]), sd_rw_R);
  sd_rw_R ~ normal(0, 0.2);
  // Zone deviations: stationary sd sd_ar_zone / sqrt(1 - rho^2), about 0.5 when
  // rho = 0.8 and sd_ar_zone = 0.3. A loose sum-to-zero within each province at each
  // knot keeps the province random walk near the zones' average (a tight one leaves a
  // narrow ridge that saturates the sampler's tree depth).
  rho_zone ~ beta(4, 2);
  sd_ar_zone ~ normal(0, 0.3);
  to_vector(eps_zone) ~ std_normal();
  to_vector(prov_mean * a_zone) ~ normal(0, 0.2);
  phi_move ~ beta(2, 58);
  log_seed ~ normal(0, 1.5);
  inv_sqrt_nb ~ normal(0, 1);

  // Delay: about half of a week's cases published in the first week after onset
  // (rate 0.1 a day); each later week's rate is centred on the week before's.
  log_haz[1] ~ normal(log(0.1), 1);
  log_haz[2:6] ~ normal(log_haz[1:5], 1);
  sd_prov_haz ~ normal(0, 0.5);
  to_vector(z_prov_haz) ~ std_normal();
  sd_step_haz ~ normal(0, 0.3);
  to_vector(eps_haz) ~ std_normal();
  inv_sqrt_kappa ~ normal(0, 0.3);
}

generated quantities {
  // Eventual counts for every zone and week: observed weeks get the reported count plus
  // a draw of the cases still to be published; forecast weeks get a full draw. Province
  // log R continues its random walk into the forecast weeks; zone deviations follow the
  // spline, whose last coefficients come from the AR(1).
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
        lr += normal_rng(0, sd_rw_R);
        log_R_all[p, W + h] = lr;
      }
    }
    matrix[Z, W + H] zone_dev_all = a_zone * Bs';
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
    matrix[P, W + H] trend_all = c_haz * Bs';
    for (z in 1:Z) {
      int p = prov[z];
      for (w in 1:(W + H)) {
        real mu = expected_onsets[z, w] + eps;
        real F_now = w <= W ? F_week[p, w] : 0;
        real F_truth = fmax(week_completeness(truth_age[w],
                                              exp(log_haz_prov[p]' + trend_all[p, w]),
                                              d_max),
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
