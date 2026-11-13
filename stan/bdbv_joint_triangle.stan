// Joint nowcast/forecast model for BDBV confirmed cases by week of symptom onset.
//
// Transmission: deterministic daily SEIR in each health zone, coupled by a gravity
// mixing kernel (people spend a fraction phi_move of their time in other zones).
// Log R is a weekly random walk per province (starting from a national mean) plus a
// zone offset. With zone_rw = 1 the zone offset also follows its own weekly random walk
// (non-centered: most zones have few cases), so a zone can surge while its province
// declines.
//
// Reporting: a case appears in the published data after its onset-to-report delay D
// (discretized lognormal; log-median is a weekly random walk per province, starting
// from a national mean) plus, for a fraction pi (by province) of cases, an extra
// late-entry lag L (discretized lognormal). L covers cases entered into the data after
// their recorded report date.
// The target ("final" count) is the number published within final_age days of the
// onset week's end, and F(v) is the fraction of that published by vintage v. Without
// the cutoff, the slow trickle into months-old weeks leaves the size of the
// never-yet-published tail unidentified (fits put most cases there). Three likelihood
// terms:
//   1. Counts in the fitted (latest) vintage ~ NegBin(mean = onsets * F(latest)).
//   2. Reporting triangle: each earlier vintage's count given the next vintage's
//      ~ BetaBinomial(later, F(earlier) / F(later)), overdispersed for reporting
//      batches.
//   3. Individual delays D by province and onset week, conditioned on the case being
//      published by the latest vintage (right truncation).

functions {
  // P(delay <= x) for integer x, delay = floor(continuous lognormal delay) in days.
  real delay_cumprob(int x, real m, real s) {
    if (x < 0) return 0;
    return Phi((log(x + 1.0) - m) / s);
  }

  // log P(delay == d), computed on the side of the distribution that keeps precision.
  real delay_lpmf(int d, real m, real s) {
    real hi = (log(d + 1.0) - m) / s;
    if (d == 0) return std_normal_lcdf(hi);
    real lo = (log(d * 1.0) - m) / s;
    if (lo > 0) return log_diff_exp(std_normal_lccdf(lo), std_normal_lccdf(hi));
    return log_diff_exp(std_normal_lcdf(hi), std_normal_lcdf(lo));
  }

  // P(D = d) for d = 0..max_delay.
  vector delay_pmf(int max_delay, real m, real s) {
    vector[max_delay + 1] pmf;
    real prev = 0;
    for (d in 0:max_delay) {
      real cur = delay_cumprob(d, m, s);
      pmf[d + 1] = cur - prev;
      prev = cur;
    }
    return pmf;
  }

  // P(D + L <= x) for a late-entry case: sum over d of P(D = d) P(L <= x - d), with
  // G[k + 1] = P(L <= k).
  real late_cumprob(int x, vector pmf, vector G) {
    if (x < 0) return 0;
    int n = x;
    if (n > num_elements(pmf) - 1) n = num_elements(pmf) - 1;
    return dot_product(pmf[1:(n + 1)], reverse(G[(x - n + 1):(x + 1)]));
  }

  // Fraction of a week's onsets published within `age` days of the week's last day,
  // with onsets spread evenly over the week's 7 days.
  real week_completeness(int age, real m, real s, real pi, vector pmf, vector G) {
    real total = 0;
    for (j in 0:6)
      total += (1 - pi) * delay_cumprob(age + j, m, s) + pi * late_cumprob(age + j, pmf, G);
    return total / 7;
  }

  // P(published within `age` days of the week's end | D = d), averaged over onset day.
  real published_given_delay(int age, int d, real pi, vector G) {
    real total = 0;
    for (j in 0:6) {
      int x = age + j - d;
      if (x >= 0) total += (1 - pi) + pi * G[x + 1];
    }
    return total / 7;
  }

  // Zone log R deviation from its province [Z, n_weeks]: a static offset, plus a weekly
  // random walk when steps are given (columns of eps = standardized increments).
  matrix zone_deviation(int n_weeks, vector zone_offset, matrix eps, real sd_step) {
    int Z = rows(zone_offset);
    matrix[Z, n_weeks] dev = rep_matrix(zone_offset, n_weeks);
    if (cols(eps) > 0 && rows(eps) == Z) {
      for (w in 2:n_weeks) dev[, w] = dev[, w - 1] + sd_step * eps[, w - 1];
    }
    return dev;
  }

  // Run the SEIR for n_weeks from state (S, E, I) and return weekly onsets [Z, n_weeks]
  // and the final state [Z, 3] (S, E, I).
  tuple(matrix, matrix) seir_weekly(int n_weeks, matrix beta_zw, vector S0, vector E0,
                                    vector I0, vector N_star, matrix K, matrix Kt,
                                    real phi, real p_onset, real p_rec) {
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
    return (onsets, append_col(append_col(S, E), I));
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

  array[Z, W] int<lower=0> cases;       // counts in the fitted (latest) vintage

  int<lower=1> V;                       // vintages up to and including the fitted one
  array[V, W] int age_v;                // vintage date - week end date, in days
  int<lower=1> max_delay;               // longest delay D allowed, days
  int<lower=1> final_age;               // counts are final this many days after week end
  int<lower=0> N_tri;                   // triangle pairs (zone, week, earlier vintage)
  array[N_tri] int<lower=1, upper=Z> tri_zone;
  array[N_tri] int<lower=1, upper=W> tri_week;
  array[N_tri] int<lower=1, upper=V - 1> tri_v;
  array[N_tri] int<lower=0> tri_early;  // count in vintage tri_v
  array[N_tri] int<lower=1> tri_late;   // count in vintage tri_v + 1

  int<lower=0> S_n;                     // zones seeded at model start
  array[S_n] int<lower=1, upper=Z> seed_idx;

  int<lower=0> N_del;                   // delay cells (province, onset week, delay)
  array[N_del] int<lower=1, upper=P> del_prov;
  array[N_del] int<lower=1, upper=W> del_week;
  array[N_del] int<lower=0> del_d;
  array[N_del] int<lower=1> del_n;

  real<lower=0> mean_incubation;        // days, onset at end of E
  real<lower=0> mean_infectious;        // days, onset to end of infectiousness
  real<lower=0> eps;                    // background expected cases per zone-week

  int<lower=0, upper=1> zone_rw;        // 1 = zone offsets vary weekly
  int<lower=0, upper=1> do_gq;          // 1 to generate nowcast/forecast samples
}

transformed data {
  matrix[Z, Z] Kt = K';
  vector[Z] Ktpop = Kt * pop;
  real rec_rate = 1 / mean_infectious;
  real p_onset = -expm1(-1 / mean_incubation);
  real p_rec = -expm1(-rec_rate);
  array[W] int age = age_v[V];
  int max_age = max(max(age), final_age) + 7;
}

parameters {
  // Transmission. Province log R is a centered weekly random walk: the counts pin down
  // each week's value directly, so the centered form avoids the strong correlations a
  // non-centered sum of increments creates. Week 1 (before the outbreak, little data)
  // is non-centered around the national mean.
  real log_R0;
  real<lower=0> sd_prov_R;
  vector[P] z_prov_R;
  matrix[P, W - 1] log_R_later;
  real<lower=0> sd_rw_R;
  real<lower=0> sd_zone_R;
  matrix[zone_rw ? Z : 0, W - 1] eps_zone;
  array[zone_rw] real<lower=0> sd_rw_zone;
  vector[Z] z_zone_R;
  real<lower=0, upper=1> phi_move;
  vector[S_n] log_seed;
  real<lower=0> inv_sqrt_nb;

  // Reporting delay: province log-median is a centered weekly random walk (thousands of
  // delays identify each province-week directly). Week 1, before the outbreak, has no
  // delays, so it is non-centered around the national mean.
  real delay_mu;
  real<lower=0.1> delay_sigma;
  real<lower=0> sd_prov_d;
  vector[P] z_prov_d;
  matrix[P, W - 1] delay_m_later;
  real<lower=0> sd_rw_d;

  // Late entries: fraction pi by province (partially pooled), lag L lognormal.
  real logit_pi0;
  real<lower=0> sd_pi;
  vector[P] z_pi;
  real late_m;                          // log-median of L, days
  real<lower=0.1> late_s;
  real<lower=0> inv_sqrt_kappa;         // triangle overdispersion (0 = binomial)
}

transformed parameters {
  matrix[P, W] log_R_prov = append_col(log_R0 + sd_prov_R * z_prov_R, log_R_later);
  matrix[P, W] delay_m = append_col(delay_mu + sd_prov_d * z_prov_d, delay_m_later);
  vector[P] pi_late = inv_logit(logit_pi0 + sd_pi * z_pi);
  // G[k + 1] = P(L <= k), L = floor of a lognormal.
  vector[max_age + 1] G;
  for (k in 0:max_age) G[k + 1] = delay_cumprob(k, late_m, late_s);
  // F_v[v][p, w]: fraction of week w's final count (province p) published by vintage
  // v. F_pub_latest: fraction of all of week w's cases published by the latest vintage,
  // with no cutoff, for the right-truncation term of the delay likelihood.
  array[V] matrix[P, W] F_v;
  matrix[P, W] F_pub_latest;
  for (p in 1:P) {
    for (w in 1:W) {
      vector[max_delay + 1] pmf = delay_pmf(max_delay, delay_m[p, w], delay_sigma);
      real final = week_completeness(final_age, delay_m[p, w], delay_sigma, pi_late[p],
                                     pmf, G);
      for (v in 1:V) {
        if (age_v[v, w] >= final_age) {
          F_v[v][p, w] = 1;
        } else {
          F_v[v][p, w] = week_completeness(age_v[v, w], delay_m[p, w], delay_sigma,
                                           pi_late[p], pmf, G) / final;
        }
      }
      F_pub_latest[p, w] = week_completeness(age[w], delay_m[p, w], delay_sigma,
                                             pi_late[p], pmf, G);
    }
  }
  matrix[P, W] F_week = F_v[V];
}

model {
  matrix[Z, W] zone_dev = zone_deviation(W, sd_zone_R * z_zone_R, eps_zone,
                                         zone_rw ? sd_rw_zone[1] : 0);
  matrix[Z, W] beta_zw;
  for (z in 1:Z)
    beta_zw[z] = exp(log_R_prov[prov[z]] + zone_dev[z]) * rec_rate;

  vector[Z] E0 = rep_vector(0, Z);
  E0[seed_idx] = exp(log_seed);
  vector[Z] I0 = 0.5 * E0;
  vector[Z] N_star = (1 - phi_move) * pop + phi_move * Ktpop;
  matrix[Z, W] onsets = seir_weekly(W, beta_zw, pop - E0 - I0, E0, I0, N_star, K, Kt,
                                    phi_move, p_onset, p_rec).1;

  real nb_phi = inv_square(inv_sqrt_nb);
  for (z in 1:Z)
    cases[z] ~ neg_binomial_2((onsets[z] + eps) .* F_week[prov[z]], nb_phi);

  real kappa = inv_square(inv_sqrt_kappa);
  for (i in 1:N_tri) {
    int p = prov[tri_zone[i]];
    int w = tri_week[i];
    int v = tri_v[i];
    real ratio = fmin(fmax(F_v[v][p, w] / F_v[v + 1][p, w], 1e-9), 1 - 1e-9);
    tri_early[i] ~ beta_binomial(tri_late[i], ratio * kappa, (1 - ratio) * kappa);
  }

  for (i in 1:N_del) {
    int p = del_prov[i];
    int w = del_week[i];
    target += del_n[i] * (delay_lpmf(del_d[i] | delay_m[p, w], delay_sigma)
                          + log(published_given_delay(age[w], del_d[i], pi_late[p], G))
                          - log(F_pub_latest[p, w]));
  }

  // Priors
  log_R0 ~ normal(log(1.5), 0.5);
  sd_prov_R ~ normal(0, 0.3);
  z_prov_R ~ std_normal();
  to_vector(log_R_prov[, 2:W]) ~ normal(to_vector(log_R_prov[, 1:(W - 1)]), sd_rw_R);
  sd_rw_R ~ normal(0, 0.2);
  sd_zone_R ~ normal(0, 0.5);
  z_zone_R ~ std_normal();
  to_vector(eps_zone) ~ std_normal();
  sd_rw_zone ~ normal(0, 0.1);
  phi_move ~ beta(2, 58);
  log_seed ~ normal(0, 1.5);
  inv_sqrt_nb ~ normal(0, 1);

  delay_mu ~ normal(log(5), 0.5);
  delay_sigma ~ normal(0.8, 0.3);
  sd_prov_d ~ normal(0, 0.3);
  z_prov_d ~ std_normal();
  to_vector(delay_m[, 2:W]) ~ normal(to_vector(delay_m[, 1:(W - 1)]), sd_rw_d);
  sd_rw_d ~ normal(0, 0.1);

  logit_pi0 ~ normal(logit(0.1), 1);
  sd_pi ~ normal(0, 0.5);
  z_pi ~ std_normal();
  late_m ~ normal(log(21), 0.7);
  late_s ~ normal(1, 0.5);
  inv_sqrt_kappa ~ normal(0, 0.3);
}

generated quantities {
  // Final (fully reported) counts for every zone and week: observed weeks get the
  // reported count plus a draw of the cases still to be reported; forecast weeks get a
  // full draw. Province log R continues its random walk into the forecast weeks.
  array[do_gq ? Z : 0, do_gq ? W + H : 0] int y_final;
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
    matrix[Z, W + H - 1] eps_all = rep_matrix(0, Z, W + H - 1);
    if (zone_rw) {
      eps_all[, 1:(W - 1)] = eps_zone;
      for (z in 1:Z)
        for (h in 1:H) eps_all[z, W - 1 + h] = std_normal_rng();
    }
    matrix[Z, W + H] zone_dev_all = zone_deviation(W + H, sd_zone_R * z_zone_R,
                                                   eps_all, zone_rw ? sd_rw_zone[1] : 0);
    matrix[Z, W + H] beta_zw;
    for (z in 1:Z)
      beta_zw[z] = exp(log_R_all[prov[z]] + zone_dev_all[z]) * rec_rate;

    vector[Z] E0 = rep_vector(0, Z);
    E0[seed_idx] = exp(log_seed);
    vector[Z] I0 = 0.5 * E0;
    vector[Z] N_star = (1 - phi_move) * pop + phi_move * Ktpop;
    expected_onsets = seir_weekly(W + H, beta_zw, pop - E0 - I0, E0, I0, N_star, K, Kt,
                                  phi_move, p_onset, p_rec).1;
    R_prov = exp(log_R_all);

    for (z in 1:Z) {
      for (w in 1:W) {
        real remaining = (expected_onsets[z, w] + eps) * (1 - F_week[prov[z], w]);
        y_final[z, w] = cases[z, w] + (remaining > 1e-8 ? neg_binomial_2_rng(remaining, nb_phi) : 0);
      }
      for (h in 1:H)
        y_final[z, W + h] = neg_binomial_2_rng(expected_onsets[z, W + h] + eps, nb_phi);
    }
  }
}
