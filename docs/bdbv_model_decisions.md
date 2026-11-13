# BDBV joint nowcast/forecast model: decisions to review

A record of every judgment call made while building the first version of the joint
transmission + reporting-delay model (R + Stan). Each entry says what was decided, why,
and what to change if you disagree. Entries marked **[review]** are the ones most likely
to matter for forecast quality.

Files:

| File | Role |
|---|---|
| `src/fcast/bdbv_inputs.py` | Python: zone mapping, gravity kernel, snapshot of hub auxiliary data |
| `R/bdbv_data.R` | Builds the Stan data from one data vintage |
| `stan/bdbv_joint.stan` | The joint model |
| `R/bdbv_fit.R` | Fitting, posterior samples, hub-format output |
| `R/bdbv_plots.R` | Diagnostic plots |
| `R/bdbv_synthetic.R`, `scripts/synthetic_check_bdbv.R` | Synthetic recovery check |
| `scripts/run_bdbv.R` | End-to-end run for one vintage |

## Data

**D1. Two data sources, factored so neither is used twice.**
`time-series.csv` gives counts by onset week (one row per vintage). The hub's
`weekly-reporting-delays-by-province.csv` gives each case's onset-to-report delay by
province and onset week. The model uses (a) the latest vintage's counts and (b) the
individual delays, conditioned on being reported by the as-of date. For a line list
these two are exactly the likelihood factored into "how many reported so far" and "when
each was reported", so delay information isn't counted twice. Earlier vintages are
*not* used in the fit; they're for checking the nowcast and for backtests.
*Alternative:* fit to increments between vintages instead of the delay file. There are
only 5 vintages (one in August, four weekly from Oct 9), so the delay file is far more
informative.

**D2. [review] Delays vary by province, not by health zone.** The delay file has no
zone identifier, so zone-level delay differences aren't identifiable from it. Every
zone uses its province's delay distribution. The hierarchy is country → province, plus
national week and province-week effects. Adding zone effects would need zone-level
delays (the daily line list, if available) or the vintage increments.

**D3. Hub zone codes mapped to GRID3 by DHIS2 code, then by name.** 498 of 519 match
by DHIS2 code and 21 by normalized name within province (`NAME_OVERRIDES` in
`bdbv_inputs.py` lists hand-matched spellings). The hub mapping gives Bandalungwa and
Kinshasa the same DHIS2 code, so DHIS2 codes listed twice are ignored and those two
match by name. Manika is in Haut-Katanga in the hub but Lualaba in GRID3; the model
uses hub provinces throughout.

**D4. Data cutoff = as_of date** (`report_lag_days = 0`). The delay files are
consistent with cases reported up to the release date. If the backtest shows the
newest week is consistently over-nowcast, raise this to 1–2.

**D5. Negative delays (35 cases) and delays over 90 days (12) are dropped** as data
errors, out of 8,811. Delays the onset week couldn't produce by the as-of date are also
dropped (none currently).

**D6. Cases with imputed onset dates are treated like observed ones.** The hub imputes
onsets for ~15% of cases from its own fitted delay distribution. They appear in the
counts but not in the delay file. Since they were imputed from a delay distribution,
they should roughly follow the same completeness.

**D7. Model window starts 3 weeks before the outbreak.** Outbreak start = first week of
two consecutive weeks with cases nationally (week ending 2026-04-25), so the model
starts with the week ending 2026-04-04. Two isolated cases (weeks ending 02-14 and
03-07) are before the window and ignored.

## Geography and movement

**G1. [review] Only provinces with ≥10 cumulative cases are modeled.** Currently
Ituri, Nord-Kivu, Haut-Uele, Tshopo, Bas-Uele: 117 zones, 10,408 of 10,417 cases.
Zones elsewhere (7 cases total across Sud-Kivu, Sud-Ubangi, Nord-Ubangi) appear in the
hub output with their reported counts for past weeks and 0 for future weeks. Change
`min_province_cases` to widen. Fitting all 519 zones should work but is slower (the
mixing step scales with zones²).

**G2. Gravity kernel restricted to the modeled zones and renormalized.** The national
kernel is too diffuse to sparsify (keeping 95% of each zone's travel takes ~260
destinations), so the model uses the dense kernel over modeled zones. On average 77%
of a modeled zone's travel stays within the modeled provinces (min 32%); the rest is
dropped. Kernel exponents match the epymorph `Gravity` defaults (population^1, (travel
time + 1 h)^−2).

**G3. Movement as mixing, not explicit commuting.** People spend a fraction `phi_move`
of their time in other zones, distributed by the kernel. Prior Beta(2, 58), mean ≈
0.033, which matches the epymorph `Gravity` default (10% of people travel for 1/3 of a
day).

## Transmission

**T1. [review] Deterministic SEIR, noise only in the observations.** HMC can't sample
discrete counts of people. Uncertainty comes from the R random walk and a negative
binomial on counts. This is weakest for zones with a handful of cases, where chance
extinction and sporadic spread matter.

**T2. [review] Fixed natural history: incubation 7.5 days, infectious period 8 days after onset, both exponential.**
Together they set the mean generation interval, 15.5 days, which links the growth rate to R. The values were 9 and 7 days before 2026-11-04, and every fit so far used those.

*Incubation: use the Bundibugyo estimates.*

| Source | Outbreak | n | Estimate |
|---|---|---|---|
| [Wamala et al. 2010](https://doi.org/10.3201/eid1607.091525) | Uganda 2007 | 116 | median 7 days (range 2–20); mean 7.68 (SD 3.53) as re-reported by a [2026 airport-screening preprint](https://www.medrxiv.org/content/10.64898/2026.06.11.26355442v1.full) |
| [MacNeil et al. 2010](https://doi.org/10.3201/eid1612.100627) | Uganda 2007 | 24 | mean 6.3 days (95% CI 5.2–7.3) |
| [Kratz et al. 2015](https://doi.org/10.1371/journal.pone.0129333) | Isiro, DRC 2012 | 3 | mean 11.3 days (2, 5 and 27 days) |

- I used 7.5 days, about the case-weighted mean of the two Uganda studies. The Isiro estimate is too small to use.
- It is a little shorter than the pooled Ebola estimate: 8.5 days (95% CI 7.7–9.2), mostly Zaire, from [Nash et al. 2024](https://doi.org/10.1016/S1473-3099(24)00374-8). That review found incubation differs significantly between species (p = 0.02).

*Infectious period: no Bundibugyo estimate exists, so it is set from the Zaire serial interval.*
- No Bundibugyo serial interval, generation time or infectious period has been published. [Bongiovanni et al. 2026](https://doi.org/10.1016/j.ijid.2026.109037) confirms this.
- The most comprehensive estimates are from the Nash et al. 2024 meta-analysis (522 papers, almost all Zaire):
  - serial interval 15.4 days (95% CI 13.2–17.5), pooling 19 Zaire and 2 Sudan estimates;
  - this matches the WHO Ebola Response Team's 15.3 ± 9.3 days from West Africa ([NEJM 2014](https://doi.org/10.1056/NEJMoa1411100)).
- In this SEIR the mean generation interval is the mean incubation plus the mean infectious period, so 15.4 − 7.5 ≈ 8 days.
- Nash et al. also pooled an infectious period of 5.0 days (3.7–6.3). I didn't use it:
  - it rests on only three low-quality estimates;
  - with 7.5 days of incubation it would give a 12.5-day generation interval, which disagrees with the serial interval.
- An effective 8 days after onset fits other evidence:
  - Bundibugyo onset to death is 9–11 days. Examples: Wamala, a median of 10 days; [Imperial 2026 report](https://www.imperial.ac.uk/media/imperial-college/medicine/mrc-gida/Report-ebola-update-20-05-2026.pdf), mean 11.4 (SD 5.4).
  - Unsafe burials continue transmission after death.
  - Onset to admission in Isiro 2012 was about 4 days (Funk & Abbott 2026, cited in the preprint). That bounds how quickly a case is removed when the response works well.

*Caveats.*
- Exponential stages are far more dispersed than the data. An exponential with mean 7.5 has SD 7.5, while the observed SD is about 3.5. Splitting E (and I) into two or more stages each would fix this.
- The values are fixed. Turning them into informative priors would put uncertainty in the generation interval into R and the forecasts:
  - incubation ~ N(7.5, 0.7);
  - generation interval ~ N(15.4, 1.1).

**T3. Onset = leaving E.** The target is cases by week of symptom onset, so the
modeled observable is the E → I flow. Ascertainment is fixed at 1: the model describes
the confirmed-case epidemic, since ascertainment and transmission can't be separated
from case counts alone.

**T4. [review] R structure: log R(zone, week) = province weekly random walk (week 1
drawn around a national mean) + static zone offset.** The random walk is
parameterized *centered* (each province-week's log R is a parameter); see I5. Weekly
step SD prior is half-normal(0.2): the first fit wanted ~0.2, double the original
0.1 prior. No zone-level time variation, so zone
epidemics differ in timing only through seeding and movement. If zone fits look poor
(e.g., one zone surging while its province declines), the next step is a zone-level
random walk partially pooled toward the province's.

**T5. Seeding.** Zones with cases in the outbreak's first 3 weeks (10 zones: 6 in
Ituri, 4 in Nord-Kivu) get an estimated number of initial exposed people (prior
lognormal, median 1) plus half that infectious. All other zones start fully susceptible
and are reached only by movement. No external importation term, because with a
deterministic model it would seed small epidemics in every zone.

**T6. Background rate `eps = 0.01` expected cases per zone-week.** Keeps the
likelihood finite where the deterministic model predicts essentially zero but a case
was reported. Adds about 1 expected case per week across all modeled zones.

## Reporting

**R1. [review] Delay distribution: discretized lognormal.** Delay in days = floor of
a continuous lognormal; onsets spread evenly within the week. The hub fitted a
generalized gamma for its imputation; lognormal is simpler. The empirical 99th
percentile (35 days) is consistent with a lognormal, though a long tail (~5–8% of cases
added more than 6 weeks later in the August → October revisions) may be
under-represented.

**R2. Delay log-median = province weekly random walk, week 1 drawn around a national
mean;** one shared log-scale SD. The first version also had a national weekly random
walk shared by all provinces and an iid province-week term. Those were removed (I5):
with ~9,000 delays each province-week is identified directly, and the extra layers
were redundant and made sampling very slow. Consequence: a nationwide slow-reporting
week isn't pooled across provinces explicitly; each province's random walk picks it up
from its own data.

**R3. Completeness F for a week of age A = average over the week's 7 onset days of
P(delay ≤ A + j).** The count likelihood is NegBin(onsets × F); the delay likelihood
divides by F for right truncation.

## Inference and outputs

**I1. NUTS: 4 chains × 500 warmup + 500 samples, adapt_delta 0.9,** started from
Pathfinder draws (falling back to near-prior-median inits). Random inits can start the
SEIR in explosive or extinct regimes.

**I2. Forecasts:** province log R continues its random walk; the SEIR runs forward
from the fitted state. Final counts for observed weeks = reported so far + a negative
binomial draw of the remainder (same dispersion as the main likelihood — an
approximation).

**I3. Hub output.** 100 posterior draws → 100 trajectories. Weeks: reference date − 5
weeks to + 3 weeks (9 weeks, matching current submissions). The reference date is the
Saturday ending the week of the first Sunday on or after as_of, since submissions are
made Sunday–Tuesday (as_of 2026-11-01, a Sunday → 2026-11-07; as_of 2026-10-24, a
Saturday → 2026-10-31). Provinces
and country are sums over zones within each sample, so all levels agree. Team/model
name defaults to `NAU-jointseir` (set `BDBV_TEAM_MODEL`); it's a placeholder.

**I4. Outputs go to `output/bdbv/<as_of>/` (gitignored).** Hub auxiliary files are
snapshotted to `data/raw/bdbv/` (also gitignored) by `fcast.bdbv_inputs --hub
<path>`, defaulting to `../BDBV-Modeling-Hub`.

**I5. Centered random walks (sampling efficiency).** The first version wrote province
log R and delay log-medians as non-centered sums (mean + offset + cumulative sum of
scaled increments). The data identify each week's value tightly, which made the pieces
strongly correlated: every draw hit the maximum tree depth (1,023 gradient
evaluations) with step size ~0.004, so a fit took ~2.5 hours with poor mixing. The
centered parameterization raised the step size ~9× and cut tree depth to ~7 (~130
gradient evaluations) in a short test. Week 1 (pre-outbreak, little data) stays
non-centered for both.

## Results of the first run

Outputs are in `output/bdbv/<as_of>/`. That directory holds the hub CSV, the plots, `diagnostics.txt`, and `fit.rds`.

**Runs and sampling.** Three vintages were fitted, each with 4 chains × 500 warmup + 500 samples. The two jobs shared 8 cores, so wall time would be lower on an idle machine.

| as_of | Reference date | Wall time | Divergences | Treedepth hits | Worst R-hat (min bulk ESS) |
|---|---|---|---|---|---|
| 2026-11-01 | 2026-11-07 | 57 min | 0 | 2 / 2000 | 1.06 (55, on `sd_rw_R` and `lp__`) |
| 2026-10-24 | 2026-10-31 | 53 min | 0 | 0 | 1.07 (58, on `z_zone_R`) |
| 2026-10-17 | 2026-10-24 | 39 min | 0 | 0 | 1.06 (70, on `sd_zone_R` and `z_zone_R`) |

These fits are usable but short. The slow directions are zone R offsets ↔ `sd_zone_R`, and `sd_rw_R`. For production, run 1000 warmup and 1000 samples, or more chains. A non-centered/centered switch for zones with many cases is another option.

**Posterior estimates (11-01).** The three vintages agree closely.
- Movement is far below the gravity default. `phi_move` ≈ 0.0037 (90% CI 0.0025–0.0051) against a prior mean of 0.033. Zones behave nearly independently once seeded, and spread to new zones is attributed to large positive zone offsets.
- Zone heterogeneity is large: `sd_zone_R` ≈ 0.53, so a zone at +1 SD has R × 1.7.
- Province R moves about 25% per week (`sd_rw_R` ≈ 0.25).
- Overdispersion: the NB size is ≈ 2.5 (`inv_sqrt_nb` ≈ 0.64).
- Delay: median ≈ e^2.28 ≈ 10 days at week 1, and the province random walk moves it about 18% per week. `delay_sigma` ≈ 0.88. The completeness of the newest week is ≈ 0.4–0.5 in the main provinces.
- Bas-Uele's delays are barely identified, because it has very few cases.

**Backtest.** Script: `scripts/backtest_bdbv.R`, output `output/bdbv/backtest.csv`. The 10-17 and 10-24 runs were scored against the 11-01 vintage, for country and provinces, on weeks at least 14 days old.
- Absolute error of the median vs the naive "reported as-is" nowcast:
  - Weeks −1 to −3 before the reference: the model is 2–4× better (e.g. week −1 of 10-17: MAE 9.4 vs 24.7).
  - Weeks −4 and −5: about even.
- Coverage pooled over provinces is about right: 90% intervals cover 85–100%, and 50% intervals cover ~85% (too wide for small provinces).
- **[review] The country nowcast is biased low and its intervals are too narrow.** The 90% interval covered the country truth in only 3 of 9 weeks. The model under-predicts even weeks 4–5 weeks old: on 10-17, week 09-19 had 453 reported, the model predicted 460, and the truth is 486. The data show why:
  - Late additions after 4 weeks are larger than the delay file implies. Week 09-19 grew 7% between ages 28 and 42 days, and August weeks grew 15–20% between the 08-18 and 10-09 vintages.
  - In the delay file, only 2.4% of delays exceed 28 days.
  - So cases are being added late without showing up as long delays. Likely sources are imputed-onset cases (D6, ~15%, absent from the delay file) and backfilled records.
  - Country intervals are too narrow because the remaining-case draws are independent across zones (I2), so they average out when summed.

Plots: `output/bdbv/backtest_provinces.png` and `backtest_zones.png` (the six zones with the most recent cases). Two things show up at zone level.
- **Backfill into old weeks.** In CD-IT-26, the week ending 09-19 went from 23 cases in the 10-17 vintage to 41 in 11-01, at an age of 4–6 weeks.
- **Surges the model misses.** CD-NK-03 and CD-IT-26 grew in the weeks after the as-of date while their provinces declined. A zone offset that is fixed over time (T4) can't follow that.

**[review] Why the nowcasts are overconfident: the delay file's report date is not when a case appears in the data.** This undermines D1.

The check: compare the 11-01 delay file with the 10-17 one. Many cases whose recorded report date is on or before 10-17 were not yet in the 10-17 release.

| Province | Onset week ending | Delay ≤ 21 d in the 11-01 file | In the 10-17 file |
|---|---|---|---|
| Nord-Kivu | 09-26 | 118 | 93 |
| Ituri | 09-26 | 214 | 187 |

That is 5–20% of cases across recent weeks, including weeks already 4 weeks old. Cases are entered late and back-dated. Two related facts:
- About 12% of time-series cases have no row in the delay file at all (imputed onsets, D6).

The model therefore learns the delay from report to lab or surveillance, not to publication, and completeness F comes out too high.

| Onset week ending | Model F (90% CI) | Observed fraction against the incomplete 11-01 data |
|---|---|---|
| Ituri, 09-26 | 0.92–0.95 | 0.85 |
| Nord-Kivu, 10-10 | 0.69–0.83 | 0.47 |

F is also very precise: about 9,000 delays pin down a parametric curve, but they measure the wrong delay.

The right fix is to learn completeness from the reporting triangle, i.e. counts at each `as_of`. The delay file would still supply the distribution's shape where vintages are sparse. A fallback is to estimate a "publication lag" on top of the delay-file delay from the vintages.

**Synthetic recovery.** Script: `scripts/synthetic_check_bdbv.R`. Results are in `output/bdbv/2026-11-01/synthetic/recovery.txt`.

The truth was the 11-01 posterior median. From it I simulated counts and individual delays case by case, then refitted, which took 35 minutes with no divergences.
- **Nowcasts:** the 90% interval covers the true final province counts in all 30 province-weeks tested (last 6 weeks). Intervals are narrow, as they should be, once the true delay distribution is lognormal.
- **Parameter coverage:**

  | Parameter | 50% coverage | 90% coverage |
  |---|---|---|
  | province log R | 0.52 | 0.96 |
  | completeness F | 0.73 | 0.97 |
  | delay medians | 0.73 | 0.95 |

  All scalar parameters except one fall inside their 90% intervals.
- **The exception is `sd_zone_R`:** estimated 0.35 (90% CI 0.28–0.44) against a "truth" of 0.52. Much of this comes from how the truth was built. The posterior-median zone offsets are shrunk, so the spread actually simulated was 0.52 × sd(median z) = 0.42, which is inside the interval. Zones with no cases also carry no information about their offset.

So the machinery is self-consistent. The backtest problems come from the data not matching the model's assumptions (the late-addition tail), not from a bug in the fit.

**Next steps, in rough priority order.**
1. **Late-addition tail.** This is the biggest source of error. Options:
   - a mixture delay: lognormal plus a small long-tail component, with weight estimated from the vintage revisions;
   - fit the revisions between vintages directly (D1 alternative) alongside the delay file;
   - add a correlated (country- or province-level) multiplicative factor on the remaining counts, so the uncertainty doesn't average out across zones.
2. Longer chains (see above).
3. A zone-level R random walk (T4), if zone fits matter for scoring.
4. All 519 zones (G1). The cost scales with zones².
5. Informative natural-history priors (T2).

## Second iteration: reporting triangle and zone R random walk

### What changed

**R1′. Completeness comes from the reporting triangle.** [review]
- **Triangle term:** for each zone-week, each earlier vintage's count is modeled given the next vintage's count:
  `early ~ beta_binomial(later, ratio·κ, (1 − ratio)·κ)`, where ratio = F(age at the earlier vintage) / F(age at the later vintage).
- **Delay file:** still sets the shape of the onset-to-report delay D.
  - Its right-truncation term now uses publication completeness.
  - The delay file's report dates come before publication, so they can't set completeness on their own.

**R2′. Late entry.** [review]
- A fraction `pi_late[p]` of each province's cases gets an extra lag L before it is published.
  - L is a discretized lognormal, with parameters `late_m` and `late_s`.
  - `pi_late` is partially pooled across provinces.
- **Why:** the triangle shows a slow trickle of cases into weeks that are already months old.

**R3′. A week's count is final at `final_age = 90` days.** [review]
- The final count is defined as the cases published within 90 days of the week's end.
- Triangle pairs where the earlier vintage is already ≥ 90 days old are dropped.
- **Why:** without a cutoff, the size of the tail is not identified. Earlier fits drifted to `pi ≈ 0.75` with a 230-day median lag.

**R4′. Overdispersion between vintages.**
- The beta-binomial concentration κ allows extra noise between vintages beyond binomial sampling.
- It is estimated at about 5–6 (`inv_sqrt_kappa` ≈ 0.41–0.45), which is a lot of noise.

**T4. Zone R random walk** (`zone_rw`, env `BDBV_ZONE_RW=1`). [review]
- Each zone's log-R offset follows a non-centered weekly random walk on top of the province walk.
- It continues into the forecast weeks.
- **Default:** off.

### Fits

All fits used 4 chains × 500/500 iterations.

| Run | Divergences | Max-treedepth hits | Worst R-hat | Time |
|---|---|---|---|---|
| triangle, 10-17 | 0 | 5 | 1.06 (zone offsets) | ~2.4 h |
| triangle, 10-24 | 0 | 0 | — | ~2.4 h |
| + zone RW, 10-17 | 14 | 1,495 / 2,000 | 1.03 | ~5 h |
| + zone RW, 10-24 | 0 | 2,000 / 2,000 | 1.02 | ~5 h |

**Zone RW sampling:**
- **Tree depth:** these fits sit at the tree-depth limit almost every iteration. The chains mix, but slowly.
- **Variance:** most zone variation moves out of the static offset (`sd_zone_R` 0.53 → 0.13–0.15) and into the walk (`sd_rw_zone` ≈ 0.21–0.23 per week).

**Late-entry estimates:** still large in Nord-Kivu.
- `pi_late`: about 0.70–0.78 in Nord-Kivu, 0.26–0.34 in Ituri.
- Late-lag median: exp(`late_m`) ≈ 60–100 days, with a wide spread.

### Backtest

Truth is the 11-01 vintage, for weeks at least 14 days old. CRPS is averaged over the country and provinces, and summed over zones.

| Level | as_of | v1 (delay file only) | triangle | triangle + zone RW | naive |
|---|---|---|---|---|---|
| country CRPS | 10-17 | **37** | 81 | 78 | — |
| country CRPS | 10-24 | **19** | 128 | 126 | — |
| province CRPS | 10-17 | **10.6** | 16.0 | 15.3 | — |
| province CRPS | 10-24 | **4.8** | 23.7 | 23.7 | — |
| zone CRPS (sum) | 10-17 | **470** | 654 | 553 | — |
| zone CRPS (sum) | 10-24 | **173** | 417 | 438 | — |
| province MAE | 10-24 | **6.5** | 29.1 | 27.9 | 11.1 |

**The triangle model is a regression, and it overpredicts.**
- **Nord-Kivu:** at 10-24, for weeks 2–5 weeks old, the median is 1.5–1.8× the 11-01 count. For example, the week ending 10-03 has a median of 308 against a truth of 207.
- **Ituri:** overpredicted by about 15%.
- **The scoring explains only part of the gap.**
  - The model targets the count at 90 days, while the 11-01 truth for these weeks is only 2–6 weeks old.
  - But the vintages give no support for 50–80% growth after 6 weeks:

    | Onset weeks | Growth from about 6 weeks to about 90 days |
    |---|---|
    | Ituri, June–July | +2–6% |
    | Nord-Kivu, June–July | +32–39% |

  - The Nord-Kivu weeks had only 25–38 cases each, early in the outbreak.
- **Root cause:**
  - The only information about growth after 6 weeks comes from the 08-18 → 10-09 vintage pair, i.e. from June–August weeks.
  - With a time-constant `pi_late`, Nord-Kivu's early reporting backlog is applied to every current week.
  - Current Nord-Kivu weeks show modest growth from age ~4 to ~6 weeks (week ending 09-19: 156 → 167).

**Zone RW vs the triangle model.**
- **Mixed result:**
  - Zone CRPS is 15% better at 10-17 (553 vs 654) and 5% worse at 10-24 (438 vs 417).
  - Country and province scores are about the same.
- **Not a fair test:** the reporting bias above swamps it, so the comparison doesn't show whether the zone RW helps.
- **Cost:** it doubles the run time and saturates tree depth.
- **Decision:** keep it off by default.

### Next steps
1. **Make late entry time-varying or recent-only.** Options:
   - a weekly random walk on logit `pi_late`;
   - fitting the triangle only on pairs from the last few vintages;
   - a cutoff at roughly 6–8 weeks instead of 90 days.

   Until then, v1 is the better model to submit.
2. **Score on the same target.** Score each model's predicted count *as published in the truth vintage*, not its 90-day final. That separates errors in the reporting model from differences in what is being scored.
3. **Revisit the zone RW** once the reporting bias is fixed. It needs a reparameterization or a tighter prior to get off the tree-depth limit.

## Third iteration (v3): releases only, spline delays, zone AR(1) spline

### What changed (2026-11-04)

**R1″. Reporting from the releases only.** [review]
- A case is published T days after onset, T a discretized lognormal with no upper limit.
- Its log-median and log-spread are cubic B-splines by province, with random-walk coefficients across knots (non-centered, with a tight step prior on the spread).
- F(age) = P(T ≤ age) is the share of the eventual count. The hub scores against eventual values.
- Likelihood:
  - the latest release's counts ~ NegBin(onsets × F);
  - the reporting triangle (beta-binomial), using every release pair.
- Removed:
  - the delay file (it measures onset to report, not to publication);
  - the late-entry mixture and `final_age`.
- Previous model: `stan/bdbv_joint_triangle.stan`.

**T4″. Zone transmission.**
- Zone log-R deviation = cubic B-spline. Knots every 4 weeks, counted back from the last data week, so the newest segment has the same length in every fit.
- Coefficients follow a stationary AR(1) toward 0, i.e. toward the province average, with shared `rho_zone` and `sd_ar_zone`.
- A loose sum-to-zero (sd 0.2) applies within each province at each knot. A tight one (0.05) saturated tree depth.
- The spline basis covers the forecast weeks, so their coefficients come from the AR(1).
- This replaces the static zone offsets and the zone random walk.

**T2′. Stage durations.**
- The daily exit probability is now 1/mean, so the mean durations are exactly 7.5 and 8 days.
- With 1 − exp(−1/mean), the means came out about half a day longer.

**Backtest target.** Runs fit with `BDBV_TRUTH` also predict each count *as published in the truth release* (`samples_truth.csv.gz`). `backtest_published.csv/.png` scores those, which is fair for every week.

### Fits

| Setting | Value |
|---|---|
| Chains, iterations | 3 × 300 warmup + 300 sampling |
| Maximum tree depth | 9 (reduced to meet a deadline) |
| Fitting time | ~45 min per as_of, two runs in parallel |
| Divergences | 0 |
| Tree depth | every iteration hit the limit of 9 |
| Worst R-hat | 1.08 (`sd_rw_R` at 10-17) |

| Estimate | Value | What it means |
|---|---|---|
| `rho_zone` | 0.33–0.36 | a zone departure halves in ~2.6 weeks |
| `sd_ar_zone` | ≈ 0.5–0.56 | |
| NegBin overdispersion `inv_sqrt_nb` | 0.53 → 0.20 | the zone spline now absorbs much of the noise in counts |

### Backtest

Truth is the 11-01 release.

**(a) Counts as published in the 11-01 release.** All weeks from 5 weeks before the reference date. v3 only.

| as_of | Country CRPS | Province CRPS | Province 90% coverage | Zone CRPS (sum) |
|---|---|---|---|---|
| 10-17 | 22.7 | 14.5 | 0.66 | 684 |
| 10-24 | 10.6 | 2.9 | 0.90 | 320 |

- **At 10-24:** what will be published is nowcast well and well calibrated.
- **At 10-17** (only three releases): Ituri's latest weeks are overpredicted, and Nord-Kivu's underpredicted from the 10-03 week on.

**(b) Eventual counts against the 11-01 counts.** Weeks at least 14 days old. Same table as earlier.

| Level | as_of | v1 | triangle | v3 |
|---|---|---|---|---|
| country | 10-17 / 10-24 | **37 / 19** | 81 / 128 | 65 / 76 |
| province | 10-17 / 10-24 | **10.6 / 4.8** | 16.0 / 23.7 | 18.0 / 12.6 |
| zone (sum) | 10-17 / 10-24 | **470 / 173** | 654 / 417 | 557 / 270 |

**Problems found.** v3 is not ready to submit; v1 is still the best model to submit.
1. **The lognormal tail is too heavy.**
   - With a spread of ~1 and a median of 10–15 days, about 10–15% of cases are still to come after 5 weeks.
   - The releases show a few percent: e.g. the week ending 09-26 grew 379 → 423 between age 21 and 36 days, and about 2% between 90 and 120 days.
   - The fitted spread is set by the fast early completion, and with one lognormal that same spread forces a long tail.
   - Eventual counts are therefore inflated by roughly 10–40% in recent weeks.
   - **Options:** a lighter-tailed distribution (gamma or Weibull), a tighter prior on the spread, or a check of the fitted tail against the May–July weeks.
2. **Forecast weeks explode.**

   | Run | National forecast 1 week ahead | 4 weeks ahead |
   |---|---|---|
   | v3, 10-24 | median 945 | median 3,111, 95th percentile 37,397 |
   | v1 | median ~300, 95th percentile ~500 | 95th percentile ~1,300 |

   - **Causes:**
     - the inflated eventual counts in recent weeks;
     - low `rho_zone` with large `sd_ar_zone`, so forecast-week zone deviations are drawn almost fresh with sd ≈ 0.6;
     - the summed tail is dominated by zones that happen to draw high R.
   - **Options:** a tighter `sd_ar_zone` prior or a higher `rho_zone` prior. The low rho / high sd combination looks like the spline absorbing count noise, so the two may be linked to the overdispersion change.
3. **Sampling:** every iteration hit tree depth 9. Production runs need a reparameterization of the zone spline, or tree depth 10 and longer runs.

**Files:**
- `output/bdbv_v3/`: `backtest_published.png`, `backtest_provinces.png`, `backtest_zones.png`
- `output/backtest_comparison.png` (all three models; bands cut at 3× the largest observed count)

## Pseudo-releases from the delay files (2026-11-05)

**P1. The 10-09 delay release is recovered from the hub's git history.** The hub only
started archiving delay releases at 10-17, so the 10-09 one (7,113 cases) survives only
in commit `3e6a3c5`. `fcast.bdbv_inputs` now writes any delay release that is in the
history but not archived into `data/raw/bdbv/archived-delay-data/`. It also snapshots
`daily-confirmed-case.csv` and `daily-missing-symptom-by-report-date.csv`.

**P2. Each delay release is a line list of the observed-onset cases.** For every
release and province-week, its rows equal the `imputed = FALSE` counts in
`daily-confirmed-case.csv` exactly (389 of 389 cells). So onset week + delay gives each
case's report date, to within the 7 days of its onset week. Counting the cases
reported by t gives a "pseudo-release" at any t back to May.

**P3. [review] Report date ≠ publication, and the gap varies by release.**
`fcast.bdbv_pseudo` builds the pseudo-releases two ways:
- naive: a case counts once it is reported;
- lag-adjusted: a publication lag G (days from report to publication, a point mass at
  day 0 plus a lognormal) is fit to the real releases by binomial likelihood. A case
  reported on day r counts as published by t with probability G(t − r) / G(a − r),
  where a is the source release.

Check: the source for each real release is the next delay release after it. Values are
reconstructed / real, observed-onset cases, country, with G fit leaving that release out
(`data/processed/bdbv_pseudo_validation.csv`).

| Release | Age 0–6 d | 7–13 | 14–27 | 28–55 | 56+ |
|---|---|---|---|---|---|
| 08-18, naive / lag-adjusted | 1.19 / 0.72 | 1.07 / 0.77 | 1.07 / 0.88 | 1.05 / 0.95 | 1.03 / 0.98 |
| 10-09 | 1.24 / 1.06 | 1.08 / 0.99 | 1.01 / 0.98 | 1.00 / 0.99 | 1.01 / 1.01 |
| 10-17 | 1.87 / 1.60 | 1.18 / 1.09 | 1.16 / 1.13 | 1.07 / 1.05 | 1.02 / 1.01 |
| 10-24 | 1.55 / 1.22 | 1.23 / 1.08 | 1.02 / 0.96 | 1.00 / 0.98 | 1.00 / 0.99 |

- Naive pseudo-releases always run ahead of the real ones, by 20–90% in the newest
  week. This is the same back-dating found in the first run.
- The lag helps, but it isn't stable. G(7 d) fit to single releases ranges from 0.40
  (10-17) to 0.89 (10-24), and the 08-18 release alone can't be fit.
- **By province:** Ituri is reconstructed within 2–7% at every age from 10-09 on.
  Nord-Kivu is the problem. In the 10-17 and 10-24 releases its newest two weeks were
  published 40–46% below the lag-adjusted reconstruction, and at 10-17 weeks ≥ 14 days
  old were 14% below too. That's a backlog specific to that release and province, which
  no single lag can predict.
- Weeks at least 2 weeks old are reconstructed within ~5% (worst: 08-18, and Nord-Kivu
  at 10-17). Weeks under 2 weeks old are off by 20–60%.

**P4. Decision (2026-11-05): the model uses the real releases only.** The pseudo-releases
are not used in any fit; `fcast.bdbv_pseudo` is kept as a diagnostic. The reasoning:
- Before 08-18 they're the only reporting-triangle data at all (June–August onset
  weeks, all ages), and they fill the seven-week gap between the 08-18 and 10-09
  releases with weekly points.
- From 10-09 on they mostly duplicate the real releases, which are weekly anyway.
  Using both would count the same increments twice.
- They're expected values, so they carry none of the release-to-release noise that κ
  estimates (P3).

So the clean option would be province-level pseudo pairs only before the first real
release. That informs the delay spline for June–August, not the recent weeks the
nowcast depends on. The more direct use is the check v3 problem 1 asks for: the
completeness tail for May–July weeks. Here is that curve from the 11-01 source, June–July
onset weeks, as a share of the 11-01 count. Naive is an upper bound; lag-adjusted
depends on the unstable G.

| Age (days) | 0 | 7 | 14 | 21 | 28 | 42 | 56 | 84 |
|---|---|---|---|---|---|---|---|---|
| Ituri, naive / lag-adj. | 0.31 / 0.19 | 0.82 / 0.59 | 0.94 / 0.76 | 0.97 / 0.83 | 0.99 / 0.87 | 0.99 / 0.92 | 1.00 / 0.95 | 1.00 / 0.98 |
| Nord-Kivu | 0.42 / 0.27 | 0.75 / 0.57 | 0.87 / 0.71 | 0.93 / 0.79 | 0.95 / 0.85 | 0.98 / 0.91 | 0.99 / 0.94 | 1.00 / 0.98 |

Outputs: `data/processed/bdbv_pseudo_vintages.csv` (one row per source release,
pseudo-release Saturday, province and onset week; `observation` lag-adjusted, `naive`,
`source_count`), `bdbv_pseudo_validation.csv`, `bdbv_pseudo_lag.json`.

## Fourth iteration (v4): reporting triangle on province totals

### What changed (2026-11-05)

**R5. The triangle uses province totals, not zone counts.** [review]
- Each province-week's count (summed over modeled zones) in one release is modeled given the next release's, as before. Zone counts still enter through the latest release.
- **Why:** fit to zone counts, the reporting part alone reproduces v3's heavy Nord-Kivu tail (11% unpublished at 91 days, 19% at 119). Fit to province totals it leaves 5% and 12%, and it predicts the 10-24 → 11-01 province growth much better (log predictive density −113 vs −134). Two reasons:
  - Cases move between zones of a province from one release to the next. Zone counts are much noisier between releases (concentration κ ≈ 7, against ≈ 36 for province totals), so they barely constrain the tail.
  - Capping zone counts that drop turns those moves into apparent late additions: weeks 5–9 weeks old grow 3% at zone level and 2% at province level in Nord-Kivu (3% and 0% in Tshopo).
- **Also tried, not adopted:** a Weibull delay (same median/spread parameterization, light right tail). On the triangle alone it predicts the newest week's growth worse (×3.0 against an observed ×4.2; lognormal ×3.7) and scores worse on the 11-01 release (−116 vs −113). Over the 2–4-week spans that matter for final counts, it predicts almost the same growth as the lognormal.

### Fits

Same settings as v3: 3 chains × 300 + 300 iterations, tree depth 9, 0 divergences. Tree-depth hits: 600/900 at 10-17, 900/900 at 10-24. `inv_sqrt_kappa` ≈ 0.12–0.17 (κ ≈ 35–70). Other estimates are close to v3's.

### Backtest

Final counts scored against the 11-01 counts (weeks at least 14 days old), CRPS:

| Level | as_of | v1 | v3 | v4 |
|---|---|---|---|---|
| country | 10-17 / 10-24 | **37 / 19** | 65 / 76 | 88 / 94 |
| province | 10-17 / 10-24 | **10.6 / 4.8** | 18.0 / 12.6 | 24.3 / 16.6 |
| zone (sum) | 10-17 / 10-24 | **470 / 173** | 557 / 270 | 697 / 328 |

Counts as published in the 11-01 release: country CRPS 30.5 / 12.3, province 17.4 / 4.0 (v3: 22.7 / 10.6 and 14.5 / 2.9).

**v4 is worse than v3, although its far tail is fixed.** The heavy tail beyond about two months was not what inflated final counts.
- **Final counts are observed / F.** Decomposing the 10-24 fit, the final-count median is within ~2% of reported / F for every recent week. SEIR onsets are not the driver.
- **F is too low at 3–6 weeks.** At 10-24, the week ending 09-26 (28 days old) has country F = 0.89 in v4 and 0.92 in v3: 11% and 8% still to come. The releases point to about 3–5%: that week went from 419 to 423 by 11-01, and weeks 5–9 weeks old grow 1–2% over the following 3 weeks.
- **Every delay shape tried overpredicts mid-age growth.** Joint-model lognormal, reporting-only lognormal and reporting-only Weibull all predict 4–8% growth from 10-09 to 11-01 for weeks 3–4 weeks old, against 1–2% observed. The very steep first two weeks (×4, then ×1.5, then ×1.13) force a spread that leaves too much for weeks 3–6.
- **Release-to-release growth is lumpy.** The 10-17 → 10-24 release added 8% to weeks 3–5 weeks old, while other releases added 1–2%. A delay that is smooth in onset week can't represent a release that catches up a backlog. It averages over releases instead.

### Next steps
1. **A non-parametric delay.** Use a discrete hazard per week of age (0, 1, 2, 3, 4, 5+), varying by province and smoothly in time, so the early and late parts of the curve aren't tied together by one spread parameter.
2. **Release effects.** Add a random effect per release on the hazard, so a backlog catch-up is absorbed by that release instead of the curve.
3. Keep R5: it is better founded and fixes the far tail, but it doesn't fix final counts on its own.

## Fifth iteration (v5): weekly-hazard delay

### What changed (2026-11-05)

**R6. Delay = a constant daily publication rate within each week since onset.** [review]
- Weeks since onset: 0, 1, 2, 3, 4, and 5+.
- Log rate = national rate for that week + a province × week deviation + a province trend over onset weeks. The trend is a cubic B-spline with random-walk coefficients starting at 0, and it scales all of a province's rates together.
- **Why:** with a lognormal, one spread parameter ties the fast first two weeks to the tail. Here every week since onset has its own rate.
- **The eventual count is defined as cases published within `d_max = 120` days of onset**, so F(age) = P(T ≤ age | T ≤ 120).
  - Without a cap the tail isn't identified: a small late rate on a large leftover pool fits the releases as well as a large rate on a small one. The uncapped reporting-only fit left 10–35% of cases unpublished at 28 days.
  - Triangle pairs whose earlier release is already ≥ 120 days old are dropped, as in v2 (R3′).
- **Files and settings:**
  - Model: `stan/bdbv_joint_hazard.stan`.
  - Select it with `BDBV_DELAY=hazard`, or `delay = "hazard"` in `bdbv_config()`.
  - `stan/bdbv_joint.stan` is still the lognormal (v4).
- **Reporting-only test** (triangle alone, scored on the 11-01 release, weeks under 91 days old):

  | Delay | Log score, 10-24 | Log score, 10-17 |
  |---|---|---|
  | Lognormal | −82.6 | −112.2 |
  | Hazard, no cap | −82.7 | −107.9 |
  | Hazard, 120-day cap | −82.4 | −106.7 |

  All three predict the same 10-09 → 11-01 growth for weeks 3–6 weeks old: 1.06, against 1.03 observed.

### Fits

Same settings as v3 and v4. 0 divergences, and every iteration hit tree depth 9.
- Daily rates (national, posterior median) by week since onset: 0.033, 0.08, 0.03, 0.02, 0.011, then 0.007 from week 5. The first week is slower than the second, which looks like a publication lag of a few days.
- `sd_prov_haz` ≈ 0.35–0.6 and `sd_step_haz` ≈ 0.45–0.58, so provinces differ substantially and drift over time.

### Backtest

Final counts against the 11-01 counts (weeks at least 14 days old), CRPS:

| Level | as_of | v1 | v3 | v4 | v5 |
|---|---|---|---|---|---|
| country | 10-17 / 10-24 | **37 / 19** | 65 / 76 | 88 / 94 | 53 / 87 |
| province | 10-17 / 10-24 | **10.6 / 4.8** | 18.0 / 12.6 | 24.3 / 16.6 | 16.9 / 15.2 |
| zone (sum) | 10-17 / 10-24 | **470 / 173** | 557 / 270 | 697 / 328 | 557 / 302 |

Counts as they appear in the 11-01 release (the next-release nowcast): v5 is the best of v3–v5.

| | v3 | v4 | v5 |
|---|---|---|---|
| Country CRPS, 10-17 / 10-24 | 22.7 / 10.6 | 30.5 / 12.3 | **18.2 / 11.4** |
| Province CRPS, 10-17 / 10-24 | 14.5 / 2.9 | 17.4 / 4.0 | 13.9 / 4.0 |

Country 90% coverage for v5 is 0.71 at 10-17 and 0.83 at 10-24.

Per zone (`scripts/compare_zones.R`, all 117 zones, weeks from 8 weeks before as_of; figures in `output/zone_comparison_<as_of>.pdf` and `output/zone_crps.png`):

| as_of | v5 better (by >10%) | v1 better (by >10%) |
|---|---|---|
| 10-17 | 21 zones | 80 zones |
| 10-24 | 15 zones | 96 zones |

**Conclusion: don't submit v5; v1 stays the submission model.**
- The hazard predicts the *next* releases best, but its final counts are still high.
- At 10-24 the country median for the week ending 09-26 is 507. 420 had been reported, and 11-01 shows 423.
- Final counts sit a few percent above the reported counts even for weeks 6–8 weeks old. The model expects 5–15% of a week to arrive between weeks 4 and 17, and the releases so far show only 1–3% of that.
- **Still unresolved:** what eventual counts will be. Every release-based model (v3–v5) says weeks 4–8 weeks old will grow by 5–15%. v1 says about 1–3%. The 11-01 release can't settle it.
  - If the release-based models are right, the backtest's truth is too low and v1's advantage is partly a scoring artifact.
  - If v1 is right, the releases' early catch-up batches (e.g. 10-17 → 10-24) are being extrapolated into a tail that doesn't exist.
  - **The 11-08 and 11-15 releases will decide it.** Score the saved 10-17 and 10-24 runs against them before changing models.

## Forecast R damping, from finalized-data backtests (2026-11-06)

Backtests on finalized data (`scripts/backtest_final.R`, `scripts/score_final.R`): 13 weekly origins 2026-05-09 to 2026-08-01, perfect reporting, truth from the 11-01 release (weeks at least 60 days old), WIS on 1–4 week horizons.

**Decision: forecasts decay each zone's log R by φ = 0.5 per forecast week** (province R and the zone deviation alike). This is now the default (`forecast_R = "damped"`, `damp_R = 0.5` in `bdbv_config()`), so it applies to `scripts/run_bdbv.R`.

WIS relative to the flatline baseline (all horizons):

| Model | Country | Province | Zone |
|---|---|---|---|
| **E0_damp050** (φ = 0.5) | **0.71** | **0.78** | **0.93** |
| E0_flat (province R held, zone spline continued) | 7.37 | 5.99 | 3.37 |
| E3_ar_mean (province AR(1), φ_R fitted) | 1.51 | 1.43 | 1.24 |
| F1_mean (E3 plus weekly AR(1) zone deviations) | 1.38 | 1.33 | 1.18 |

- **Why the undamped model overpredicts:** R at the last fitted weeks is poorly identified (it reflects infections about a generation interval earlier). The forecast carries that estimate forward and compounds it.
- **Why 0.5:** it is roughly true persistence (about 0.6–0.75 per week) times the reliability of the last-week R estimate (about 0.65). Regressing later R on the real-time estimate gives a slope of 0.46, against 0.62 from well-identified retrospective estimates. Forward-chained selection (choosing φ from earlier origins only) also picks 0.5.
- **The optimum is flat:** country relative WIS is 0.77 at φ = 0.3, 0.71 at 0.5 and 0.77 at 0.6. Damping helps most in the growth phase. In the July–August plateau it is about even with the flatline or slightly worse.
- **Estimating φ inside the fit doesn't work.** The fitted persistence (φ_R ≈ 0.73–0.79; weekly zone ρ ≈ 0.8–0.86 from late June) describes how persistent true R is. It doesn't discount the noisy last-week estimate. The AR(1) options stay in the model (`ar_R`, `zone_weekly`) but are off by default.
- **Not tried:** an ARMA(1,1) or transient weekly shock term on R, which would let the fit separate lasting changes from one-off bumps.

**Revisit when:**
- **The reporting model is in place.** Rerun the φ grid on real-time data snapshots. Incomplete recent weeks should call for more damping, so expect a value below 0.5.
- **The 2026-08-08 and 2026-08-15 origins can be scored** (the 11-08 and 11-15 releases, 60-day rule). The backtests should cover weekly forecast dates from mid-May to mid-August. We kept the 60-day truth rule rather than relaxing it, so when a new release lands, extend the origins in `scripts/backtest_final.R` and `truth_end` in `scripts/plot_forecasts_final.R` (`scripts/score_final.R` applies the 60-day rule to the latest release itself) to the latest forecast date whose 4-week horizon is at least 60 days old.
- **The damped model loses to the flatline for several weeks running** in live scoring.
