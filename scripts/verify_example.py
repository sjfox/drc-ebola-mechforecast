"""
Verify the epymorph metapopulation model against independent calculations.

Runs one example (SEIR with Centroids movement on the Ituri and Nord-Kivu health
zones) many times and checks it six ways, each computed without epymorph's own code:

1. Transitions: each tick's events match Poisson draws at the rates in the equations.
2. Movement: daily visitor counts match the distance kernel.
3. Reference model: the average of the runs matches our own implementation of the
   same model, in discrete steps and in continuous time.
4. Growth rate: early exponential growth matches SEIR theory.
5. Final size: the share ever infected matches the final-size equation.
6. Conservation: no one is created or lost.

Parameter values are placeholders, not fitted to the outbreak.

Run with:  uv run python -u scripts/verify_example.py
"""

import json
from pathlib import Path

import numpy as np
from epymorph.kit import BasicSimulator, SingleStrataRUME, TimeFrame, default_rng, init, ipm, mm

from fcast.geo import load_geo

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports" / "verify_example.json"

PROVINCES = ["Ituri", "Nord-Kivu"]
SEED_ZONES = ["Bunia", "Mongbwalu", "Rwampara"]
SEED_PER_ZONE = 100  # infectious people placed in each seed zone on day 0
START = "2026-04-01"
DAYS = 600
RUNS = 20

BETA = 0.20
SIGMA = 1 / 9
GAMMA = 1 / 10
PHI = 40.0  # miles: epymorph's haversine distances are in miles
COMMUTER_PROPORTION = 0.05
R0 = BETA / GAMMA

EARTH_RADIUS_MILES = 3958.8
S, E, I, R = range(4)


# ---------------------------------------------------------------------------
# Independent pieces (no epymorph code)
# ---------------------------------------------------------------------------


def haversine_miles(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    lon, lat = np.radians(lon), np.radians(lat)
    dlon = lon[None, :] - lon[:, None]
    dlat = lat[None, :] - lat[:, None]
    a = np.sin(dlat / 2) ** 2 + np.cos(lat[:, None]) * np.cos(lat[None, :]) * np.sin(dlon / 2) ** 2
    return 2 * EARTH_RADIUS_MILES * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def travel_matrix(pop: np.ndarray, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    """
    Expected travelers from home i to destination j on each day's first step.

    floor(pop_i * commuter_proportion) people pick destinations with probability
    proportional to exp(-distance / phi), including their own zone; anyone who
    "picks" home stays home.
    """
    kernel = np.exp(-haversine_miles(lon, lat) / PHI)
    kernel /= kernel.sum(axis=1, keepdims=True)
    n = np.floor(pop * COMMUTER_PROPORTION)
    trips = n[:, None] * kernel
    np.fill_diagonal(trips, 0.0)
    return trips


def growth_rate_theory(r0: float, sigma: float, gamma: float) -> float:
    """Exponential growth rate of SEIR: solves (1 + r/sigma)(1 + r/gamma) = R0."""
    a = 1 / (sigma * gamma)
    b = 1 / sigma + 1 / gamma
    return float((-b + np.sqrt(b * b + 4 * a * (r0 - 1))) / (2 * a))


def final_size_theory(r0: float) -> float:
    """Final share ever infected: solves z = 1 - exp(-R0 z)."""
    z = 0.9
    for _ in range(500):
        z = 1 - np.exp(-r0 * z)
    return float(z)


def reference_model(
    initial: np.ndarray, pop: np.ndarray, trips: np.ndarray, steps: tuple[float, float], substeps: int | None
) -> np.ndarray:
    """
    Deterministic version of the same metapopulation SEIR.

    Each day: a share trips[i, j] / pop_i of every compartment of zone i visits zone j
    for the first step, while infection happens wherever people are (force of
    infection beta * I / N among everyone present); then everyone goes home for the
    second step.

    substeps=None advances each step in one go (the expected value of epymorph's
    discrete steps). substeps=k integrates each step in continuous time with k RK4
    substeps per day-fraction unit (true ODE solution).

    Returns daily new infectious cases (E -> I) by home zone, shape (DAYS, N).
    """
    n_zones = len(pop)
    share = trips / pop[:, None]
    stay = 1 - share.sum(axis=1)
    mix = share.copy()
    mix[np.arange(n_zones), np.arange(n_zones)] = stay  # mix[i, j]: share of i at j

    def rates(x):  # x: (home, location, compartment)
        present = x.sum(axis=0)  # (location, compartment)
        n_present = np.maximum(present.sum(axis=1), 1.0)
        foi = BETA * present[:, I] / n_present  # per location
        se = foi[None, :] * x[:, :, S]
        ei = SIGMA * x[:, :, E]
        ir = GAMMA * x[:, :, I]
        d = np.empty_like(x)
        d[:, :, S] = -se
        d[:, :, E] = se - ei
        d[:, :, I] = ei - ir
        d[:, :, R] = ir
        return d, ei

    def advance(x, tau):
        if substeps is None:
            d, ei = rates(x)
            return x + d * tau, ei.sum(axis=1) * tau
        k = max(1, round(substeps * tau))
        h = tau / k
        cases = np.zeros(n_zones)
        for _ in range(k):
            # Track cumulative E -> I flow as an extra state integrated with RK4.
            d1, f1 = rates(x)
            d2, f2 = rates(x + h / 2 * d1)
            d3, f3 = rates(x + h / 2 * d2)
            d4, f4 = rates(x + h * d3)
            x = x + h / 6 * (d1 + 2 * d2 + 2 * d3 + d4)
            cases += h / 6 * (f1 + 2 * f2 + 2 * f3 + f4).sum(axis=1)
        return x, cases

    home = initial.astype(np.float64)
    daily = np.zeros((DAYS, n_zones))
    for day in range(DAYS):
        x = mix[:, :, None] * home[:, None, :]  # everyone splits per the travel matrix
        x, c0 = advance(x, steps[0])
        home = x.sum(axis=1)  # everyone returns home
        x = np.zeros((n_zones, n_zones, 4))
        x[np.arange(n_zones), np.arange(n_zones)] = home
        x, c1 = advance(x, steps[1])
        home = x.sum(axis=1)
        daily[day] = c0 + c1
    return daily


def fit_growth(series: np.ndarray, lo_day: int, hi_day: int) -> float:
    days = np.arange(lo_day, hi_day + 1)
    return float(np.polyfit(days, np.log(series[lo_day : hi_day + 1]), 1)[0])


# ---------------------------------------------------------------------------
# Checks on epymorph output
# ---------------------------------------------------------------------------


def check_transitions(out, taus: np.ndarray) -> dict:
    """Compare each tick's events with Poisson(rate x tau) at the visited location."""
    post = out.visit_compartments.astype(np.float64)  # after this tick's events
    ev = out.visit_events.astype(np.float64)  # (ticks, N, 3): S->E, E->I, I->R
    pre = post.copy()  # undo this tick's events to get the state the rates saw
    pre[..., S] += ev[..., 0]
    pre[..., E] += -ev[..., 0] + ev[..., 1]
    pre[..., I] += -ev[..., 1] + ev[..., 2]
    pre[..., R] -= ev[..., 2]
    n = np.maximum(pre.sum(axis=2), 1.0)
    tau = taus[:, None]
    expected = np.stack(
        [BETA * pre[..., S] * pre[..., I] / n * tau, SIGMA * pre[..., E] * tau, GAMMA * pre[..., I] * tau],
        axis=2,
    )
    n_steps = len(taus) // DAYS
    names = ["S → E (infection)", "E → I (becoming infectious)", "I → R (recovery)"]
    bins = np.arange(-4, 4.01, 0.25)
    result = []
    for k, name in enumerate(names):
        o, x = ev[..., k], expected[..., k]
        mask = x >= 10  # normal approximation to Poisson is good here
        z = (o[mask] - x[mask]) / np.sqrt(x[mask])
        hist, _ = np.histogram(np.clip(z, -3.999, 3.999), bins=bins)
        result.append(
            {
                "name": name,
                "observed_total": float(o.sum()),
                "expected_total": float(x.sum()),
                "ratio": float(o.sum() / x.sum()),
                "n_cells": int(mask.sum()),
                "z_mean": float(z.mean()),
                "z_sd": float(z.std()),
                "dispersion": float((((o - x) ** 2)[x > 0] / x[x > 0]).mean()),
                "daily_observed": o.reshape(DAYS, n_steps, -1).sum(axis=(1, 2)).tolist(),
                "daily_expected": x.reshape(DAYS, n_steps, -1).sum(axis=(1, 2)).round(2).tolist(),
                "z_hist": hist.tolist(),
            }
        )
    return {"transitions": result, "z_bins": bins.tolist()}


def check_movement(out, pop: np.ndarray, trips: np.ndarray, n_steps: int) -> dict:
    """Visitor counts on each day's first step versus the travel matrix."""
    present = out.visit_compartments.sum(axis=2).astype(np.float64)  # (ticks, N)
    away = present[0::n_steps]  # first step of each day
    home_again = present[1::n_steps]
    out_trips, in_trips = trips.sum(axis=1), trips.sum(axis=0)
    expected = pop - out_trips + in_trips

    # Daily variance of the count present at each zone. epymorph samples destinations in
    # two stages: c ~ Multinomial(n_i, kernel_i) gives requested trips (the stay-home
    # share is dropped), then the m_i = n_i - c_ii actual travelers are re-assigned with
    # Multinomial(m_i, c_i / m_i). Leaving is unaffected by stage two; arriving gains
    # E[c_ij (1 - c_ij / m_i)] of extra variance.
    n = np.floor(pop * COMMUTER_PROPORTION)
    k = trips / n[:, None]  # kernel with the diagonal removed
    k_stay = 1 - k.sum(axis=1)
    var_leave = n * k_stay * (1 - k_stay)
    var_c = n[:, None] * k * (1 - k)
    m = n * (1 - k_stay)
    second_stage = n[:, None] * k - (var_c + (n[:, None] * k) ** 2) / m[:, None]
    var_one_stage = var_leave + var_c.sum(axis=0)
    var_two_stage = var_leave + (var_c + second_stage).sum(axis=0)

    mean_obs = away.mean(axis=0)
    days = away.shape[0]
    z1 = (mean_obs - expected) / np.sqrt(var_one_stage / days)
    z = (mean_obs - expected) / np.sqrt(var_two_stage / days)
    # Day-to-day spread, observed vs predicted, as a direct check on the variance.
    sd_ratio = away.std(axis=0, ddof=1) / np.sqrt(var_two_stage)
    return {
        "expected_net": (expected - pop).round(2).tolist(),
        "observed_net_mean": (mean_obs - pop).round(2).tolist(),
        "z": z.round(3).tolist(),
        "z_mean": float(z.mean()),
        "z_sd": float(z.std()),
        "max_abs_z": float(np.abs(z).max()),
        "z_sd_one_stage": float(z1.std()),
        "daily_sd_ratio": sd_ratio.round(3).tolist(),
        "daily_sd_ratio_median": float(np.median(sd_ratio)),
        "expected_travelers_per_day": float(out_trips.sum()),
        "all_home_on_second_step": bool(np.array_equal(home_again, np.broadcast_to(pop, home_again.shape))),
    }


# ---------------------------------------------------------------------------


def main() -> None:
    geo = load_geo(provinces=PROVINCES)
    t = geo.table
    pop = geo.params["population"].astype(np.float64)
    n_zones = len(pop)
    seed_ids = geo.ids_for(SEED_ZONES, province="Ituri")
    seed_idx = [int(np.flatnonzero(t["node_id"] == s)[0]) for s in seed_ids]

    initial = np.zeros((n_zones, 4), dtype=np.int64)
    initial[:, S] = geo.params["population"]
    initial[seed_idx, S] -= SEED_PER_ZONE
    initial[seed_idx, I] = SEED_PER_ZONE

    rume = SingleStrataRUME.build(
        ipm=ipm.SEIRS(),
        mm=mm.Centroids(),
        init=init.Explicit(initial),
        scope=geo.scope,
        time_frame=TimeFrame.of(START, DAYS),
        params={
            **geo.params,
            "beta": BETA,
            "sigma": SIGMA,
            "gamma": GAMMA,
            "xi": 0.0,
            "phi": PHI,
            "commuter_proportion": COMMUTER_PROPORTION,
        },
    )
    steps = tuple(rume.tau_step_lengths)
    n_steps = len(steps)
    taus = np.tile(np.array(steps), DAYS)
    sim = BasicSimulator(rume)

    cases = np.zeros((RUNS, DAYS, n_zones))
    final = np.zeros((RUNS, n_zones, 4))
    conserved_home = conserved_total = nonnegative = True
    for r in range(RUNS):
        out = sim.run(rng_factory=default_rng(1000 + r))
        home = out.data_by_home
        cases[r] = home.events_per_day[:, :, 1]
        final[r] = out.home_compartments[-1]
        conserved_home &= bool(np.all(out.home_compartments.sum(axis=2) == geo.params["population"]))
        conserved_total &= bool(np.all(out.visit_compartments.sum(axis=(1, 2)) == pop.sum()))
        nonnegative &= bool(out.home_compartments.min() >= 0 and out.visit_compartments.min() >= 0)
        if r == 0:
            featured = out
        print(f"run {r + 1}/{RUNS} done", flush=True)

    trips = travel_matrix(pop, t["centroid_lon"].to_numpy(), t["centroid_lat"].to_numpy())

    print("checking transitions and movement ...", flush=True)
    transitions = check_transitions(featured, taus)
    movement = check_movement(featured, pop, trips, n_steps)

    print("running reference models ...", flush=True)
    ref_discrete = reference_model(initial, pop, trips, steps, substeps=None)
    ref_continuous = reference_model(initial, pop, trips, steps, substeps=48)

    total = cases.sum(axis=2)  # (RUNS, DAYS)
    mean_total = total.mean(axis=0)
    ref_d_total, ref_c_total = ref_discrete.sum(axis=1), ref_continuous.sum(axis=1)

    # Growth: fit from day 30 (past the start-up transient) until cases reach 2% of peak.
    peak = mean_total.max()
    hi_day = int(np.argmax(mean_total >= 0.02 * peak))
    lo_day = 30
    r_theory = growth_rate_theory(R0, SIGMA, GAMMA)
    growth = {
        "theory": r_theory,
        "runs_mean": fit_growth(mean_total, lo_day, hi_day),
        "per_run": [fit_growth(total[k], lo_day, hi_day) for k in range(RUNS)],
        "reference_discrete": fit_growth(ref_d_total, lo_day, hi_day),
        "reference_continuous": fit_growth(ref_c_total, lo_day, hi_day),
        "window": [lo_day, hi_day],
        # Earlier windows have less immunity in the seeded zones, so the continuous ODE
        # should approach theory; a gap that remains for the stepped versions is the
        # time-step effect.
        "window_sensitivity": [
            {
                "window": [a, b],
                "runs_mean": fit_growth(mean_total, a, b),
                "reference_discrete": fit_growth(ref_d_total, a, b),
                "reference_continuous": fit_growth(ref_c_total, a, b),
            }
            for a, b in [(20, 50), (30, 60), (30, 80), (lo_day, hi_day)]
        ],
    }

    attack = (pop - final[:, :, S]) / pop  # (RUNS, N)
    z_theory = final_size_theory(R0)
    still_active = (final[:, :, E] + final[:, :, I]).sum(axis=1)

    def peak_stats(series):
        return {"day": int(np.argmax(series)), "size": float(series.max()), "total": float(series.sum())}

    result = {
        "setup": {
            "provinces": PROVINCES,
            "n_zones": n_zones,
            "population": int(pop.sum()),
            "seed_zones": SEED_ZONES,
            "seed_ids": seed_ids,
            "seed_per_zone": SEED_PER_ZONE,
            "start": START,
            "days": DAYS,
            "runs": RUNS,
            "steps": steps,
            "beta": BETA,
            "sigma": SIGMA,
            "gamma": GAMMA,
            "r0": R0,
            "phi_miles": PHI,
            "commuter_proportion": COMMUTER_PROPORTION,
        },
        "zones": [
            {"node_id": z.node_id, "name": z.name, "province": z.province, "population": int(z.population)}
            for z in t.itertuples()
        ],
        "conservation": {
            "home_totals_match": conserved_home,
            "overall_total_matches": conserved_total,
            "no_negative_counts": nonnegative,
            "ticks_checked": int(RUNS * DAYS * n_steps),
        },
        "transitions": transitions,
        "movement": movement,
        "reference": {
            "runs_mean": mean_total.round(1).tolist(),
            "runs_lo": np.percentile(total, 5, axis=0).round(1).tolist(),
            "runs_hi": np.percentile(total, 95, axis=0).round(1).tolist(),
            "featured_run": total[0].tolist(),
            "discrete": ref_d_total.round(1).tolist(),
            "continuous": ref_c_total.round(1).tolist(),
            "summary": {
                "runs_mean": peak_stats(mean_total),
                "discrete": peak_stats(ref_d_total),
                "continuous": peak_stats(ref_c_total),
            },
            "zone_peak_day": {
                "runs_median": np.median(cases.argmax(axis=1), axis=0).tolist(),
                "discrete": ref_discrete.argmax(axis=0).tolist(),
                "continuous": ref_continuous.argmax(axis=0).tolist(),
            },
        },
        "growth": growth,
        "final_size": {
            "theory": z_theory,
            "attack_rate_median": np.median(attack, axis=0).round(4).tolist(),
            "attack_rate_lo": np.percentile(attack, 5, axis=0).round(4).tolist(),
            "attack_rate_hi": np.percentile(attack, 95, axis=0).round(4).tolist(),
            "overall_median": float(np.median((pop.sum() - final[:, :, S].sum(axis=1)) / pop.sum())),
            "still_infected_at_end_median": float(np.median(still_active)),
        },
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, separators=(",", ":")))
    print(f"Wrote {OUT} ({OUT.stat().st_size / 1e6:.2f} MB)")

    print("\nSummary")
    print(f"  conservation: home {conserved_home}, total {conserved_total}, nonneg {nonnegative}")
    for tr in transitions["transitions"]:
        print(
            f"  {tr['name']:<30} obs/exp {tr['ratio']:.4f}  z mean {tr['z_mean']:+.3f} sd {tr['z_sd']:.3f}"
            f"  dispersion {tr['dispersion']:.3f}"
        )
    print(
        f"  movement: z mean {movement['z_mean']:+.3f} sd {movement['z_sd']:.3f} max|z| {movement['max_abs_z']:.2f}"
        f"  all home on step 2: {movement['all_home_on_second_step']}"
    )
    s = result["reference"]["summary"]
    for k in ("runs_mean", "discrete", "continuous"):
        print(f"  {k:<11} peak day {s[k]['day']}, peak {s[k]['size']:,.0f}, total {s[k]['total']:,.0f}")
    print(
        f"  growth: theory {r_theory:.5f}, runs {growth['runs_mean']:.5f}, discrete {growth['reference_discrete']:.5f},"
        f" continuous {growth['reference_continuous']:.5f}, window {growth['window']}"
    )
    fs = result["final_size"]
    ar = np.array(fs["attack_rate_median"])
    print(
        f"  final size: theory {z_theory:.4f}, overall {fs['overall_median']:.4f}, zones min {ar.min():.4f} max {ar.max():.4f},"
        f" still infected at end {fs['still_infected_at_end_median']:.0f}"
    )


if __name__ == "__main__":
    main()
