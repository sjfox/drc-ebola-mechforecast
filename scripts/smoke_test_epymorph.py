"""
Smoke test: run a short SEIR metapopulation simulation on the processed geography,
once with epymorph's `Centroids` movement and once with `fcast.movement.Gravity`.

This only checks that the data plug into epymorph correctly. The parameter values are
placeholders, not calibrated to the outbreak.

Run with:  uv run python scripts/smoke_test_epymorph.py
"""

import numpy as np
from epymorph.kit import BasicSimulator, SingleStrataRUME, TimeFrame, default_rng, init, ipm, mm

from fcast.geo import load_geo
from fcast.movement import Gravity

geo = load_geo(provinces=["Ituri", "Nord-Kivu"])
seeds = geo.ids_for(["Bunia", "Mongbwalu", "Rwampara"], province="Ituri")


def run(movement):
    rume = SingleStrataRUME.build(
        ipm=ipm.SEIRS(),
        mm=movement,
        init=init.LabeledLocations(labels=seeds, seed_size=30),
        scope=geo.scope,
        time_frame=TimeFrame.of("2026-05-01", 90),
        params={
            **geo.params,
            "beta": 0.25,
            "sigma": 1 / 9,
            "gamma": 1 / 10,
            "xi": 0.0,  # no waning immunity: SEIR
            "phi": 40.0,
            "commuter_proportion": 0.05,
        },
    )

    out = BasicSimulator(rume).run(rng_factory=default_rng(42))

    # home_compartments: (ticks, nodes, compartments), counted by home zone rather than
    # by where people are during the tick.
    final = out.home_compartments[-1]
    infected = final[:, 1] + final[:, 2]
    ever = geo.params["population"] - final[:, 0]
    t = geo.table.assign(final_EI=infected, ever_infected=ever)
    print(f"{len(t)} health zones, {t['population'].sum():,} people")
    print(f"Zones ever infected by day 90: {(t['ever_infected'] > 0).sum()}")
    print(
        t.sort_values("ever_infected", ascending=False)[
            ["name", "province", "population", "ever_infected", "final_EI"]
        ]
        .head(10)
        .to_string(index=False)
    )
    assert np.array_equal(out.initial.sum(axis=1), geo.params["population"])
    assert np.array_equal(final.sum(axis=1), geo.params["population"])
    print("OK: initial and final populations match the GRID3 node populations.")


for name, movement in [("Centroids", mm.Centroids()), ("Gravity", Gravity())]:
    print(f"--- {name} movement ---")
    run(movement)
