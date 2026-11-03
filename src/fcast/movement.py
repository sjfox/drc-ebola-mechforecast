"""
Movement models for the DRC health zones.

`Gravity` is a drop-in alternative to epymorph's `Centroids`: the same daily
commute-and-return pattern, but destinations are chosen by a gravity kernel over road
travel time instead of an exponential kernel over straight-line distance.

    from fcast.movement import Gravity
    rume = SingleStrataRUME.build(..., mm=Gravity(), params={**geo.params, ...})

`load_geo` supplies the `travel_time` matrix (hours) alongside `population`.
"""

from functools import cached_property

import numpy as np
from numpy.typing import NDArray
from typing_extensions import override

from epymorph.attribute import AttributeDef
from epymorph.data_shape import Shapes
from epymorph.data_type import SimDType
from epymorph.error import DataAttributeError
from epymorph.movement_model import EveryDay, MovementClause, MovementModel
from epymorph.simulation import Tick, TickDelta, TickIndex
from epymorph.util import row_normalize


def gravity_kernel(
    population: NDArray,
    travel_time: NDArray,
    pop_exp: float = 1.0,
    time_exp: float = 2.0,
    time_offset: float = 1.0,
) -> NDArray[np.float64]:
    """
    Row-normalized probability that a traveler from zone i goes to zone j:

        K[i, j] ∝ population[j] ** pop_exp / (travel_time[i, j] + time_offset) ** time_exp

    with K[i, i] = 0, so every traveler leaves their home zone. The origin's population
    enters through the number of travelers, not the kernel. `time_offset` (hours) keeps
    neighboring urban zones, minutes apart, from dominating.
    """
    pop = np.asarray(population, dtype=np.float64)
    t = np.asarray(travel_time, dtype=np.float64)
    weight = pop[None, :] ** pop_exp / (t + time_offset) ** time_exp
    np.fill_diagonal(weight, 0.0)
    return row_normalize(weight)


class GravityClause(MovementClause):
    """The clause of the gravity model."""

    requirements = (
        AttributeDef(
            "population", int, Shapes.N, comment="The total population at each node."
        ),
        AttributeDef(
            "travel_time",
            float,
            Shapes.NxN,
            comment="Road travel time in hours from node i (row) to node j (column).",
        ),
        AttributeDef(
            "gravity_pop_exp",
            float,
            Shapes.Scalar,
            default_value=1.0,
            comment="Exponent on destination population.",
        ),
        AttributeDef(
            "gravity_time_exp",
            float,
            Shapes.Scalar,
            default_value=2.0,
            comment="Exponent on travel time: higher keeps travel more local.",
        ),
        AttributeDef(
            "gravity_time_offset",
            float,
            Shapes.Scalar,
            default_value=1.0,
            comment="Hours added to every travel time before applying the exponent.",
        ),
        AttributeDef(
            "commuter_proportion",
            float,
            Shapes.Scalar,
            default_value=0.1,
            comment="The proportion of each zone's population that travels each day.",
        ),
    )

    predicate = EveryDay()
    leaves = TickIndex(step=0)
    returns = TickDelta(step=1, days=0)

    @cached_property
    def dispersal_kernel(self) -> NDArray[np.float64]:
        if self.data("gravity_time_offset") <= 0:
            err = "Gravity movement model parameter 'gravity_time_offset' must be positive."
            raise DataAttributeError(err)
        return gravity_kernel(
            self.data("population"),
            self.data("travel_time"),
            pop_exp=self.data("gravity_pop_exp"),
            time_exp=self.data("gravity_time_exp"),
            time_offset=self.data("gravity_time_offset"),
        )

    @override
    def evaluate(self, tick: Tick, available: NDArray[SimDType]) -> NDArray[np.int64]:
        comm_prop = self.data("commuter_proportion")
        if comm_prop < 0:
            err = (
                "Gravity movement model parameter 'commuter_proportion' must be "
                "greater than or equal to zero."
            )
            raise DataAttributeError(err)
        n_commuters = np.floor(available * comm_prop).astype(SimDType)
        return self.rng.multinomial(n_commuters, self.dispersal_kernel)


class Gravity(MovementModel):
    """
    A fixed proportion of each zone's population travels every day to another zone,
    spends 1/3 of a day there, and returns home for the remaining 2/3. Destinations
    follow a gravity kernel: more likely when the destination is populous and close by
    road.
    """

    steps = (1 / 3, 2 / 3)
    clauses = (GravityClause(),)
