"""Strang-split evolution driver for the independent laboratory."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .conduction import conductive_luminosity, implicit_conduction_step
from .config import LaboratoryConfig
from .diagnostics import EnergyBudget, energy_budget, potential_centers
from .equilibrium import HydrostaticNFW
from .grid import SphericalGrid
from .hydrodynamics import HydroOperator, NumericalFailure, primitive


@dataclass(frozen=True)
class Snapshot:
    time: float
    state: np.ndarray
    potential: np.ndarray
    energy: EnergyBudget
    conductive_luminosity: np.ndarray


class CuspCoreSolver:
    """Float64 solver containing no Reotransductor-specific fields or imports."""

    def __init__(
        self,
        cells: int,
        conductivity_hat: float,
        config: LaboratoryConfig | None = None,
    ) -> None:
        self.config = config or LaboratoryConfig()
        if conductivity_hat < 0.0:
            raise ValueError("Conductivity cannot be negative")
        self.conductivity_hat = float(conductivity_hat)
        self.grid = SphericalGrid.uniform(cells, self.config.radius)
        self.equilibrium = HydrostaticNFW.build(self.grid)
        self.hydrodynamics = HydroOperator(
            self.grid, self.equilibrium, self.config.gamma
        )
        self.state = self.equilibrium.conserved(self.config.gamma)
        self.time = 0.0
        primitive(self.state, self.config.gamma)

    def snapshot(self) -> Snapshot:
        return Snapshot(
            time=self.time,
            state=self.state.copy(),
            potential=potential_centers(self.state[:, 0], self.grid),
            energy=energy_budget(self.state, self.grid, self.config.gamma),
            conductive_luminosity=conductive_luminosity(
                self.state,
                self.grid,
                self.config.gamma,
                self.conductivity_hat,
            ),
        )

    def step(self, dt: float) -> None:
        if dt <= 0.0 or not np.isfinite(dt):
            raise ValueError("A finite positive timestep is required")
        half = 0.5 * dt
        first = implicit_conduction_step(
            self.state,
            self.grid,
            self.config.gamma,
            self.conductivity_hat,
            half,
        )
        hydro = self.hydrodynamics.step_rk2(first, dt)
        self.state = implicit_conduction_step(
            hydro,
            self.grid,
            self.config.gamma,
            self.conductivity_hat,
            half,
        )
        primitive(self.state, self.config.gamma)
        self.time += dt

    def evolve_to(self, target_time: float) -> Snapshot:
        if target_time < self.time:
            raise ValueError("Cannot evolve backward in time")
        while self.time < target_time:
            dt = min(
                self.hydrodynamics.maximum_timestep(self.state, self.config.cfl),
                target_time - self.time,
            )
            self.step(dt)
        self.time = float(target_time)
        return self.snapshot()


__all__ = ["CuspCoreSolver", "NumericalFailure", "Snapshot"]
