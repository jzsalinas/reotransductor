"""Frozen physical and explicitly classified numerical configuration."""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt


G_SI = 6.67430e-11
MSUN_KG = 1.98847e30
PC_M = 3.085677581491367e16
KPC_M = 1.0e3 * PC_M


@dataclass(frozen=True)
class PhysicalScales:
    """Reference realization; these are scale definitions, not fitted data."""

    r_s_kpc: float = 1.0
    rho_s_msun_pc3: float = 1.0

    @property
    def r_s_m(self) -> float:
        return self.r_s_kpc * KPC_M

    @property
    def rho_s_kg_m3(self) -> float:
        return self.rho_s_msun_pc3 * MSUN_KG / PC_M**3

    @property
    def mass_scale_kg(self) -> float:
        return self.rho_s_kg_m3 * self.r_s_m**3

    @property
    def velocity_scale_m_s(self) -> float:
        return sqrt(G_SI * self.rho_s_kg_m3 * self.r_s_m**2)

    @property
    def time_scale_s(self) -> float:
        return 1.0 / sqrt(G_SI * self.rho_s_kg_m3)


@dataclass(frozen=True)
class LaboratoryConfig:
    """Frozen model values plus declared numerical-method parameters."""

    gamma: float = 5.0 / 3.0
    radius: float = 10.0
    conductivity_hat_star: float = 1.5
    cfl: float = 0.35
    minimum_resolved_cells: int = 8
    required_core_cells: int = 16
    mass_tolerance: float = 1.0e-12
    energy_tolerance: float = 1.0e-6
    control_velocity_tolerance: float = 1.0e-10
    control_density_l1_tolerance: float = 1.0e-10
    control_slope_tolerance: float = 1.0e-8

    def reference_conduction_time(self, conductivity_hat: float) -> float:
        """Return t_cond,0/t0, not the local time at r_s."""
        if conductivity_hat <= 0.0:
            return float("inf")
        return 1.0 / ((self.gamma - 1.0) * conductivity_hat)
