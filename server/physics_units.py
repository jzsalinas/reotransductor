"""
Physical Constants, Planck Scale Derivations, and Cosmological Unit Transformations.
Provides rigorous dimensional calibration between dimensionless computational lattice units
and fundamental physical units (SI, Planck, and Astrophysical: Mpc, Gyr, M_sun).
"""

from dataclasses import dataclass
from typing import Dict, Any
import numpy as np


@dataclass(frozen=True)
class FundamentalConstants:
    """
    CODATA 2022 Recommended Values for Fundamental Physical Constants.
    """
    # Speed of light in vacuum (m / s)
    C: float = 299792458.0

    # Newtonian constant of gravitation (m^3 / (kg * s^2))
    G: float = 6.67430e-11

    # Reduced Planck constant (J * s = kg * m^2 / s)
    HBAR: float = 1.054571817e-34

    # Boltzmann constant (J / K = kg * m^2 / (s^2 * K))
    K_B: float = 1.380649e-23

    # Astronomical & Astrophysical Constants
    # Solar mass (kg)
    M_SUN: float = 1.98847e30

    # Megaparsec in meters (m)
    MPC: float = 3.08567758128e22

    # Year in seconds (Julian year = 365.25 days)
    YEAR_S: float = 31557600.0

    # Million years (Myr in s)
    MYR_S: float = 3.15576e13

    # Billion years (Gyr in s)
    GYR_S: float = 3.15576e16

    # Baseline Hubble Constant H0 (km / s / Mpc)
    H0_BASELINE: float = 70.0


@dataclass(frozen=True)
class PlanckScales:
    """
    Derived Planck Scale Units from Fundamental Constants.
    """
    # Planck length: l_P = sqrt(hbar * G / c^3) (m)
    LENGTH: float = np.sqrt(FundamentalConstants.HBAR * FundamentalConstants.G / (FundamentalConstants.C ** 3))

    # Planck time: t_P = sqrt(hbar * G / c^5) (s)
    TIME: float = np.sqrt(FundamentalConstants.HBAR * FundamentalConstants.G / (FundamentalConstants.C ** 5))

    # Planck mass: m_P = sqrt(hbar * c / G) (kg)
    MASS: float = np.sqrt(FundamentalConstants.HBAR * FundamentalConstants.C / FundamentalConstants.G)

    # Planck density: rho_P = c^5 / (hbar * G^2) (kg / m^3)
    DENSITY: float = (FundamentalConstants.C ** 5) / (FundamentalConstants.HBAR * (FundamentalConstants.G ** 2))

    # Planck temperature: T_P = sqrt(hbar * c^5 / (G * k_B^2)) (K)
    TEMPERATURE: float = np.sqrt(
        (FundamentalConstants.HBAR * (FundamentalConstants.C ** 5)) /
        (FundamentalConstants.G * (FundamentalConstants.K_B ** 2))
    )

    # Fundamental Rheotransducer Coupling Constant:
    # kappa_0 = (l_P^3 * t_P) / k_B = (hbar^2 * G^2) / (c^7 * k_B) (m * s^3 * K / kg)
    KAPPA_0: float = (
        (FundamentalConstants.HBAR ** 2) * (FundamentalConstants.G ** 2)
    ) / (
        (FundamentalConstants.C ** 7) * FundamentalConstants.K_B
    )


class CosmologicalUnits:
    """
    Bidirectional dimensional converter between computational grid units and physical units.
    Calibrates a 3D simulation lattice of box size L_box (default 100 Mpc) and resolution N (default 32).
    """

    def __init__(
        self,
        box_size_mpc: float = 100.0,
        grid_resolution: int = 32,
        c_code: float = 2.5,
        h0_km_s_mpc: float = 70.0
    ):
        self.constants = FundamentalConstants()
        self.planck = PlanckScales()

        self.box_size_mpc = float(box_size_mpc)
        self.grid_resolution = int(grid_resolution)
        self.c_code = float(c_code)
        self.h0_kms_mpc = float(h0_km_s_mpc)

        # Spatial step per cell in meters: Delta x = (L_box / N) in meters
        self.cell_size_mpc = self.box_size_mpc / self.grid_resolution
        self.delta_x_m = self.cell_size_mpc * self.constants.MPC

        # Temporal scale per code time unit Delta t = 1.0 (in seconds)
        # Defined such that c_code cells / Delta t_code equals real speed of light c:
        # c = (c_code * Delta x_m) / Delta t_s  ==>  Delta t_s = (c_code * Delta x_m) / c
        self.time_unit_s = (self.c_code * self.delta_x_m) / self.constants.C
        self.time_unit_myr = self.time_unit_s / self.constants.MYR_S

        # Critical Cosmic Density: rho_crit = 3 * H0^2 / (8 * pi * G) in kg / m^3
        # H0 in SI (1 / s): H0_si = H0 * 1000 / (1 Mpc in m)
        self.h0_si = (self.h0_kms_mpc * 1000.0) / self.constants.MPC
        self.rho_crit_si = (3.0 * (self.h0_si ** 2)) / (8.0 * np.pi * self.constants.G)
        self.rho_crit_msun_mpc3 = self.rho_crit_si * (self.constants.MPC ** 3) / self.constants.M_SUN

        # Velocity scale: 1 cell / code_time_unit in km / s
        self.velocity_unit_km_s = (self.delta_x_m / self.time_unit_s) / 1000.0

        # Gravitational potential scale: (Delta x / Delta t)^2 in (m/s)^2
        self.potential_unit_si = (self.delta_x_m / self.time_unit_s) ** 2

    # =========================================================================
    # DENSITY CONVERSIONS (rho_code == 1.0 represents cosmic mean critical density)
    # =========================================================================

    def density_code_to_si(self, rho_code: np.ndarray | float) -> np.ndarray | float:
        """Converts dimensionless code density to physical density in kg / m^3."""
        return rho_code * self.rho_crit_si

    def density_si_to_code(self, rho_si: np.ndarray | float) -> np.ndarray | float:
        """Converts physical density in kg / m^3 to dimensionless code density."""
        return rho_si / self.rho_crit_si

    def density_code_to_astrophysical(self, rho_code: np.ndarray | float) -> np.ndarray | float:
        """Converts dimensionless code density to M_sun / Mpc^3."""
        return rho_code * self.rho_crit_msun_mpc3

    def density_astrophysical_to_code(self, rho_msun_mpc3: np.ndarray | float) -> np.ndarray | float:
        """Converts M_sun / Mpc^3 to dimensionless code density."""
        return rho_msun_mpc3 / self.rho_crit_msun_mpc3

    # =========================================================================
    # VELOCITY CONVERSIONS
    # =========================================================================

    def velocity_code_to_km_s(self, v_code: np.ndarray | float) -> np.ndarray | float:
        """Converts velocity from code cells / Delta t to km / s."""
        return v_code * self.velocity_unit_km_s

    def velocity_km_s_to_code(self, v_km_s: np.ndarray | float) -> np.ndarray | float:
        """Converts velocity from km / s to code cells / Delta t."""
        return v_km_s / self.velocity_unit_km_s

    # =========================================================================
    # TIME CONVERSIONS
    # =========================================================================

    def time_code_to_seconds(self, t_code: np.ndarray | float) -> np.ndarray | float:
        """Converts code time units to physical seconds."""
        return t_code * self.time_unit_s

    def time_code_to_myr(self, t_code: np.ndarray | float) -> np.ndarray | float:
        """Converts code time units to Million Years (Myr)."""
        return t_code * self.time_unit_myr

    def time_code_to_gyr(self, t_code: np.ndarray | float) -> np.ndarray | float:
        """Converts code time units to Billion Years (Gyr)."""
        return (t_code * self.time_unit_myr) / 1000.0

    def time_myr_to_code(self, t_myr: np.ndarray | float) -> np.ndarray | float:
        """Converts Million Years (Myr) to code time units."""
        return t_myr / self.time_unit_myr

    # =========================================================================
    # TEMPERATURE CONVERSIONS
    # =========================================================================

    @staticmethod
    def temperature_code_to_astrophysical(t_code: np.ndarray | float) -> np.ndarray | float:
        """
        Maps normalized code temperature T in [2.73, 50.0] to physical plasma Kelvins.
        2.73 code units corresponds to CMB base (2.73 K).
        High-energy shock zones reach up to ~6,000 K astrophysical.
        """
        return t_code * 120.0

    @staticmethod
    def temperature_astrophysical_to_code(t_kelvin: np.ndarray | float) -> np.ndarray | float:
        """Maps physical plasma Kelvins to code temperature units."""
        return t_kelvin / 120.0

    # =========================================================================
    # KAPPA & DISSIPATIVE COUPLING
    # =========================================================================

    def get_fundamental_kappa_0(self) -> float:
        """
        Returns the derived fundamental Planck-Boltzmann coupling constant kappa_0 in SI units:
        kappa_0 = (hbar^2 * G^2) / (c^7 * k_B) [m * s^3 * K / kg]
        """
        return self.planck.KAPPA_0

    def summary(self) -> Dict[str, Any]:
        """Returns structured dictionary of the cosmological dimensional calibration."""
        return {
            "box_size_mpc": self.box_size_mpc,
            "grid_resolution": f"{self.grid_resolution}^3 ({self.grid_resolution**3:,} cells)",
            "cell_size_mpc": round(self.cell_size_mpc, 4),
            "cell_size_meters": f"{self.delta_x_m:.4e} m",
            "time_unit_myr": round(self.time_unit_myr, 2),
            "time_unit_seconds": f"{self.time_unit_s:.4e} s",
            "rho_crit_si": f"{self.rho_crit_si:.4e} kg/m^3",
            "rho_crit_astrophysical": f"{self.rho_crit_msun_mpc3:.4e} M_sun/Mpc^3",
            "velocity_unit_km_s": round(self.velocity_unit_km_s, 2),
            "planck_kappa_0_si": f"{self.planck.KAPPA_0:.4e} m*s^3*K/kg",
            "planck_length_m": f"{self.planck.LENGTH:.4e} m",
            "planck_time_s": f"{self.planck.TIME:.4e} s",
            "planck_density_si": f"{self.planck.DENSITY:.4e} kg/m^3",
            "planck_temperature_k": f"{self.planck.TEMPERATURE:.4e} K"
        }

    def print_summary(self) -> None:
        """Prints a human-readable table of physical units calibration."""
        s = self.summary()
        print("=" * 70)
        print("  REOTRANSDUCTOR: COSMOLOGICAL PHYSICAL UNITS & PLANCK CALIBRATION")
        print("=" * 70)
        for k, v in s.items():
            print(f"  * {k.replace('_', ' ').title():<32} : {v}")
        print("=" * 70)
