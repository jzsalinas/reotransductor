"""Hydrostatic NFW initial state without point evaluation at the cusp."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss

from .grid import SphericalGrid


FOUR_PI = 4.0 * np.pi
_QUAD_X, _QUAD_W = leggauss(32)


def nfw_mass(radius: np.ndarray | float) -> np.ndarray | float:
    """Exact dimensionless NFW mass M(<r) for G=r_s=rho_s=1."""
    x = np.asarray(radius, dtype=np.float64)
    if np.any(x < 0.0):
        raise ValueError("Radius cannot be negative")
    regular = np.log1p(x) - x / (1.0 + x)
    small = x < 1.0e-4
    if np.any(small):
        xs = x[small] if x.ndim else x
        series = (
            0.5 * xs**2
            - (2.0 / 3.0) * xs**3
            + 0.75 * xs**4
            - 0.8 * xs**5
            + (5.0 / 6.0) * xs**6
        )
        if x.ndim:
            regular = regular.copy()
            regular[small] = series
        else:
            regular = np.asarray(series)
    result = FOUR_PI * regular
    return float(result) if np.ndim(radius) == 0 else result


def nfw_density(radius: np.ndarray | float) -> np.ndarray | float:
    """Pointwise NFW density for strictly positive radii only."""
    x = np.asarray(radius, dtype=np.float64)
    if np.any(x <= 0.0):
        raise ValueError("The NFW density is not evaluated at r <= 0")
    result = 1.0 / (x * (1.0 + x) ** 2)
    return float(result) if np.ndim(radius) == 0 else result


def nfw_acceleration(radius: np.ndarray | float) -> np.ndarray | float:
    """Inward acceleration magnitude G M(<r)/r^2 with G=1."""
    x = np.asarray(radius, dtype=np.float64)
    if np.any(x <= 0.0):
        raise ValueError("Acceleration is only evaluated at positive radii")
    result = np.asarray(nfw_mass(x)) / x**2
    return float(result) if np.ndim(radius) == 0 else result


def _integrate_cells(grid: SphericalGrid, function) -> np.ndarray:
    """High-order deterministic quadrature on every finite-volume interval."""
    a = grid.faces[:-1, None]
    b = grid.faces[1:, None]
    nodes = 0.5 * (a + b) + 0.5 * (b - a) * _QUAD_X[None, :]
    values = function(nodes)
    return 0.5 * (b[:, 0] - a[:, 0]) * (values @ _QUAD_W)


@dataclass(frozen=True)
class HydrostaticNFW:
    """Cell averages and finite face values for the hydrostatic state."""

    density: np.ndarray
    pressure: np.ndarray
    face_pressure: np.ndarray
    enclosed_mass_faces: np.ndarray
    outer_theta: float
    hydrostatic_drop: np.ndarray

    @classmethod
    def build(cls, grid: SphericalGrid) -> "HydrostaticNFW":
        face_mass = np.asarray(nfw_mass(grid.faces), dtype=np.float64)
        density = np.diff(face_mass) / grid.volumes

        radius = grid.faces[-1]
        rho_outer = nfw_density(radius)
        g_outer = nfw_acceleration(radius)
        dlogrho_outer = -1.0 / radius - 2.0 / (1.0 + radius)
        theta_outer = -g_outer / dlogrho_outer
        pressure_outer = rho_outer * theta_outer

        def rho_g(r: np.ndarray) -> np.ndarray:
            return nfw_density(r) * nfw_acceleration(r)

        # The first pressure drop diverges logarithmically and is never formed.
        pressure_drop = np.full(grid.centers.size, np.nan, dtype=np.float64)
        if grid.centers.size > 1:
            subgrid = SphericalGrid(
                faces=grid.faces[1:].copy(),
                centers=grid.centers[1:].copy(),
                volumes=grid.volumes[1:].copy(),
                areas=grid.areas[1:].copy(),
                dr=grid.dr,
            )
            pressure_drop[1:] = _integrate_cells(subgrid, rho_g)

        face_pressure = np.full(grid.faces.size, np.nan, dtype=np.float64)
        face_pressure[-1] = pressure_outer
        for i in range(grid.centers.size - 1, 0, -1):
            face_pressure[i] = face_pressure[i + 1] + pressure_drop[i]

        # Integration by parts gives the exact volume-average definition while
        # avoiding P(0): integral(r^2 P dr) = [r^3 P]/3 + integral(r^3 rho g dr)/3.
        gravity_moment = _integrate_cells(
            grid, lambda r: r**3 * nfw_density(r) * nfw_acceleration(r)
        )
        outer_term = grid.faces[1:] ** 3 * face_pressure[1:]
        inner_term = np.zeros_like(outer_term)
        inner_term[1:] = grid.faces[1:-1] ** 3 * face_pressure[1:-1]
        pressure_integral = (outer_term - inner_term + gravity_moment) / 3.0
        pressure = FOUR_PI * pressure_integral / grid.volumes

        if not np.all(np.isfinite(density)) or not np.all(density > 0.0):
            raise ArithmeticError("Invalid NFW cell-average density")
        if not np.all(np.isfinite(pressure)) or not np.all(pressure > 0.0):
            raise ArithmeticError("Invalid hydrostatic cell-average pressure")
        return cls(
            density=density,
            pressure=pressure,
            face_pressure=face_pressure,
            enclosed_mass_faces=face_mass,
            outer_theta=float(theta_outer),
            hydrostatic_drop=pressure_drop,
        )

    def conserved(self, gamma: float) -> np.ndarray:
        state = np.zeros((self.density.size, 3), dtype=np.float64)
        state[:, 0] = self.density
        state[:, 2] = self.pressure / (gamma - 1.0)
        return state
