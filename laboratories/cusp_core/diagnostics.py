"""Conservation, profile, core, and physical-unit diagnostics."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import PhysicalScales
from .grid import SphericalGrid
from .hydrodynamics import enclosed_mass_faces, primitive


@dataclass(frozen=True)
class EnergyBudget:
    mass: float
    internal: float
    kinetic: float
    gravitational: float
    total: float


def gravitational_energy(density: np.ndarray, grid: SphericalGrid) -> float:
    """Exact W=-integral M(r)/r dM for piecewise-constant shell densities."""
    face_mass = enclosed_mass_faces(density, grid)
    a = grid.faces[:-1]
    b = grid.faces[1:]
    shell_factor = 4.0 * np.pi * density
    self_coefficient = shell_factor / 3.0
    integral = face_mass[:-1] * 0.5 * (b**2 - a**2)
    integral += self_coefficient * (
        (b**5 - a**5) / 5.0 - 0.5 * a**3 * (b**2 - a**2)
    )
    return float(-np.sum(shell_factor * integral))


def potential_faces(density: np.ndarray, grid: SphericalGrid) -> np.ndarray:
    """Potential with Phi(infinity)=0 for the truncated spherical fluid."""
    mass = enclosed_mass_faces(density, grid)
    phi = np.empty(density.size + 1, dtype=np.float64)
    phi[-1] = -mass[-1] / grid.faces[-1]
    for i in range(density.size - 1, -1, -1):
        a = grid.faces[i]
        b = grid.faces[i + 1]
        coefficient = 4.0 * np.pi * density[i] / 3.0
        if a == 0.0:
            acceleration_integral = 0.5 * coefficient * b**2
        else:
            constant = mass[i] - coefficient * a**3
            acceleration_integral = (
                constant * (1.0 / a - 1.0 / b)
                + 0.5 * coefficient * (b**2 - a**2)
            )
        phi[i] = phi[i + 1] - acceleration_integral
    return phi


def potential_centers(density: np.ndarray, grid: SphericalGrid) -> np.ndarray:
    face_phi = potential_faces(density, grid)
    face_mass = enclosed_mass_faces(density, grid)
    result = np.empty_like(density)
    for i, radius in enumerate(grid.centers):
        a = grid.faces[i]
        coefficient = 4.0 * np.pi * density[i] / 3.0
        if a == 0.0:
            acceleration_integral = 0.5 * coefficient * radius**2
        else:
            constant = face_mass[i] - coefficient * a**3
            acceleration_integral = (
                constant * (1.0 / a - 1.0 / radius)
                + 0.5 * coefficient * (radius**2 - a**2)
            )
        result[i] = face_phi[i] + acceleration_integral
    return result


def energy_budget(state: np.ndarray, grid: SphericalGrid, gamma: float) -> EnergyBudget:
    rho, velocity, pressure = primitive(state, gamma)
    internal_density = pressure / (gamma - 1.0)
    kinetic_density = 0.5 * rho * velocity**2
    internal = float(np.sum(internal_density * grid.volumes))
    kinetic = float(np.sum(kinetic_density * grid.volumes))
    gravitational = gravitational_energy(rho, grid)
    return EnergyBudget(
        mass=float(np.sum(rho * grid.volumes)),
        internal=internal,
        kinetic=kinetic,
        gravitational=gravitational,
        total=internal + kinetic + gravitational,
    )


def logarithmic_slope(radius: np.ndarray, density: np.ndarray) -> np.ndarray:
    """Second-order finite-difference d ln(rho)/d ln(r)."""
    radius = np.asarray(radius, dtype=np.float64)
    density = np.asarray(density, dtype=np.float64)
    if radius.ndim != 1 or density.shape != radius.shape or radius.size < 3:
        raise ValueError("Slope estimation requires matching 1-D arrays of length >= 3")
    if np.any(radius <= 0.0) or np.any(density <= 0.0):
        raise ValueError("Radius and density must be positive")
    return np.gradient(np.log(density), np.log(radius), edge_order=2)


def core_radius(
    radius: np.ndarray,
    density: np.ndarray,
    dx: float,
    minimum_resolved_cells: int = 8,
    required_core_cells: int = 16,
) -> float | None:
    """Outward gamma=-0.5 crossing, accepted only when r_core/dx >= 16."""
    slope = logarithmic_slope(radius, density)
    valid = radius >= minimum_resolved_cells * dx
    indices = np.flatnonzero(valid)
    if indices.size < 2 or slope[indices[0]] <= -0.5:
        return None
    crossing = None
    for left, right in zip(indices[:-1], indices[1:]):
        if slope[left] > -0.5 and slope[right] <= -0.5:
            fraction = (-0.5 - slope[left]) / (slope[right] - slope[left])
            crossing = float(
                np.exp(np.log(radius[left]) + fraction * (np.log(radius[right]) - np.log(radius[left])))
            )
            break
    if crossing is None or crossing / dx < required_core_cells:
        return None
    return crossing


def rotation_curve(
    density: np.ndarray, grid: SphericalGrid, scales: PhysicalScales | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """Predict v_c=sqrt(G M(<r)/r), with optional fixed physical units."""
    mass = enclosed_mass_faces(density, grid)[1:]
    radius = grid.faces[1:]
    velocity_hat = np.sqrt(mass / radius)
    if scales is None:
        return radius, velocity_hat
    radius_kpc = radius * scales.r_s_kpc
    velocity_km_s = velocity_hat * scales.velocity_scale_m_s / 1.0e3
    return radius_kpc, velocity_km_s


def relative_error(value: float, reference: float) -> float:
    denominator = abs(reference)
    if denominator == 0.0:
        return abs(value - reference)
    return abs(value - reference) / denominator
