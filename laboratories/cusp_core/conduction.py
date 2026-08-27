"""Conservative implicit Fourier conduction through spherical luminosities."""

from __future__ import annotations

import numpy as np

from .grid import SphericalGrid
from .hydrodynamics import NumericalFailure, primitive


def conductive_luminosity(
    state: np.ndarray, grid: SphericalGrid, gamma: float, conductivity_hat: float
) -> np.ndarray:
    """Return L=4*pi*r^2*q; boundary luminosities are exactly zero."""
    if conductivity_hat < 0.0:
        raise ValueError("Conductivity cannot be negative")
    rho, _, pressure = primitive(state, gamma)
    theta = pressure / rho
    luminosity = np.zeros(rho.size + 1, dtype=np.float64)
    if conductivity_hat == 0.0:
        return luminosity
    gradients = (theta[1:] - theta[:-1]) / np.diff(grid.centers)
    luminosity[1:-1] = -conductivity_hat * grid.areas[1:-1] * gradients
    # luminosity[0] = 0 implements L_cond(0)=0 without evaluating q(0).
    # luminosity[-1] = 0 is the insulating outer wall.
    return luminosity


def _solve_tridiagonal(
    lower: np.ndarray, diagonal: np.ndarray, upper: np.ndarray, rhs: np.ndarray
) -> np.ndarray:
    """Float64 Thomas solve for a strictly diagonally dominant system."""
    n = diagonal.size
    cprime = np.empty(max(n - 1, 0), dtype=np.float64)
    dprime = np.empty(n, dtype=np.float64)
    pivot = diagonal[0]
    if pivot <= 0.0 or not np.isfinite(pivot):
        raise NumericalFailure("Invalid implicit-conduction matrix pivot")
    if n > 1:
        cprime[0] = upper[0] / pivot
    dprime[0] = rhs[0] / pivot
    for i in range(1, n):
        pivot = diagonal[i] - lower[i - 1] * cprime[i - 1]
        if pivot <= 0.0 or not np.isfinite(pivot):
            raise NumericalFailure("Invalid implicit-conduction matrix pivot")
        if i < n - 1:
            cprime[i] = upper[i] / pivot
        dprime[i] = (rhs[i] - lower[i - 1] * dprime[i - 1]) / pivot
    solution = np.empty(n, dtype=np.float64)
    solution[-1] = dprime[-1]
    for i in range(n - 2, -1, -1):
        solution[i] = dprime[i] - cprime[i] * solution[i + 1]
    return solution


def implicit_conduction_step(
    state: np.ndarray,
    grid: SphericalGrid,
    gamma: float,
    conductivity_hat: float,
    dt: float,
) -> np.ndarray:
    """Crank--Nicolson conservative heat step at fixed rho and momentum."""
    if conductivity_hat < 0.0:
        raise ValueError("Conductivity cannot be negative")
    if dt < 0.0:
        raise ValueError("Timestep cannot be negative")
    if conductivity_hat == 0.0 or dt == 0.0:
        return state.copy()
    rho, _, pressure = primitive(state, gamma)
    theta = pressure / rho
    capacity = grid.volumes * rho / (gamma - 1.0)
    conductance = conductivity_hat * grid.areas[1:-1] / np.diff(grid.centers)

    half_dt = 0.5 * dt
    diffusion = np.zeros_like(theta)
    face_exchange = conductance * (theta[1:] - theta[:-1])
    diffusion[:-1] += face_exchange
    diffusion[1:] -= face_exchange
    diagonal = capacity.copy()
    diagonal[:-1] += half_dt * conductance
    diagonal[1:] += half_dt * conductance
    lower = -half_dt * conductance.copy()
    upper = -half_dt * conductance.copy()
    theta_new = _solve_tridiagonal(
        lower, diagonal, upper, capacity * theta + half_dt * diffusion
    )
    if np.any(theta_new <= 0.0) or not np.all(np.isfinite(theta_new)):
        raise NumericalFailure("Implicit conduction produced non-positive temperature")

    updated = state.copy()
    kinetic = 0.5 * state[:, 1] ** 2 / rho
    updated[:, 2] = rho * theta_new / (gamma - 1.0) + kinetic
    primitive(updated, gamma)
    return updated
