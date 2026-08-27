"""Second-order well-balanced spherical Euler--Poisson discretization."""

from __future__ import annotations

import numpy as np

from .equilibrium import HydrostaticNFW
from .grid import SphericalGrid


class NumericalFailure(RuntimeError):
    """Raised when a state violates a non-negotiable numerical condition."""


def primitive(state: np.ndarray, gamma: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert conserved variables (rho, rho*u, E) without state repair."""
    rho = state[:, 0]
    if not np.all(np.isfinite(state)):
        raise NumericalFailure("Non-finite conserved state")
    if np.any(rho <= 0.0):
        raise NumericalFailure("Density became non-positive")
    velocity = state[:, 1] / rho
    internal = state[:, 2] - 0.5 * state[:, 1] ** 2 / rho
    if np.any(internal <= 0.0) or not np.all(np.isfinite(internal)):
        raise NumericalFailure("Internal energy became non-positive or non-finite")
    pressure = (gamma - 1.0) * internal
    return rho, velocity, pressure


def _minmod3(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    same_sign = (np.sign(a) == np.sign(b)) & (np.sign(b) == np.sign(c))
    magnitude = np.minimum(np.abs(a), np.minimum(np.abs(b), np.abs(c)))
    return np.where(same_sign, np.sign(a) * magnitude, 0.0)


def monotonized_central(values: np.ndarray) -> np.ndarray:
    """Monotonized-central (MC) limited increments on a uniform mesh."""
    slope = np.zeros_like(values)
    backward = values[1:-1] - values[:-2]
    centered = 0.5 * (values[2:] - values[:-2])
    forward = values[2:] - values[1:-1]
    slope[1:-1] = _minmod3(2.0 * backward, centered, 2.0 * forward)
    return slope


def _conserved_from_primitive(
    rho: np.ndarray, velocity: np.ndarray, pressure: np.ndarray, gamma: float
) -> np.ndarray:
    result = np.empty((rho.size, 3), dtype=np.float64)
    result[:, 0] = rho
    result[:, 1] = rho * velocity
    result[:, 2] = pressure / (gamma - 1.0) + 0.5 * rho * velocity**2
    return result


def hllc_flux(left: np.ndarray, right: np.ndarray, gamma: float) -> np.ndarray:
    """HLLC flux for one-dimensional Euler states stored as primitive rows."""
    rho_l, u_l, p_l = left.T
    rho_r, u_r, p_r = right.T
    if np.any(rho_l <= 0.0) or np.any(rho_r <= 0.0):
        raise NumericalFailure("MUSCL reconstruction produced non-positive density")
    if np.any(p_l <= 0.0) or np.any(p_r <= 0.0):
        raise NumericalFailure("MUSCL reconstruction produced non-positive pressure")

    e_l = p_l / (gamma - 1.0) + 0.5 * rho_l * u_l**2
    e_r = p_r / (gamma - 1.0) + 0.5 * rho_r * u_r**2
    c_l = np.sqrt(gamma * p_l / rho_l)
    c_r = np.sqrt(gamma * p_r / rho_r)
    speed_l = np.minimum(u_l - c_l, u_r - c_r)
    speed_r = np.maximum(u_l + c_l, u_r + c_r)
    denominator = rho_l * (speed_l - u_l) - rho_r * (speed_r - u_r)
    if np.any(denominator == 0.0):
        raise NumericalFailure("Degenerate HLLC contact-wave denominator")
    speed_m = (
        p_r
        - p_l
        + rho_l * u_l * (speed_l - u_l)
        - rho_r * u_r * (speed_r - u_r)
    ) / denominator

    state_l = np.column_stack((rho_l, rho_l * u_l, e_l))
    state_r = np.column_stack((rho_r, rho_r * u_r, e_r))
    flux_l = np.column_stack(
        (rho_l * u_l, rho_l * u_l**2 + p_l, u_l * (e_l + p_l))
    )
    flux_r = np.column_stack(
        (rho_r * u_r, rho_r * u_r**2 + p_r, u_r * (e_r + p_r))
    )

    p_star_l = p_l + rho_l * (speed_l - u_l) * (speed_m - u_l)
    p_star_r = p_r + rho_r * (speed_r - u_r) * (speed_m - u_r)
    p_star = 0.5 * (p_star_l + p_star_r)
    rho_star_l = rho_l * (speed_l - u_l) / (speed_l - speed_m)
    rho_star_r = rho_r * (speed_r - u_r) / (speed_r - speed_m)
    e_star_l = (
        (speed_l - u_l) * e_l - p_l * u_l + p_star * speed_m
    ) / (speed_l - speed_m)
    e_star_r = (
        (speed_r - u_r) * e_r - p_r * u_r + p_star * speed_m
    ) / (speed_r - speed_m)
    star_l = np.column_stack((rho_star_l, rho_star_l * speed_m, e_star_l))
    star_r = np.column_stack((rho_star_r, rho_star_r * speed_m, e_star_r))

    result = np.empty_like(flux_l)
    mask_l = speed_l >= 0.0
    mask_star_l = (speed_l < 0.0) & (speed_m >= 0.0)
    mask_star_r = (speed_m < 0.0) & (speed_r > 0.0)
    mask_r = speed_r <= 0.0
    result[mask_l] = flux_l[mask_l]
    result[mask_star_l] = (
        flux_l[mask_star_l]
        + speed_l[mask_star_l, None] * (star_l[mask_star_l] - state_l[mask_star_l])
    )
    result[mask_star_r] = (
        flux_r[mask_star_r]
        + speed_r[mask_star_r, None] * (star_r[mask_star_r] - state_r[mask_star_r])
    )
    result[mask_r] = flux_r[mask_r]
    return result


def enclosed_mass_faces(density: np.ndarray, grid: SphericalGrid) -> np.ndarray:
    mass = np.empty(density.size + 1, dtype=np.float64)
    mass[0] = 0.0
    mass[1:] = np.cumsum(density * grid.volumes)
    return mass


def center_acceleration(density: np.ndarray, grid: SphericalGrid) -> np.ndarray:
    """Piecewise-constant finite-volume enclosed-mass acceleration at centers."""
    face_mass = enclosed_mass_faces(density, grid)
    partial_volume = (4.0 * np.pi / 3.0) * (
        grid.centers**3 - grid.faces[:-1] ** 3
    )
    center_mass = face_mass[:-1] + density * partial_volume
    return center_mass / grid.centers**2


class HydroOperator:
    """Euler--Poisson operator preserving the specified equilibrium exactly."""

    def __init__(
        self, grid: SphericalGrid, equilibrium: HydrostaticNFW, gamma: float
    ) -> None:
        self.grid = grid
        self.equilibrium = equilibrium
        self.gamma = gamma
        self.equilibrium_acceleration = center_acceleration(equilibrium.density, grid)

    def _face_fluxes(self, state: np.ndarray) -> np.ndarray:
        rho, velocity, pressure = primitive(state, self.gamma)
        perturbation = np.column_stack(
            (
                rho - self.equilibrium.density,
                velocity,
                pressure - self.equilibrium.pressure,
            )
        )
        slopes = np.column_stack(
            [monotonized_central(perturbation[:, i]) for i in range(3)]
        )
        count = rho.size
        flux = np.zeros((count + 1, 3), dtype=np.float64)
        equilibrium_flux = np.zeros_like(flux)

        if count > 1:
            face_eq = np.column_stack(
                (
                    np.asarray(self.equilibrium.enclosed_mass_faces[1:-1] * 0.0),
                    np.zeros(count - 1),
                    self.equilibrium.face_pressure[1:-1],
                )
            )
            face_eq[:, 0] = 1.0 / (
                self.grid.faces[1:-1] * (1.0 + self.grid.faces[1:-1]) ** 2
            )
            left = face_eq + perturbation[:-1] + 0.5 * slopes[:-1]
            right = face_eq + perturbation[1:] - 0.5 * slopes[1:]
            flux[1:-1] = hllc_flux(left, right, self.gamma)
            equilibrium_flux[1:-1, 1] = self.equilibrium.face_pressure[1:-1]

        # Both boundaries are impermeable. At r=0 only the integrated flux is
        # meaningful and its area is exactly zero. The outer wall transmits a
        # pressure force but no mass or energy.
        outer_pressure = (
            self.equilibrium.face_pressure[-1]
            + perturbation[-1, 2]
            + 0.5 * slopes[-1, 2]
        )
        if outer_pressure <= 0.0:
            raise NumericalFailure("Outer wall reconstruction produced non-positive pressure")
        flux[-1, 1] = outer_pressure
        equilibrium_flux[-1, 1] = self.equilibrium.face_pressure[-1]
        return flux - equilibrium_flux

    def rhs(self, state: np.ndarray) -> np.ndarray:
        rho, velocity, pressure = primitive(state, self.gamma)
        residual_flux = self._face_fluxes(state)
        integrated_flux = self.grid.areas[:, None] * residual_flux
        result = -(integrated_flux[1:] - integrated_flux[:-1]) / self.grid.volumes[:, None]

        acceleration = center_acceleration(rho, self.grid)
        pressure_delta = pressure - self.equilibrium.pressure
        geometric_integral = (
            4.0 * np.pi * pressure_delta * (self.grid.faces[1:] ** 2 - self.grid.faces[:-1] ** 2)
        )
        gravity_delta = rho * acceleration - self.equilibrium.density * self.equilibrium_acceleration
        result[:, 1] += geometric_integral / self.grid.volumes - gravity_delta
        result[:, 2] += -rho * velocity * acceleration
        return result

    def maximum_timestep(self, state: np.ndarray, cfl: float) -> float:
        rho, velocity, pressure = primitive(state, self.gamma)
        signal_speed = np.abs(velocity) + np.sqrt(self.gamma * pressure / rho)
        maximum = float(np.max(signal_speed))
        if maximum <= 0.0 or not np.isfinite(maximum):
            raise NumericalFailure("Invalid Euler signal speed")
        return cfl * self.grid.dr / maximum

    def step_rk2(self, state: np.ndarray, dt: float) -> np.ndarray:
        """Second-order SSP Runge--Kutta update."""
        stage = state + dt * self.rhs(state)
        primitive(stage, self.gamma)
        updated = 0.5 * state + 0.5 * (stage + dt * self.rhs(stage))
        primitive(updated, self.gamma)
        return updated
