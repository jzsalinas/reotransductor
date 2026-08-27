"""CCL-NUM-2 well-balanced spherical Euler--Poisson discretization."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .equilibrium import HydrostaticNFW, nfw_density
from .grid import SphericalGrid


class NumericalFailure(RuntimeError):
    """Raised when a non-negotiable numerical condition is violated."""


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
    """Monotonized-central limited cell increments on a uniform mesh."""
    slope = np.zeros_like(values)
    backward = values[1:-1] - values[:-2]
    centered = 0.5 * (values[2:] - values[:-2])
    forward = values[2:] - values[1:-1]
    slope[1:-1] = _minmod3(2.0 * backward, centered, 2.0 * forward)
    return slope


def _primitive_rows(state: np.ndarray, gamma: float) -> np.ndarray:
    rho, velocity, pressure = primitive(state, gamma)
    return np.column_stack((rho, velocity, pressure))


def hllc_flux(left: np.ndarray, right: np.ndarray, gamma: float) -> np.ndarray:
    """HLLC flux for one-dimensional Euler states stored as primitive rows."""
    rho_l, u_l, p_l = left.T
    rho_r, u_r, p_r = right.T
    if np.any(rho_l <= 0.0) or np.any(rho_r <= 0.0):
        raise NumericalFailure("Reconstruction produced non-positive density")
    if np.any(p_l <= 0.0) or np.any(p_r <= 0.0):
        raise NumericalFailure("Reconstruction produced non-positive pressure")

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
        p_r - p_l
        + rho_l * u_l * (speed_l - u_l)
        - rho_r * u_r * (speed_r - u_r)
    ) / denominator

    state_l = np.column_stack((rho_l, rho_l * u_l, e_l))
    state_r = np.column_stack((rho_r, rho_r * u_r, e_r))
    flux_l = np.column_stack((rho_l * u_l, rho_l * u_l**2 + p_l, u_l * (e_l + p_l)))
    flux_r = np.column_stack((rho_r * u_r, rho_r * u_r**2 + p_r, u_r * (e_r + p_r)))

    p_star_l = p_l + rho_l * (speed_l - u_l) * (speed_m - u_l)
    p_star_r = p_r + rho_r * (speed_r - u_r) * (speed_m - u_r)
    p_star = 0.5 * (p_star_l + p_star_r)
    rho_star_l = rho_l * (speed_l - u_l) / (speed_l - speed_m)
    rho_star_r = rho_r * (speed_r - u_r) / (speed_r - speed_m)
    e_star_l = ((speed_l - u_l) * e_l - p_l * u_l + p_star * speed_m) / (speed_l - speed_m)
    e_star_r = ((speed_r - u_r) * e_r - p_r * u_r + p_star * speed_m) / (speed_r - speed_m)
    star_l = np.column_stack((rho_star_l, rho_star_l * speed_m, e_star_l))
    star_r = np.column_stack((rho_star_r, rho_star_r * speed_m, e_star_r))

    result = np.empty_like(flux_l)
    mask_l = speed_l >= 0.0
    mask_star_l = (speed_l < 0.0) & (speed_m >= 0.0)
    mask_star_r = (speed_m < 0.0) & (speed_r > 0.0)
    mask_r = speed_r <= 0.0
    result[mask_l] = flux_l[mask_l]
    result[mask_star_l] = flux_l[mask_star_l] + speed_l[mask_star_l, None] * (
        star_l[mask_star_l] - state_l[mask_star_l]
    )
    result[mask_star_r] = flux_r[mask_star_r] + speed_r[mask_star_r, None] * (
        star_r[mask_star_r] - state_r[mask_star_r]
    )
    result[mask_r] = flux_r[mask_r]
    return result


def enclosed_mass_faces(density: np.ndarray, grid: SphericalGrid) -> np.ndarray:
    mass = np.empty(density.size + 1, dtype=np.float64)
    mass[0] = 0.0
    mass[1:] = np.cumsum(density * grid.volumes)
    return mass


def center_acceleration(density: np.ndarray, grid: SphericalGrid) -> np.ndarray:
    """Piecewise-constant enclosed-mass acceleration retained for diagnostics."""
    face_mass = enclosed_mass_faces(density, grid)
    partial_volume = (4.0 * np.pi / 3.0) * (grid.centers**3 - grid.faces[:-1] ** 3)
    center_mass = face_mass[:-1] + density * partial_volume
    return center_mass / grid.centers**2


def _admissibility_polynomial(average: np.ndarray, delta: np.ndarray) -> tuple[float, float, float]:
    rho, momentum, energy = average
    drho, dmomentum, denergy = delta
    return (
        rho * energy - 0.5 * momentum**2,
        rho * denergy + energy * drho - momentum * dmomentum,
        drho * denergy - 0.5 * dmomentum**2,
    )


def _first_internal_energy_root(average: np.ndarray, raw: np.ndarray) -> float | None:
    """First alpha in (0,1] at which rho*E-momentum^2/2 reaches zero."""
    c0, c1, c2 = _admissibility_polynomial(average, raw - average)
    if c0 <= 0.0 or not np.isfinite(c0):
        raise NumericalFailure("Cell average has non-positive internal energy")
    scale = max(abs(c0), abs(c1), abs(c2), 1.0)
    roots: list[float] = []
    if abs(c2) <= np.finfo(np.float64).eps * scale:
        if c1 < 0.0:
            roots.append(-c0 / c1)
    else:
        discriminant = c1 * c1 - 4.0 * c2 * c0
        if discriminant >= 0.0:
            square_root = np.sqrt(discriminant)
            roots.extend(((-c1 - square_root) / (2.0 * c2), (-c1 + square_root) / (2.0 * c2)))
    admissible_roots = [root for root in roots if np.isfinite(root) and 0.0 < root <= 1.0]
    return min(admissible_roots) if admissible_roots else None


def conservative_positivity_scale(
    average: np.ndarray, raw_states: np.ndarray
) -> tuple[np.ndarray, float]:
    """Use the largest common alpha that makes every supplied state admissible."""
    average = np.asarray(average, dtype=np.float64)
    raw_states = np.asarray(raw_states, dtype=np.float64)
    if average.shape != (3,) or raw_states.ndim != 2 or raw_states.shape[1] != 3:
        raise ValueError("Expected one conserved average and a row of reconstructed states")
    primitive(average[None, :], gamma=5.0 / 3.0)  # gamma does not affect admissibility.
    alpha = 1.0
    for raw in raw_states:
        delta = raw - average
        if raw[0] <= 0.0 and delta[0] < 0.0:
            alpha = min(alpha, -average[0] / delta[0])
        root = _first_internal_energy_root(average, raw)
        if root is not None:
            alpha = min(alpha, root)
    if alpha < 1.0:
        alpha = float(np.nextafter(alpha, 0.0))
    if not 0.0 < alpha <= 1.0:
        raise NumericalFailure("Positivity limiter could not find an admissible scaling")
    limited = average + alpha * (raw_states - average)
    # At a quadratic root, one nextafter may still round the evaluated state
    # onto the boundary. Move only by representable floats toward the known
    # admissible cell average; this is not a thermodynamic floor or insertion.
    for _ in range(64):
        try:
            primitive(limited, gamma=5.0 / 3.0)
            break
        except NumericalFailure:
            alpha = float(np.nextafter(alpha, 0.0))
            limited = average + alpha * (raw_states - average)
    else:
        raise NumericalFailure("Positivity scaling failed without an admissible repair")
    return limited, alpha


@dataclass(frozen=True)
class GravityWorkReport:
    delta_gas_energy: float
    delta_gravitational_energy: float
    conservation_residual: float
    local_work_density: np.ndarray
    integrated_mass_flux: np.ndarray


class HydroOperator:
    """Euler--Poisson operator with compatible W_h work and positivity scaling."""

    def __init__(self, grid: SphericalGrid, equilibrium: HydrostaticNFW, gamma: float) -> None:
        self.grid = grid
        self.equilibrium = equilibrium
        self.gamma = gamma
        self.equilibrium_state = equilibrium.conserved(gamma)
        self._equilibrium_acceleration = self.gravity_acceleration(equilibrium.density)
        self.limiter_activations = 0
        self.limiter_cell_reconstructions = 0
        self.last_gravity_work: GravityWorkReport | None = None

    @property
    def limiter_activation_frequency(self) -> float:
        if self.limiter_cell_reconstructions == 0:
            return 0.0
        return self.limiter_activations / self.limiter_cell_reconstructions

    def gravity_acceleration(self, density: np.ndarray) -> np.ndarray:
        """Gradient of the conjugate potential of the same discrete W_h."""
        from .diagnostics import conjugate_gravitational_potential

        psi = conjugate_gravitational_potential(density * self.grid.volumes, self.grid)
        return np.gradient(psi, self.grid.volume_centroids, edge_order=2)

    def _equilibrium_face_state(self, face_indices: np.ndarray) -> np.ndarray:
        radius = self.grid.faces[face_indices]
        result = np.zeros((face_indices.size, 3), dtype=np.float64)
        result[:, 0] = nfw_density(radius)
        result[:, 2] = self.equilibrium.face_pressure[face_indices] / (self.gamma - 1.0)
        return result

    def reconstruct_conserved(self, state: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Reconstruct equilibrium plus a conservative perturbation polynomial."""
        primitive(state, self.gamma)
        perturbation = state - self.equilibrium_state
        slopes = np.column_stack([monotonized_central(perturbation[:, j]) for j in range(3)])
        centroids = self.grid.volume_centroids
        count = state.shape[0]
        raw_left = np.empty_like(state)
        raw_right = np.empty_like(state)
        # The zero-area origin is not evaluated for the weak NFW cusp.
        raw_left[0] = state[0]
        if count > 1:
            indices = np.arange(1, count)
            eq_left = self._equilibrium_face_state(indices)
            raw_left[1:] = eq_left + perturbation[1:] + slopes[1:] * (
                (self.grid.faces[1:-1] - centroids[1:]) / self.grid.dr
            )[:, None]
        indices = np.arange(1, count + 1)
        eq_right = self._equilibrium_face_state(indices)
        raw_right = eq_right + perturbation + slopes * (
            (self.grid.faces[1:] - centroids) / self.grid.dr
        )[:, None]

        limited_left = raw_left.copy()
        limited_right = raw_right.copy()
        alpha = np.ones(count, dtype=np.float64)
        for i in range(count):
            raw = raw_right[i : i + 1] if i == 0 else np.vstack((raw_left[i], raw_right[i]))
            limited, alpha[i] = conservative_positivity_scale(state[i], raw)
            if i == 0:
                limited_right[i] = limited[0]
            else:
                limited_left[i], limited_right[i] = limited
        self.limiter_cell_reconstructions += count
        self.limiter_activations += int(np.count_nonzero(alpha < 1.0))
        return limited_left, limited_right, alpha

    def _face_fluxes(self, state: np.ndarray) -> np.ndarray:
        left_cell, right_cell, _ = self.reconstruct_conserved(state)
        count = state.shape[0]
        flux = np.zeros((count + 1, 3), dtype=np.float64)
        equilibrium_flux = np.zeros_like(flux)
        if count > 1:
            left = _primitive_rows(right_cell[:-1], self.gamma)
            right = _primitive_rows(left_cell[1:], self.gamma)
            flux[1:-1] = hllc_flux(left, right, self.gamma)
            equilibrium_flux[1:-1, 1] = self.equilibrium.face_pressure[1:-1]
        outer = _primitive_rows(right_cell[-1:], self.gamma)[0]
        flux[-1, 1] = outer[2]
        equilibrium_flux[-1, 1] = self.equilibrium.face_pressure[-1]
        return flux - equilibrium_flux

    def spatial_operator(self, state: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Return the non-gravitational-energy RHS and integrated face fluxes."""
        rho, _, pressure = primitive(state, self.gamma)
        residual_flux = self._face_fluxes(state)
        integrated_flux = self.grid.areas[:, None] * residual_flux
        result = -(integrated_flux[1:] - integrated_flux[:-1]) / self.grid.volumes[:, None]
        acceleration = self.gravity_acceleration(rho)
        pressure_delta = pressure - self.equilibrium.pressure
        geometric_integral = 4.0 * np.pi * pressure_delta * (
            self.grid.faces[1:] ** 2 - self.grid.faces[:-1] ** 2
        )
        gravity_delta = rho * acceleration - self.equilibrium.density * self._equilibrium_acceleration
        result[:, 1] += geometric_integral / self.grid.volumes - gravity_delta
        return result, integrated_flux

    def rhs(self, state: np.ndarray) -> np.ndarray:
        """Semidiscrete diagnostic RHS with the continuum-form energy source."""
        result, _ = self.spatial_operator(state)
        rho, velocity, _ = primitive(state, self.gamma)
        result[:, 2] -= rho * velocity * self.gravity_acceleration(rho)
        return result

    def _gravity_work(
        self,
        mass_before: np.ndarray,
        mass_after: np.ndarray,
        integrated_mass_flux: np.ndarray,
        dt: float,
    ) -> tuple[np.ndarray, GravityWorkReport]:
        from .diagnostics import conjugate_gravitational_potential, gravitational_energy_from_mass

        midpoint_mass = 0.5 * (mass_before + mass_after)
        psi = conjugate_gravitational_potential(midpoint_mass, self.grid)
        face_work = -integrated_mass_flux[1:-1] * np.diff(psi)
        cell_work = np.zeros_like(mass_before)
        cell_work[:-1] += 0.5 * face_work
        cell_work[1:] += 0.5 * face_work
        delta_w = gravitational_energy_from_mass(mass_after, self.grid) - gravitational_energy_from_mass(
            mass_before, self.grid
        )
        delta_gas = float(np.sum(cell_work))
        residual = delta_gas + delta_w
        scale = max(abs(delta_gas), abs(delta_w), 1.0)
        if abs(residual) / scale > 5.0e-13:
            raise NumericalFailure("Discrete gravity-work conservation identity failed")
        work_density = cell_work / (self.grid.volumes * dt)
        return cell_work, GravityWorkReport(
            delta_gas,
            delta_w,
            residual,
            work_density,
            integrated_mass_flux.copy(),
        )

    def maximum_timestep(self, state: np.ndarray, cfl: float) -> float:
        rho, velocity, pressure = primitive(state, self.gamma)
        signal_speed = np.abs(velocity) + np.sqrt(self.gamma * pressure / rho)
        maximum = float(np.max(signal_speed))
        if maximum <= 0.0 or not np.isfinite(maximum):
            raise NumericalFailure("Invalid Euler signal speed")
        return cfl * self.grid.dr / maximum

    def step_rk2(self, state: np.ndarray, dt: float) -> np.ndarray:
        """SSP-RK2 using its integrated mass flux for exact discrete gravity work."""
        if dt <= 0.0 or not np.isfinite(dt):
            raise ValueError("A finite positive timestep is required")
        primitive(state, self.gamma)
        mass_before = state[:, 0] * self.grid.volumes

        rhs_first, flux_first = self.spatial_operator(state)
        stage = state + dt * rhs_first
        mass_stage = stage[:, 0] * self.grid.volumes
        stage_mass_flux = dt * flux_first[:, 0]
        stage_work, _ = self._gravity_work(mass_before, mass_stage, stage_mass_flux, dt)
        stage[:, 2] += stage_work / self.grid.volumes
        primitive(stage, self.gamma)

        rhs_second, flux_second = self.spatial_operator(stage)
        updated = state + 0.5 * dt * (rhs_first + rhs_second)
        integrated_mass_flux = 0.5 * dt * (flux_first[:, 0] + flux_second[:, 0])
        mass_after = updated[:, 0] * self.grid.volumes
        gravity_work, report = self._gravity_work(
            mass_before, mass_after, integrated_mass_flux, dt
        )
        updated[:, 2] += gravity_work / self.grid.volumes
        primitive(updated, self.gamma)
        self.last_gravity_work = report
        return updated


__all__ = [
    "GravityWorkReport",
    "HydroOperator",
    "NumericalFailure",
    "center_acceleration",
    "conservative_positivity_scale",
    "enclosed_mass_faces",
    "hllc_flux",
    "monotonized_central",
    "primitive",
]
