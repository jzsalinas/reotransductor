"""Certified matrix-exponential Fourier conduction in spherical volumes."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import ExponentialActionContract
from .grid import SphericalGrid
from .hydrodynamics import NumericalFailure, primitive


@dataclass(frozen=True)
class ExponentialActionReport:
    """A-posteriori numerical certification for one exponential action."""

    krylov_dimension: int
    action_error_estimate: float
    residual_estimate: float
    weighted_invariant_error: float
    constant_mode_error: float
    minimum_principle_violation: float
    maximum_principle_violation: float
    minimum_temperature: float


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
    # The origin condition is the finite-volume statement L_cond(0)=0.
    # The outer wall is insulating.
    return luminosity


def conduction_operators(
    density: np.ndarray,
    grid: SphericalGrid,
    gamma: float,
    conductivity_hat: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return H and the diagonal/off-diagonal representation of symmetric B."""
    density = np.asarray(density, dtype=np.float64)
    if density.shape != grid.centers.shape or np.any(density <= 0.0):
        raise NumericalFailure("Conduction received a non-positive cell density")
    if conductivity_hat < 0.0:
        raise ValueError("Conductivity cannot be negative")
    capacity = grid.volumes * density / (gamma - 1.0)
    conductance = conductivity_hat * grid.areas[1:-1] / np.diff(grid.centers)
    diagonal = np.zeros_like(capacity)
    diagonal[:-1] -= conductance / capacity[:-1]
    diagonal[1:] -= conductance / capacity[1:]
    off_diagonal = conductance / np.sqrt(capacity[:-1] * capacity[1:])
    return capacity, diagonal, off_diagonal


def _symmetric_tridiagonal_matvec(
    diagonal: np.ndarray, off_diagonal: np.ndarray, vector: np.ndarray
) -> np.ndarray:
    result = diagonal * vector
    result[:-1] += off_diagonal * vector[1:]
    result[1:] += off_diagonal * vector[:-1]
    return result


def _projected_exponential(
    alpha: np.ndarray, beta: np.ndarray, time: float
) -> tuple[np.ndarray, float]:
    """Compute exp(time*T)e1 for the small symmetric Lanczos projection."""
    dimension = alpha.size
    projected = np.diag(alpha)
    if dimension > 1:
        projected += np.diag(beta[: dimension - 1], 1)
        projected += np.diag(beta[: dimension - 1], -1)
    eigenvalues, eigenvectors = np.linalg.eigh(projected)
    coefficients = eigenvectors @ (np.exp(time * eigenvalues) * eigenvectors[0, :])
    return coefficients, float(coefficients[-1])


def _lanczos_exponential_action(
    diagonal: np.ndarray,
    off_diagonal: np.ndarray,
    vector: np.ndarray,
    time: float,
    contract: ExponentialActionContract,
) -> tuple[np.ndarray, int, float, float]:
    """Adaptive Lanczos exp(B*time) action with an a-posteriori contract."""
    norm = float(np.linalg.norm(vector))
    if norm == 0.0:
        return vector.copy(), 0, 0.0, 0.0

    size = vector.size
    maximum = min(size, contract.maximum_krylov_dimension)
    basis = np.empty((size, maximum), dtype=np.float64)
    alpha = np.empty(maximum, dtype=np.float64)
    beta = np.empty(max(maximum - 1, 0), dtype=np.float64)
    basis[:, 0] = vector / norm
    previous_vector = np.zeros(size, dtype=np.float64)
    previous_beta = 0.0
    previous_action: np.ndarray | None = None
    action = vector.copy()
    action_error = float("inf")
    residual = float("inf")
    breakdown = False

    for j in range(maximum):
        work = _symmetric_tridiagonal_matvec(diagonal, off_diagonal, basis[:, j])
        if j:
            work -= previous_beta * previous_vector
        alpha[j] = np.dot(basis[:, j], work)
        work -= alpha[j] * basis[:, j]
        # Full reorthogonalization protects the explicitly split null mode.
        projection = basis[:, : j + 1].T @ work
        work -= basis[:, : j + 1] @ projection
        next_beta = float(np.linalg.norm(work))
        if j < maximum - 1:
            beta[j] = next_beta

        dimension = j + 1
        checkpoint = (
            dimension >= contract.initial_krylov_dimension
            and dimension % contract.krylov_dimension_increment == 0
        ) or next_beta <= np.finfo(np.float64).eps * max(1.0, abs(alpha[j]))
        checkpoint = checkpoint or dimension == maximum
        if checkpoint:
            coefficients, last_coefficient = _projected_exponential(
                alpha[:dimension], beta[: max(dimension - 1, 0)], time
            )
            action = norm * (basis[:, :dimension] @ coefficients)
            residual = abs(time * next_beta * norm * last_coefficient)
            if previous_action is not None:
                action_error = float(np.linalg.norm(action - previous_action))
            scale = max(float(np.linalg.norm(action)), norm)
            tolerance = contract.absolute_action_tolerance + contract.relative_action_tolerance * scale
            breakdown = next_beta <= np.finfo(np.float64).eps * max(1.0, abs(alpha[j]))
            if breakdown and residual <= tolerance:
                # An invariant Krylov subspace is exact in exact arithmetic;
                # the smaller preceding projection is not its error estimate.
                action_error = 0.0
                return action, dimension, action_error, residual
            if previous_action is not None and action_error <= tolerance and residual <= tolerance:
                return action, dimension, action_error, residual
            previous_action = action.copy()

        if dimension == maximum:
            break
        previous_vector = basis[:, j].copy()
        previous_beta = next_beta
        if next_beta == 0.0:
            break
        basis[:, j + 1] = work / next_beta

    scale = max(float(np.linalg.norm(action)), norm)
    tolerance = contract.absolute_action_tolerance + contract.relative_action_tolerance * scale
    if not breakdown and (action_error > tolerance or residual > tolerance):
        raise NumericalFailure(
            "Matrix-exponential action did not satisfy the pre-specified "
            f"Krylov contract (m={maximum}, error={action_error:.3e}, "
            f"residual={residual:.3e}, tolerance={tolerance:.3e})"
        )
    return action, maximum, action_error, residual


def exponential_conduction_step(
    state: np.ndarray,
    grid: SphericalGrid,
    gamma: float,
    conductivity_hat: float,
    dt: float,
    contract: ExponentialActionContract | None = None,
    *,
    return_report: bool = False,
) -> np.ndarray | tuple[np.ndarray, ExponentialActionReport]:
    """Apply exp(dt H^-1 A) without clipping or thermodynamic repair."""
    contract = contract or ExponentialActionContract()
    if conductivity_hat < 0.0:
        raise ValueError("Conductivity cannot be negative")
    if dt < 0.0:
        raise ValueError("Timestep cannot be negative")
    rho, _, pressure = primitive(state, gamma)
    theta = pressure / rho
    if conductivity_hat == 0.0 or dt == 0.0:
        report = ExponentialActionReport(0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, float(np.min(theta)))
        result = state.copy()
        return (result, report) if return_report else result

    capacity, diagonal, off_diagonal = conduction_operators(rho, grid, gamma, conductivity_hat)
    sqrt_capacity = np.sqrt(capacity)
    invariant_before = float(np.dot(capacity, theta))
    theta_constant = invariant_before / float(np.sum(capacity))
    constant_vector = sqrt_capacity * theta_constant
    deviation = sqrt_capacity * theta - constant_vector
    deviation -= sqrt_capacity * (
        float(np.dot(sqrt_capacity, deviation)) / float(np.sum(capacity))
    )
    evolved, dimension, action_error, residual = _lanczos_exponential_action(
        diagonal, off_diagonal, deviation, dt, contract
    )
    evolved -= sqrt_capacity * (
        float(np.dot(sqrt_capacity, evolved)) / float(np.sum(capacity))
    )
    theta_new = (constant_vector + evolved) / sqrt_capacity

    invariant_after = float(np.dot(capacity, theta_new))
    invariant_scale = max(abs(invariant_before), np.finfo(np.float64).tiny)
    invariant_error = abs(invariant_after - invariant_before) / invariant_scale
    # The weighted constant mode is split before Lanczos and reinserted
    # algebraically. Certify the action, not cancellation in B@sqrt(H).
    restored_constant = constant_vector / sqrt_capacity
    constant_error = float(np.max(np.abs(restored_constant - theta_constant))) / max(
        abs(theta_constant), 1.0
    )
    theta_scale = max(float(np.max(np.abs(theta))), 1.0)
    minimum_violation = max(0.0, float(np.min(theta) - np.min(theta_new)))
    maximum_violation = max(0.0, float(np.max(theta_new) - np.max(theta)))
    dmp_allowance = contract.maximum_principle_allowance * theta_scale
    positivity_allowance = contract.positivity_allowance * theta_scale

    if invariant_error > contract.weighted_invariant_tolerance:
        raise NumericalFailure("Matrix-exponential thermal invariant check failed")
    if constant_error > contract.constant_mode_tolerance:
        raise NumericalFailure("Matrix-exponential constant-mode check failed")
    if minimum_violation > dmp_allowance or maximum_violation > dmp_allowance:
        raise NumericalFailure("Matrix-exponential discrete maximum-principle check failed")
    if float(np.min(theta_new)) < -positivity_allowance:
        raise NumericalFailure("Matrix-exponential positivity certification failed")
    if np.any(theta_new <= 0.0) or not np.all(np.isfinite(theta_new)):
        raise NumericalFailure("Matrix-exponential action produced a non-positive temperature")

    updated = state.copy()
    kinetic = 0.5 * state[:, 1] ** 2 / rho
    updated[:, 2] = rho * theta_new / (gamma - 1.0) + kinetic
    primitive(updated, gamma)
    report = ExponentialActionReport(
        dimension,
        action_error,
        residual,
        invariant_error,
        constant_error,
        minimum_violation,
        maximum_violation,
        float(np.min(theta_new)),
    )
    return (updated, report) if return_report else updated


# Compatibility name: this now executes CCL-NUM-2, not Crank--Nicolson.
implicit_conduction_step = exponential_conduction_step


__all__ = [
    "ExponentialActionReport",
    "conduction_operators",
    "conductive_luminosity",
    "exponential_conduction_step",
]
