"""Pre-registered CCL-NUM-2 verification matrix (tests 1--16)."""

from __future__ import annotations

import unittest

import numpy as np

from laboratories.cusp_core.conduction import (
    conduction_operators,
    conductive_luminosity,
    exponential_conduction_step,
)
from laboratories.cusp_core.diagnostics import energy_budget
from laboratories.cusp_core.equilibrium import HydrostaticNFW
from laboratories.cusp_core.grid import SphericalGrid
from laboratories.cusp_core.hydrodynamics import (
    HydroOperator,
    NumericalFailure,
    conservative_positivity_scale,
    primitive,
)
from laboratories.cusp_core.solver import CuspCoreSolver


GAMMA = 5.0 / 3.0
CONDUCTIVITY = 1.5


def state_from_theta(grid: SphericalGrid, theta: np.ndarray, density: np.ndarray | None = None) -> np.ndarray:
    rho = np.ones(grid.centers.size) if density is None else np.asarray(density)
    return np.column_stack((rho, np.zeros_like(rho), rho * theta / (GAMMA - 1.0)))


def spherical_mode_average(grid: SphericalGrid, wave_number: float) -> np.ndarray:
    a = grid.faces[:-1]
    b = grid.faces[1:]

    def antiderivative(radius: np.ndarray) -> np.ndarray:
        return (
            -radius * np.cos(wave_number * radius) / wave_number**2
            + np.sin(wave_number * radius) / wave_number**3
        )

    return 3.0 * (antiderivative(b) - antiderivative(a)) / (b**3 - a**3)


class CCLNUM2Verification(unittest.TestCase):
    """Software verification only; all low-resolution runs are DEVELOPMENT."""

    def test_01_constant_theta(self) -> None:
        grid = SphericalGrid.uniform(96, 1.0)
        state = state_from_theta(grid, np.full(96, 2.25))
        updated, report = exponential_conduction_step(
            state, grid, GAMMA, CONDUCTIVITY, 7.0, return_report=True
        )
        np.testing.assert_allclose(updated, state, rtol=0.0, atol=2.0e-14)
        self.assertLessEqual(report.constant_mode_error, 5.0e-14)

    def test_02_positive_gaussian_diffusion(self) -> None:
        grid = SphericalGrid.uniform(128, 1.0)
        theta = 1.0 + 0.4 * np.exp(-((grid.centers - 0.45) / 0.08) ** 2)
        state = state_from_theta(grid, theta)
        updated, report = exponential_conduction_step(
            state, grid, GAMMA, CONDUCTIVITY, 2.0e-4, return_report=True
        )
        theta_new = (GAMMA - 1.0) * updated[:, 2]
        self.assertGreater(float(np.min(theta_new)), 0.0)
        self.assertLess(float(np.max(theta_new)), float(np.max(theta)))
        self.assertGreaterEqual(float(np.min(theta_new)), float(np.min(theta)) - 5.0e-13)
        self.assertEqual(report.minimum_principle_violation, 0.0)

    def test_03_spherical_analytic_diffusion_mode(self) -> None:
        # tan(kR)=kR is the insulating spherical l=0 eigencondition.
        wave_number = 4.493409457909064
        grid = SphericalGrid.uniform(128, 1.0)
        mode = spherical_mode_average(grid, wave_number)
        theta = 1.0 + 0.1 * mode
        time = 0.003
        updated = exponential_conduction_step(
            state_from_theta(grid, theta), grid, GAMMA, CONDUCTIVITY, time
        )
        theta_numerical = (GAMMA - 1.0) * updated[:, 2]
        theta_exact = 1.0 + 0.1 * mode * np.exp(
            -(GAMMA - 1.0) * CONDUCTIVITY * wave_number**2 * time
        )
        error = np.sqrt(
            np.sum((theta_numerical - theta_exact) ** 2 * grid.volumes)
            / np.sum(grid.volumes)
        )
        self.assertLess(error, 2.0e-7)

    def test_04_thermal_conservation(self) -> None:
        grid = SphericalGrid.uniform(128, 1.0)
        density = 0.7 + grid.centers**2
        theta = 1.2 + 0.2 * np.sin(2.0 * np.pi * grid.centers)
        state = state_from_theta(grid, theta, density)
        before = float(np.sum(state[:, 2] * grid.volumes))
        updated, report = exponential_conduction_step(
            state, grid, GAMMA, CONDUCTIVITY, 0.02, return_report=True
        )
        after = float(np.sum(updated[:, 2] * grid.volumes))
        self.assertLess(abs(after - before) / abs(before), 5.0e-13)
        self.assertLessEqual(report.weighted_invariant_error, 5.0e-13)

    def test_05_origin_luminosity_is_zero(self) -> None:
        grid = SphericalGrid.uniform(64, 1.0)
        state = state_from_theta(grid, 1.0 + grid.centers)
        luminosity = conductive_luminosity(state, grid, GAMMA, CONDUCTIVITY)
        self.assertEqual(luminosity[0], 0.0)
        self.assertNotEqual(luminosity[1], 0.0)

    def test_06_stiff_diffusion_mu_through_1e5(self) -> None:
        grid = SphericalGrid.uniform(32, 1.0)
        theta = 1.0 + 0.2 * np.cos(np.pi * grid.centers)
        state = state_from_theta(grid, theta)
        diffusivity = (GAMMA - 1.0) * CONDUCTIVITY
        for mu in (1.0e2, 1.0e3, 1.0e4, 1.0e5):
            dt = mu * grid.dr**2 / diffusivity
            updated, report = exponential_conduction_step(
                state, grid, GAMMA, CONDUCTIVITY, dt, return_report=True
            )
            theta_new = (GAMMA - 1.0) * updated[:, 2]
            self.assertGreater(float(np.min(theta_new)), 0.0, mu)
            self.assertLessEqual(report.weighted_invariant_error, 5.0e-13, mu)
            self.assertLessEqual(report.minimum_principle_violation, 5.0e-13, mu)
            self.assertLessEqual(report.maximum_principle_violation, 5.0e-13, mu)

    def test_07_nfw_central_weak_solution(self) -> None:
        for cells in (64, 128, 256):
            grid = SphericalGrid.uniform(cells, 10.0)
            equilibrium = HydrostaticNFW.build(grid)
            state = equilibrium.conserved(GAMMA)
            before = float(np.sum(state[:, 2] * grid.volumes))
            updated = exponential_conduction_step(
                state, grid, GAMMA, CONDUCTIVITY, 1.0e-4
            )
            rho, _, pressure = primitive(updated, GAMMA)
            self.assertTrue(np.all(pressure / rho > 0.0))
            after = float(np.sum(updated[:, 2] * grid.volumes))
            self.assertLess(abs(after - before) / abs(before), 5.0e-13)

    def test_08_temporal_convergence(self) -> None:
        def run(steps: int) -> np.ndarray:
            solver = CuspCoreSolver(64, conductivity_hat=0.0)
            radius = solver.grid.centers
            rho, _, pressure = primitive(solver.state, GAMMA)
            velocity = 0.02 * np.sin(np.pi * radius / 10.0)
            pressure = pressure * (1.0 + 0.05 * np.sin(2.0 * np.pi * radius / 10.0))
            solver.state[:, 1] = rho * velocity
            solver.state[:, 2] = pressure / (GAMMA - 1.0) + 0.5 * rho * velocity**2
            for _ in range(steps):
                solver.step(0.02 / steps)
            return solver.state

        solutions = [run(steps) for steps in (4, 8, 16, 32)]
        differences = [
            np.linalg.norm(coarse - fine)
            for coarse, fine in zip(solutions[:-1], solutions[1:])
        ]
        orders = np.log2(np.asarray(differences[:-1]) / differences[1:])
        self.assertTrue(np.all((orders >= 1.8) & (orders <= 2.2)), orders)

    def test_09_spatial_convergence(self) -> None:
        wave_number = 4.493409457909064
        errors = []
        for cells in (32, 64, 128):
            grid = SphericalGrid.uniform(cells, 1.0)
            mode = spherical_mode_average(grid, wave_number)
            time = 0.003
            state = state_from_theta(grid, 1.0 + 0.1 * mode)
            updated = exponential_conduction_step(
                state, grid, GAMMA, CONDUCTIVITY, time
            )
            exact = 1.0 + 0.1 * mode * np.exp(
                -(GAMMA - 1.0) * CONDUCTIVITY * wave_number**2 * time
            )
            numerical = (GAMMA - 1.0) * updated[:, 2]
            errors.append(
                np.sqrt(
                    np.sum((numerical - exact) ** 2 * grid.volumes)
                    / np.sum(grid.volumes)
                )
            )
        orders = np.log2(np.asarray(errors[:-1]) / errors[1:])
        self.assertTrue(np.all((orders >= 1.8) & (orders <= 2.2)), orders)

    def test_10_hydro_gravity_energy(self) -> None:
        solver = CuspCoreSolver(128, conductivity_hat=0.0)
        rho, _, pressure = primitive(solver.state, GAMMA)
        velocity = 0.02 * np.sin(np.pi * solver.grid.centers / 10.0)
        solver.state[:, 1] = rho * velocity
        solver.state[:, 2] = pressure / (GAMMA - 1.0) + 0.5 * rho * velocity**2
        before = energy_budget(solver.state, solver.grid, GAMMA).total
        solver.step(1.0e-3)
        after = energy_budget(solver.state, solver.grid, GAMMA).total
        self.assertLess(abs(after - before) / abs(before), 5.0e-13)
        self.assertLess(abs(solver.hydrodynamics.last_gravity_work.conservation_residual), 5.0e-13)

    def test_11_extreme_positive_muscl_reconstruction(self) -> None:
        average = np.array([1.0, 0.4, 2.0])
        raw = np.array([[-100.0, 30.0, -4.0], [0.2, -20.0, 0.1]])
        limited, alpha = conservative_positivity_scale(average, raw)
        self.assertGreater(alpha, 0.0)
        self.assertLess(alpha, 1.0)
        rho, _, pressure = primitive(limited, GAMMA)
        self.assertTrue(np.all(rho > 0.0))
        self.assertTrue(np.all(pressure > 0.0))
        ratios = (limited - average) / (raw - average)
        np.testing.assert_allclose(ratios[np.isfinite(ratios)], alpha, rtol=2.0e-14)

    def test_12_exact_nfw_equilibrium(self) -> None:
        grid = SphericalGrid.uniform(256, 10.0)
        equilibrium = HydrostaticNFW.build(grid)
        operator = HydroOperator(grid, equilibrium, GAMMA)
        state = equilibrium.conserved(GAMMA)
        _, _, alpha = operator.reconstruct_conserved(state)
        np.testing.assert_array_equal(alpha, np.ones_like(alpha))
        self.assertLess(float(np.max(np.abs(operator.rhs(state)))), 5.0e-12)
        updated = operator.step_rk2(state, 1.0e-3)
        np.testing.assert_allclose(updated, state, rtol=0.0, atol=5.0e-14)
        self.assertEqual(operator.limiter_activation_frequency, 0.0)

    def test_13_local_collapse_and_expansion_gravity_work(self) -> None:
        for sign in (-1.0, 1.0):
            errors = []
            for cells in (64, 128, 256):
                grid = SphericalGrid.uniform(cells, 10.0)
                equilibrium = HydrostaticNFW.build(grid)
                operator = HydroOperator(grid, equilibrium, GAMMA)
                state = equilibrium.conserved(GAMMA)
                rho, _, pressure = primitive(state, GAMMA)
                velocity = sign * 0.02 * np.sin(np.pi * grid.volume_centroids / 10.0)
                state[:, 1] = rho * velocity
                state[:, 2] = pressure / (GAMMA - 1.0) + 0.5 * rho * velocity**2
                operator.step_rk2(state, 1.0e-6)
                report = operator.last_gravity_work
                expected = -rho * velocity * operator.gravity_acceleration(rho)
                resolved = (grid.centers > 0.5) & (grid.centers < 9.5)
                errors.append(
                    np.sum(np.abs(report.local_work_density[resolved] - expected[resolved]) * grid.volumes[resolved])
                    / np.sum(np.abs(expected[resolved]) * grid.volumes[resolved])
                )
            orders = np.log2(np.asarray(errors[:-1]) / errors[1:])
            self.assertTrue(np.all(orders > 1.8), (sign, errors, orders))

    def test_14_no_silent_repair(self) -> None:
        average = np.array([1.0, 0.0, -1.0])
        with self.assertRaises(NumericalFailure):
            conservative_positivity_scale(average, np.array([[1.0, 0.0, 1.0]]))
        grid = SphericalGrid.uniform(16, 1.0)
        state = state_from_theta(grid, np.ones(16))
        state[3, 2] = -1.0
        with self.assertRaises(NumericalFailure):
            exponential_conduction_step(state, grid, GAMMA, CONDUCTIVITY, 0.1)

    def test_15_exponential_action_against_dense_reference(self) -> None:
        grid = SphericalGrid.uniform(20, 1.0)
        density = 0.8 + grid.centers
        theta = 1.0 + 0.15 * np.cos(2.0 * np.pi * grid.centers)
        state = state_from_theta(grid, theta, density)
        time = 0.013
        updated = exponential_conduction_step(state, grid, GAMMA, CONDUCTIVITY, time)
        capacity, diagonal, off = conduction_operators(density, grid, GAMMA, CONDUCTIVITY)
        dense = np.diag(diagonal) + np.diag(off, 1) + np.diag(off, -1)
        eigenvalues, eigenvectors = np.linalg.eigh(dense)
        transformed = np.sqrt(capacity) * theta
        exact_transformed = eigenvectors @ (
            np.exp(time * eigenvalues) * (eigenvectors.T @ transformed)
        )
        exact = exact_transformed / np.sqrt(capacity)
        numerical = (GAMMA - 1.0) * updated[:, 2] / density
        np.testing.assert_allclose(numerical, exact, rtol=2.0e-12, atol=2.0e-13)

    def test_16_thermal_invariant_and_constant_mode(self) -> None:
        grid = SphericalGrid.uniform(80, 1.0)
        density = 0.5 + 2.0 * grid.centers**2
        theta = 1.0 + 0.2 * np.sin(np.pi * grid.centers)
        state = state_from_theta(grid, theta, density)
        _, report = exponential_conduction_step(
            state, grid, GAMMA, CONDUCTIVITY, 0.1, return_report=True
        )
        self.assertLessEqual(report.weighted_invariant_error, 5.0e-13)
        self.assertLessEqual(report.constant_mode_error, 5.0e-14)
        self.assertLessEqual(report.action_error_estimate, 5.0e-12)
        self.assertLessEqual(report.residual_estimate, 5.0e-12)


if __name__ == "__main__":
    unittest.main()
