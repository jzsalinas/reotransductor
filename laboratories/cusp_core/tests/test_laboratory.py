"""Required verification suite; no observational or project-physics inputs."""

from __future__ import annotations

import ast
import unittest
from pathlib import Path

import numpy as np

from laboratories.cusp_core.conduction import (
    conductive_luminosity,
    implicit_conduction_step,
)
from laboratories.cusp_core.config import LaboratoryConfig, PhysicalScales
from laboratories.cusp_core.diagnostics import (
    core_radius,
    energy_budget,
    gravitational_energy,
    potential_faces,
    rotation_curve,
)
from laboratories.cusp_core.equilibrium import (
    HydrostaticNFW,
    nfw_acceleration,
    nfw_density,
    nfw_mass,
)
from laboratories.cusp_core.grid import SphericalGrid
from laboratories.cusp_core.hydrodynamics import (
    HydroOperator,
    center_acceleration,
    enclosed_mass_faces,
)


class LaboratoryTests(unittest.TestCase):
    def setUp(self) -> None:
        self.gamma = 5.0 / 3.0
        self.grid = SphericalGrid.uniform(128, 10.0)
        self.equilibrium = HydrostaticNFW.build(self.grid)
        self.state = self.equilibrium.conserved(self.gamma)

    def test_exact_nfw_enclosed_mass_formula(self) -> None:
        radii = np.array([1.0e-5, 0.1, 1.0, 10.0])
        expected = 4.0 * np.pi * (np.log1p(radii) - radii / (1.0 + radii))
        np.testing.assert_allclose(nfw_mass(radii), expected, rtol=2.0e-11, atol=1.0e-20)
        derivative = 4.0 * np.pi * radii**2 * nfw_density(radii)
        step = 1.0e-6 * radii
        numerical = (nfw_mass(radii + step) - nfw_mass(radii - step)) / (2.0 * step)
        np.testing.assert_allclose(numerical, derivative, rtol=2.0e-7)

    def test_finite_volume_average_of_central_cusp(self) -> None:
        grid = SphericalGrid.uniform(100_000, 1.0)
        equilibrium = HydrostaticNFW.build(grid)
        dx = grid.dr
        exact_average = nfw_mass(dx) / ((4.0 * np.pi / 3.0) * dx**3)
        self.assertEqual(equilibrium.density[0], exact_average)
        self.assertAlmostEqual(equilibrium.density[0] * dx, 1.5, delta=3.0e-5)
        with self.assertRaises(ValueError):
            nfw_density(0.0)

    def test_poisson_enclosed_mass_and_acceleration(self) -> None:
        mass = enclosed_mass_faces(self.equilibrium.density, self.grid)
        np.testing.assert_allclose(
            mass, self.equilibrium.enclosed_mass_faces, rtol=3.0e-15, atol=2.0e-14
        )
        acceleration = nfw_acceleration(self.grid.faces[1:])
        np.testing.assert_allclose(
            mass[1:] / self.grid.faces[1:] ** 2, acceleration, rtol=3.0e-15
        )
        self.assertTrue(np.all(center_acceleration(self.equilibrium.density, self.grid) > 0.0))
        potential = potential_faces(self.equilibrium.density, self.grid)
        self.assertTrue(np.all(np.diff(potential) > 0.0))
        self.assertAlmostEqual(
            potential[-1], -mass[-1] / self.grid.faces[-1], places=14
        )

    def test_hydrostatic_initial_pressure_and_outer_condition(self) -> None:
        face_pressure = self.equilibrium.face_pressure
        np.testing.assert_allclose(
            face_pressure[1:-1] - face_pressure[2:],
            self.equilibrium.hydrostatic_drop[1:],
            rtol=2.0e-14,
            atol=2.0e-14,
        )
        radius = self.grid.faces[-1]
        residual = -nfw_acceleration(radius) - self.equilibrium.outer_theta * (
            -1.0 / radius - 2.0 / (1.0 + radius)
        )
        self.assertAlmostEqual(residual, 0.0, places=14)
        operator = HydroOperator(self.grid, self.equilibrium, self.gamma)
        self.assertLess(np.max(np.abs(operator.rhs(self.state))), 2.0e-13)

    def test_central_zero_luminosity_boundary(self) -> None:
        state = self.state.copy()
        state[:, 2] *= np.linspace(1.0, 2.0, state.shape[0])
        luminosity = conductive_luminosity(state, self.grid, self.gamma, 1.5)
        self.assertEqual(luminosity[0], 0.0)
        self.assertNotEqual(luminosity[1], 0.0)

    def test_constant_theta_has_zero_conductive_flux(self) -> None:
        state = self.state.copy()
        rho = state[:, 0]
        theta = 2.0
        state[:, 2] = rho * theta / (self.gamma - 1.0)
        luminosity = conductive_luminosity(state, self.grid, self.gamma, 1.5)
        np.testing.assert_allclose(luminosity, 0.0, rtol=0.0, atol=1.0e-11)

    def test_conductive_energy_conservation(self) -> None:
        state = self.state.copy()
        state[:, 2] *= 1.0 + 0.1 * np.sin(self.grid.centers)
        before = float(np.sum(state[:, 2] * self.grid.volumes))
        updated = implicit_conduction_step(state, self.grid, self.gamma, 1.5, 0.01)
        after = float(np.sum(updated[:, 2] * self.grid.volumes))
        self.assertLess(abs(after - before) / abs(before), 3.0e-15)

    def test_mass_conservation(self) -> None:
        operator = HydroOperator(self.grid, self.equilibrium, self.gamma)
        rhs = operator.rhs(self.state)
        mass_rate = float(np.sum(rhs[:, 0] * self.grid.volumes))
        self.assertLess(abs(mass_rate), 1.0e-13)

    def test_total_energy_accounting(self) -> None:
        uniform_grid = SphericalGrid.uniform(64, 2.0)
        density = np.ones(64)
        mass = 4.0 * np.pi * 2.0**3 / 3.0
        expected_w = -3.0 * mass**2 / (5.0 * 2.0)
        self.assertAlmostEqual(gravitational_energy(density, uniform_grid), expected_w, places=12)
        budget = energy_budget(self.state, self.grid, self.gamma)
        self.assertEqual(budget.total, budget.internal + budget.kinetic + budget.gravitational)

    def test_core_slope_estimator_on_analytic_profile(self) -> None:
        dx = 0.0025
        radius = (np.arange(4000) + 0.5) * dx
        scale = 2.0
        density = 1.0 / (1.0 + (radius / scale) ** 2)
        estimate = core_radius(radius, density, dx)
        self.assertIsNotNone(estimate)
        self.assertAlmostEqual(estimate, scale / np.sqrt(3.0), delta=2.0e-4)
        cusp = 1.0 / (radius * (1.0 + radius) ** 2)
        self.assertIsNone(core_radius(radius, cusp, dx))

    def test_rotation_curve_on_uniform_sphere(self) -> None:
        grid = SphericalGrid.uniform(64, 2.0)
        radius, velocity = rotation_curve(np.ones(64), grid)
        expected = np.sqrt(4.0 * np.pi / 3.0) * radius
        np.testing.assert_allclose(velocity, expected, rtol=2.0e-15)
        physical_radius, physical_velocity = rotation_curve(
            np.ones(64), grid, PhysicalScales()
        )
        np.testing.assert_allclose(physical_radius, radius)
        np.testing.assert_allclose(
            physical_velocity,
            velocity * PhysicalScales().velocity_scale_m_s / 1.0e3,
        )

    def test_no_reotransductor_physics_imports(self) -> None:
        package = Path(__file__).resolve().parents[1]
        forbidden = {"server", "observational", "experiments"}
        for source_path in package.rglob("*.py"):
            tree = ast.parse(source_path.read_text(encoding="utf-8"), filename=str(source_path))
            for node in ast.walk(tree):
                names: list[str] = []
                if isinstance(node, ast.Import):
                    names = [alias.name for alias in node.names]
                elif isinstance(node, ast.ImportFrom) and node.module:
                    names = [node.module]
                for name in names:
                    self.assertNotIn(name.split(".")[0], forbidden, source_path)

    def test_frozen_conductivity_and_timescale(self) -> None:
        config = LaboratoryConfig()
        self.assertEqual(config.conductivity_hat_star, 1.5)
        self.assertAlmostEqual(config.reference_conduction_time(1.5), 1.0)
        self.assertAlmostEqual(config.reference_conduction_time(1.5) / 4.0, 0.25)


if __name__ == "__main__":
    unittest.main()
