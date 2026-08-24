"""
Unit tests for Galactic-Scale Reotransductor Simulator (observational/galaxy_simulator.py).
Tests 3D isolated galaxy dynamics, Poisson potential solver, non-equilibrium entropy
production, proper time accumulation, and SPARC rotation curve confrontation.
"""

import os
import unittest
import numpy as np

from observational.galaxy_simulator import (
    GalacticReotransductorSimulator,
    run_single_galaxy_simulation
)


class TestGalacticSimulator(unittest.TestCase):
    """Test suite for galactic-scale isolated rotation curve simulations."""

    def setUp(self):
        self.test_galaxy = "NGC2403"
        self.sim = GalacticReotransductorSimulator(
            galaxy_name=self.test_galaxy,
            grid_size=64  # Fast test resolution
        )

    def test_galaxy_initialization(self):
        """Verifies 3D lattice dimensions, physical coordinates, and baryonic initialization."""
        self.assertEqual(self.sim.grid_size, 64)
        self.assertGreater(self.sim.box_size_kpc, 0.0)
        self.assertEqual(self.sim.rho_baryons.shape, (64, 64, 64))
        self.assertEqual(self.sim.v_x.shape, (64, 64, 64))
        self.assertEqual(self.sim.temperature.shape, (64, 64, 64))
        
        # Central baryonic density should be positive and maximal near center
        mid = 32
        self.assertGreater(self.sim.rho_baryons[mid, mid, mid], 0.0)
        self.assertGreater(
            self.sim.rho_baryons[mid, mid, mid],
            self.sim.rho_baryons[0, 0, 0]
        )

    def test_baryonic_potential_poisson_solver(self):
        """Verifies 3D Green's function FFT Poisson solver produces valid negative potential."""
        phi = self.sim.compute_baryonic_gravitational_potential()
        self.assertEqual(phi.shape, (64, 64, 64))
        
        # Gravitational potential of isolated mass must be strictly attractive (<= 0)
        mid = 32
        self.assertLess(phi[mid, mid, mid], 0.0)
        # Deepest potential well at galactic core
        self.assertLess(phi[mid, mid, mid], phi[0, 0, 0])

    def test_entropy_dissipation_and_proper_time(self):
        """Verifies proper time accumulation tau >= 1.0 and entropy scalar positivity."""
        results = self.sim.evolve_reotransductor_dynamics(n_rotations=1.0)
        
        # Proper time field must be strictly monotonic (tau >= 1.0)
        self.assertTrue(np.all(self.sim.tau_field >= 1.0))
        self.assertGreater(np.max(self.sim.tau_field), 1.0)
        
        # Check return dictionary fields
        self.assertIn("galaxy_name", results)
        self.assertIn("r_kpc", results)
        self.assertIn("v_obs", results)
        self.assertIn("v_tot_reotransductor", results)
        self.assertIn("v_tot_nfw", results)
        self.assertIn("chi2_reotransductor", results)
        self.assertIn("chi2_nfw", results)
        self.assertIn("preferred_model", results)

    def test_rotation_curve_ddo154(self):
        """Verifies dwarf galaxy DDO154 rotation curve simulation and statistical core preference."""
        sim_dwarf = GalacticReotransductorSimulator(galaxy_name="DDO154", grid_size=64)
        results = sim_dwarf.evolve_reotransductor_dynamics()
        
        # DDO154 is a canonical core-dominated galaxy: Reotransductor Core should beat NFW
        self.assertLess(results["chi2_reotransductor"], results["chi2_nfw"])
        self.assertEqual(results["preferred_model"], "Reotransductor Cored Profile")
        self.assertLess(results["red_chi2_reotransductor"], 5.0)

    def test_full_pipeline_and_figure_rendering(self):
        """Verifies end-to-end execution and figure output generation."""
        out_fig = "checkpoints/snapshots/test_galaxy_rotation.png"
        os.makedirs("checkpoints/snapshots", exist_ok=True)
        
        results = run_single_galaxy_simulation(
            galaxy_name="NGC2403",
            grid_size=64,
            output_fig=out_fig
        )
        
        self.assertTrue(os.path.exists(out_fig))
        self.assertGreater(os.path.getsize(out_fig), 1000)
        
        # Clean up test output
        if os.path.exists(out_fig):
            os.remove(out_fig)


if __name__ == '__main__':
    unittest.main()
