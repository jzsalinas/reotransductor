"""
Unit tests for the SPARC 2020 Halo Radial Profile & Cusp-Core Problem Pipeline.
Tests analytical NFW and Burkert models, 3D radial profile extraction,
logarithmic density slopes, enclosed mass, and circular velocity.
"""

import unittest
import numpy as np
from observational.halo_data import SPARCHaloData
from observational.halo_analyzer import HaloRadialProfileAnalyzer


class TestHaloPipeline(unittest.TestCase):

    def setUp(self):
        self.sparc = SPARCHaloData()
        self.analyzer = HaloRadialProfileAnalyzer(grid_size=32, box_size_mpc=100.0, n_shells=16)

    def test_sparc_dataset_loading(self):
        """Verify that SPARC galaxy database loads correctly with positive radii and velocities."""
        gals = self.sparc.list_galaxies()
        self.assertIn("DDO_154", gals)
        ddo = self.sparc.get_galaxy("DDO_154")
        self.assertGreater(len(ddo["r_kpc"]), 5)
        self.assertTrue((ddo["r_kpc"] > 0.0).all())
        self.assertTrue((ddo["v_obs_kms"] > 0.0).all())

    def test_nfw_and_burkert_analytical_slopes(self):
        """Verify analytical inner slopes: NFW -> -1.0 (Cusp), Burkert -> 0.0 (Flat Core)."""
        r_inner = np.array([1e-3, 1e-2])
        gamma_nfw = self.sparc.nfw_log_slope(r_inner, r_s=5.0)
        gamma_burkert = self.sparc.burkert_log_slope(r_inner, r_0=5.0)

        # NFW inner slope must approach -1.0
        self.assertAlmostEqual(float(gamma_nfw[0]), -1.0, places=2)
        # Burkert inner slope must approach 0.0
        self.assertAlmostEqual(float(gamma_burkert[0]), 0.0, places=2)

    def test_logarithmic_slope_computation(self):
        """Verify numerical slope gamma(r) on a pure synthetic power law rho(r) = r^(-2.0)."""
        r_arr = np.linspace(1.0, 50.0, 50)
        rho_arr = 100.0 * (r_arr**(-2.0))
        gamma_calc = self.analyzer.compute_logarithmic_slope(r_arr, rho_arr)

        # Interior slope points should equal -2.0
        self.assertAlmostEqual(float(np.mean(gamma_calc[5:-5])), -2.0, places=1)

    def test_halo_center_detection(self):
        """Verify that locate_halo_center finds the exact coordinates of an inserted core."""
        grid = 32
        rho = np.ones((grid, grid, grid), dtype=np.float32)
        target_center = (12, 18, 24)
        rho[target_center] = 50.0

        detected_center = self.analyzer.locate_halo_center(rho)
        self.assertEqual(detected_center, target_center)

    def test_halo_diagnostics_execution(self):
        """Verify end-to-end halo diagnostics on a synthetic 3D cored halo."""
        grid = 32
        cx, cy, cz = 16, 16, 16
        idx = np.arange(grid)
        dx = (idx[:, None, None] - cx + grid // 2) % grid - grid // 2
        dy = (idx[None, :, None] - cy + grid // 2) % grid - grid // 2
        dz = (idx[None, None, :] - cz + grid // 2) % grid - grid // 2
        r = np.sqrt(dx**2 + dy**2 + dz**2) * (100.0 / grid)

        # Synthetic cored Burkert profile
        r_core = 12.0
        rho = (10.0 / ((1.0 + r / r_core) * (1.0 + (r / r_core)**2))).astype(np.float32)

        res = self.analyzer.evaluate_halo_diagnostics(rho)
        self.assertIn("inner_slope_gamma0", res)
        self.assertIn("is_cored", res)
        self.assertIn("v_circular", res)
        self.assertTrue(res["is_cored"])
        self.assertGreater(res["inner_slope_gamma0"], -0.5)


if __name__ == '__main__':
    unittest.main()
