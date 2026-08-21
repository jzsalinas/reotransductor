"""
Unit tests for the NANOGrav 15-Year Pulsar Timing & Proper Time Micro-Drift Pipeline.
Tests Hellings-Downs analytical limits, characteristic strain spectrum,
line-of-sight delay integration, and celestial cross-correlations.
"""

import unittest
import numpy as np
from observational.nanograv_data import NANOGravPulsarData
from observational.pulsar_analyzer import PulsarTimingAnalyzer


class TestNANOGravPipeline(unittest.TestCase):

    def setUp(self):
        self.nanograv = NANOGravPulsarData()
        self.analyzer = PulsarTimingAnalyzer(grid_size=16, box_size_mpc=100.0, n_pulsars=20, n_bins=10)

    def test_nanograv_dataset_loading(self):
        """Verify that NANOGrav 15-Year dataset loads with valid angular separation bins."""
        data = self.nanograv.get_dataset()
        self.assertGreater(len(data["zeta_deg"]), 5)
        self.assertTrue((data["zeta_deg"] >= 0.0).all())
        self.assertTrue((data["zeta_deg"] <= 180.0).all())
        self.assertTrue((data["err_gamma"] > 0.0).all())

    def test_hellings_downs_analytical_limits(self):
        """Verify fundamental Hellings-Downs (1983) analytical properties."""
        # 1. Autocorrelation limit at zeta -> 0
        gamma_0 = self.nanograv.hellings_downs(np.array([0.0]))
        self.assertAlmostEqual(float(gamma_0[0]), 0.50, places=2)

        # 2. Antipodal limit at zeta = 180 deg
        gamma_180 = self.nanograv.hellings_downs(np.array([180.0]))
        self.assertAlmostEqual(float(gamma_180[0]), 0.25, places=2)

        # 3. Quadrupolar minimum near zeta = 90 deg
        gamma_90 = self.nanograv.hellings_downs(np.array([90.0]))
        self.assertAlmostEqual(float(gamma_90[0]), -0.15, delta=0.03)

    def test_characteristic_strain_scaling(self):
        """Verify characteristic strain h_c(f) decreases with frequency as f^(-2/3)."""
        f1 = 1e-9  # 1 nHz
        f2 = 1e-7  # 100 nHz
        hc1 = self.nanograv.characteristic_strain(np.array([f1]))
        hc2 = self.nanograv.characteristic_strain(np.array([f2]))

        self.assertGreater(float(hc1[0]), float(hc2[0]))
        # Ratio (100)^(-2/3) ~ 0.0464
        ratio = float(hc2[0] / hc1[0])
        self.assertAlmostEqual(ratio, 100.0**(-2.0 / 3.0), places=2)

    def test_pulsar_network_generation(self):
        """Verify pulsar generation produces unit direction vectors and valid positions."""
        pulsars = self.analyzer.generate_galactic_pulsar_network(center=(8, 8, 8), seed=123)
        self.assertEqual(len(pulsars), 20)

        for p in pulsars:
            norm = np.linalg.norm(p["n_vec"])
            self.assertAlmostEqual(norm, 1.0, places=4)
            self.assertGreater(p["dist"], 0.0)

    def test_line_of_sight_and_cross_correlation(self):
        """Verify line-of-sight delay integration and angular cross-correlation calculation."""
        grid = 16
        center = (8, 8, 8)
        tau = np.ones((grid, grid, grid), dtype=np.float32)
        # Quadrupole-like perturbation
        idx = np.arange(grid)
        dx = idx[:, None, None] - 8
        dy = idx[None, :, None] - 8
        dz = idx[None, None, :] - 8
        r = np.sqrt(dx**2 + dy**2 + dz**2)
        tau += (dx**2 - dy**2) / (r**2 + 1.0)

        pulsars = self.analyzer.generate_galactic_pulsar_network(center=center, seed=42)
        delays = self.analyzer.compute_transverse_traceless_timing_response(tau, pulsars, center=center, seed=42)
        self.assertEqual(len(delays), 20)

        zeta_arr, gamma_sim, counts = self.analyzer.compute_angular_cross_correlation(pulsars, delays)
        self.assertEqual(len(zeta_arr), 10)
        self.assertEqual(len(gamma_sim), 10)
        self.assertTrue(np.all(np.isfinite(gamma_sim)))

    def test_chi2_calculation(self):
        """Verify Chi-squared calculation gives 0.0 on identical models."""
        gamma_obs = np.array([0.4, 0.1, -0.1, 0.2])
        err_gamma = np.array([0.05, 0.05, 0.05, 0.05])
        chi2 = self.nanograv.compute_chi2(gamma_obs, gamma_obs, err_gamma)
        self.assertAlmostEqual(chi2, 0.0, places=4)


if __name__ == '__main__':
    unittest.main()
