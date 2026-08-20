"""
Unit tests for the DESI 2024 BAO Observational Pipeline.
Tests 3D Wiener-Khinchin autocorrelation, spherical radial averaging,
DESI 2024 dataset integrity, and acoustic peak detection.
"""

import unittest
import numpy as np
from observational.bao_data import DESI2024BAOData
from observational.bao_analyzer import BAOSpatialCorrelationAnalyzer


class TestBAOPipeline(unittest.TestCase):

    def setUp(self):
        self.desi = DESI2024BAOData()
        self.analyzer = BAOSpatialCorrelationAnalyzer(grid_size=32, box_size_mpc=100.0, n_bins=16)

    def test_desi_dataset_shape_and_errors(self):
        """Verify DESI 2024 dataset shape, positive errors, and radial span."""
        r, xi, err = self.desi.get_desi_dataset()
        self.assertGreater(len(r), 10)
        self.assertTrue((err > 0.0).all())
        self.assertTrue((np.diff(r) > 0.0).all())
        self.assertGreater(float(np.max(r)), 120.0)

    def test_wiener_khinchin_identity(self):
        """Verify that spatial autocorrelation at zero separation equals field variance: xi(0) == var(delta)."""
        np.random.seed(42)
        grid = 32
        # Synthetic random density field
        rho = np.random.lognormal(mean=0.0, sigma=0.5, size=(grid, grid, grid))
        delta = (rho - np.mean(rho)) / np.mean(rho)
        var_expected = float(np.var(delta))

        xi_3d = self.analyzer.compute_3d_autocorrelation(rho)
        xi_0 = float(xi_3d[0, 0, 0])

        self.assertAlmostEqual(xi_0, var_expected, places=4)

    def test_isotropic_spherical_averaging(self):
        """Verify that spherical radial averaging bins correctly and covers the box radius."""
        r_centers, xi_1d, counts = self.analyzer.spherically_average_correlation(
            np.ones((32, 32, 32))
        )
        self.assertEqual(len(r_centers), 16)
        self.assertTrue((counts > 0).all())
        self.assertAlmostEqual(float(r_centers[0]), 100.0 / (32 * 2), delta=2.0)

    def test_acoustic_peak_detection(self):
        """Verify that evaluate_cosmological_fields detects acoustic peak radius and returns valid metrics."""
        np.random.seed(123)
        grid = 32
        rho = np.ones((grid, grid, grid), dtype=np.float64)
        # Add acoustic shell around center
        cx, cy, cz = grid // 2, grid // 2, grid // 2
        idx = np.arange(grid)
        x = (idx[:, None, None] - cx) % grid
        y = (idx[None, :, None] - cy) % grid
        z = (idx[None, None, :] - cz) % grid
        r = np.sqrt(x**2 + y**2 + z**2) * (100.0 / grid)
        rho += 0.5 * np.exp(-0.5 * ((r - 30.0) / 4.0)**2)

        res = self.analyzer.evaluate_cosmological_fields(rho, scale_factor=1.0)
        self.assertIn("bao_peak_radius_mpc", res)
        self.assertIn("correlation_length_r0_mpc", res)
        self.assertGreater(res["bao_peak_radius_mpc"], 0.0)

    def test_chi2_calculation(self):
        """Verify Chi-squared calculation against DESI 2024 points."""
        r_obs, xi_obs, _ = self.desi.get_desi_dataset()
        stats = self.desi.compute_chi2(r_obs, xi_obs)
        self.assertAlmostEqual(stats["chi2"], 0.0, places=2)


if __name__ == '__main__':
    unittest.main()
