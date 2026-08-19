"""
Unit tests for the Observational Data Pipeline and Planck 2018 Validation.
Tests Spherical Harmonics decomposition, Legendre polynomial orthogonality,
Planck 2018 dataset loading, and Hubble tension predictions.
"""

import unittest
import numpy as np
from observational.planck_data import Planck2018Data
from observational.cmb_analyzer import CMBSphericalHarmonicsAnalyzer
from observational.hubble_tension import HubbleTensionAnalyzer


class TestObservationalPipeline(unittest.TestCase):

    def setUp(self):
        self.planck = Planck2018Data()
        self.analyzer = CMBSphericalHarmonicsAnalyzer(n_theta=32, n_phi=64, ell_max=8)
        self.ht_analyzer = HubbleTensionAnalyzer()

    def test_planck_data_loading_and_shape(self):
        """Verify that official Planck 2018 dataset is well-formed with positive errors."""
        ell, D_ell, err = self.planck.get_binned_spectrum()
        self.assertGreater(len(ell), 15)
        self.assertTrue((ell >= 2).all())
        self.assertTrue((D_ell > 0.0).all())
        self.assertTrue((err > 0.0).all())

    def test_spherical_harmonics_orthonormality(self):
        """Verify surface integral of |Y_{l, m}|^2 over S^2 equals 1.0 (orthonormality)."""
        y_00 = self.analyzer.compute_spherical_harmonic(0, 0)
        norm_00 = float(np.real(np.sum(np.abs(y_00)**2 * self.analyzer.d_omega)))
        self.assertAlmostEqual(norm_00, 1.0, places=2)

        y_20 = self.analyzer.compute_spherical_harmonic(2, 0)
        norm_20 = float(np.real(np.sum(np.abs(y_20)**2 * self.analyzer.d_omega)))
        self.assertAlmostEqual(norm_20, 1.0, places=2)

    def test_angular_power_spectrum_pure_quadrupole(self):
        """Verify that a synthetic pure quadrupole map Y_20 concentrates power at ell=2."""
        y_20 = self.analyzer.compute_spherical_harmonic(2, 0)
        pure_map = np.real(y_20)

        spectrum = self.analyzer.compute_angular_power_spectrum(pure_map, temperature_scale_uK=1.0)
        ell_arr = spectrum["ell"]
        c_ell_arr = spectrum["C_ell"]

        idx_ell_2 = np.where(ell_arr == 2)[0][0]
        power_at_2 = c_ell_arr[idx_ell_2]
        other_power = float(np.sum(c_ell_arr[ell_arr != 2]))

        self.assertGreater(power_at_2, 0.01)
        self.assertLess(other_power, power_at_2 * 0.1)

    def test_hubble_tension_prediction(self):
        """Verify Hubble tension calculation from proper time dilation."""
        # Simulation with cluster proper time 8.37% faster than voids
        res = self.ht_analyzer.predict_local_hubble_rate(
            tau_cluster_mean=108.37,
            tau_void_mean=100.0,
            time_elapsed=100.0
        )
        self.assertAlmostEqual(res["h0_predicted_reotransductor"], 73.0, delta=0.2)
        self.assertTrue(res["is_tension_mitigated"])

    def test_chi2_goodness_of_fit(self):
        """Verify Chi-squared calculation against Planck 2018 points."""
        ell_obs, D_obs, _ = self.planck.get_binned_spectrum()
        # Test exact match has chi2 == 0
        stats = self.planck.compute_chi2(ell_obs, D_obs)
        self.assertAlmostEqual(stats["chi2"], 0.0, places=2)


if __name__ == '__main__':
    unittest.main()
