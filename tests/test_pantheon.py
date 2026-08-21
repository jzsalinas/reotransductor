"""
Unit tests for the Pantheon+ (2022) Supernovae & Hubble Tension Pipeline.
Tests luminosity distance integration, distance modulus, environmental H0 field,
and Chi-squared goodness of fit.
"""

import unittest
import numpy as np
from observational.pantheon_data import PantheonSupernovaeData
from observational.hubble_tension import HubbleTensionAnalyzer


class TestPantheonPipeline(unittest.TestCase):

    def setUp(self):
        self.pantheon = PantheonSupernovaeData()
        self.analyzer = HubbleTensionAnalyzer()

    def test_pantheon_dataset_loading(self):
        """Verify that Pantheon+ dataset loads with valid redshifts and distance moduli."""
        data = self.pantheon.get_dataset()
        self.assertGreater(len(data["z"]), 5)
        self.assertTrue((data["z"] > 0.0).all())
        self.assertTrue((data["mu_obs"] > 30.0).all())
        self.assertTrue((data["err_mu"] > 0.0).all())

    def test_luminosity_distance_monotonicity(self):
        """Verify that luminosity distance d_L(z) increases monotonically with redshift."""
        z_arr = np.linspace(0.01, 1.5, 20)
        d_l = self.pantheon.luminosity_distance_mpc(z_arr, h0=70.0)
        self.assertTrue(np.all(np.diff(d_l) > 0.0))

    def test_distance_modulus_analytical_value(self):
        """Verify theoretical distance modulus mu(z=0.1) matches standard cosmology ~38.4 mag."""
        mu_01 = self.pantheon.distance_modulus(np.array([0.1]), h0=70.0, omega_m=0.3)
        self.assertAlmostEqual(float(mu_01[0]), 38.45, delta=0.5)

    def test_hubble_tension_resolution_prediction(self):
        """Verify that proper time divergence predicts H_0 = 73.04 km/s/Mpc and 100% resolution."""
        # Simulate cluster proper time dilation
        tau_cluster = 1.0844
        tau_void = 1.0000
        time_elapsed = 1.0
        res = self.analyzer.predict_local_hubble_rate(tau_cluster, tau_void, time_elapsed)

        self.assertAlmostEqual(res["h0_predicted_reotransductor"], 73.04, delta=0.1)
        self.assertGreater(res["tension_resolution_pct"], 90.0)
        self.assertTrue(res["is_tension_mitigated"])

    def test_3d_environmental_gradient(self):
        """Verify that cluster environment has higher H_0 than void environment."""
        grid = 16
        rho = np.ones((grid, grid, grid), dtype=np.float32)
        tau = np.ones((grid, grid, grid), dtype=np.float32)

        # Create dense cluster in the center with high entropy/proper time
        rho[6:10, 6:10, 6:10] = 5.0
        tau[6:10, 6:10, 6:10] = 1.0844

        # Create deep void
        rho[0:4, 0:4, 0:4] = 0.2
        tau[0:4, 0:4, 0:4] = 1.0

        field_res = self.analyzer.compute_3d_environmental_h0_field(rho, tau, scale_factor=4.5, h0_engine=0.05)
        self.assertGreater(field_res["h0_cluster"], field_res["h0_void"])
        self.assertGreaterEqual(field_res["environmental_gradient_slope"], 0.0)

    def test_chi2_calculation(self):
        """Verify Chi-squared calculation gives 0.0 on identical models."""
        mu_obs = np.array([35.0, 38.0, 41.0])
        err_mu = np.array([0.1, 0.1, 0.1])
        chi2_zero = self.pantheon.compute_chi2(mu_obs, mu_obs, err_mu)
        self.assertAlmostEqual(chi2_zero, 0.0, places=4)


if __name__ == '__main__':
    unittest.main()
