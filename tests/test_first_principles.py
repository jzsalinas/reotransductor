"""
Unit tests for First-Principles Physics & Plasma Thermodynamics in Reotransductor.
Tests dynamic adiabatic sound speed, Spitzer-Braginskii thermal conductivity,
Landauer thermal erasure rate, and Bekenstein-Hawking quadratic mass scaling.
"""

import unittest
import numpy as np
from server.physics_units import CosmologicalUnits
from server.engine import CosmologicalEngine


class TestFirstPrinciplesPhysics(unittest.TestCase):

    def setUp(self):
        self.units = CosmologicalUnits(box_size_mpc=100.0, grid_resolution=32, c_code=2.5, h0_km_s_mpc=70.0)
        self.engine = CosmologicalEngine(auto_resume=False)

    def test_adiabatic_sound_speed_monotonicity_and_relativistic_bound(self):
        """Verify that sound speed increases with temperature and stays strictly subluminal."""
        T_cold = 2.73
        T_warm = 12.0
        T_hot = 2000.0

        cs2_cold = self.units.compute_sound_speed_sq(T_cold, base_cs2=0.18)
        cs2_warm = self.units.compute_sound_speed_sq(T_warm, base_cs2=0.18)
        cs2_hot = self.units.compute_sound_speed_sq(T_hot, base_cs2=0.18)

        self.assertGreater(cs2_warm, cs2_cold)
        self.assertGreater(cs2_hot, cs2_warm)

        # Relativistic limit: c_s <= c / sqrt(3) ==> c_s^2 <= c^2 / 3
        cs2_max_allowed = (self.units.c_code ** 2) / 3.0
        self.assertAlmostEqual(cs2_hot, cs2_max_allowed, places=5)

    def test_spitzer_conductivity_scaling(self):
        """Verify that Spitzer conductivity increases with temperature and is regularized by density."""
        T_array = np.array([2.73, 6.0, 15.0], dtype=np.float32)
        rho_array = np.ones_like(T_array)

        kappa = self.units.compute_spitzer_conductivity(T_array, rho_array, base_k=0.3)
        self.assertTrue((kappa[1:] > kappa[:-1]).all())
        self.assertTrue((kappa >= 0.05).all())
        self.assertTrue((kappa <= 2.5).all())

        # Test extreme thermal plasma clamp at 2000 K
        kappa_extreme = self.units.compute_spitzer_conductivity(2000.0, 1.0, base_k=0.3)
        self.assertEqual(kappa_extreme, 2.5)

    def test_landauer_decay_thermal_proportionality(self):
        """Verify that Landauer informational erasure rate scales linearly with temperature."""
        T_1 = 2.73
        T_2 = 5.46

        gamma_1 = self.units.compute_landauer_decay(T_1, base_decay=0.015)
        gamma_2 = self.units.compute_landauer_decay(T_2, base_decay=0.015)

        self.assertAlmostEqual(gamma_2 / gamma_1, 2.0, places=4)

    def test_bekenstein_entropy_quadratic_mass_scaling(self):
        """Verify that Bekenstein-Hawking entropy limit scales quadratically with mass."""
        M_1 = 5000.0
        M_2 = 10000.0

        S_1 = self.units.compute_bekenstein_entropy_limit(M_1, m0_ref=5000.0, zeta_base=3500.0)
        S_2 = self.units.compute_bekenstein_entropy_limit(M_2, m0_ref=5000.0, zeta_base=3500.0)

        self.assertEqual(S_1, 3500.0)
        self.assertEqual(S_2, 3500.0 * 4.0)

    def test_continuous_differential_integration_stability(self):
        """Verify that continuous 50-step integration with first-principles fields is completely stable."""
        for step_idx in range(50):
            self.engine.step()

        self.assertFalse(np.isnan(self.engine.rho).any())
        self.assertFalse(np.isnan(self.engine.T).any())
        self.assertFalse(np.isnan(self.engine.I).any())
        self.assertFalse(np.isnan(self.engine.tau).any())
        self.assertTrue((self.engine.rho >= 0.02).all())
        self.assertTrue((self.engine.T >= 2.73).all())
        self.assertTrue((self.engine.I >= 0.0).all())
        self.assertTrue((self.engine.I <= 1.0).all())


if __name__ == '__main__':
    unittest.main()
