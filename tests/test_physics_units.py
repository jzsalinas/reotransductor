"""
Unit test suite for server/physics_units.py.
Tests physical constant derivations, Planck scales, dimensional validity of kappa_0,
and exact invertibility of all code-to-physical unit transformations.
"""

import unittest
import numpy as np
from server.physics_units import FundamentalConstants, PlanckScales, CosmologicalUnits


class TestPhysicsUnits(unittest.TestCase):

    def setUp(self):
        self.const = FundamentalConstants()
        self.planck = PlanckScales()
        self.units = CosmologicalUnits(box_size_mpc=100.0, grid_resolution=32, c_code=2.5, h0_km_s_mpc=70.0)

    def test_planck_scales_magnitude(self):
        """Verify derived Planck scales match standard physical orders of magnitude."""
        self.assertAlmostEqual(np.log10(self.planck.LENGTH), -34.79, places=1)
        self.assertAlmostEqual(np.log10(self.planck.TIME), -43.26, places=1)
        self.assertAlmostEqual(np.log10(self.planck.MASS), -7.66, places=1)
        self.assertAlmostEqual(np.log10(self.planck.DENSITY), 96.71, places=1)
        self.assertAlmostEqual(np.log10(self.planck.TEMPERATURE), 32.15, places=1)

    def test_kappa_0_derivation(self):
        """Verify kappa_0 formula and order of magnitude."""
        # kappa_0 = (hbar^2 * G^2) / (c^7 * k_B)
        expected_kappa_0 = (
            (self.const.HBAR ** 2) * (self.const.G ** 2)
        ) / (
            (self.const.C ** 7) * self.const.K_B
        )
        self.assertEqual(self.planck.KAPPA_0, expected_kappa_0)
        self.assertAlmostEqual(np.log10(self.planck.KAPPA_0), -124.78, delta=0.5)

    def test_density_inversion(self):
        """Verify round-trip density conversions."""
        rho_code_test = np.array([0.05, 1.0, 5.5, 12.0], dtype=np.float64)
        
        # SI Round-trip
        rho_si = self.units.density_code_to_si(rho_code_test)
        rho_code_recovered_si = self.units.density_si_to_code(rho_si)
        np.testing.assert_allclose(rho_code_test, rho_code_recovered_si, rtol=1e-12)

        # Astrophysical Round-trip (M_sun / Mpc^3)
        rho_astro = self.units.density_code_to_astrophysical(rho_code_test)
        rho_code_recovered_astro = self.units.density_astrophysical_to_code(rho_astro)
        np.testing.assert_allclose(rho_code_test, rho_code_recovered_astro, rtol=1e-12)

    def test_velocity_inversion(self):
        """Verify round-trip velocity conversions."""
        v_code_test = np.array([0.0, 0.5, 1.5, 2.5], dtype=np.float64)
        v_km_s = self.units.velocity_code_to_km_s(v_code_test)
        v_code_rec = self.units.velocity_km_s_to_code(v_km_s)
        np.testing.assert_allclose(v_code_test, v_code_rec, rtol=1e-12)

        # Ensure c_code (2.5) strictly matches physical c in km/s (approx 299,792 km/s)
        c_km_s = self.units.velocity_code_to_km_s(2.5)
        self.assertAlmostEqual(c_km_s, self.const.C / 1000.0, places=3)

    def test_time_inversion(self):
        """Verify round-trip time conversions."""
        t_code_test = np.array([1.0, 10.0, 100.0, 1000.0], dtype=np.float64)
        
        # Myr Round-trip
        t_myr = self.units.time_code_to_myr(t_code_test)
        t_code_rec = self.units.time_myr_to_code(t_myr)
        np.testing.assert_allclose(t_code_test, t_code_rec, rtol=1e-12)

    def test_temperature_inversion(self):
        """Verify round-trip temperature conversions."""
        t_code_test = np.array([2.73, 10.0, 25.0, 50.0], dtype=np.float64)
        t_astro = self.units.temperature_code_to_astrophysical(t_code_test)
        t_code_rec = self.units.temperature_astrophysical_to_code(t_astro)
        np.testing.assert_allclose(t_code_test, t_code_rec, rtol=1e-12)

    def test_critical_cosmic_density(self):
        """Verify standard cosmological critical density value at H0 = 70 km/s/Mpc."""
        # Standard rho_crit is approx 9.2e-27 kg/m^3 or 1.36e11 M_sun/Mpc^3
        self.assertAlmostEqual(self.units.rho_crit_si * 1e27, 9.205, delta=0.1)
        self.assertAlmostEqual(self.units.rho_crit_msun_mpc3 / 1e11, 1.36, delta=0.1)

    def test_hubble_code_unit_roundtrip(self):
        """Verify that H_0 is dynamically computed from physical h0_si * time_unit_s."""
        expected = self.units.h0_si * self.units.time_unit_s
        self.assertAlmostEqual(self.units.H_0, expected, places=12)
        self.assertAlmostEqual(self.units.get_hubble_code_unit(), expected, places=12)

    def test_spitzer_power_law_scaling(self):
        """Verify that Spitzer thermal conductivity scales as T^(5/2) within the unclamped range."""
        base = 0.1
        # T1 = 1.5 * 2.73 ==> spitzer_k = 0.1 * 1.5^2.5 / 1.0 approx 0.275 (within [0.05, 2.5])
        # T2 = 3.0 * 2.73 ==> spitzer_k = 0.1 * 3.0^2.5 / 1.0 approx 1.559 (within [0.05, 2.5])
        k1 = self.units.compute_spitzer_conductivity(1.5 * 2.73, 0.0, base_k=base)
        k2 = self.units.compute_spitzer_conductivity(3.0 * 2.73, 0.0, base_k=base)
        ratio = k2 / k1
        expected_ratio = 2.0 ** 2.5
        self.assertAlmostEqual(ratio, expected_ratio, delta=0.05)

    def test_kappa_0_exact_value(self):
        """Verify that kappa_0 evaluated from CODATA constants equals 1.6487e-125."""
        self.assertAlmostEqual(self.planck.KAPPA_0 / 1.6487e-125, 1.0, delta=0.01)


if __name__ == '__main__':
    unittest.main()
