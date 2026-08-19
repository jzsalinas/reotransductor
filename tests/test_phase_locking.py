"""
Unit tests for Holographic Phase-Locking in Fourier Space.
Verifies real-valued inversion, phase coherence scaling with alpha_mem,
and power spectrum conservation across eons.
"""

import unittest
import numpy as np
from server.engine import CosmologicalEngine


class TestPhaseLocking(unittest.TestCase):

    def setUp(self):
        self.engine = CosmologicalEngine(auto_resume=False)
        self.grid_size = self.engine.grid_size

        # Create a synthetic structured fossil field with clear peaks
        X, Y, Z = self.engine.X, self.engine.Y, self.engine.Z
        self.synthetic_tau = (
            500.0 * np.exp(-((X - 16.0)**2 + (Y - 16.0)**2 + (Z - 16.0)**2) / 25.0) +
            300.0 * np.exp(-((X - 24.0)**2 + (Y - 8.0)**2 + (Z - 24.0)**2) / 20.0)
        ).astype(np.float32)

    def test_real_valued_inversion(self):
        """Verify that the phase-locked fluctuation field is strictly real and non-divergent."""
        fluct = self.engine._generate_phase_locked_fluctuations(self.synthetic_tau, alpha_mem=0.35)
        self.assertEqual(fluct.shape, (self.grid_size, self.grid_size, self.grid_size))
        self.assertEqual(fluct.dtype, np.float32)
        self.assertFalse(np.isnan(fluct).any())
        self.assertFalse(np.isinf(fluct).any())

    def test_alpha_coherence_correlation(self):
        """
        Verify that higher alpha_mem yields higher cross-correlation
        with the underlying fossil field structure.
        """
        np.random.seed(123)
        # alpha = 0.0 (pure uncorrelated quantum noise)
        fluct_quantum = self.engine._generate_phase_locked_fluctuations(self.synthetic_tau, alpha_mem=0.0)
        
        np.random.seed(123)
        # alpha = 0.95 (strong holographic phase-locking)
        fluct_fossil_locked = self.engine._generate_phase_locked_fluctuations(self.synthetic_tau, alpha_mem=0.95)

        # Compute spatial Pearson correlation with smoothed fossil field
        tau_flat = self.synthetic_tau.ravel() - np.mean(self.synthetic_tau)
        tau_norm = tau_flat / np.linalg.norm(tau_flat)

        corr_quantum = float(np.dot(fluct_quantum.ravel(), tau_norm))
        corr_locked = float(np.dot(fluct_fossil_locked.ravel(), tau_norm))

        # Locked field must show significantly higher correlation with the fossil field
        self.assertGreater(corr_locked, corr_quantum)

    def test_power_spectrum_conservation(self):
        """Verify that the variance/power is well-conditioned and non-zero."""
        fluct = self.engine._generate_phase_locked_fluctuations(self.synthetic_tau, alpha_mem=0.35)
        std_val = float(np.std(fluct))
        # Scaled to approx 0.35 amplitude
        self.assertAlmostEqual(std_val, 0.35, places=2)

    def test_white_hole_bounce_execution(self):
        """Verify that _trigger_white_hole_eon_3d executes cleanly without numerical errors."""
        self.engine.rho = np.ones((self.grid_size, self.grid_size, self.grid_size), dtype=np.float32)
        self.engine.tau = self.synthetic_tau.copy()
        
        rho_new, vx, vy, vz, T_new = self.engine._trigger_white_hole_eon_3d()
        
        self.assertTrue((rho_new >= 0.05).all())
        self.assertTrue((rho_new <= 12.0).all())
        self.assertTrue((T_new >= 2.73).all())
        self.assertEqual(vx.shape, (self.grid_size, self.grid_size, self.grid_size))


if __name__ == '__main__':
    unittest.main()
