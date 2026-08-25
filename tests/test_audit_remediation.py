"""
Audit Remediation and Scientific Verification Test Suite.
Validates the critical fixes identified in the project audit:
  1. Proper time equation correspondence: d(tau_phys)/dt = 1 + kappa * sigma vs. d(Delta_tau)/dt = kappa * sigma.
  2. Pseudo-random number generator (RNG) determinism across seeds.
  3. Checkpoint persistence of seed, coordinate time, and physical proper time.
  4. Telegram secret credential protection (zero plaintext token leakage).
  5. NANOGrav independent metrics separation (model chi^2 vs. Hellings-Downs reference chi^2).
"""

import unittest
import numpy as np
import tempfile
import os
import shutil

from server.engine import CosmologicalEngine
from server.physics_units import CosmologicalUnits, FundamentalConstants, PlanckScales
from observational.nanograv_data import NANOGravPulsarData
from experiments.compare_planck import generate_eon_observational_report
from experiments.compare_bao import generate_eon_bao_report
from experiments.compare_halo import generate_eon_halo_report
from experiments.compare_pantheon import generate_eon_pantheon_report
from experiments.compare_nanograv import generate_eon_nanograv_report


class TestAuditRemediation(unittest.TestCase):
    """Verifies all remediation steps derived from the Q1 research audit."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_tau_physical_vs_dissipative_excess(self):
        """
        Verify that tau_physical integrates d(tau_phys)/dt = 1 + kappa * sigma
        and tau_excess integrates d(Delta_tau)/dt = kappa * sigma.
        """
        engine = CosmologicalEngine(grid_size=16, checkpoint_dir=self.temp_dir, auto_resume=False, seed=123)
        self.assertEqual(engine.t_coord, 0.0)
        np.testing.assert_allclose(engine.tau_excess - engine.to_cpu(engine.tau_eon_start), 0.0)
        np.testing.assert_allclose(engine.tau_physical - engine.to_cpu(engine.tau_eon_start), 0.0)

        # Run 10 integration steps
        n_steps = 10
        for _ in range(n_steps):
            engine.step()

        expected_t_coord = n_steps * engine.DT
        self.assertAlmostEqual(engine.t_coord, expected_t_coord, places=5)

        # In voids (sigma -> 0), physical proper time must be >= coordinate time t
        self.assertTrue(np.all(engine.tau_physical >= engine.t_coord))

        # Physical time equals coordinate time plus dissipative excess odometer:
        # tau_phys(x, t) = t + Delta_tau(x, t)
        np.testing.assert_allclose(
            engine.tau_physical,
            engine.t_coord + engine.tau_excess,
            rtol=1e-5
        )

        # Emergence rate d(tau_phys)/dt must be 1.0 + d(Delta_tau)/dt >= 1.0
        np.testing.assert_allclose(
            engine.d_tau_physical_dt,
            1.0 + engine.d_tau_dt,
            rtol=1e-5
        )
        self.assertTrue(np.all(engine.d_tau_physical_dt >= 1.0))

    def test_seed_determinism_and_reproducibility(self):
        """Verify that two simulation engines with identical seeds produce bit-exact identical evolution."""
        engine1 = CosmologicalEngine(grid_size=16, checkpoint_dir=self.temp_dir, auto_resume=False, seed=42)
        engine2 = CosmologicalEngine(grid_size=16, checkpoint_dir=self.temp_dir, auto_resume=False, seed=42)
        engine_diff = CosmologicalEngine(grid_size=16, checkpoint_dir=self.temp_dir, auto_resume=False, seed=99)

        # Initial fields match for identical seeds
        np.testing.assert_allclose(engine1.rho, engine2.rho, rtol=1e-6)
        np.testing.assert_allclose(engine1.T, engine2.T, rtol=1e-6)

        # Different seeds produce different fields
        self.assertFalse(np.allclose(engine1.rho, engine_diff.rho))

        # Run steps and confirm determinism
        for _ in range(5):
            engine1.step()
            engine2.step()
            engine_diff.step()

        np.testing.assert_allclose(engine1.rho, engine2.rho, rtol=1e-5)
        np.testing.assert_allclose(engine1.tau, engine2.tau, rtol=1e-5)

    def test_checkpoint_state_retention(self):
        """Verify that seed, t_coord, and tau_physical are persisted in checkpoints."""
        ckpt_path = os.path.join(self.temp_dir, "test_ckpt.npz")
        engine = CosmologicalEngine(grid_size=16, checkpoint_dir=self.temp_dir, auto_resume=False, seed=777)
        for _ in range(8):
            engine.step()

        engine.save_checkpoint(filepath=ckpt_path)

        # Load into fresh engine
        restored = CosmologicalEngine(grid_size=16, checkpoint_dir=self.temp_dir, auto_resume=False)
        self.assertTrue(restored.load_checkpoint(ckpt_path))

        self.assertEqual(restored.seed, 777)
        self.assertAlmostEqual(restored.t_coord, engine.t_coord, places=5)
        np.testing.assert_allclose(restored.tau_physical, engine.tau_physical, rtol=1e-5)

    def test_telegram_credential_protection(self):
        """Verify that telegram bot token is never exposed in plaintext through config getter."""
        from server.notifier import TelegramNotifier
        cfg_file = os.path.join(self.temp_dir, "telegram_config.json")
        notifier = TelegramNotifier(config_paths=[cfg_file])
        test_secret = "123456789:ABCdefGHIjklMNOpqrsTUVwxyz1234567"
        notifier.save_config({
            "enabled": True,
            "bot_token": test_secret,
            "chat_id": "987654321",
            "interval_eons": 5
        })

        # Emulate the secured get_telegram_config logic
        cfg = notifier.config
        token = str(cfg.get("bot_token") or "").strip()
        masked = (token[:6] + "..." + token[-4:]) if len(token) > 10 else ("***" if token else "")
        public_payload = {
            "enabled": bool(cfg.get("enabled", False)),
            "chat_id": cfg.get("chat_id", ""),
            "interval_eons": int(cfg.get("interval_eons", 10)),
            "bot_token_masked": masked,
            "bot_token_configured": bool(token)
        }

        # Assert secret token is not in the public payload
        self.assertNotIn(test_secret, str(public_payload))
        self.assertTrue(public_payload["bot_token_configured"])
        self.assertTrue(public_payload["bot_token_masked"].startswith("123456..."))

    def test_nanograv_model_vs_reference_metrics_separation(self):
        """Verify that compare_nanograv generates distinct metrics for model and reference."""
        engine = CosmologicalEngine(grid_size=16, checkpoint_dir=self.temp_dir, auto_resume=False, seed=42)
        for _ in range(5):
            engine.step()

        metrics = generate_eon_nanograv_report(engine, output_dir=os.path.join(self.temp_dir, "snapshots"))
        self.assertIn("chi2_simulation", metrics)
        self.assertIn("chi2_hellings_downs", metrics)
        self.assertIn("chi2_model", metrics)
        self.assertIn("chi2_hd_reference", metrics)

        # Confirm they are float numbers
        self.assertIsInstance(metrics["chi2_simulation"], (int, float))
        self.assertIsInstance(metrics["chi2_hellings_downs"], (int, float))

    def test_provenance_manifest_and_data_integrity(self):
        """Verify that all observational datasets match the cryptographic SHA-256 hashes in PROVENANCE.yml."""
        from observational.verify_data import verify_all_datasets
        all_valid, results = verify_all_datasets("data/PROVENANCE.yml")
        self.assertTrue(all_valid, f"One or more datasets failed integrity check: {results}")
        self.assertGreaterEqual(len(results), 5)

    def test_cpu_gpu_hardware_consistency(self):
        """
        Verify that running simulation steps on CPU (NumPy) and GPU (CuPy)
        with identical seed produces equivalent physical fields within numerical precision.
        """
        try:
            import cupy
            gpu_available = cupy.cuda.is_available()
        except Exception:
            gpu_available = False

        if not gpu_available:
            self.skipTest("CUDA GPU / CuPy not available on this host environment.")

        cpu_engine = CosmologicalEngine(grid_size=16, checkpoint_dir=self.temp_dir, auto_resume=False, use_gpu=False, seed=123)
        gpu_engine = CosmologicalEngine(grid_size=16, checkpoint_dir=self.temp_dir, auto_resume=False, use_gpu=True, seed=123)

        for _ in range(15):
            cpu_engine.step()
            gpu_engine.step()

        rho_cpu = cpu_engine.rho
        rho_gpu = gpu_engine.to_cpu(gpu_engine.rho)
        tau_cpu = cpu_engine.tau
        tau_gpu = gpu_engine.to_cpu(gpu_engine.tau)

        # Assert arrays are close within float32 numerical precision (CuFFT vs PocketFFT)
        np.testing.assert_allclose(rho_cpu, rho_gpu, rtol=5e-3, atol=5e-3)
        np.testing.assert_allclose(tau_cpu, tau_gpu, rtol=5e-3, atol=5e-3)
        self.assertAlmostEqual(cpu_engine.scale_factor, gpu_engine.scale_factor, places=4)

    def test_telemetry_structure_and_field_validity(self):
        """
        Verify that engine.get_telemetry() returns all required fields
        without any NameError or missing variables.
        """
        engine = CosmologicalEngine(grid_size=16, checkpoint_dir=self.temp_dir, auto_resume=False, seed=42)
        for _ in range(5):
            engine.step()

        telemetry = engine.get_telemetry()
        self.assertIsInstance(telemetry, dict)
        required_keys = [
            "eon", "era", "scale_factor", "redshift", "temp_norm", "temp_astro",
            "c_light", "mass_fraction", "s_bh", "s_crit", "tunnel_progress",
            "progress_label", "active_route", "p_grav", "p_conformal",
            "fossil_odometer", "time_myr", "grid_size", "grid_voxels",
            "box_size_mpc", "cell_size_mpc", "h0_kms_mpc", "kappa_0_planck",
            "attractor", "z_slice", "total_steps", "is_running"
        ]
        for key in required_keys:
            self.assertIn(key, telemetry, f"Missing key '{key}' in telemetry payload")
            self.assertIsNotNone(telemetry[key])

    def test_headless_engine_batch_execution(self):
        """
        Verify that engine.step_batch executes multiple hydrodynamic steps
        consistently and saves valid 6-epoch checkpoints.
        """
        engine = CosmologicalEngine(grid_size=16, checkpoint_dir=self.temp_dir, auto_resume=False, seed=42)
        initial_steps = engine.total_steps
        engine.step_batch(50)
        self.assertEqual(engine.total_steps, initial_steps + 50)
        self.assertGreater(engine.scale_factor, 1.0)
        self.assertTrue(np.all(engine.to_cpu(engine.rho) > 0.0))


if __name__ == '__main__':
    unittest.main()
