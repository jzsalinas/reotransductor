import unittest
import numpy as np
import os
import tempfile
from server.engine import CosmologicalEngine

class TestRestartReproducibility(unittest.TestCase):
    def test_resume_matches_uninterrupted_run(self):
        """
        Verify that N+M continuous steps perfectly match N steps -> checkpoint -> M steps.
        This tests the proper conservation of RNG state and state arrays.
        """
        # Create continuous engine
        engine_continuous = CosmologicalEngine(auto_resume=False)
        engine_continuous.seed = 42
        engine_continuous.reset_simulation(archive_existing=False)

        # Run 3 steps continuously
        for _ in range(3):
            engine_continuous.step()

        # Create checkpoint engine
        engine_ckpt = CosmologicalEngine(auto_resume=False)
        engine_ckpt.seed = 42
        engine_ckpt.reset_simulation(archive_existing=False)

        # Run 2 steps and save checkpoint
        for _ in range(2):
            engine_ckpt.step()
            
        with tempfile.NamedTemporaryFile(suffix='.npz', delete=False) as tmp:
            ckpt_path = tmp.name

        try:
            engine_ckpt.save_checkpoint(filepath=ckpt_path)

            # Create resumed engine and load checkpoint
            engine_resumed = CosmologicalEngine(auto_resume=False)
            loaded = engine_resumed.load_checkpoint(ckpt_path)
            self.assertTrue(loaded)

            # Run 1 more step
            engine_resumed.step()

            # Compare arrays
            np.testing.assert_allclose(engine_continuous.to_cpu(engine_continuous.rho), engine_resumed.to_cpu(engine_resumed.rho), rtol=1e-5, atol=1e-6)
            np.testing.assert_allclose(engine_continuous.to_cpu(engine_continuous.T), engine_resumed.to_cpu(engine_resumed.T), rtol=1e-5, atol=1e-6)
            np.testing.assert_allclose(engine_continuous.to_cpu(engine_continuous.v_x), engine_resumed.to_cpu(engine_resumed.v_x), rtol=1e-5, atol=1e-6)
            np.testing.assert_allclose(engine_continuous.to_cpu(engine_continuous.tau), engine_resumed.to_cpu(engine_resumed.tau), rtol=1e-5, atol=1e-6)
            np.testing.assert_allclose(engine_continuous.to_cpu(engine_continuous.I), engine_resumed.to_cpu(engine_resumed.I), rtol=1e-5, atol=1e-6)
            
            # Check RNG state
            self.assertEqual(engine_continuous.rng.bit_generator.state, engine_resumed.rng.bit_generator.state)
        finally:
            if os.path.exists(ckpt_path):
                os.remove(ckpt_path)

if __name__ == '__main__':
    unittest.main()
