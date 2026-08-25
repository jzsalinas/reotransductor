#!/usr/bin/env python3
"""
Reotransductor Headless Cosmological Simulator (High-Performance CLI).
Executes pure 3D hydrodynamic & proper time field integration at maximum hardware speed
without FastAPI web server, WebSocket streaming, or JSON serialization overhead.
Ideal for large grid resolutions (64^3, 128^3, 256^3) and long multi-eon production runs.
"""

import os
import sys
import time
import argparse
import numpy as np

# Add project root to Python module path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from server.engine import CosmologicalEngine
from server.physics_units import CosmologicalUnits


def parse_args():
    parser = argparse.ArgumentParser(
        description="Reotransductor Headless Production Simulator (Ultra-Fast 3D Engine)"
    )
    parser.add_argument(
        "--grid",
        type=int,
        default=32,
        choices=[16, 32, 64, 128, 256],
        help="Spatial lattice grid resolution N (default: 32 for 32x32x32)",
    )
    parser.add_argument(
        "--box",
        type=float,
        default=500.0,
        help="Physical comoving cosmological box size in Mpc (default: 500.0)",
    )
    parser.add_argument(
        "--eons",
        type=int,
        default=1,
        help="Number of cosmological eons to simulate (default: 1)",
    )
    parser.add_argument(
        "--gpu",
        action="store_true",
        help="Enable NVIDIA CUDA GPU acceleration with CuPy",
    )
    parser.add_argument(
        "--cpu",
        action="store_true",
        help="Force CPU multi-core execution with NumPy / OpenBLAS",
    )
    parser.add_argument(
        "--reset",
        action="store_true",
        help="Force clean reset to Eon 1 upon startup",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducible initial perturbations (default: 42)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=200,
        help="Number of hydrodynamic integration steps per compute batch (default: 200)",
    )
    parser.add_argument(
        "--checkpoint-dir",
        type=str,
        default="checkpoints",
        help="Directory to persist epoch and snapshot checkpoints (default: checkpoints)",
    )
    parser.add_argument(
        "--max-steps",
        type=int,
        default=0,
        help="Optional maximum total steps (0 = run until target eons complete)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Hardware selection
    use_gpu = False
    gpu_name = "N/A"
    if args.gpu and not args.cpu:
        try:
            import cupy as cp
            if cp.cuda.is_available() and cp.cuda.runtime.getDeviceCount() > 0:
                use_gpu = True
                gpu_name = cp.cuda.runtime.getDeviceProperties(0)['name'].decode()
        except Exception as ex:
            print(f"⚠️  GPU requested but CuPy/CUDA initialization failed: {ex}. Falling back to CPU.")
            use_gpu = False

    # Set environment variables for CosmologicalEngine
    os.environ["REOTRANSDUCTOR_GRID_SIZE"] = str(args.grid)
    os.environ["REOTRANSDUCTOR_BOX_SIZE_MPC"] = str(args.box)
    os.environ["REOTRANSDUCTOR_USE_GPU"] = "1" if use_gpu else "0"

    print("=" * 75)
    print("  🌌 REOTRANSDUCTOR 3D HEADLESS SIMULATION RUNNER (MAXIMUM THROUGHPUT)")
    print("=" * 75)
    print(f"• Grid Resolution:     {args.grid}³ ({args.grid**3:,} voxels)")
    print(f"• Comoving Box Size:   {args.box:.1f} Mpc (dx = {args.box / args.grid:.2f} Mpc/cell)")
    print(f"• Target Eons:         {args.eons}")
    print(f"• Compute Hardware:    {'🟢 GPU: ' + gpu_name + ' (CuPy / CUDA)' if use_gpu else '🔵 CPU: Multi-Core (NumPy / OpenBLAS)'}")
    print(f"• Batch Step Size:     {args.batch_size} steps/iteration")
    print(f"• Random Seed:         {args.seed}")
    print(f"• Mode:                {'🔴 CLEAN RESET (Eon 1)' if args.reset else '🟢 RESUME/CONTINUOUS'}")
    print(f"• Checkpoint Dir:      {args.checkpoint_dir}")
    print("=" * 75)

    engine = CosmologicalEngine(
        grid_size=args.grid,
        checkpoint_dir=args.checkpoint_dir,
        auto_resume=not args.reset,
        use_gpu=use_gpu,
        seed=args.seed
    )

    initial_eon = engine.eon
    target_end_eon = initial_eon + args.eons
    start_time = time.time()
    last_log_time = start_time
    last_log_steps = engine.total_steps
    step_count = 0

    print(f"🚀 Starting simulation at Eon {engine.eon}, Scale Factor a = {engine.scale_factor:.3f} ...\n")

    try:
        while engine.eon < target_end_eon:
            engine.step_batch(args.batch_size)
            step_count += args.batch_size

            now = time.time()
            if now - last_log_time >= 3.0:
                elapsed = now - last_log_time
                steps_done = engine.total_steps - last_log_steps
                sps = steps_done / max(1e-5, elapsed)

                rho_max = float(np.max(engine.to_cpu(engine.rho)))
                rho_mean = float(np.mean(engine.to_cpu(engine.rho)))
                t_mean = float(np.mean(engine.to_cpu(engine.T)))
                tau_cpu = engine.to_cpu(engine.tau)
                tau_max = float(np.max(tau_cpu))
                bh_count = int(np.sum(tau_cpu >= 0.99e7))

                print(
                    f"[{time.strftime('%H:%M:%S')}] Eon {engine.eon:<2} | "
                    f"a = {engine.scale_factor:6.3f} | "
                    f"Prog: {engine.progress * 100.0:5.1f}% | "
                    f"ρ̄={rho_mean:.3f} (max={rho_max:5.2f}) | "
                    f"T̄={t_mean:5.1f}K | "
                    f"τ_max={tau_max:6.1f} (BHs: {bh_count}) | "
                    f"⚡ {sps:6.1f} steps/s"
                )
                last_log_time = now
                last_log_steps = engine.total_steps

            if args.max_steps > 0 and step_count >= args.max_steps:
                print(f"\nReached requested maximum step limit ({args.max_steps} steps). Halting.")
                break

    except KeyboardInterrupt:
        print("\n\n🛑 Simulation interrupted by user (Ctrl+C). Saving current state...")

    finally:
        engine.save_checkpoint()
        total_time = time.time() - start_time
        total_steps = engine.total_steps
        avg_sps = step_count / max(1e-5, total_time)

        print("\n" + "=" * 75)
        print("  🏁 HEADLESS SIMULATION COMPLETED")
        print("=" * 75)
        print(f"• Final Eon:           {engine.eon}")
        print(f"• Final Scale Factor:  a = {engine.scale_factor:.4f}")
        print(f"• Total Steps Run:     {step_count:,} steps ({total_steps:,} cumulative)")
        print(f"• Total Wall Time:     {total_time:.2f} s ({total_time / 60.0:.2f} min)")
        print(f"• Average Throughput:  {avg_sps:.1f} steps/second")
        print(f"• Checkpoint Saved:    {os.path.join(args.checkpoint_dir, f'latest_g{args.grid}.npz')}")
        print("=" * 75 + "\n")


if __name__ == "__main__":
    main()
