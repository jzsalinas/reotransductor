#!/usr/bin/env python3
"""
Top-level entrypoint for the Reotransductor 24/7 Cosmological Server.
Optimized for deployment on Dell PowerEdge R820 multi-core servers.
"""

import argparse
import sys
import os
import uvicorn

def main():
    parser = argparse.ArgumentParser(description="Reotransductor 3D 24/7 Cosmological Server")
    default_host = os.getenv("REOTRANSDUCTOR_HOST", "127.0.0.1")
    parser.add_argument("--host", type=str, default=default_host, help=f"Host address to bind (default: {default_host})")
    parser.add_argument("--port", type=int, default=8000, help="Port to bind (default: 8000)")
    parser.add_argument("--grid", type=int, default=32, choices=[16, 32, 64, 128, 256], help="Spatial lattice grid resolution N (default: 32 for 32x32x32)")
    parser.add_argument("--box", type=float, default=100.0, help="Comoving cosmological box size in Mpc (default: 100.0, use 300.0 or 500.0 for BAO acoustic horizon sampling)")
    parser.add_argument("--speed", type=int, default=20, help="Initial steps per frame (default: 20)")
    parser.add_argument("--gpu", action="store_true", help="Enable NVIDIA CUDA GPU acceleration with CuPy")
    parser.add_argument("--cpu", action="store_true", help="Force CPU multi-core execution with NumPy / OpenBLAS (default)")
    parser.add_argument("--reset", action="store_true", help="Force clean reset to Eon 1 upon startup (archives previous checkpoints)")
    parser.add_argument("--headless", action="store_true", help="Execute in ultra-fast headless background mode without web server or UI overhead")
    parser.add_argument("--log-level", type=str, default="info", help="Uvicorn log level (default: info)")
    
    args = parser.parse_args()

    if args.box <= 0:
        parser.error("--box must be > 0")

    if args.headless:
        # Delegate directly to headless high-speed production runner
        from scripts.run_headless_simulation import main as headless_main
        sys.argv = [sys.argv[0], "--grid", str(args.grid), "--box", str(args.box)]
        if args.gpu:
            sys.argv.append("--gpu")
        if args.cpu:
            sys.argv.append("--cpu")
        if args.reset:
            sys.argv.append("--reset")
        headless_main()
        return

    # Determine hardware preference (default to CPU unless --gpu is explicitly requested)
    gpu_requested = bool(args.gpu)
    if args.cpu:
        gpu_requested = False

    # Test GPU availability
    gpu_available = False
    gpu_name = "N/A"
    if gpu_requested:
        try:
            import cupy as cp
            if cp.cuda.is_available() and cp.cuda.runtime.getDeviceCount() > 0:
                gpu_available = True
                gpu_name = cp.cuda.runtime.getDeviceProperties(0)['name'].decode()
        except Exception:
            gpu_available = False

    use_gpu = gpu_requested and gpu_available
    os.environ["REOTRANSDUCTOR_USE_GPU"] = "1" if use_gpu else "0"
    os.environ["REOTRANSDUCTOR_GRID_SIZE"] = str(args.grid)
    os.environ["REOTRANSDUCTOR_BOX_SIZE_MPC"] = str(args.box)

    if args.reset:
        os.environ["REOTRANSDUCTOR_FORCE_RESET"] = "1"
    if args.speed:
        os.environ["REOTRANSDUCTOR_INITIAL_SPEED"] = str(args.speed)

    print("=" * 70)
    print("  🌌 REOTRANSDUCTOR 3D COSMOLOGICAL SERVER (24/7 AUTONOMOUS ENGINE)")
    print("=" * 70)
    print(f"• Host Binding:       http://{args.host}:{args.port}")
    print(f"• Grid Resolution:    {args.grid}³ ({args.grid**3:,} voxels)")
    print(f"• Comoving Box Size:  {args.box:.1f} Mpc (dx = {args.box/args.grid:.2f} Mpc/cell)")
    print(f"• Initial Speed:      {args.speed} steps/frame")
    print(f"• Compute Hardware:   {'🟢 GPU: ' + gpu_name + ' (CuPy / CUDA)' if use_gpu else '🔵 CPU: Multi-Core (NumPy / OpenBLAS)'}")
    print(f"• Mode:               {'🔴 CLEAN RESET (Eon 1 Primordial)' if args.reset else '🟢 AUTO-RESUME (Continuous)'}")
    print(f"• Checkpoint Path:    ./checkpoints/latest_g{args.grid}.npz")
    print("=" * 70)
    print("Press Ctrl+C to stop the server gracefully (auto-saves state).")
    print("=" * 70)

    uvicorn.run(
        "server.app:app",
        host=args.host,
        port=args.port,
        log_level=args.log_level,
        access_log=False
    )

if __name__ == "__main__":
    main()
