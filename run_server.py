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
    parser.add_argument("--host", type=str, default="0.0.0.0", help="Host address to bind (default: 0.0.0.0)")
    parser.add_argument("--port", type=int, default=8000, help="Port to bind (default: 8000)")
    parser.add_argument("--speed", type=int, default=20, help="Initial steps per frame (default: 20)")
    parser.add_argument("--reset", action="store_true", help="Force clean reset to Eon 1 upon startup (archives previous checkpoints)")
    parser.add_argument("--log-level", type=str, default="info", help="Uvicorn log level (default: info)")
    
    args = parser.parse_args()

    if args.reset:
        os.environ["REOTRANSDUCTOR_FORCE_RESET"] = "1"

    print("=" * 70)
    print("  🌌 REOTRANSDUCTOR 3D COSMOLOGICAL SERVER (24/7 AUTONOMOUS ENGINE)")
    print("=" * 70)
    print(f"• Host Binding:       http://{args.host}:{args.port}")
    print(f"• Initial Speed:      {args.speed} steps/frame")
    print(f"• Mode:               {'🔴 CLEAN RESET (Eon 1 Primordial)' if args.reset else '🟢 AUTO-RESUME (Continuous)'}")
    print(f"• Checkpoint Path:    ./checkpoints/latest.npz")
    print(f"• Multi-Core CPU:     Enabled (NumPy / OpenBLAS multi-threading)")
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
