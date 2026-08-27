"""Command-line entry point for the protected K_hat=0 validation."""

from __future__ import annotations

import argparse
import json

from .validation import CONTROL_RESOLUTIONS, CONTROL_TIMES, run_control_validation


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", nargs="+", type=int, default=CONTROL_RESOLUTIONS)
    parser.add_argument("--times", nargs="+", type=float, default=CONTROL_TIMES)
    arguments = parser.parse_args()
    passed, records = run_control_validation(
        tuple(arguments.resolutions), tuple(arguments.times)
    )
    print("CONTROL VALIDATION — NOT A CONDUCTIVE SCIENTIFIC RUN")
    print(json.dumps(records, indent=2, sort_keys=True))
    print("CONTROL VALIDATED" if passed else "CONTROL VALIDATION FAILED")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
