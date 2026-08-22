"""
Observational Dataset Verification and Integrity Check Tool.
Validates SHA-256 cryptographic hashes of all local datasets against data/PROVENANCE.yml.
"""

import os
import sys
import hashlib
from typing import Dict, Any, Tuple, List


def compute_sha256(filepath: str) -> str:
    """Computes SHA-256 hash of a file in 64KB chunks."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_simple_yaml_provenance(provenance_path: str = "data/PROVENANCE.yml") -> Dict[str, Dict[str, str]]:
    """
    Lightweight YAML parser for PROVENANCE.yml (works without external PyYAML dependency).
    """
    if not os.path.exists(provenance_path):
        raise FileNotFoundError(f"Provenance manifest not found at: {provenance_path}")

    datasets = {}
    current_key = None

    with open(provenance_path, "r", encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            if line.startswith("  ") and not line.startswith("    ") and stripped.endswith(":"):
                current_key = stripped[:-1].strip()
                datasets[current_key] = {}
            elif line.startswith("    ") and current_key is not None and ":" in stripped:
                k, v = stripped.split(":", 1)
                k = k.strip()
                v = v.strip().strip('"').strip("'")
                datasets[current_key][k] = v

    return datasets


def verify_all_datasets(provenance_path: str = "data/PROVENANCE.yml") -> Tuple[bool, List[Dict[str, Any]]]:
    """
    Verifies that all registered datasets exist and have valid SHA-256 hashes.
    Returns (all_valid, results_list).
    """
    datasets = parse_simple_yaml_provenance(provenance_path)
    results = []
    all_valid = True

    for key, meta in datasets.items():
        path = meta.get("path", "")
        expected_sha = meta.get("sha256", "")
        name = meta.get("name", key)

        if not os.path.exists(path):
            all_valid = False
            results.append({
                "key": key,
                "name": name,
                "path": path,
                "status": "MISSING",
                "expected_sha256": expected_sha,
                "actual_sha256": None,
                "valid": False
            })
            continue

        actual_sha = compute_sha256(path)
        is_match = (actual_sha.lower() == expected_sha.lower())
        if not is_match:
            all_valid = False

        results.append({
            "key": key,
            "name": name,
            "path": path,
            "status": "VALID" if is_match else "CHECKSUM_MISMATCH",
            "expected_sha256": expected_sha,
            "actual_sha256": actual_sha,
            "valid": is_match
        })

    return all_valid, results


def run_verification_cli():
    """CLI runner for dataset integrity verification."""
    print("=" * 80)
    print("  🔭 REOTRANSDUCTOR: OBSERVATIONAL DATA INTEGRITY & PROVENANCE CHECK")
    print("=" * 80)

    try:
        all_valid, results = verify_all_datasets()
        for r in results:
            symbol = "✅" if r["valid"] else "❌"
            print(f"{symbol} [{r['key']}] {r['name']}")
            print(f"   • Path: {r['path']}")
            print(f"   • SHA-256: {r['actual_sha256'] or 'N/A'}")
            if not r["valid"]:
                print(f"   ⚠️ Expected SHA-256: {r['expected_sha256']}")
            print("-" * 80)

        print("=" * 80)
        if all_valid:
            print("  ✨ ALL OBSERVATIONAL DATASETS VERIFIED (100% INTEGRITY)")
            print("=" * 80)
            return 0
        else:
            print("  ⚠️ INTEGRITY CHECK FAILED: One or more datasets are missing or corrupt.")
            print("=" * 80)
            return 1
    except Exception as e:
        print(f"❌ Error verifying datasets: {e}")
        return 1


if __name__ == '__main__':
    sys.exit(run_verification_cli())
