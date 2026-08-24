#!/usr/bin/env bash
# Automated Evaluation Pipeline for Paper III (Hubble Tension & NANOGrav)
# This script executes the latest version of the validation scripts on the generated checkpoints
# to reproduce the scientific artifacts located in results/paper_III/

set -e

echo "========================================================="
echo " RUNNING REOTRANSDUCTOR PAPER III VALIDATION PIPELINE    "
echo "========================================================="

echo "\n[1] Evaluating Cosmological Hubble Tension Environmental Gradient..."
# Searches for the latest Pantheon epoch checkpoint
LATEST_PANTHEON=$(ls -t checkpoints/pantheon_eon_*_g*.npz 2>/dev/null | head -n 1)

if [ -z "$LATEST_PANTHEON" ]; then
    echo "No pantheon_eon_N_gX.npz checkpoint found! Please run the engine up to a=4.50 first."
else
    echo "Found latest checkpoint: $LATEST_PANTHEON"
    .venv/bin/python experiments/compare_pantheon.py --checkpoint "$LATEST_PANTHEON"
    echo "✓ Hubble Tension validation complete."
fi

echo "\n[2] Evaluating Galactic Proper Time Micro-Drifts (NANOGrav 15-Yr)..."
LATEST_NANOGRAV=$(ls -t checkpoints/eon_*_g*.npz 2>/dev/null | head -n 1)

if [ -z "$LATEST_NANOGRAV" ]; then
    echo "No eon_N_gX.npz final eon checkpoint found! Please run a full eon first."
else
    echo "Found latest checkpoint: $LATEST_NANOGRAV"
    .venv/bin/python experiments/compare_nanograv.py --checkpoint "$LATEST_NANOGRAV"
    echo "✓ NANOGrav 15-Year Hellings-Downs correlation validation complete."
fi

echo "\n========================================================="
echo " PIPELINE COMPLETE. Artifacts saved in results/paper_III/ "
echo "========================================================="
