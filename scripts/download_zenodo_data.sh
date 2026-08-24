#!/usr/bin/env bash
#
# download_zenodo_data.sh
# 
# Script to download the high-resolution Cosmological Checkpoints (.npz files) 
# from the official Zenodo repository for the Reotransductor Paper III.
#
# Usage:
#   ./scripts/download_zenodo_data.sh

set -e

# TODO: Replace with the actual Zenodo record ID once the dataset is published.
ZENODO_RECORD_ID="XXXXXXX"
ZENODO_FILE="reotransductor_checkpoints_high_res.zip"

echo "========================================================="
echo " DOWNLOADING REOTRANSDUCTOR CHECKPOINTS FROM ZENODO      "
echo "========================================================="

if [ "$ZENODO_RECORD_ID" = "XXXXXXX" ]; then
    echo "[!] ERROR: The dataset has not been published yet."
    echo "[!] Please update the ZENODO_RECORD_ID variable in this script with the official ID."
    echo "    (e.g., if your DOI is 10.5281/zenodo.1234567, your RECORD_ID is 1234567)"
    exit 1
fi

echo "Target Record ID: $ZENODO_RECORD_ID"
echo "Downloading $ZENODO_FILE..."

# Download the zip file directly from the Zenodo API
wget -q --show-progress "https://zenodo.org/record/${ZENODO_RECORD_ID}/files/${ZENODO_FILE}?download=1" -O "${ZENODO_FILE}"

echo "Extracting checkpoints to checkpoints/ directory..."
# -n prevents overwriting existing files if the user already has some
unzip -n "${ZENODO_FILE}" -d checkpoints/

echo "Cleaning up..."
rm "${ZENODO_FILE}"

echo "✓ Download and extraction complete! The high-resolution matrices are now available."
echo "You can now run the validation pipelines (e.g., ./scripts/pipeline_paper_III.sh)."
