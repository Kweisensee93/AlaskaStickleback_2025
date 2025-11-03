#!/bin/bash
#SBATCH --job-name=patch_GRIDSS
#SBATCH --output=/storage/homefs/kw23y068/logfiles/patch_GRIDSS_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/patch_GRIDSS_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=epyc2

set -euo pipefail

# Enable verbose output for debugging
export APPTAINER_DEBUG=true
export APPTAINER_VERBOSE=true

# --- Paths ---
GRIDSS_DIR=/storage/homefs/kw23y068/software/gridss/
ORIG_IMAGE=${GRIDSS_DIR}/GRIDSS.sif
PATCHED_IMAGE=${GRIDSS_DIR}/Patched_GRIDSS.sif
TMP_DIR=${GRIDSS_DIR}/gridss_patch_tmp
SANDBOX_DIR=${TMP_DIR}/sandbox

echo "Starting GRIDSS patching process at $(date)"
echo "============================================"

# --- Prepare working directory ---
echo "Cleaning up any existing temporary directory..."
rm -rf "${TMP_DIR}"
mkdir -p "${TMP_DIR}"

# --- Build a writable sandbox copy ---
echo "Step 1: Creating writable sandbox..."
apptainer build --sandbox "${SANDBOX_DIR}" "${ORIG_IMAGE}"
echo "Sandbox created at: ${SANDBOX_DIR}"
echo ""

# --- Apply the Rscript path fix inside the sandbox ---
echo "Step 2: Patching *.R scripts inside sandbox..."
apptainer exec --writable --no-mount bind-paths "${SANDBOX_DIR}" bash -c '
    echo "Finding and patching R scripts..."
    find /opt/gridss/ -name "*.R" -exec sed -i "s/Rscript/\/usr\/bin\/Rscript/g" {} \;
    echo "Rscript paths patched."
'
echo ""

# --- Verify the patch ---
echo "Step 3: Verifying the patch..."
apptainer exec --no-mount bind-paths "${SANDBOX_DIR}" bash -c '
    echo "Checking for patched Rscript paths..."
    if grep -r "/usr/bin/Rscript" /opt/gridss/*.R > /dev/null 2>&1; then
        echo "✓ Found /usr/bin/Rscript in patched files"
    else
        echo "⚠ Warning: /usr/bin/Rscript not found in R files"
    fi
    
    echo ""
    echo "Verifying Rscript availability..."
    which Rscript || echo "⚠ Rscript not found in PATH"
    
    echo ""
    echo "Rscript version:"
    Rscript --version 2>&1 || echo "⚠ Could not get Rscript version"
'
echo ""

# --- Rebuild immutable SIF ---
echo "Step 4: Repacking patched image to ${PATCHED_IMAGE}..."
apptainer build "${PATCHED_IMAGE}" "${SANDBOX_DIR}"
echo "Patched image created successfully."
echo ""

# --- Clean up ---
echo "Step 5: Cleaning up temporary files..."
rm -rf "${TMP_DIR}"
echo "Temporary directory removed."
echo ""

echo "============================================"
echo "Patched image created at: ${PATCHED_IMAGE}"
echo "Patching complete at $(date)"