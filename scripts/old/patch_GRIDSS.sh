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
PATCHED_IMAGE=${GRIDSS_DIR}/Patched_GRIDSS_v03.sif
TMP_DIR=${GRIDSS_DIR}/gridss_patch_tmp
SANDBOX_DIR=${TMP_DIR}/sandbox


echo "Starting GRIDSS patching process at $(date)"
echo "============================================"

# --- Prepare working directory ---
echo "Cleaning up any existing temporary directory..."
rm -rf "${TMP_DIR}"
mkdir -p "${TMP_DIR}"
mkdir -p "${SANDBOX_DIR}"

# --- Build a writable sandbox copy ---
echo "Step 1: Creating writable sandbox..."
apptainer build --sandbox --force "${SANDBOX_DIR}" "${ORIG_IMAGE}"
echo "Sandbox created at: ${SANDBOX_DIR}"
echo ""

# Fix the storage directory issue(it needs to exist for apptainer))
mkdir -p "${SANDBOX_DIR}/storage"

# --- Apply the Rscript path fix inside the sandbox ---

## This would be replacement of GitHub issue fix:
# echo "Step 2: Patching *.R scripts inside sandbox..."
# apptainer exec --writable --no-mount all "${SANDBOX_DIR}" bash -c '
#     echo "Finding and patching R scripts..."
#     find /opt/gridss/ -name "*.R" -exec sed -i "s/Rscript/\/usr\/bin\/Rscript/g" {} \;
#     echo "Rscript paths patched."
# '
# echo ""

# We have files that use #!/usr/bin/env Rscript ; We find these with:
# grep -R "Rscript" /storage/homefs/kw23y068/software/gridss/gridss_patch_tmp/sandbox/opt/gridss 2>/dev/null
# We find
# opt/gridss/gridss_extract_overlapping_fragments
# opt/gridss/gridss_somatic_filter
# opt/gridss/gridss
# opt/gridss/virusbreakend
# Hence we need to patch both the R scripts and the shell wrappers that call Rscript
echo "Step 2: Patching Rscript references in GRIDSS shell wrappers..."
apptainer exec --writable --no-mount all "${SANDBOX_DIR}" bash -c '
    echo "Finding and patching Rscript calls in wrapper scripts..."
    find /opt/gridss -type f \
        \( -name "gridss*" -o -name "virusbreakend*" \) \
        -exec sed -i "s/\bRscript\b/\/usr\/bin\/Rscript/g" {} +
    echo "Rscript paths patched."
'


# --- Verify the patch ---
echo "Step 3: Verifying the patch..."
apptainer exec --no-mount all "${SANDBOX_DIR}" bash -c '
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

apptainer exec --no-mount all "${SANDBOX_DIR}" /usr/bin/Rscript --version

cd ${GRIDSS_DIR}
echo "Patched files:"
#find /opt /usr/local /gridss -type f -name "*" -exec grep -l "/usr/bin/Rscript" {} + 2>/dev/null || true

# --- Clean up ---
#If you need the tmp files for debugging, comment out:
##echo "Step 5: Cleaning up temporary files..."
#rm -rf "${TMP_DIR}"
#echo "Temporary directory removed."
#echo ""



echo "============================================"
echo "Patched image created at: ${PATCHED_IMAGE}"
echo "Patching complete at $(date)"