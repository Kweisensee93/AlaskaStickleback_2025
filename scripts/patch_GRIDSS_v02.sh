#!/bin/bash
#SBATCH --job-name=patch_GRIDSS
#SBATCH --output=/storage/homefs/kw23y068/logfiles/patch_GRIDSS_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/patch_GRIDSS_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=epyc2

set -euo pipefail

# --- Paths ---
GRIDSS_DIR=/storage/homefs/kw23y068/software/gridss/
ORIG_IMAGE=${GRIDSS_DIR}/GRIDSS.sif
PATCHED_IMAGE=${GRIDSS_DIR}/Patched_GRIDSS_v02.sif
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

# --- Create necessary directories ---
echo "Step 2: Creating required directories in sandbox..."
mkdir -p "${SANDBOX_DIR}/storage"
mkdir -p "${SANDBOX_DIR}/usr/bin"
echo "Directories created."
echo ""

# --- Create a wrapper script for Rscript ---
echo "Step 3: Creating Rscript wrapper script..."
cat > "${SANDBOX_DIR}/usr/bin/Rscript" << 'EOF'
#!/bin/bash
# Wrapper script to find and use Rscript from host
# This script will be called by GRIDSS

# Try to find Rscript in common locations
for rscript_path in \
    /usr/local/bin/Rscript \
    /usr/bin/Rscript.real \
    $(which Rscript 2>/dev/null) \
    /opt/R/*/bin/Rscript \
    ; do
    if [ -x "$rscript_path" ] 2>/dev/null; then
        exec "$rscript_path" "$@"
    fi
done

# If we get here, try using 'Rscript' from PATH
if command -v Rscript >/dev/null 2>&1; then
    exec Rscript "$@"
fi

echo "ERROR: Rscript not found. Please bind mount R into the container." >&2
echo "Example: apptainer exec --bind /path/to/R/bin/Rscript:/usr/bin/Rscript.real ..." >&2
exit 1
EOF

chmod +x "${SANDBOX_DIR}/usr/bin/Rscript"
echo "Rscript wrapper created."
echo ""

# --- Verify the wrapper was created ---
echo "Step 4: Verifying wrapper creation..."
if [ -x "${SANDBOX_DIR}/usr/bin/Rscript" ]; then
    echo "✓ Wrapper script created successfully at /usr/bin/Rscript"
    echo "Wrapper contents:"
    head -n 5 "${SANDBOX_DIR}/usr/bin/Rscript"
else
    echo "✗ ERROR: Failed to create wrapper script"
    exit 1
fi
echo ""

# --- Optional: Patch any hardcoded paths (from original GitHub fix) ---
#echo "Step 5: Patching any hardcoded Rscript paths (optional)..."
#find "${SANDBOX_DIR}/opt/gridss/" -name "*.R" -type f -exec sed -i 's|#!/usr/bin/env Rscript|#!/usr/bin/Rscript|g' {} \; 2>/dev/null || true
#echo "Patching complete."
#echo ""

# --- Rebuild immutable SIF ---
echo "Step 6: Repacking patched image to ${PATCHED_IMAGE}..."
apptainer build "${PATCHED_IMAGE}" "${SANDBOX_DIR}"
echo "Patched image created successfully."
echo ""

# --- Clean up ---
echo "Step 7: Cleaning up temporary files..."
rm -rf "${TMP_DIR}"
echo "Temporary directory removed."
echo ""

echo "============================================"
echo "Patched image created at: ${PATCHED_IMAGE}"
echo "Patching complete at $(date)"
echo ""
echo "USAGE INSTRUCTIONS:"
echo "When running GRIDSS, bind your host Rscript like this:"
echo "  apptainer exec --bind \$(which Rscript):/usr/bin/Rscript.real \\"
echo "    ${PATCHED_IMAGE} /opt/gridss/gridss ..."