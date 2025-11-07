#!/bin/bash
#SBATCH --job-name=GRIDSS
#SBATCH --output=/storage/homefs/kw23y068/logfiles/GRIDSS_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/GRIDSS_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=9
#SBATCH --mem=32G
#SBATCH --partition=epyc2

set -euo pipefail

# Load R module
module load R/4.4.2-gfbf-2024a

# Get R installation paths
RSCRIPT_PATH=$(which Rscript)
R_BIN_DIR=$(dirname "${RSCRIPT_PATH}")
R_BASE=$(dirname "${R_BIN_DIR}")

echo "R base directory: ${R_BASE}"
echo "R bin directory: ${R_BIN_DIR}"
echo "Rscript path: ${RSCRIPT_PATH}"

SAMPLE_NAME="FG_CC_19T_031"

# Define paths
GRIDSS_IMAGE=/storage/homefs/kw23y068/software/gridss/Patched_GRIDSS.sif
PROJECT_DIR=/storage/research/iee_evol/Korbi
REFERENCE=${PROJECT_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna
RUN_DIR=/storage/scratch/iee_evol/kw23y068/Gridss/${SAMPLE_NAME}

# Create run directory if it doesn't exist
if [ ! -d "${RUN_DIR}" ]; then
    mkdir -p "${RUN_DIR}"
fi

# Input files
BAM=${PROJECT_DIR}/bams_real/${SAMPLE_NAME}.fixmate.coordsorted.bam

# Verify inputs exist
echo "Checking input files..."
[ -f "${GRIDSS_IMAGE}" ] || { echo "ERROR: GRIDSS image not found: ${GRIDSS_IMAGE}"; exit 1; }
[ -f "${REFERENCE}" ] || { echo "ERROR: Reference not found: ${REFERENCE}"; exit 1; }
[ -f "${BAM}" ] || { echo "ERROR: BAM file not found: ${BAM}"; exit 1; }
echo "All input files found."
echo ""

# Run GRIDSS with R bound into the container
echo "Starting GRIDSS..."
apptainer exec \
    --bind ${RSCRIPT_PATH}:/usr/bin/Rscript \
    --bind ${R_BASE}/lib/R:/usr/lib/R \
    --bind ${PROJECT_DIR} \
    --bind ${RUN_DIR} \
    ${GRIDSS_IMAGE} \
    /opt/gridss/gridss \
        --reference ${REFERENCE} \
        --output ${RUN_DIR}/${SAMPLE_NAME}.vcf.gz \
        --threads ${SLURM_CPUS_PER_TASK} \
        ${BAM}

echo "GRIDSS complete at $(date)"