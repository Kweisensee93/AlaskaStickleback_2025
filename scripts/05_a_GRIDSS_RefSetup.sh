#!/bin/bash
#SBATCH --job-name=GRIDSS_refsetup
#SBATCH --output=/storage/homefs/kw23y068/logfiles/GRIDSS_refsetup_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/GRIDSS_refsetup_%j.err
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=epyc2

# This script sets up the reference genome for GRIDSS.
# Only needs to be run ONCE before any preprocessing or calling.

module load Anaconda3/2024.02-1
source $(conda info --base)/etc/profile.d/conda.sh
conda activate gridss

# Paths
PROJECT_DIR=/storage/research/iee_evol/Korbi
REFERENCE=${PROJECT_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna
WORK_DIR=/storage/scratch/iee_evol/kw23y068/Gridss_joint_reference_setup

mkdir -p "${WORK_DIR}"

# Clean up any potential leftover lock files
rm -f "${REFERENCE}"*.tmp.gridsslock

# Check if reference is already set up, otherwise build it
if [ -f "${REFERENCE}.gridsscache" ]; then
    echo "Reference already set up. Skipping setupreference step."
else
    echo "Setting up GRIDSS reference..."
    gridss -s setupreference -r "${REFERENCE}" -w "${WORK_DIR}"
fi

SETUP_EXIT_CODE=$?

echo ""
echo "========================================"
echo "Finished reference setup"
echo "Date: $(date)"
echo "========================================"

conda deactivate

if [ ${SETUP_EXIT_CODE} -eq 0 ]; then
    echo "SUCCESS: Reference setup completed"
    echo "Reference is stored as: ${REFERENCE}.gridsscache"
else
    echo "ERROR: Reference setup failed"
    exit ${SETUP_EXIT_CODE}
fi