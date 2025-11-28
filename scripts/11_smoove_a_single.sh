#!/bin/bash
#SBATCH --job-name=smoove_single
#SBATCH --array=1-80
#SBATCH --output=/storage/homefs/kw23y068/logfiles/smoove_single_%A_%a.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/smoove_single_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=epyc2

# The puffin reference could use joint-calling in one go, since they only have 18 puffins to analyze
# We need to split this to several substeps, in accordance with the smoove GitHub:
# single-call, merge, genotype (we leave out the steps afterwards)
# https://github.com/brentp/smoove

echo "==========================================  "
echo "Single-sample smoove calling"
echo "Sample: ${SLURM_ARRAY_TASK_ID}"
echo "Date: $(date)"
echo "=========================================="

# Define paths
SAMPLE_LIST=/storage/research/iee_evol/Korbi/output/bams_full.csv
SAMPLE_PATH=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")
SAMPLE_FILE=$(basename "${SAMPLE_PATH}")
SAMPLE_NAME="${SAMPLE_FILE%.fixmate.coordsorted.bam}"

PROJECT_DIR=/storage/research/iee_evol/Korbi
REFERENCE=${PROJECT_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna
OUTPUT_DIR=/storage/scratch/iee_evol/kw23y068/smoove_single/${SAMPLE_NAME}
mkdir -p ${OUTPUT_DIR}

# Check if BAM exists
if [ ! -f "${SAMPLE_PATH}" ]; then
    echo "ERROR: BAM file not found: ${SAMPLE_PATH}"
    exit 1
fi

# Load conda
module load Anaconda3/2024.02-1
source $(conda info --base)/etc/profile.d/conda.sh
conda activate smoove

# Setup temp directory
TMP_DIR=${OUTPUT_DIR}/tmp_${SAMPLE_NAME}
mkdir -p ${TMP_DIR}
export TMPDIR=${TMP_DIR}

# Change to output directory
cd ${OUTPUT_DIR}

echo "Processing ${SAMPLE_NAME}..."

# Run smoove call for single sample
# Remove -x flag to keep all chromosomes initially
smoove call \
    --name ${SAMPLE_NAME} \
    --fasta ${REFERENCE} \
    -p ${SLURM_CPUS_PER_TASK} \
    --genotype \
    --duphold \
    ${SAMPLE_PATH}

EXIT_CODE=$?

if [ ${EXIT_CODE} -ne 0 ]; then
    echo "ERROR: smoove failed for ${SAMPLE_NAME} with exit code ${EXIT_CODE}"
fi

# Cleanup
rm -rf ${TMP_DIR}
conda deactivate

echo "Complete for ${SAMPLE_NAME}"
echo "Date: $(date)"
