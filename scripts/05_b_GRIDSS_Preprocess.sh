#!/bin/bash
#SBATCH --job-name=GRIDSS_prep
#SBATCH --output=/storage/homefs/kw23y068/logfiles/GRIDSS_prep_%A_%a.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/GRIDSS_prep_%A_%a.err
#SBATCH --array=1-80
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=epyc2

# According to GRIDSS documentation, it is optimized for 8 CPU threads and 32GB RAM - lets see

####################
# PARAMETERS TO SET
####################
# May be altered for feeding in arguments to script
# Define which batch of samples to process
BATCH_START=1    # Change to 81, 161, 241, 321, 401 for other batches
BATCH_END=80     # Change to 160, 240, 320, 400, 480 for other batches


# This script is seting up the reference and preprocessing BAM files for GRIDSS SV calling.
# This is a sub-step of 05_GRIDSS.sh

echo "========================================"
echo "GRIDSS Preprocessing"
echo "Batch range: ${BATCH_START}-${BATCH_END}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Date: $(date)"
echo "========================================"

# Load GRIDSS environment
module load Anaconda3/2024.02-1
source $(conda info --base)/etc/profile.d/conda.sh
conda activate gridss


# Paths
PROJECT_DIR=/storage/research/iee_evol/Korbi
REFERENCE=${PROJECT_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna
SAMPLE_LIST=/storage/research/iee_evol/Korbi/output/bams_full.csv

# Set up run directory based on batch; with 0 padding to 3 digits
START_PAD=$(printf "%03d" "${BATCH_START}")
END_PAD=$(printf "%03d" "${BATCH_END}")
RUN_DIR="/storage/scratch/iee_evol/kw23y068/Gridss_joint_${START_PAD}_${END_PAD}/preprocessed"
mkdir -p "${RUN_DIR}"

#RUN_DIR=/storage/scratch/iee_evol/kw23y068/Gridss_joint_001_080/preprocessed
#RUN_DIR=/storage/scratch/iee_evol/kw23y068/Gridss_joint_081_160/preprocessed
#RUN_DIR=/storage/scratch/iee_evol/kw23y068/Gridss_joint_161_240/preprocessed
#RUN_DIR=/storage/scratch/iee_evol/kw23y068/Gridss_joint_241_320/preprocessed
#RUN_DIR=/storage/scratch/iee_evol/kw23y068/Gridss_joint_321_400/preprocessed
#RUN_DIR=/storage/scratch/iee_evol/kw23y068/Gridss_joint_401_480/preprocessed

# Get sample
# Calculate actual sample index (offset by batch start)
SAMPLE_INDEX=$((BATCH_START + SLURM_ARRAY_TASK_ID - 1))

# Sample path only works, if we have only one sample per line in the csv!
SAMPLE_PATH=$(sed -n "${SAMPLE_INDEX}p" "${SAMPLE_LIST}")
SAMPLE_FILE=$(basename "${SAMPLE_PATH}")
SAMPLE_NAME="${SAMPLE_FILE%.fixmate.coordsorted.bam}"
BAM=${PROJECT_DIR}/bams_real/${SAMPLE_NAME}.fixmate.coordsorted.bam

# Verify BAM exists
if [ ! -f "${BAM}" ]; then
    echo "ERROR: BAM file not found: ${BAM}"
    exit 1
fi

# Step 2: Preprocess
gridss -s preprocess -r ${REFERENCE} -w "${RUN_DIR}" "${BAM}"

PREPROCESS_EXIT_CODE=$?

echo ""
echo "========================================"
echo "Finished preprocessing"
echo "Sample: ${SAMPLE_NAME}"
echo "Exit code: ${PREPROCESS_EXIT_CODE}"
echo "Date: $(date)"
echo "========================================"

conda deactivate

if [ ${PREPROCESS_EXIT_CODE} -eq 0 ]; then
    echo "SUCCESS: Preprocessing completed for ${SAMPLE_NAME}"
else
    echo "ERROR: Preprocessing failed for ${SAMPLE_NAME}"
    exit ${PREPROCESS_EXIT_CODE}
fi