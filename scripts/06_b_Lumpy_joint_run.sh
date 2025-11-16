#!/bin/bash
#SBATCH --job-name=Lumpy_joint_run
#SBATCH --output=/storage/homefs/kw23y068/logfiles/Lumpy_joint_run_%A_%a.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/Lumpy_joint_run_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G
#SBATCH --partition=epyc2


####################
# PARAMETERS TO SET
####################
# May be altered for feeding in arguments to script
# Define which batch of samples to process
BATCH_START=1    # Change to 81, 161, 241, 321, 401 for other batches
BATCH_END=81     # Change to 160, 240, 320, 400, 480 for other batches

echo "========================================"
echo "Lumpy Joint Preprocessing"
echo "Batch range: ${BATCH_START}-${BATCH_END}"
echo "Date: $(date)"
echo "========================================"

# Activate Lumpy environment
module load Anaconda3/2024.02-1
source $(conda info --base)/etc/profile.d/conda.sh
conda activate lumpy

#paths
PROJECT_DIR=/storage/research/iee_evol/Korbi
#REFERENCE=${PROJECT_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna

START_PAD=$(printf "%03d" "${BATCH_START}")
END_PAD=$(printf "%03d" "${BATCH_END}")
OUT_DIR="/storage/scratch/iee_evol/kw23y068/lumpy_joint_${START_PAD}_${END_PAD}/"
if [ ! -d "${OUT_DIR}" ]; then
    echo "${OUT_DIR} does not exist. Please run preprocessing first."
    exit 1
fi

# Build input files:
SAMPLE_LIST=/storage/research/iee_evol/Korbi/output/bams_full.csv
# Initialize empty strings for BAMs, splitters, discordants:
BAM_FILES=""
SPLIT_FILES=""
DISC_FILES=""
# Loop through samples in the batch to collect files
for (( IDX=${BATCH_START}; IDX<=${BATCH_END}; IDX++ )); do

    SAMPLE_PATH=$(sed -n "${IDX}p" "${SAMPLE_LIST}")
    SAMPLE_FILE=$(basename "${SAMPLE_PATH}")
    SAMPLE_NAME="${SAMPLE_FILE%.fixmate.coordsorted.bam}"

    RUN_DIR="${OUT_DIR}/${SAMPLE_NAME}"
    if [ ! -d "$RUN_DIR" ]; then
        echo "ERROR: Missing preprocessed run directory for ${SAMPLE_NAME}"
        echo "Expected: $RUN_DIR"
        exit 1
    fi

    BAM="${PROJECT_DIR}/bams_real/${SAMPLE_NAME}.fixmate.coordsorted.bam"
    SPLIT="${RUN_DIR}/${SAMPLE_NAME}.splitters.bam"
    DISC="${RUN_DIR}/${SAMPLE_NAME}.discordants.bam"

    if [[ ! -f "$SPLIT" || ! -f "$DISC" ]]; then
        echo "ERROR: Missing preprocessed file for ${SAMPLE_NAME}"
        echo "Expected:"
        echo "  $SPLIT"
        echo "  $DISC"
        exit 1
    fi

    BAM_FILES+="${BAM},"
    SPLIT_FILES+="${SPLIT},"
    DISC_FILES+="${DISC},"

done

# Set output VCF path
OUTPUT_VCF="${OUT_DIR}/lumpy_joint_${START_PAD}_${END_PAD}.vcf"

# Remove trailing commas
BAM_FILES="${BAM_FILES%,}"
SPLIT_FILES="${SPLIT_FILES%,}"
DISC_FILES="${DISC_FILES%,}"

echo "Running LUMPY joint call"

lumpyexpress \
    -B "${BAM_FILES}" \
    -S "${SPLIT_FILES}" \
    -D "${DISC_FILES}" \
    -o "${OUTPUT_VCF}"

echo "Done. Output written to: ${OUTPUT_VCF}"

# Compress VCF
# the conda environment has bgzip installed
bgzip -c "${OUTPUT_VCF}" > "${OUTPUT_VCF}.gz"

conda deactivate
# No BFCtools in lumpy env, so load it separately, afte conda deactivate

module load BCFtools/1.12-GCC-10.3.0

bcftools index "${OUTPUT_VCF}.gz"

# Get job statistics
echo ""
echo "========================================"
echo "Resource Usage Summary:"
echo "========================================"
sacct -j ${SLURM_JOB_ID} --format=JobID,MaxRSS,MaxVMSize,Elapsed,CPUTime,TotalCPU -P
