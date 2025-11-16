#!/bin/bash
#SBATCH --job-name=Lumpy_joint_preprocess
#SBATCH --output=/storage/homefs/kw23y068/logfiles/Lumpy_joint_preprocess_%A_%a.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/Lumpy_joint_preprocess_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --array=1-80
#SBATCH --partition=epyc2


####################
# PARAMETERS TO SET
####################
# May be altered for feeding in arguments to script
# Define which batch of samples to process
BATCH_START=1    # Change to 81, 161, 241, 321, 401 for other batches
BATCH_END=80     # Change to 160, 240, 320, 400, 480 for other batches

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
REFERENCE=${PROJECT_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna

START_PAD=$(printf "%03d" "${BATCH_START}")
END_PAD=$(printf "%03d" "${BATCH_END}")
OUT_DIR="/storage/scratch/iee_evol/kw23y068/lumpy_joint_${START_PAD}_${END_PAD}/"
if [ ! -d "${OUT_DIR}" ]; then
    mkdir -p "${OUT_DIR}"
fi


SAMPLE_LIST=/storage/research/iee_evol/Korbi/output/bams_full.csv
# Get sample
# Calculate actual sample index (offset by batch start)
SAMPLE_INDEX=$((BATCH_START + SLURM_ARRAY_TASK_ID - 1))
# Sample path only works, if we have only one sample per line in the csv!
SAMPLE_PATH=$(sed -n "${SAMPLE_INDEX}p" "${SAMPLE_LIST}")
SAMPLE_FILE=$(basename "${SAMPLE_PATH}")
SAMPLE_NAME="${SAMPLE_FILE%.fixmate.coordsorted.bam}"
BAM=${PROJECT_DIR}/bams_real/${SAMPLE_NAME}.fixmate.coordsorted.bam

# Set up run directory for sample
RUN_DIR="${OUT_DIR}/${SAMPLE_NAME}"
if [ ! -d "$RUN_DIR" ]; then
    mkdir -p "$RUN_DIR"
fi

# Sample preparation loop
# Extract discordant & split-read alignments
DISCORDANT="${RUN_DIR}/${SAMPLE_NAME}.discordants.unsorted.bam"
SPLITREAD="${RUN_DIR}/${SAMPLE_NAME}.splitters.unsorted.bam"

# Extract the discordant paired-end alignments.
samtools view -b -F 1294 ${BAM} > ${DISCORDANT}
# Extract the split-read alignments
samtools view -h ${BAM} \
    | extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > ${SPLITREAD}


# Sort both alignments
#samtools sort ${DISCORDANT} ${SAMPLE_NAME}.discordants
#samtools sort ${SPLITREAD} ${SAMPLE_NAME}.splitters

# Sort both alignments
samtools sort -o "${RUN_DIR}/${SAMPLE_NAME}.discordants.bam" "${DISCORDANT}"
samtools sort -o "${RUN_DIR}/${SAMPLE_NAME}.splitters.bam" "${SPLITREAD}"


CALL_EXIT_CODE=$?

echo ""
echo "========================================"
echo "Finished  Lumpy preprocessing for sample: ${SAMPLE_NAME}"
echo "Exit code: ${CALL_EXIT_CODE}"
echo "Date: $(date)"
echo "========================================"

conda deactivate

# Get job statistics
echo ""
echo "========================================"
echo "Resource Usage Summary:"
echo "========================================"
sacct -j ${SLURM_JOB_ID} --format=JobID,MaxRSS,MaxVMSize,Elapsed,CPUTime,TotalCPU -P
