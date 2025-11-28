#!/bin/bash
#SBATCH --job-name=smoove
#SBATCH --output=/storage/homefs/kw23y068/logfiles/smoove_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/smoove_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=19
#SBATCH --mem=128G
#SBATCH --partition=epyc2

####################
# PARAMETERS TO SET
####################
# May be altered for feeding in arguments to script
# Define which batch of samples to process
BATCH_START=1    # Change to 81, 161, 241, 321, 401 for other batches
BATCH_END=80     # Change to 160, 240, 320, 400, 480 for other batches

START_PAD=$(printf "%03d" "${BATCH_START}")
END_PAD=$(printf "%03d" "${BATCH_END}")

echo "SV filtering with smoove start"
echo "Processing samples ${BATCH_START} to ${BATCH_END}"
echo "Date: $(date)"

# Set paths
PROJECT_DIR=/storage/research/iee_evol/Korbi
REFERENCE=${PROJECT_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna
OUTPUT_DIR=/storage/scratch/iee_evol/kw23y068/smoove
OUT_DIR=${OUTPUT_DIR}/smoove_${START_PAD}_${END_PAD}
mkdir -p ${OUT_DIR}

# Set files
SAMPLE_LIST=/storage/research/iee_evol/Korbi/output/bams_full.csv
# Loop through samples in the batch
BAM_FILES=()
SAMPLE_NAMES=()

for i in $(seq ${BATCH_START} ${BATCH_END}); do
    SAMPLE_PATH=$(sed -n "${i}p" "${SAMPLE_LIST}")
    
    if [ -z "${SAMPLE_PATH}" ]; then
        echo "WARNING: Could not read sample at line ${i}"
        continue
    fi
    
    SAMPLE_FILE=$(basename "${SAMPLE_PATH}")
    SAMPLE_NAME="${SAMPLE_FILE%.fixmate.coordsorted.bam}"
    
    BAM="${PROJECT_DIR}/bams_real/${SAMPLE_NAME}.fixmate.coordsorted.bam"
    
    BAM_FILES+=("${BAM}")
    SAMPLE_NAMES+=("${SAMPLE_NAME}")
    echo "  [${i}] ${SAMPLE_NAME}"
done

module load Anaconda3/2024.02-1
# Ensure conda commands work in batch mode
source $(conda info --base)/etc/profile.d/conda.sh
conda activate smoove

# Setup temporary directory
TMP_DIR=${OUT_DIR}/tmp
mkdir -p ${TMP_DIR}
export TMPDIR=${TMP_DIR}

# smoove writes to current directory
cd ${OUT_DIR}

# They use the -x option to filter out sex chromosomes; We leave it out for now
smoove call \
    --name smoove_${START_PAD}_${END_PAD} \
    --fasta ${REFERENCE} \
    -p 18 \
    --genotype \
    --duphold \
    "${BAM_FILES[@]}"

rm -rf ${TMP_DIR}

conda deactivate
