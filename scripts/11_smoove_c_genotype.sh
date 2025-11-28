#!/bin/bash
#SBATCH --job-name=smoove_genotype
#SBATCH --output=/storage/homefs/kw23y068/logfiles/smoove_genotype_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/smoove_genotype_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
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

# Define paths
PROJECT_DIR=/storage/research/iee_evol/Korbi
REFERENCE=${PROJECT_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna
INPUT_DIR=/storage/scratch/iee_evol/kw23y068/smoove_single
OUTPUT_DIR=/storage/scratch/iee_evol/kw23y068/smoove_${START_PAD}_${END_PAD}
MERGED_VCF=${OUTPUT_DIR}/smoove_${START_PAD}_${END_PAD}_merged.sites.vcf.gz

# Get samples
SAMPLE_LIST=/storage/research/iee_evol/Korbi/output/bams_full.csv
BAM_FILES=()

# Loop through samples in the batch
for i in $(seq ${BATCH_START} ${BATCH_END}); do
    SAMPLE_PATH=$(sed -n "${i}p" "${SAMPLE_LIST}")
    SAMPLE_FILE=$(basename "${SAMPLE_PATH}")
    SAMPLE_NAME="${SAMPLE_FILE%.fixmate.coordsorted.bam}"

    # Check if VCF was successfully created (should match what was merged)
    # Only add BAMS where smoove ran successfully
    VCF="${INPUT_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}-smoove.genotyped.vcf.gz"
    if [[ -f "${VCF}" ]]; then
        BAM_FILES+=("${SAMPLE_PATH}")
        echo "  [${i}] ${SAMPLE_NAME} (VCF found — added)"
    else
        echo "  [${i}] ${SAMPLE_NAME} (VCF missing — skipped)"
    fi

done


# Load conda
module load Anaconda3/2024.02-1
source $(conda info --base)/etc/profile.d/conda.sh
conda activate smoove

# Setup temp directory
TMP_DIR=${OUTPUT_DIR}/tmp
mkdir -p ${TMP_DIR}
export TMPDIR=${TMP_DIR}

# Change to output directory
cd ${OUTPUT_DIR}

# We leave out the following options from the smoove workflow
# -d -x 

smoove genotype \
    -p ${SLURM_CPUS_PER_TASK} \
    --name smoove_merged_genotyped_${START_PAD}_${END_PAD} \
    --outdir ${OUTPUT_DIR} \
    --fasta ${REFERENCE} \
    --vcf ${MERGED_VCF} \
    --duphold \
    "${BAM_FILES[@]}"

conda deactivate
