#!/bin/bash
#SBATCH --job-name=smoove_merge
#SBATCH --output=/storage/homefs/kw23y068/logfiles/smoove_merge_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/smoove_merge_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
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

# We follow population call from the smoove GitHub up to step 2:  https://github.com/brentp/smoove
# So we first run per sample, then merge # one may genotype

# Define paths
INPUT_DIR=/storage/scratch/iee_evol/kw23y068/smoove_single
OUTPUT_DIR=/storage/scratch/iee_evol/kw23y068/smoove_${START_PAD}_${END_PAD}
PROJECT_DIR=/storage/research/iee_evol/Korbi
REFERENCE=${PROJECT_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna

# Get samples
SAMPLE_LIST=/storage/research/iee_evol/Korbi/output/bams_full.csv
VCF_FILES=()
MISSING_VCF=0
#/storage/scratch/iee_evol/kw23y068/smoove_single/FG_CC_19T_031/FG_CC_19T_031-smoove.genotyped.vcf.gz
# Loop through samples in the batch
for i in $(seq ${BATCH_START} ${BATCH_END}); do
    SAMPLE_PATH=$(sed -n "${i}p" "${SAMPLE_LIST}")
    SAMPLE_FILE=$(basename "${SAMPLE_PATH}")
    SAMPLE_NAME="${SAMPLE_FILE%.fixmate.coordsorted.bam}"

    # Add warning if VCF file is missing
    VCF="${INPUT_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}-smoove.genotyped.vcf.gz"

        if [ ! -f "${VCF}" ]; then
        echo "WARNING: VCF file not found: ${VCF}"
        MISSING_VCF=$((MISSING_VCF + 1))
        continue
    fi

    # Append to VCF file list
        VCF_FILES+=("${VCF}")
    echo "  [${i}] ${SAMPLE_NAME}"
done

#Activate smoove
module load Anaconda3/2024.02-1
source $(conda info --base)/etc/profile.d/conda.sh
conda activate smoove

TMP_DIR=${OUTPUT_DIR}/tmp
mkdir -p ${TMP_DIR}
export TMPDIR=${TMP_DIR}

# Change to output directory
cd ${OUTPUT_DIR}

smoove merge \
    --name smoove_${START_PAD}_${END_PAD}_merged \
    --fasta ${REFERENCE} \
    --outdir ${OUTPUT_DIR} \
    "${VCF_FILES[@]}"

conda deactivate
