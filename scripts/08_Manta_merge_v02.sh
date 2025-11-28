#!/bin/bash
#SBATCH --job-name=manta_merge_bcftools
#SBATCH --output=/storage/homefs/kw23y068/logfiles/manta_merge_bcftools_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/manta_merge_bcftools_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=epyc2

####################
# PARAMETERS
####################
BATCH_START=1
BATCH_END=80
START_PAD=$(printf "%03d" "${BATCH_START}")
END_PAD=$(printf "%03d" "${BATCH_END}")

# Paths - adjust to your filtered Manta directory
INPUT_DIR=/storage/scratch/iee_evol/kw23y068/Manta_filtered_${START_PAD}_${END_PAD}
OUTPUT_DIR=/storage/scratch/iee_evol/kw23y068/manta_${START_PAD}_${END_PAD}_merged
SAMPLE_LIST=/storage/research/iee_evol/Korbi/output/bams_full.csv
MERGED_VCF=${OUTPUT_DIR}/manta_${START_PAD}_${END_PAD}_merged.vcf.gz

mkdir -p ${OUTPUT_DIR}

# Load bcftools
module load HTSlib/1.12-GCC-10.3.0
module load BCFtools/1.12-GCC-10.3.0

VCF_FILES=()
MISSING_VCF=0
EMPTY_VCF=0

echo "Collecting filtered Manta VCF files for batch ${BATCH_START}-${BATCH_END}..."
echo ""

# Loop through samples
for i in $(seq ${BATCH_START} ${BATCH_END}); do
    SAMPLE_PATH=$(sed -n "${i}p" "${SAMPLE_LIST}")
    SAMPLE_FILE=$(basename "${SAMPLE_PATH}")
    SAMPLE_NAME="${SAMPLE_FILE%.fixmate.coordsorted.bam}"
    
    VCF="${INPUT_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_manta_filtered.vcf.gz"
    
    if [[ ! -f "${VCF}" ]]; then
        echo "  [${i}] ${SAMPLE_NAME} - VCF not found"
        MISSING_VCF=$((MISSING_VCF + 1))
        continue
    fi
    
    # Check if VCF has variants
    NUM_VARIANTS=$(zcat ${VCF} 2>/dev/null | grep -v "^#" | wc -l)
    
    if [ $? -ne 0 ]; then
        echo " [${i}] ${SAMPLE_NAME} - VCF corrupted/unreadable"
        MISSING_VCF=$((MISSING_VCF + 1))
        continue
    fi
    
    if [ ${NUM_VARIANTS} -eq 0 ]; then
        echo " [${i}] ${SAMPLE_NAME} - VCF empty (0 variants)"
        EMPTY_VCF=$((EMPTY_VCF + 1))
        continue
    fi
    
    # Ensure VCF is indexed
    if [[ ! -f "${VCF}.tbi" ]]; then
        echo "  [${i}] ${SAMPLE_NAME} - Indexing failed"
        MISSING_VCF=$((MISSING_VCF + 1))
        continue
    fi
    
    VCF_FILES+=("${VCF}")
    echo " [${i}] ${SAMPLE_NAME} - OK (${NUM_VARIANTS} variants)"
done

echo ""
echo "=========================================="
echo "Summary:"
echo "  Total samples: $((BATCH_END - BATCH_START + 1))"
echo "  Valid VCFs: ${#VCF_FILES[@]}"
echo "  Empty VCFs: ${EMPTY_VCF}"
echo "  Missing/Corrupted: ${MISSING_VCF}"
echo "=========================================="
echo ""

# Check if we have VCFs to merge
if [ ${#VCF_FILES[@]} -eq 0 ]; then
    echo "ERROR: No valid VCF files found!"
    exit 1
fi

echo "Merging ${#VCF_FILES[@]} VCFs with bcftools..."
echo ""

# bcftools merge options:
# --merge none: creates separate records for variants at same position (conservative)
# --merge all: merges variants at same position (more aggressive)
# For SVs, 'none' is usually safer

bcftools merge \
    --threads ${SLURM_CPUS_PER_TASK} \
    --merge none \
    --output-type z \
    --output ${MERGED_VCF} \
    "${VCF_FILES[@]}"

MERGE_EXIT=$?

if [ ${MERGE_EXIT} -ne 0 ]; then
    echo ""
    echo "ERROR: bcftools merge failed with exit code ${MERGE_EXIT}"
    exit ${MERGE_EXIT}
fi

echo ""
echo "Indexing merged VCF..."
tabix -p vcf ${MERGED_VCF}

if [ $? -ne 0 ]; then
    echo "ERROR: Indexing failed!"
    exit 1
fi


echo ""
echo "Input statistics:"
echo "  Samples processed: $((BATCH_END - BATCH_START + 1))"
echo "  VCFs merged: ${#VCF_FILES[@]}"
echo "  Empty VCFs: ${EMPTY_VCF}"
echo "  Missing/Corrupted: ${MISSING_VCF}"
echo ""
echo "Output files:"
echo "  Merged VCF: ${MERGED_VCF}"
echo "  Index: ${MERGED_VCF}.tbi"
echo ""
