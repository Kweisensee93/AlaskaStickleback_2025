#!/bin/bash
#SBATCH --job-name=manta_merge_survivor
#SBATCH --output=/storage/homefs/kw23y068/logfiles/manta_merge_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/manta_merge_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=epyc2

####################
# PARAMETERS
####################
BATCH_START=1
BATCH_END=80
START_PAD=$(printf "%03d" "${BATCH_START}")
END_PAD=$(printf "%03d" "${BATCH_END}")

# Paths
INPUT_DIR=/storage/scratch/iee_evol/kw23y068/Manta
OUTPUT_DIR=/storage/scratch/iee_evol/kw23y068/manta_${START_PAD}_${END_PAD}_merged
SAMPLE_LIST=/storage/research/iee_evol/Korbi/output/bams_full.csv
MERGED_VCF=${OUTPUT_DIR}/manta_${START_PAD}_${END_PAD}_merged.vcf

mkdir -p ${OUTPUT_DIR}

# Create file to store VCF paths (SURVIVOR needs this as input)
VCF_LIST=${OUTPUT_DIR}/vcf_files.txt
> ${VCF_LIST}  # Clear/create empty file

MISSING_VCF=0
FOUND_VCF=0
#/storage/scratch/iee_evol/kw23y068/Manta/FG_CC_19T_031/results/variants/diploidSV.vcf.gz
# Loop through samples in the batch
for i in $(seq ${BATCH_START} ${BATCH_END}); do
    SAMPLE_PATH=$(sed -n "${i}p" "${SAMPLE_LIST}")
    SAMPLE_FILE=$(basename "${SAMPLE_PATH}")
    SAMPLE_NAME="${SAMPLE_FILE%.fixmate.coordsorted.bam}"

    # Check if VCF was successfully created (should match what was merged)
    # Only add VCFs where manta ran successfully
    VCF="${INPUT_DIR}/${SAMPLE_NAME}/results/variants/diploidSV.vcf.gz"
    if [[ -f "${VCF}" ]]; then
        NUM_VARIANTS=$(zcat ${VCF} 2>/dev/null | grep -v "^#" | wc -l)
        if [[ "${NUM_VARIANTS}" -eq 0 ]]; then
            echo "  [${i}] ${SAMPLE_NAME} (VCF empty — skipped)"
            MISSING_VCF=$((MISSING_VCF + 1))
            continue
        else
            echo "${VCF}" >> ${VCF_LIST}  # Add to file
            echo "  [${i}] ${SAMPLE_NAME} (VCF found — added)"
            FOUND_VCF=$((FOUND_VCF + 1))
        fi

    else
        echo "  [${i}] ${SAMPLE_NAME} (VCF missing — skipped)"
        MISSING_VCF=$((MISSING_VCF + 1))
    fi

done

echo "Total samples in batch: $((BATCH_END - BATCH_START + 1))"
echo "VCFs found: ${FOUND_VCF}"
echo "VCFs missing: ${MISSING_VCF}"

if [ ${FOUND_VCF} -eq 0 ]; then
    echo "ERROR: No VCF files found to merge!"
    exit 1
fi

module load Anaconda3/2024.02-1
source $(conda info --base)/etc/profile.d/conda.sh
conda activate survivor

# SURVIVOR merge parameters:
# vcf_list: file with list of VCFs
# max_distance: maximum distance between breakpoints (1000 = 1kb)
# min_support: minimum number of supporting caller (1 = keep all)
# use_type: take SV type into account (1=yes, 0=no)
# use_strand: take strand into account (1=yes, 0=no)
# estimate_distance: estimate distance (0=no, 1=yes)
# min_size: minimum SV size (-1 = all sizes)
# output_vcf: output filename

SURVIVOR merge \
    ${VCF_LIST} \
    1000 \
    1 \
    1 \
    1 \
    0 \
    50 \
    ${MERGED_VCF} 2>&1 | tee ${OUTPUT_DIR}/survivor_merge.log

conda deactivate

module -q reset
module load HTSlib/1.12-GCC-10.3.0
module load BCFtools/1.12-GCC-10.3.0
#module load VCFtools/0.1.16-GCC-10.3.0

# Compress and index
bgzip -c ${MERGED_VCF} > ${MERGED_VCF}.gz
tabix -p vcf ${MERGED_VCF}.gz

echo "Samples input: $((BATCH_END - BATCH_START + 1))"
echo "Missing samples: ${MISSING_VCF}"
echo "Merged VCFs: ${FOUND_VCF}"
