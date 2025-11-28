#!/bin/bash
#SBATCH --job-name=manta_vcf_check
#SBATCH --output=/storage/homefs/kw23y068/logfiles/manta_vcf_check_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/manta_vcf_check_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=epyc2

BATCH_START=1
BATCH_END=80

INPUT_DIR=/storage/scratch/iee_evol/kw23y068/Manta
SAMPLE_LIST=/storage/research/iee_evol/Korbi/output/bams_full.csv
OUTPUT_FILE=/storage/homefs/kw23y068/logfiles/manta_vcf_validation.txt

module load BCFtools/1.12-GCC-10.3.0

echo "Checking Manta VCF files..." | tee ${OUTPUT_FILE}
echo "" | tee -a ${OUTPUT_FILE}

VALID=0
INVALID=0
EMPTY=0
MISSING=0

for i in $(seq ${BATCH_START} ${BATCH_END}); do
    SAMPLE_PATH=$(sed -n "${i}p" "${SAMPLE_LIST}")
    SAMPLE_FILE=$(basename "${SAMPLE_PATH}")
    SAMPLE_NAME="${SAMPLE_FILE%.fixmate.coordsorted.bam}"
    
    VCF="${INPUT_DIR}/${SAMPLE_NAME}/results/variants/diploidSV.vcf.gz"
    
    if [[ ! -f "${VCF}" ]]; then
        echo "[${i}] ${SAMPLE_NAME}: MISSING" | tee -a ${OUTPUT_FILE}
        MISSING=$((MISSING + 1))
        continue
    fi
    
    # Check if file is empty or corrupted
    NUM_VARIANTS=$(zcat ${VCF} 2>/dev/null | grep -v "^#" | wc -l)
    
    if [ $? -ne 0 ]; then
        echo "[${i}] ${SAMPLE_NAME}: CORRUPTED (can't read)" | tee -a ${OUTPUT_FILE}
        INVALID=$((INVALID + 1))
        continue
    fi
    
    if [ ${NUM_VARIANTS} -eq 0 ]; then
        echo "[${i}] ${SAMPLE_NAME}: EMPTY (0 variants)" | tee -a ${OUTPUT_FILE}
        EMPTY=$((EMPTY + 1))
        continue
    fi
    
    # Quick validation with bcftools
    bcftools view -h ${VCF} > /dev/null 2>&1
    
    if [ $? -ne 0 ]; then
        echo "[${i}] ${SAMPLE_NAME}: INVALID VCF format" | tee -a ${OUTPUT_FILE}
        INVALID=$((INVALID + 1))
        continue
    fi
    
    echo "[${i}] ${SAMPLE_NAME}: OK (${NUM_VARIANTS} variants)" | tee -a ${OUTPUT_FILE}
    VALID=$((VALID + 1))
done

echo "" | tee -a ${OUTPUT_FILE}
echo "========================================" | tee -a ${OUTPUT_FILE}
echo "Summary:" | tee -a ${OUTPUT_FILE}
echo "  Valid VCFs: ${VALID}" | tee -a ${OUTPUT_FILE}
echo "  Empty VCFs: ${EMPTY}" | tee -a ${OUTPUT_FILE}
echo "  Invalid/Corrupted: ${INVALID}" | tee -a ${OUTPUT_FILE}
echo "  Missing: ${MISSING}" | tee -a ${OUTPUT_FILE}
echo "========================================" | tee -a ${OUTPUT_FILE}

if [ ${INVALID} -gt 0 ]; then
    echo "" | tee -a ${OUTPUT_FILE}
    echo "WARNING: Found ${INVALID} invalid/corrupted VCF files!" | tee -a ${OUTPUT_FILE}
    echo "These must be excluded from merging." | tee -a ${OUTPUT_FILE}
fi

if [ ${EMPTY} -gt 0 ]; then
    echo "" | tee -a ${OUTPUT_FILE}
    echo "Note: Found ${EMPTY} empty VCF files (no variants called)." | tee -a ${OUTPUT_FILE}
    echo "Consider excluding these from merging." | tee -a ${OUTPUT_FILE}
fi