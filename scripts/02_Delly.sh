#!/bin/bash
#SBATCH --job-name=delly_call
#SBATCH --output=/storage/homefs/kw23y068/logfiles/delly_call_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/delly_call_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=epyc2

# Delly uses a Threshold for "Pass"
#Pass/LowQual is indeed related to paired end mapping. SVs are flagged as "PASS" if >=3 paired-ends support the variant (for translocations >=5). In addition, the mean mapping quality has to be >=20
# https://groups.google.com/g/delly-users/c/6Mq2juBraRY
# https://groups.google.com/g/delly-users/c/x1YAvpAhOAI

echo "=== Starting Delly SV calling ==="
echo "Date: $(date)"
echo "Host: $(hostname)"

# Load necessary modules
module load BCFtools/1.12-GCC-10.3.0

# Define paths
DELLY_SIF=/storage/homefs/kw23y068/software/Delly/delly_v1.5.0.sif
PROJECT_DIR=/storage/research/iee_evol/Korbi
REFERENCE=${PROJECT_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna
OUTPUT_DIR=/storage/scratch/iee_evol/kw23y068/Delly/FG_CC_19T_031_delly
# files
BAM=${PROJECT_DIR}/bams_real/FG_CC_19T_031.fixmate.coordsorted.bam
DELLY_VCF=${OUTPUT_DIR}/FG_CC_19T_031.vcf.gz
FILTERED_VCF=${OUTPUT_DIR}/FG_CC_19T_031.filtered.vcf.gz

LENGTH_THRESHOLD=1000000

if [ ! -d "${OUTPUT_DIR}" ]
then
    mkdir ${OUTPUT_DIR}
fi


# Run Delly for different SV types
echo ">>> Running Delly call for all SV types"
apptainer exec \
    --bind ${PROJECT_DIR} \
    --bind ${OUTPUT_DIR} \
    ${DELLY_SIF} delly call \
    -g ${REFERENCE} \
    -o ${OUTPUT_DIR}/FG_CC_19T_031.bcf \
    ${BAM}

# Convert BCF to VCF for easier viewing
echo ">>> Converting BCF to VCF"
bcftools view \
    ${OUTPUT_DIR}/FG_CC_19T_031.bcf \
    -O v \
    -o ${OUTPUT_DIR}/FG_CC_19T_031.vcf

# Compress and index VCF
echo ">>> Compressing and indexing VCF"
bgzip -c ${OUTPUT_DIR}/FG_CC_19T_031.vcf > ${OUTPUT_DIR}/FG_CC_19T_031.vcf.gz
bcftools index ${OUTPUT_DIR}/FG_CC_19T_031.vcf.gz

echo ">>> Filtering for PASS variants only and applying length threshold: ${LENGTH_THRESHOLD} bp"
FILTER_EXPR="FILTER=\"PASS\" && (abs(POS - INFO/END) <= ${LENGTH_THRESHOLD})"

bcftools view \
    -i "${FILTER_EXPR}" \
    -O z \
    -o "${FILTERED_VCF}" \
    "${DELLY_VCF}"

# Index the filtered VCF
echo ">>> Indexing filtered VCF"
bcftools index ${FILTERED_VCF}

# Generate basic statistics
# echo ">>> Generating statistics"
# echo "Total SVs called:" > ${OUTPUT_DIR}/stats.txt
# bcftools view -H ${OUTPUT_DIR}/FG_CC_19T_031.vcf.gz | wc -l >> ${OUTPUT_DIR}/stats.txt

# echo "" >> ${OUTPUT_DIR}/stats.txt
# echo "SV types:" >> ${OUTPUT_DIR}/stats.txt
# bcftools query -f '%SVTYPE\n' ${OUTPUT_DIR}/FG_CC_19T_031.vcf.gz | sort | uniq -c >> ${OUTPUT_DIR}/stats.txt

# echo "=== Delly analysis complete ==="
# echo "Results in: ${OUTPUT_DIR}"
# cat ${OUTPUT_DIR}/stats.txt
# echo "Date: $(date)"

# Generate statistics for filtered variants
# echo ">>> Generating statistics for PASS variants"
# echo "PASS SVs:" > ${OUTPUT_DIR}/stats_filtered.txt
# bcftools view -H ${FILTERED_VCF} | wc -l >> ${OUTPUT_DIR}/stats_filtered.txt

# echo "" >> ${OUTPUT_DIR}/stats_filtered.txt
# echo "PASS SV types:" >> ${OUTPUT_DIR}/stats_filtered.txt
# bcftools query -f '%SVTYPE\n' ${FILTERED_VCF} | sort | uniq -c >> ${OUTPUT_DIR}/stats_filtered.txt

# # Compare with original
# echo "" >> ${OUTPUT_DIR}/stats_filtered.txt
# echo "Comparison:" >> ${OUTPUT_DIR}/stats_filtered.txt
# echo "Total SVs: $(bcftools view -H ${INPUT_VCF} | wc -l)" >> ${OUTPUT_DIR}/stats_filtered.txt
# echo "PASS SVs: $(bcftools view -H ${FILTERED_VCF} | wc -l)" >> ${OUTPUT_DIR}/stats_filtered.txt

# echo "=== Filtering complete ==="
# echo "Filtered results in: ${FILTERED_VCF}"
# cat ${OUTPUT_DIR}/stats_filtered.txt
# echo "Date: $(date)"
