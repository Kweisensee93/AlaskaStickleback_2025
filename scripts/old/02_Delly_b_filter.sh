#!/bin/bash
#SBATCH --job-name=delly_filter
#SBATCH --output=/storage/homefs/kw23y068/logfiles/delly_filter_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/delly_filter_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=epyc2

# Delly uses a Threshold for "Pass"
#Pass/LowQual is indeed related to paired end mapping. SVs are flagged as "PASS" if >=3 paired-ends support the variant (for translocations >=5). In addition, the mean mapping quality has to be >=20
# https://groups.google.com/g/delly-users/c/6Mq2juBraRY
# https://groups.google.com/g/delly-users/c/x1YAvpAhOAI

echo "=== Starting Delly SV filtering ==="
echo "Date: $(date)"
echo "Host: $(hostname)"

# Load necessary modules
module load BCFtools/1.12-GCC-10.3.0

# Define paths
OUTPUT_DIR=/storage/scratch/iee_evol/kw23y068/Delly/FG_CC_19T_031_delly
INPUT_VCF=${OUTPUT_DIR}/FG_CC_19T_031.vcf.gz
FILTERED_VCF=${OUTPUT_DIR}/FG_CC_19T_031.PASS.ENDsize.vcf.gz
# Maybe use a variable for length threshold later
LENGTH_THRESHOLD=1000000
# Delly also creates rather large SVs. Hence an upper size filter is used.
# as for its upper limit i set 1,000,000 bp (1Mb) as variants larger than this are often spurious.
# This upper boundary is an arbitrary choice and may be adjusted if deemed necessary.

# Filter for PASS variants only
echo ">>> Filtering for PASS variants only and applying length threshold: ${LENGTH_THRESHOLD} bp"
# difference -f PASS and -i 'FILTER="PASS"'?
#FILTER_EXPR='FILTER="PASS" && SVLEN <= '${LENGTH_THRESHOLD}' && SVLEN >= -'${LENGTH_THRESHOLD}
# The filter above kicks out variants with no value for SVLEN
#FILTER_EXPR='FILTER="PASS" && (abs(POS - INFO/END )<=${LENGTH_THRESHOLD})'
FILTER_EXPR="FILTER=\"PASS\" && (abs(POS - INFO/END) <= ${LENGTH_THRESHOLD})"


bcftools view \
    -i "${FILTER_EXPR}" \
    -O z \
    -o "${FILTERED_VCF}" \
    "${INPUT_VCF}"

# # Delly also creates rather large SVs. Hence an upper size filter is used.
# # as for its upper limit i set 1,000,000 bp (1Mb) as variants larger than this are often spurious.
# # This upper boundary is an arbitrary choice and may be adjusted if deemed necessary.
# echo ">>> Applying size filter: keeping variants <= 1,000,000 bp"
# # if we also exclude below 50bp, use:
# #-i 'SVLEN <= 1000000 && SVLEN >= -1000000 && (SVLEN >= 50 || SVLEN <= -50)' \
# bcftools filter \
#     -i 'SVLEN <= 1000000 && SVLEN >= -1000000' \
#     -O z \
#     -o ${OUTPUT_DIR}/FG_CC_19T_031.PASS.size_filtered.vcf.gz \
#     ${FILTERED_VCF}

# Index the filtered VCF
echo ">>> Indexing filtered VCF"
bcftools index ${FILTERED_VCF}

# Generate statistics for filtered variants
echo ">>> Generating statistics for PASS variants"
echo "PASS SVs:" > ${OUTPUT_DIR}/stats_filtered.txt
bcftools view -H ${FILTERED_VCF} | wc -l >> ${OUTPUT_DIR}/stats_filtered.txt

echo "" >> ${OUTPUT_DIR}/stats_filtered.txt
echo "PASS SV types:" >> ${OUTPUT_DIR}/stats_filtered.txt
bcftools query -f '%SVTYPE\n' ${FILTERED_VCF} | sort | uniq -c >> ${OUTPUT_DIR}/stats_filtered.txt

# Compare with original
echo "" >> ${OUTPUT_DIR}/stats_filtered.txt
echo "Comparison:" >> ${OUTPUT_DIR}/stats_filtered.txt
echo "Total SVs: $(bcftools view -H ${INPUT_VCF} | wc -l)" >> ${OUTPUT_DIR}/stats_filtered.txt
echo "PASS SVs: $(bcftools view -H ${FILTERED_VCF} | wc -l)" >> ${OUTPUT_DIR}/stats_filtered.txt

echo "=== Filtering complete ==="
echo "Filtered results in: ${FILTERED_VCF}"
cat ${OUTPUT_DIR}/stats_filtered.txt
echo "Date: $(date)"