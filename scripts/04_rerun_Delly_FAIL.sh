#!/bin/bash
#SBATCH --job-name=delly_rerun
#SBATCH --array=1-6
#SBATCH --output=/storage/homefs/kw23y068/logfiles/delly_rerun_%A_%a.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/delly_rerun_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=epyc2


# Adjust array for number of failed samples. Use:
#grep '^FAILED\|^PARTIAL' /storage/scratch/iee_evol/kw23y068/Delly/*/status.txt | \
#    sed -E 's#.*/Delly/([^/]+)/status\.txt.*#\1#' | sort -u | wc -l
# Use number to adjust --array=1-X

# Define base directory
DELLY_BASE_DIR=/storage/scratch/iee_evol/kw23y068/Delly

# retrieve sample names of failed samples using
# grep -h '^SUCCESS\\|^FAILED\\|^PARTIAL' /storage/scratch/iee_evol/kw23y068/Delly/*/status.txt | sort | uniq -c | awk -F'/' '{print $(NF-1)}'
#grep '^FAILED\|^PARTIAL' /storage/scratch/iee_evol/kw23y068/Delly/*/status.txt | sed -E 's#.*/Delly/([^/]+)/status\.txt.*#\1#' | sort -u
# Get list of failed sample names
# remove -h for grep command:
FAILED_SAMPLES=$(grep '^FAILED\|^PARTIAL' ${DELLY_BASE_DIR}/*/status.txt | \
    sed -E 's#.*/Delly/([^/]+)/status\.txt.*#\1#' | sort -u)

# Convert to array; Retrieve Samplename
FAILED_ARRAY=($FAILED_SAMPLES)
SAMPLE_NAME="${FAILED_ARRAY[$((SLURM_ARRAY_TASK_ID - 1))]}"


#Define Samples from input list
SAMPLE_LIST=/storage/research/iee_evol/Korbi/output/bams_full.csv
SAMPLE_PATH=$(grep "${SAMPLE_NAME}.fixmate.coordsorted.bam" "${SAMPLE_LIST}")
SAMPLE_FILE=$(basename "${SAMPLE_PATH}")
#SAMPLE_NAME="${SAMPLE_FILE%.fixmate.coordsorted.bam}" #not needed

# Load necessary modules
module load BCFtools/1.12-GCC-10.3.0

echo "Processing sample: ${SAMPLE_NAME}"

# Define paths
DELLY_SIF=/storage/homefs/kw23y068/software/Delly/delly_v1.5.0.sif
PROJECT_DIR=/storage/research/iee_evol/Korbi
REFERENCE=${PROJECT_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna
OUTPUT_DIR=/storage/scratch/iee_evol/kw23y068/Delly/${SAMPLE_NAME}
# files
#SAMPLE_PATH replaces:
#BAM=${PROJECT_DIR}/bams_real/FG_CC_19T_031.fixmate.coordsorted.bam
DELLY_VCF=${OUTPUT_DIR}/${SAMPLE_NAME}.vcf.gz
FILTERED_VCF=${OUTPUT_DIR}/${SAMPLE_NAME}.filtered.vcf.gz

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
    -o ${OUTPUT_DIR}/${SAMPLE_NAME}.bcf \
    ${SAMPLE_PATH}

sleep 5s

# Convert BCF to VCF for easier viewing
echo ">>> Converting BCF to VCF"
bcftools view \
    ${OUTPUT_DIR}/${SAMPLE_NAME}.bcf \
    -O v \
    -o ${OUTPUT_DIR}/${SAMPLE_NAME}.vcf

sleep 5s

# Compress and index VCF
echo ">>> Compressing and indexing VCF"
bgzip -c ${OUTPUT_DIR}/${SAMPLE_NAME}.vcf > ${OUTPUT_DIR}/${SAMPLE_NAME}.vcf.gz
sleep 5s
bcftools index ${OUTPUT_DIR}/${SAMPLE_NAME}.vcf.gz
sleep 5s

echo ">>> Filtering for PASS variants only and applying length threshold: ${LENGTH_THRESHOLD} bp"
FILTER_EXPR="FILTER=\"PASS\" && (abs(POS - INFO/END) <= ${LENGTH_THRESHOLD})"

bcftools view \
    -i "${FILTER_EXPR}" \
    -O z \
    -o "${FILTERED_VCF}" \
    "${DELLY_VCF}"

sleep 5s
# Index the filtered VCF
echo ">>> Indexing filtered VCF"
bcftools index ${FILTERED_VCF}
