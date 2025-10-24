#!/bin/bash
#SBATCH --job-name=delly_call
#SBATCH --output=/storage/homefs/kw23y068/logfiles/delly_call_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/delly_call_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=epyc2

echo "=== Starting Delly SV calling ==="
echo "Date: $(date)"
echo "Host: $(hostname)"

# Load necessary modules
module load BCFtools/1.12-GCC-10.3.0

# Define paths
DELLY_SIF=/storage/homefs/kw23y068/software/Delly/delly_v1.5.0.sif
PROJECT_DIR=/storage/research/iee_evol/Korbi
BAM=${PROJECT_DIR}/bams_real/FG_CC_19T_031.fixmate.coordsorted.bam
REFERENCE=${PROJECT_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna
OUTPUT_DIR=/storage/scratch/iee_evol/kw23y068/Delly/FG_CC_19T_031_delly

if [ ! -d "${OUTPUT_DIR}" ]
then
    mkdir ${OUTPUT_DIR}
fi

# # BAM and Reference are already indexed! Use Samtools if needed otherwise :)
# # Check if BAM index exists
# if [ ! -f "${BAM}.bai" ]; then
#     echo ">>> BAM index not found, creating it..."
#     singularity exec ${DELLY_SIF} samtools index ${BAM}
# fi

# # Check if reference index exists
# if [ ! -f "${REFERENCE}.fai" ]; then
#     echo ">>> Reference index not found, creating it..."
#     singularity exec ${DELLY_SIF} samtools faidx ${REFERENCE}
# fi

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

# Generate basic statistics
echo ">>> Generating statistics"
echo "Total SVs called:" > ${OUTPUT_DIR}/stats.txt
bcftools view -H ${OUTPUT_DIR}/FG_CC_19T_031.vcf.gz | wc -l >> ${OUTPUT_DIR}/stats.txt

echo "" >> ${OUTPUT_DIR}/stats.txt
echo "SV types:" >> ${OUTPUT_DIR}/stats.txt
bcftools query -f '%SVTYPE\n' ${OUTPUT_DIR}/FG_CC_19T_031.vcf.gz | sort | uniq -c >> ${OUTPUT_DIR}/stats.txt

echo "=== Delly analysis complete ==="
echo "Results in: ${OUTPUT_DIR}"
cat ${OUTPUT_DIR}/stats.txt
echo "Date: $(date)"