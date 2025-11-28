#!/bin/bash
#SBATCH --job-name=Lumpy_postprocess
#SBATCH --output=/storage/homefs/kw23y068/logfiles/Lumpy_postprocess_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/Lumpy_postprocess_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=epyc2

## We may need to mimic smoove steps within lumpy ???

###########################
# PARAMETERS TO SET
###########################
BATCH_START=1
BATCH_END=80

PROJECT_DIR=/storage/research/iee_evol/Korbi
REFERENCE=${PROJECT_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna

START_PAD=$(printf "%03d" "${BATCH_START}")
END_PAD=$(printf "%03d" "${BATCH_END}")
OUT_DIR="/storage/scratch/iee_evol/kw23y068/lumpy_joint_${START_PAD}_${END_PAD}"

VCF_RAW="${OUT_DIR}/lumpy_joint_${START_PAD}_${END_PAD}.vcf.gz"
VCF_GT="${OUT_DIR}/lumpy_joint_${START_PAD}_${END_PAD}.genotyped.vcf.gz"
VCF_DH="${OUT_DIR}/lumpy_joint_${START_PAD}_${END_PAD}.genotyped.duphold.vcf.gz"
VCF_FILTERED="${OUT_DIR}/lumpy_joint_${START_PAD}_${END_PAD}.final.filtered.vcf.gz"

echo "========================================"
echo " LUMPY POSTPROCESSING PIPELINE"
echo " Batch: ${BATCH_START}-${BATCH_END}"
echo " Start: $(date)"
echo "========================================"

module load Anaconda3/2024.02-1
source $(conda info --base)/etc/profile.d/conda.sh
conda activate lumpy

##########################################
# 1. GENOTYPE SVs (svtyper)
##########################################

echo "[1/3] Genotyping with svtyper..."
svtyper \
    -i "${VCF_RAW}" \
    -B $(echo ${PROJECT_DIR}/bams_real/*.fixmate.coordsorted.bam | tr ' ' ',') \
    -l \
| bgzip -c > "${VCF_GT}"

bcftools index "${VCF_GT}"

##########################################
# 2. RUN DUPHOLD (depth annotation)
##########################################

echo "[2/3] Running duphold..."

# duphold needs a directory of BAMs and the reference index
duphold \
    -f "${REFERENCE}" \
    -o "${VCF_DH}" \
    -b "${PROJECT_DIR}/bams_real" \
    "${VCF_GT}"

bcftools index "${VCF_DH}"

##########################################
# 3. SIMPLE SV FILTERING (smoove-like)
##########################################

echo "[3/3] Filtering low-quality SVs..."

# Very similar to `smoove` default filters:
bcftools filter \
    -e 'INFO:SU<4 || INFO:DHBFC<0.7 || INFO:DHFFC<0.7' \
    -O z \
    -o "${VCF_FILTERED}" \
    "${VCF_DH}"

bcftools index "${VCF_FILTERED}"

conda deactivate

echo ""
echo "========================================"
echo "DONE!"
echo " Final VCF: ${VCF_FILTERED}"
echo " End: $(date)"
echo "========================================"


