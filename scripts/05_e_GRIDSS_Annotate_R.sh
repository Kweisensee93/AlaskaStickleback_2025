#!/bin/bash
#SBATCH --job-name=GRIDSS_R_annotate
#SBATCH --output=/storage/homefs/kw23y068/logfiles/GRIDSS_R_annotate_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/GRIDSS_R_annotate_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=300G
#SBATCH --partition=epyc2

# OOM kill at 128G

####################
# PARAMETERS TO SET
####################
# May be altered for feeding in arguments to script
# Define which batch of samples to process
BATCH_START=1    # Change to 81, 161, 241, 321, 401 for other batches
BATCH_END=80     # Change to 160, 240, 320, 400, 480 for other batches

START_PAD=$(printf "%03d" "${BATCH_START}")
END_PAD=$(printf "%03d" "${BATCH_END}")

# Paths
PROJECT_DIR=/storage/research/iee_evol/Korbi
RUN_DIR="/storage/scratch/iee_evol/kw23y068/Gridss_joint_${START_PAD}_${END_PAD}/preprocessed"

VCF_OUT=${RUN_DIR}/assembly/joint.vcf.gz
VCF_ANNOTATED=${RUN_DIR}/assembly/joint.annotated.vcf
SIMPLE_BED=${RUN_DIR}/assembly/joint.simple.bed

module load Anaconda3/2024.02-1
source $(conda info --base)/etc/profile.d/conda.sh
conda activate R_gridss_annotation

# Annotate with StructuralVariantAnnotation
echo "Starting SV annotation..."
echo "Date: $(date)"

Rscript /storage/homefs/kw23y068/software/scripts/Annotate_GRIDSS_v02.R \
    ${VCF_OUT} \
    "GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna" \
    ${VCF_ANNOTATED} \
    ${SIMPLE_BED}

conda deactivate
# We use the same versions of samtools/bcftools in both environments, so we use gridss env again
conda activate gridss
bgzip ${RUN_DIR}/${SAMPLE_NAME}.annotated.vcf
tabix -p vcf ${RUN_DIR}/${SAMPLE_NAME}.annotated.vcf.gz
conda deactivate
