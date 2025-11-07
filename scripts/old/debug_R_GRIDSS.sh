#!/bin/bash
#SBATCH --job-name=R_GRIDSS
#SBATCH --output=/storage/homefs/kw23y068/logfiles/debug_R_GRIDSS_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/debug_R_GRIDSS_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=9
#SBATCH --mem=32G
#SBATCH --partition=epyc2

module load Anaconda3/2024.02-1

SAMPLE_NAME="FG_CC_19T_031"

GRIDSS_CONDA=/storage/homefs/kw23y068/.conda/envs/gridss/bin/gridss
PROJECT_DIR=/storage/research/iee_evol/Korbi
RUN_DIR=/storage/scratch/iee_evol/kw23y068/Gridss/${SAMPLE_NAME}

# Output files
VCF_OUT=${RUN_DIR}/${SAMPLE_NAME}.vcf.gz
ASSEMBLY_BAM=${RUN_DIR}/${SAMPLE_NAME}.assembly.bam
VCF_ANNOTATED=${RUN_DIR}/${SAMPLE_NAME}.annotated.vcf
SIMPLE_BED=${RUN_DIR}/${SAMPLE_NAME}.simple.bed

source $(conda info --base)/etc/profile.d/conda.sh
conda activate R_gridss_annotation

Rscript ${RUN_DIR}/annotate_sv.R \
    ${VCF_OUT} \
    "GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna" \
    ${VCF_ANNOTATED} \
    ${SIMPLE_BED}
