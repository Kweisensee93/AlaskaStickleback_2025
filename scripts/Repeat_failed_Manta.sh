#!/bin/bash
#SBATCH --job-name=manta_repeat
#SBATCH --array=1-9
#SBATCH --output=/storage/homefs/kw23y068/logfiles/manta_repeat_%A_%a.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/manta_repeat_%A_%a.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=epyc2

# Array job is defined by 1 to
# wc -l </storage/scratch/iee_evol/kw23y068/Manta/failed_bams.csv

#Python version must be 2.6< Version <3.0
# We use latest Python 2 on cluster
module load Python/2.7.18-GCCcore-11.3.0-bare

# failed sample list created as follows:
# > failed_runs.txt
# for d in */; do
#     if [ ! -f "$d/results/variants/diploidSV.vcf.gz" ]; then
#         echo "${d%/}" >> failed_runs.txt
#     fi
# done
# awk -v base="/storage/research/iee_evol/Korbi/bams_real" '{gsub("/",""); \
#   print base"/"$1".fixmate.coordsorted.bam"}' failed_runs.txt > failed_bams.csv


## ==========================================================
##  ORIGINAL SAMPLE INPUT (commented out)
## ==========================================================
#Define Samples from input list
# SAMPLE_LIST=/storage/research/iee_evol/Korbi/output/bams_full.csv
# SAMPLE_PATH=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")
# SAMPLE_FILE=$(basename "${SAMPLE_PATH}")
# SAMPLE_NAME="${SAMPLE_FILE%.fixmate.coordsorted.bam}"

## ==========================================================
##  UPDATED INPUT: failed samples only
## ==========================================================
SAMPLE_LIST=/storage/scratch/iee_evol/kw23y068/Manta/failed_bams.csv
SAMPLE_PATH=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")
SAMPLE_FILE=$(basename "${SAMPLE_PATH}")
SAMPLE_NAME="${SAMPLE_FILE%.fixmate.coordsorted.bam}"


# Define paths
MANTA_DIR=/storage/homefs/kw23y068/software
CONFIG_SCRIPT=${MANTA_DIR}/manta-1.6.0.centos6_x86_64/bin/configManta.py
PROJECT_DIR=/storage/research/iee_evol/Korbi

# Input files
#SAMPLE_PATH replaces:
#BAM=${PROJECT_DIR}/bams_real/FG_CC_19T_031.fixmate.coordsorted.bam
# Manta script adds .fai automatically, so I remove it here
REFERENCE=${PROJECT_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna
RUN_DIR=/storage/scratch/iee_evol/kw23y068/Manta/${SAMPLE_NAME}

# Clean up previous run directory if it exists
if [ -d "${RUN_DIR}" ]; then
    rm -rf "${RUN_DIR}"
fi

# Create run directory (after cleanup)
if [ ! -d "${RUN_DIR}" ]
then
    mkdir ${RUN_DIR}
fi


# Step 1: Configure Manta workflow
python $CONFIG_SCRIPT \
    --bam $SAMPLE_PATH \
    --referenceFasta $REFERENCE \
    --runDir $RUN_DIR

# Step 2: Run Manta
python $RUN_DIR/runWorkflow.py -m local -j ${SLURM_CPUS_PER_TASK}
