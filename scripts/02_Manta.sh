#!/bin/bash
#SBATCH --job-name=manta_config
#SBATCH --output=/storage/homefs/kw23y068/logfiles/manta_config_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/manta_config_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=epyc2

#Python version must be 2.6< Version <3.0
# We use latest Python 2 on cluster
module load Python/2.7.18-GCCcore-11.3.0-bare

# Define paths
MANTA_DIR=/storage/homefs/kw23y068/software
CONFIG_SCRIPT=${MANTA_DIR}/manta-1.6.0.centos6_x86_64/bin/configManta.py
PROJECT_DIR=/storage/research/iee_evol/Korbi

# Input files
BAM=${PROJECT_DIR}/bams_real/FG_CC_19T_031.fixmate.coordsorted.bam
# Manta script adds .fai automatically, so I remove it here
REFERENCE=${PROJECT_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna
RUN_DIR=/storage/scratch/iee_evol/kw23y068/Manta/FG_CC_19T_031_manta

if [ ! -d "${RUN_DIR}" ]
then
    mkdir ${RUN_DIR}
fi


# Step 1: Configure Manta workflow
python $CONFIG_SCRIPT \
    --bam $BAM \
    --referenceFasta $REFERENCE \
    --runDir $RUN_DIR

# Step 2: Run Manta
python $RUN_DIR/runWorkflow.py -m local -j ${SLURM_CPUS_PER_TASK}
