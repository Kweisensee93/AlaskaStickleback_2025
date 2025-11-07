#!/bin/bash
#SBATCH --job-name=GRIDSS
#SBATCH --output=/storage/homefs/kw23y068/logfiles/GRIDSS_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/GRIDSS_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=9
#SBATCH --mem=32G
#SBATCH --partition=epyc2

echo "starting GRIDSS SV calling"
echo "Date: $(date)"

# Load Anaconda module to access conda environments
module load Anaconda3/2024.02-1

SAMPLE_NAME="FG_CC_19T_031"

# Define paths
GRIDSS_CONDA=/storage/homefs/kw23y068/.conda/envs/gridss_test/bin/gridss
PROJECT_DIR=/storage/research/iee_evol/Korbi
RUN_DIR=/storage/scratch/iee_evol/kw23y068/Gridss/${SAMPLE_NAME}

# Create run directory if it doesn't exist
if [ ! -d "${RUN_DIR}" ]; then
    mkdir -p "${RUN_DIR}"
fi

# Input files
REFERENCE=${PROJECT_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna
BAM=${PROJECT_DIR}/bams_real/${SAMPLE_NAME}.fixmate.coordsorted.bam
# Output files
VCF_OUT=${RUN_DIR}/${SAMPLE_NAME}.vcf.gz
ASSEMBLY_BAM=${RUN_DIR}/${SAMPLE_NAME}.assembly.bam

# Verify inputs exist
[ -f "${REFERENCE}" ] || { echo "ERROR: Reference not found: ${REFERENCE}"; exit 1; }
[ -f "${BAM}" ] || { echo "ERROR: BAM file not found: ${BAM}"; exit 1; }

# from GRIDSS documentation:
# Usage: gridss [options] -r <reference.fa> -o <output.vcf.gz> -a <assembly.bam> input1.bam [input2.bam [...]]
# Run GRIDSS
# Activate GRIDSS environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate gridss

# Gridss uses 8 threads per default, so we set threads to SLURM_CPUS_PER_TASK-1
# I gave 9 CPUS to have a bit of overhead, needs to be adapted for multiple samples!
${GRIDSS_CONDA} \
    -r ${REFERENCE} \
    -o ${VCF_OUT} \
    -a ${ASSEMBLY_BAM} \
    -w ${RUN_DIR} \
    ${BAM}

echo "Finished GRIDSS SV calling"
echo "Date: $(date)"