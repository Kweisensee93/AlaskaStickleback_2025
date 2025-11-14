#!/bin/bash
#SBATCH --job-name=GRIDSS
#SBATCH --output=/storage/homefs/kw23y068/logfiles/GRIDSS_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/GRIDSS_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=9
#SBATCH --mem=32G
#SBATCH --partition=epyc2

set -euo pipefail

module load R/4.4.2-gfbf-2024a
New_Rscript=$(which Rscript)

SAMPLE_NAME="FG_CC_19T_031"

# Define paths
GRIDSS_IMAGE=/storage/homefs/kw23y068/software/gridss/Patched_GRIDSS.sif
PROJECT_DIR=/storage/research/iee_evol/Korbi
REFERENCE=${PROJECT_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna
RUN_DIR=/storage/scratch/iee_evol/kw23y068/Gridss/${SAMPLE_NAME}

if [ ! -d "${RUN_DIR}" ]
then
    mkdir ${RUN_DIR}
fi

# Input files
BAM=${PROJECT_DIR}/bams_real/${SAMPLE_NAME}.fixmate.coordsorted.bam

# We need to get R inside the container
#RSCRIPT_PATH=$(which Rscript)
#R_BASE=$(dirname $(dirname "${RSCRIPT_PATH}"))

# for degub: Check if binding of R from module to image works
# apptainer exec --bind ${R_BASE} \
#     ${GRIDSS_IMAGE} Rscript --version
# --bind ${R_BASE}/bin:/usr/local/bin \

# apptainer exec \
#     --bind ${PROJECT_DIR} \
#     --bind ${RUN_DIR} \
#     ${GRIDSS_IMAGE} \
#     env PATH=/usr/bin:$PATH /opt/gridss/gridss \
#     --reference ${REFERENCE} \
#     --output ${RUN_DIR}/${SAMPLE_NAME}.vcf.gz \
#     ${BAM}

apptainer exec \
    --bind ${New_Rscript}:/usr/bin/Rscript \
    --bind ${PROJECT_DIR} \
    --bind ${RUN_DIR} \
    ${GRIDSS_IMAGE} \
    /opt/gridss/gridss \
        --reference ${REFERENCE} \
        --output ${RUN_DIR}/${SAMPLE_NAME}.vcf.gz \
        ${BAM}
