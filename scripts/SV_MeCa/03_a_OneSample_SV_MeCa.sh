#!/bin/bash
#SBATCH --job-name=svmeca_OneSample
#SBATCH --cpus-per-task=16
#SBATCH --mem=70G
#SBATCH --time=04:00:00
#SBATCH --output=/storage/research/iee_evol/Korbi/logfiles/%j_SvMeCa_OneSample.log
#SBATCH --error=/storage/research/iee_evol/Korbi/logfiles/%j_SvMeCa_OneSample.err
#SBATCH --partition=epyc2

# make safer pipeline
set -Eeuo pipefail
set -o verbose

# Directories
PROJECT_DIR=/storage/research/iee_evol/Korbi
SOFTWARE_DIR=/storage/homefs/kw23y068
BASE_DIR=/storage/research/iee_evol/Korbi/output/SV-MeCa
TMP_DIR=${BASE_DIR}/tmp
BAMS_DIR="${PROJECT_DIR}/bams_real"

# ressources
CPU=$SLURM_CPUS_PER_TASK
RAM=64G     # use slightly less than allocated memory to avoid OOM kills
REFERENCE="${PROJECT_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna"

# get the container (workaround for apptainer)
CONTAINER="docker://wembasop/sv-meca:latest"

# create a local writable temp directory for Java/Nextflow
if [ ! -d "$TMP_DIR" ]; then
    mkdir $TMP_DIR
fi

SAMPLE=$(basename "$(head -n1 "${BASE_DIR}/first_sample.csv")")
SAMPLE_NAME=${SAMPLE:0:13}  # first 13 characters of the sample name

# run SV-MeCa (approx 2.5 hours)
#docker run -v $(pwd)/input:/input -v $(pwd)/output:/workspace/SV-MeCa/results wembasop/sv-meca:latest "/workspace/SV-MeCa/run_svmeca.sh bam -bam /input/HG00514.down.bam -ref /input/$REFERENCE -sample HG00514 -build hg38 -has_chr true -bed /input/hg38_centromer.bed" 
apptainer exec \
    --bind ${PROJECT_DIR} \
    --bind ${BAMS_DIR}:/input \
    --bind ${BASE_DIR}:/workspace/SV-MeCa/results \
    --bind ${TMP_DIR}:/tmp \
    $CONTAINER \
    bash -c "export TMPDIR=/tmp; export NXF_TMP=/tmp; export JAVA_CMD=/opt/conda/envs/align/bin/java; \
    /workspace/SV-MeCa/run_svmeca.sh \
    bam -bam /input/${SAMPLE} \
    -ref $REFERENCE -sample $SAMPLE_NAME -build hg38 -has_chr true"


# clean up temp directory
rm -rf ${TMP_DIR}/*

