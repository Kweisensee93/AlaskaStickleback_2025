#!/bin/bash
#SBATCH --job-name=consensusv_OneSample
#SBATCH --output=/storage/research/iee_evol/Korbi/logfiles/%j_consensusv_OneSample.out
#SBATCH --error=/storage/research/iee_evol/Korbi/logfiles/%j_consensusv_OneSample.err
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=epyc2

echo "Job started on $(date)"

# Directories
PROJECT_DIR=/storage/research/iee_evol/Korbi
SOFTWARE_DIR=/storage/homefs/kw23y068
BASE_DIR=/storage/research/iee_evol/Korbi/output/ConsensuSV
WORKING_DIR=${BASE_DIR}/working_dir
LOG_DIR=${BASE_DIR}/logfiles_luigi
TMP_DIR=${BASE_DIR}/working_dir/tmp

# Create directories if they don't exist
if [ ! -d "$WORKING_DIR" ]; then
    mkdir $WORKING_DIR
fi
if [ ! -d "$LOG_DIR" ]; then
    mkdir $LOG_DIR
fi
if [ ! -d "$TMP_DIR" ]; then
    mkdir $TMP_DIR
fi

# we need the following binds (mainly due to apptainer limitations instead of docker):
# bind whole project directory - just in case
# we use the ConsensuSV-core outside of the container, because it creates and needs writeable files/paths
# The same goes for the ConsensuSV-pipeline
# we need to bind the working dir, because we want to write there (including output files of the pipeline)
# we need to bind the luigi logdir, because luigi needs to write there
# bind the softlink for rawdata, because we want to read the fastq files from there
apptainer exec \
    --bind ${PROJECT_DIR} \
    --bind ${SOFTWARE_DIR}/software/ConsensuSV-core:/tools/ConsensuSV-core \
    --bind ${SOFTWARE_DIR}/software/ConsensuSV-pipeline:/tools/ConsensuSV-pipeline \
    --bind ${WORKING_DIR}:/test_working_dir \
    --bind ${LOG_DIR}:/var/log/luigi \
    --bind /storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/raw_BACKUP/ \
    ${SOFTWARE_DIR}/software/consensusv-pipeline.sif \
    bash -c "
        luigid --background --port 8082 --logdir /var/log/luigi && \
        sleep 5 && \
        python /workspace/run_consensusv.py RunCSVFile \
            --csv-file ${BASE_DIR}/first_sample.csv \
            --workers 4 \
            --working-dir /test_working_dir/ \
    "

echo "Job completed at $(date)"