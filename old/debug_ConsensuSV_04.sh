#!/bin/bash
#SBATCH --job-name=debug4_consensusv
#SBATCH --output=../logfiles/consensusv_debug4_%j.out
#SBATCH --error=../logfiles/consensusv_debug4_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=epyc2

# Base directory
BASE_DIR=/storage/homefs/kw23y068

# Ensure working dirs exist
mkdir -p ${BASE_DIR}/test_working_dir
mkdir -p ${BASE_DIR}/logfiles/luigi
mkdir -p ${BASE_DIR}/test_working_dir/tmp

export TMPDIR=${BASE_DIR}/test_working_dir/tmp

cd ${BASE_DIR}/software/ConsensuSV-pipeline

# Run with verbose logging
apptainer exec \
    --bind ${BASE_DIR}/software/ConsensuSV-pipeline:/tools/ConsensuSV-pipeline \
    --bind ${BASE_DIR}/test_working_dir:/test_working_dir \
    --bind ${BASE_DIR}/logfiles/luigi:/var/log/luigi \
    --pwd /tools/ConsensuSV-pipeline \
    ${BASE_DIR}/software/consensusv-pipeline.sif \
    bash -c "
        set -x
        luigid --background --port 8082 --logdir /var/log/luigi && \
        sleep 5 && \
        echo '=== Checking ConsensuSV-core ===' && \
        ls -la /tools/ConsensuSV-core/ && \
        echo '=== Checking model file ===' && \
        ls -la /tools/ConsensuSV-core/pretrained_1000g_illumina.model && \
        echo '=== Running test ===' && \
        ./test_run_csv.sh
    "

echo "Job completed"