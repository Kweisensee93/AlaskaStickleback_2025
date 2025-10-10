#!/bin/bash
#SBATCH --job-name=debug5_consensusv
#SBATCH --output=../logfiles/consensusv_debug5_%j.out
#SBATCH --error=../logfiles/consensusv_debug5_%j.err
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=epyc2

BASE_DIR=/storage/homefs/kw23y068

mkdir -p ${BASE_DIR}/test_working_dir
mkdir -p ${BASE_DIR}/logfiles/luigi
mkdir -p ${BASE_DIR}/test_working_dir/tmp

export TMPDIR=${BASE_DIR}/test_working_dir/tmp

cd ${BASE_DIR}/software/ConsensuSV-pipeline

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
        echo '=== Running with verbose Luigi logging ===' && \
        python run_consensusv.py RunCSVFile \
            --csv-file samples.csv \
            --workers 4 \
            --working-dir /test_working_dir/ \
            --log-level DEBUG 2>&1 | tee /test_working_dir/consensusv_run.log
    "

echo "Job completed"