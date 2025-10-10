#!/bin/bash
#SBATCH --job-name=consensusv_debug8
#SBATCH --output=../logfiles/consensusv_debug8_%j.out
#SBATCH --error=../logfiles/consensusv_debug8_%j.err
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
    --bind ${BASE_DIR}/software/ConsensuSV-core:/tools/ConsensuSV-core \
    --pwd /test_working_dir \
    ${BASE_DIR}/software/consensusv-pipeline.sif \
    bash -c "
        export TMPDIR=/test_working_dir/tmp
        cd /test_working_dir && \
        luigid --background --port 8082 --logdir /var/log/luigi && \
        sleep 5 && \
        python /tools/ConsensuSV-pipeline/run_consensusv.py RunCSVFile \
            --csv-file /tools/ConsensuSV-pipeline/samples.csv \
            --workers 4 \
            --working-dir /test_working_dir/ \
            --log-level DEBUG 2>&1 | tee /test_working_dir/debug_run.log
    "

echo "Job completed"