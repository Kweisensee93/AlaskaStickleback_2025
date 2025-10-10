#!/bin/bash
#SBATCH --job-name=consensusv_test2
#SBATCH --output=/storage/homefs/kw23y068/logfiles/%j_consensusv_test2.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/%j_consensusv_test2.err
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=epyc2

BASE_DIR=/storage/homefs/kw23y068

mkdir -p ${BASE_DIR}/output/consensusv/test/test_working_dir
mkdir -p ${BASE_DIR}/output/consensusv/test/logfiles/luigi
mkdir -p ${BASE_DIR}/output/consensusv/test/test_working_dir/tmp


cd ${BASE_DIR}/software/ConsensuSV-pipeline

# we need the following binds (mainly due to apptainer limitations instead of docker):
# bind whole home directory - just in case
# we use the ConsensuSV-core outside of the container, because it creates and needs writeable files/paths
# The same goes for the ConsensuSV-pipeline
# we need to bind the working dir, because we want to write there (including output files of the pipeline)
# we need to bind the luigi logdir, because luigi needs to write there
apptainer exec \
    --bind ${BASE_DIR} \
    --bind ${BASE_DIR}/software/ConsensuSV-core:/tools/ConsensuSV-core \
    --bind ${BASE_DIR}/software/ConsensuSV-pipeline:/tools/ConsensuSV-pipeline \
    --bind ${BASE_DIR}/output/consensusv/test/test_working_dir:/test_working_dir \
    --bind ${BASE_DIR}/output/consensusv/test/logfiles/luigi:/var/log/luigi \
    ${BASE_DIR}/software/consensusv-pipeline.sif \
    bash -c "
        luigid --background --port 8082 --logdir /var/log/luigi && \
        sleep 5 && \
        python /workspace/run_consensusv.py RunCSVFile \
            --csv-file /tools/ConsensuSV-pipeline/samples.csv \
            --workers 4 \
            --working-dir /test_working_dir/ \
    "

echo "Job completed"