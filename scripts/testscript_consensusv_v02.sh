#!/bin/bash
#SBATCH --job-name=consensusv_test1
#SBATCH --output=/storage/homefs/kw23y068/logfiles/%j_scratch_consensusv_test1.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/%j_scratch_consensusv_test1.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=epyc2


# Ensure working dir exists
mkdir -p /storage/homefs/kw23y068/output/consensusv/test/test_working_dir
mkdir -p /storage/homefs/kw23y068/output/consensusv/test/logfiles/luigi

cd /storage/homefs/kw23y068/software/ConsensuSV-pipeline

# Run the test
# We need a test woring directory with write permissions because we use apptainer insted of docker
apptainer exec \
    --bind /storage/homefs/kw23y068/software/ConsensuSV-pipeline:/ConsensuSV-pipeline \
    --bind /storage/homefs/kw23y068/output/consensusv/test/test_working_dir:/test_working_dir \
    /storage/homefs/kw23y068/software/consensusv-pipeline.sif \
    bash -c "
        luigid --background --port 8082 --logdir /storage/homefs/kw23y068/output/consensusv/test/logfiles/luigi && \
        sleep 5 && \
        bash /ConsensuSV-pipeline/test_run_fast.sh
    "
