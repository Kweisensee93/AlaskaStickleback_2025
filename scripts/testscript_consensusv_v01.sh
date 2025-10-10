#!/bin/bash
#SBATCH --job-name=consensusv_test
#SBATCH --output=../logfiles/consensusv_test_%j.out
#SBATCH --error=../logfiles/consensusv_test_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=epyc2

# Ensure working dir exists
mkdir -p /storage/homefs/kw23y068/test_working_dir
mkdir -p /storage/homefs/kw23y068/logfiles/luigi

cd /storage/homefs/kw23y068/software/ConsensuSV-pipeline

# Start luigid with a writable logdir
apptainer exec /storage/homefs/kw23y068/software/consensusv-pipeline.sif \
    luigid --background --port 8082 --logdir /storage/homefs/kw23y068/logfiles/luigi

# Run the test
# We need a test woring directory with write permissions because we use apptainer insted of docker
apptainer exec \
    --bind /storage/homefs/kw23y068/software/ConsensuSV-pipeline:/ConsensuSV-pipeline \
    --bind /storage/homefs/kw23y068/test_working_dir:/test_working_dir \
    /storage/homefs/kw23y068/software/consensusv-pipeline.sif \
    bash /ConsensuSV-pipeline/test_run_fast.sh
