#!/bin/bash
#SBATCH --job-name=consensusv_debug
#SBATCH --output=../logfiles/consensusv_debug_%j.out
#SBATCH --error=../logfiles/consensusv_debug_%j.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=epyc2

cd /storage/homefs/kw23y068/software/ConsensuSV-pipeline


# Run the debug command
apptainer exec /storage/homefs/kw23y068/software/consensusv-pipeline.sif \
    ls -lh /tools/ConsensuSV-core/
