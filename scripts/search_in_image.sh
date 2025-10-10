#!/bin/bash
#SBATCH --job-name=check_consensusv_path
#SBATCH --output=../../logfiles/check_consensusv_path_%j.out
#SBATCH --error=../../logfiles/check_consensusv_path_%j.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --partition=epyc2

BASE_DIR=/storage/homefs/kw23y068

echo "Searching for run_consensusv.py inside the container..."
apptainer exec ${BASE_DIR}/software/consensusv-pipeline.sif \
    bash -c "find / -type f -name run_consensusv.py 2>/dev/null || echo 'Not found.'"

echo "Done."
