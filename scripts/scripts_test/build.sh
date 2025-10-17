#!/bin/bash
#SBATCH --job-name=build_sif
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=80G
#SBATCH --time=10:00:00
#SBATCH --partition=bdw

module load apptainer
apptainer build consensusv-pipeline.sif docker://mateuszchilinski/consensusv-pipeline

