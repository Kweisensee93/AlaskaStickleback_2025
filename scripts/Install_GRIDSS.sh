#!/bin/bash
#SBATCH --job-name=build_GRIDSS
#SBATCH --output=/storage/homefs/kw23y068/logfiles/install_GRIDSS_%j.log
#SBATCH --error=/storage/homefs/kw23y068/logfiles/install_GRIDSS_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=80G
#SBATCH --time=10:00:00
#SBATCH --partition=epyc2


mkdir -p /storage/homefs/kw23y068/software/gridss
cd /storage/homefs/kw23y068/software/gridss
apptainer build GRIDSS.sif docker://gridss/gridss:latest