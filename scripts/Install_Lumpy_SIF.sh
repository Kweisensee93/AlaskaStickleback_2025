#!/bin/bash
#SBATCH --job-name=build_Lumpy
#SBATCH --output=/storage/homefs/kw23y068/logfiles/install_Lumpy_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/install_Lumpy_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=80G
#SBATCH --time=01:00:00
#SBATCH --partition=epyc2


mkdir -p /storage/homefs/kw23y068/software/lumpy
cd /storage/homefs/kw23y068/software/lumpy
apptainer build Lumpy.sif docker://quay.io/biocontainers/lumpy-sv:0.3.1--3

echo "check lumpy versions Latest GitHub is 0.3.1"
apptainer exec "/storage/homefs/kw23y068/software/lumpy/Lumpy.sif" lumpyexpress --version 2>&1
apptainer exec "/storage/homefs/kw23y068/software/lumpy/Lumpy.sif" lumpy --help 2>&1