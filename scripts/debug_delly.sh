#!/bin/bash
#SBATCH --job-name=debug_delly
#SBATCH --output=/storage/homefs/kw23y068/logfiles/debug_delly_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/debug_delly_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --partition=epyc2


DELLY_SIF=/storage/homefs/kw23y068/software/Delly/delly_v1.5.0.sif

apptainer exec \
    --bind /storage/homefs/kw23y068/ \
    ${DELLY_SIF} delly --help
