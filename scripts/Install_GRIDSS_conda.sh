#!/bin/bash
#SBATCH --job-name=build_GRIDSS_conda
#SBATCH --output=/storage/homefs/kw23y068/logfiles/install_GRIDSS_conda_%j.log
#SBATCH --error=/storage/homefs/kw23y068/logfiles/install_GRIDSS_conda_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=80G
#SBATCH --time=02:00:00
#SBATCH --partition=epyc2

module load Anaconda3/2024.02-1

# Ensure conda commands work in batch mode
source $(conda info --base)/etc/profile.d/conda.sh

# Set conda to use a writable location for environments
export CONDA_ENVS_PATH=$HOME/.conda/envs
export CONDA_PKGS_DIRS=$HOME/.conda/pkgs
mkdir -p $CONDA_ENVS_PATH $CONDA_PKGS_DIRS

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

ENV_NAME="gridss"
conda create -y -n ${ENV_NAME} gridss
conda activate ${ENV_NAME}

# Just install htslib to avoid issues (see GitHub)
conda install htslib

gridss version
echo "Installed at"
conda env list
echo "gridds shall be under ~/.conda/envs/gridss/bin/gridss"


