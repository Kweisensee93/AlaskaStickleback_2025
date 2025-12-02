#!/bin/bash
#SBATCH --job-name=build_Jasmine_conda
#SBATCH --output=/storage/homefs/kw23y068/logfiles/install_Jasmine_conda_%j.log
#SBATCH --error=/storage/homefs/kw23y068/logfiles/install_Jasmine_conda_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=80G
#SBATCH --time=02:00:00
#SBATCH --partition=epyc2
#SBATCH --account=gratis
#SBATCH --wckey=noop

module load Anaconda3/2024.02-1

# Ensure conda commands work in batch mode
source $(conda info --base)/etc/profile.d/conda.sh

# Set conda to use a writable location for environments
export CONDA_ENVS_PATH=$HOME/.conda/envs
export CONDA_PKGS_DIRS=$HOME/.conda/pkgs
mkdir -p $CONDA_ENVS_PATH $CONDA_PKGS_DIRS

# Configure channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

ENV_NAME="jasmine"
conda create -y -n ${ENV_NAME} jasminesv
conda activate ${ENV_NAME}

# Check installation
jasmine --version 2>&1 || echo "jasmine does not provide a --version flag"

echo "Installed at:"
conda env list

echo "jasmine should be under ~/.conda/envs/${ENV_NAME}/bin/"
conda deactivate
