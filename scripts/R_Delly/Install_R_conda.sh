#!/bin/bash
#SBATCH --job-name=build_R_conda
#SBATCH --output=/storage/homefs/kw23y068/logfiles/install_R_conda_%j.log
#SBATCH --error=/storage/homefs/kw23y068/logfiles/install_R_conda_%j.err
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

# Configure conda channels
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels defaults


# Create R environment that may be used for several pipelines
ENV_R_all="R_all"
# Use R to latest version of 4.4.X (Bioconductor works up to 4.5.0 as of 07.Nov.25)
# In our case we use 4.4.3 retrieved from
# conda search -c conda-forge r-base
conda create -y -n ${ENV_R_all} r-base=4.4.3
conda activate ${ENV_R_all}
conda deactivate
