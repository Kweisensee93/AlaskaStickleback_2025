#!/bin/bash
#SBATCH --job-name=add_R_packages
#SBATCH --output=/storage/homefs/kw23y068/logfiles/add_R_packages_%j.log
#SBATCH --error=/storage/homefs/kw23y068/logfiles/add_R_packages_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=32G
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

conda activate R_all

# Packages added so far:
#conda install -y -c bioconda bioconductor-genomicranges bioconductor-rsamtools
# No need for iranges, it is installed as dependency of GenomicRanges by bioconductor
#conda install -y -c bioconda bioconductor-iranges

conda deactivate
