#!/bin/bash
#SBATCH --job-name=fix_smoove_conda
#SBATCH --output=/storage/homefs/kw23y068/logfiles/fix_smoove_conda_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/fix_smoove_conda_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=80G
#SBATCH --time=02:00:00
#SBATCH --partition=epyc2

## As per GitHub smoove issue 230:
# https://github.com/brentp/smoove/issues/230
# The error **could not import: bam_hdr_destroy**
# comes from htslib version issue, that can't be fixed by updating htslib
# potential fix may be a newer version of duphold

module load Anaconda3/2024.02-1

# Ensure conda commands work in batch mode
source $(conda info --base)/etc/profile.d/conda.sh

# Set conda to use a writable location for environments
export CONDA_ENVS_PATH=$HOME/.conda/envs
export CONDA_PKGS_DIRS=$HOME/.conda/pkgs
mkdir -p $CONDA_ENVS_PATH $CONDA_PKGS_DIRS

SMOOVE_TMP="/storage/homefs/kw23y068/software/smoove_tmp"
mkdir -p ${SMOOVE_TMP}

# Configure channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

ENV_NAME="smoove"
#conda create -y -n ${ENV_NAME} smoove
conda activate ${ENV_NAME}
# as of 27.Nov.2025 Version 0.2.3 is the latest one
# Bioconda only provides 0.2.1 , but the fix needs latest!

# Check installation
echo "Check duphold before update"
which duphold
duphold --version
duphold --help


cd ${SMOOVE_TMP}

echo "Downloading duphold v0.2.3..."
wget -q https://github.com/brentp/duphold/releases/download/v0.2.3/duphold \
    -O duphold

chmod +x duphold

echo "Overwriting duphold binary in conda environment..."
#cp duphold $(conda info --base)/envs/${ENV_NAME}/bin/duphold
cp duphold $CONDA_ENVS_PATH/${ENV_NAME}/bin/duphold


# Check installation
echo "check duphold after update"
which duphold
duphold --version
ldd $(which duphold)

rm -rf ${SMOOVE_TMP}

echo "smoove should be under ~/.conda/envs/${ENV_NAME}/bin/"
conda deactivate
