#!/bin/bash
#SBATCH --job-name=install_samplot
#SBATCH --output=/storage/homefs/kw23y068/logfiles/install_samplot_%j.log
#SBATCH --error=/storage/homefs/kw23y068/logfiles/install_samplot_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

# Exit on error
set -e

# We don't have recommended Bioconda per default. Let's make it happen with Anaconda3 module.
module load Anaconda3/2024.02-1

# Make sure conda is initialized (if eval is not working use source line)
# source $(conda info --base)/etc/profile.d/conda.sh
# Initialize conda for bash shell
eval "$(conda shell.bash hook)"

echo "=== Creating Samplot environment ==="
# create enrinment with channels in this order, then install samplot
conda create -y -n samplot-env -c conda-forge -c bioconda -c defaults samplot

echo "=== Activating environment ==="
conda activate samplot-env

echo "=== Checking Samplot installation ==="
samplot --version || samplot --help

echo "=== Conda environment summary ==="
conda info --envs

echo "=== Samplot executable location ==="
which samplot

echo "=== Samplot installation completed successfully! ==="
echo " To use samplot in future jobs, add these lines (hopefully):"
echo "  module load Anaconda3/2024.02-1"
echo "  eval \"\$(conda shell.bash hook)\""
echo "  conda activate samplot-env"
