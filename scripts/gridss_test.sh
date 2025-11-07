#!/bin/bash
#SBATCH --job-name=test_gridss
#SBATCH --output=test_gridss_%j.log
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=epyc2

module load Anaconda3/2024.02-1
source $(conda info --base)/etc/profile.d/conda.sh
conda activate gridss

which gridss
gridss --version
