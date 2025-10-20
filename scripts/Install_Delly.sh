#!/bin/bash
#SBATCH --job-name=install_delly
#SBATCH --output=/storage/homefs/kw23y068/logfiles/install_delly_%j.log
#SBATCH --error=/storage/homefs/kw23y068/logfiles/install_delly_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

# download Delly from:
# https://github.com/dellytools/delly/releases/
# Version 1.5.0 used here

DELLY_LOCATION=/storage/homefs/kw23y068/software/Delly/
if [ ! -d "$DELLY_LOCATION" ]
then
    mkdir $DELLY_LOCATION
fi

# Go to installation directory
cd $DELLY_LOCATION

wget https://github.com/dellytools/delly/releases/download/v1.5.0/delly_v1.5.0.sif

echo "SIF Image downloaded to $DELLY_LOCATION/delly_v1.5.0.sif"

echo "Verifying sha256 checksum; Download (1st) and GitHub (2nd):"
echo $(sha256sum delly_v1.5.0.sif)
echo "sha256:5cead52d8c9d0bb763977523c8c9136ca091775dc3f758c026129f1b648937c9"


