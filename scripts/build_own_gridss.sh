#!/bin/bash
#SBATCH --job-name=gridss_fix
#SBATCH --output=/storage/homefs/kw23y068/logfiles/gridss_fix_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/gridss_fix_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=epyc2

# use apptainer on cluster

# Define paths
SANDBOX_DIR="/storage/homefs/kw23y068/software/gridss/gridss_patch_tmp/sandbox"
NEW_IMAGE="/storage/homefs/kw23y068/software/gridss/gridss_fixed.sif"

# Update samtools in the sandbox
apptainer exec --writable $SANDBOX_DIR bash -c '
  cd /tmp
  curl -sSL https://github.com/samtools/samtools/releases/download/1.22.1/samtools-1.22.1.tar.bz2 | tar -xj
  cd samtools-1.22.1
  ./configure --prefix=/opt/samtools
  make -j4 all && make install
  rm -rf /tmp/samtools-1.22.1*
'

# Rebuild into a new .sif image
apptainer build $NEW_IMAGE $SANDBOX_DIR
