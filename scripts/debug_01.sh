#!/bin/bash
#SBATCH --job-name=debug
#SBATCH --output=./logfiles/debug_%j.log
#SBATCH --error=./logfiles/debug_%j.err
#SBATCH --time=00:10:00
#SBATCH --partition=epyc2
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

# Directories
OUT_DIR="/storage/research/iee_evol/Korbi"

# Input files
BAM_LIST="${OUT_DIR}/samples_bamfiles.txt"
#RAW_LIST="${OUT_DIR}/samples_rawdata.txt"

# Create or clear input files
> $BAM_LIST
#> $RAW_LIST

# List BAM and RAW files

ls $OUT_DIR/bams/* > $BAM_LIST
find "${OUT_DIR}/bams" -maxdepth 1 -type f | tee "$BAM_LIST" | head
