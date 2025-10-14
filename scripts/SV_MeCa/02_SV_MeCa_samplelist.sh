#!/bin/bash
#SBATCH --job-name=svmeca_samplelist
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/storage/research/iee_evol/Korbi/logfiles/%j_SvMeCa_samplelist.log
#SBATCH --error=/storage/research/iee_evol/Korbi/logfiles/%j_SvMeCa_samplelist.err
#SBATCH --partition=epyc2

OUT_DIR=/storage/scratch/iee_evol/kw23y068
PROJECT_DIR=/storage/research/iee_evol/Korbi
BAMS_DIR="${PROJECT_DIR}/bams_real"

bams_list=${OUT_DIR}/bams_full.csv
# Create or clear the sample list file
> ${bams_list}

echo "Creating samplelist from ${DAMS_DIR}"

find "${BAMS_DIR}" -type f -name "*.bam" | sort > "${bams_list}"

num_files=$(ls "${BAMS_DIR}" | wc -l)
quarter_files=$(( num_files / 4 ))

echo "Sample list written to ${bams_list}"
echo "Total BAM files found: $(wc -l < "${bams_list}")"
echo "Files: $num_files (quarter = $quarter_files)"
