#!/bin/bash
#SBATCH --job-name=debug_sample_list_ConsensuSV
#SBATCH --output=/storage/research/iee_evol/Korbi/logfiles/%j_ConsensuSV_samplelist.log
#SBATCH --error=/storage/research/iee_evol/Korbi/logfiles/%j_ConsensuSV_samplelist.err
#SBATCH --time=00:30:00
#SBATCH --partition=epyc2
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

OUT_DIR=/storage/scratch/iee_evol/kw23y068
PROJECT_DIR=/storage/research/iee_evol/Korbi
RAW_DIR="${PROJECT_DIR}/raw_data"
sample_list=${OUT_DIR}/samples_full.csv

# Create or clear the sample list file
> ${sample_list}

echo "Creating samplelist from ${RAW_DIR}"

# Loop through all R1 fastq files using ls (works with symlinks by default)
for R1 in "${RAW_DIR}"/*_R1_*.fastq.gz; do
    [[ -e "$R1" ]] || continue  # Skip if no matches
    
    R2="${R1/_R1_/_R2_}"
    if [[ -f "$R2" ]]; then
        sample=$(basename "$R1" | sed -E 's/_S[0-9]+_L[0-9]+_R1_001\.fastq\.gz$//')
        echo "${sample},${R1},${R2}" >> "${sample_list}"
    else
        echo "Warning: No matching R2 found for $R1" >&2
        echo "Expected: ${R2}" >&2
    fi
done

echo "Sample list for ConsensuSV created at: ${sample_list}"
echo "Total samples: $(wc -l < ${sample_list})"
