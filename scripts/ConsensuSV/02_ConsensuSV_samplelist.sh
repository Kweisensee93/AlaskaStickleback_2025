#!/bin/bash
#SBATCH --job-name=sample_list_ConsensuSV
#SBATCH --output=/storage/research/iee_evol/Korbi/logfiles/%j_ConsensuSV_samplelist.log
#SBATCH --error=/storage/research/iee_evol/Korbi/logfiles/%j_ConsensuSV_samplelist.err
#SBATCH --time=00:30:00
#SBATCH --partition=epyc2
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

OUT_DIR=/storage/scratch/iee_evol/kw23y068
PROJECT_DIR=/storage/research/iee_evol/Korbi
RAW_DIR="${PROJECT_DIR}/raw_data"

# sample list for ConsensuSV - needs to be a csv without headers
# structure from https://github.com/SFGLab/ConsensuSV-pipeline needs to be:
# samplename,fastq1,fastq2
sample_list=${OUT_DIR}/samples_full.csv
# Create or clear the sample list file
> ${sample_list}

## samplenames were retrieved as follows (in sanity_check_rawdata.sh)
## ls "$PROJECT_DIR"/raw_data/* | xargs -n 1 basename > "$sample_list"
echo "Creating samplelist from ${RAW_DIR}"

# samples should look like:
#SR_CL_19T_019_2-3267191,
#/storage/research/iee_evol/Korbi/raw_data/SR_CL_19T_019_2-3267191_S32_L002_R1_001.fastq.gz,
#/storage/research/iee_evol/Korbi/raw_data/SR_CL_19T_019_2-3267191_S32_L002_R2_001.fastq.gz
# Loop through all R1 fastq files (they are fed in at the end of the while loop via <<(...) )
while IFS= read -r R1; do
    R2="${R1/_R1_/_R2_}"

    if [[ -f "$R2" ]]; then
        sample=$(basename "$R1" | sed -E 's/_S[0-9]+_L[0-9]+_R1_001\.fastq\.gz//')
        echo "${sample},${R1},${R2}" >> "${sample_list}"
    else
        echo "Warning: No matching R2 found for $R1" >&2
        echo "Expected: ${R2}" >&2
    fi
    # use find -L to follow symlinks
done < <(find -L "${RAW_DIR}" -type f -name "*_R1_*.fastq.gz" | sort)

# Verify the sample list
num_lines=$(wc -l < "${sample_list}")
# try counting with ls instead of find (should be the same, but find may get errors)
num_files=$(ls "${RAW_DIR}" | wc -l)
half_files=$(( num_files / 2 ))

echo "Samples: $num_lines"
echo "Files: $num_files (half = $half_files)"
echo "Sample list for ConsensuSV created at: ${sample_list}"
