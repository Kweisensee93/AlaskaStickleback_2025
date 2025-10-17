#!/bin/bash
#SBATCH --job-name=sanity_check_samples
#SBATCH --output=/storage/research/iee_evol/Korbi/logfiles/sanity_check_%j.log
#SBATCH --error=/storage/research/iee_evol/Korbi/logfiles/sanity_check_%j.err
#SBATCH --time=00:10:00
#SBATCH --partition=epyc2
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

# Directories
OUT_DIR="/storage/research/iee_evol/Korbi"

# Input files
BAM_LIST="${OUT_DIR}/samples_bamfiles.txt"
RAW_LIST="${OUT_DIR}/samples_rawdata.txt"

# Create or clear input files
> $BAM_LIST
> $RAW_LIST

# List BAM and RAW files

#ls $OUT_DIR/bams_real/* > $BAM_LIST
#ls $OUT_DIR/raw_data/* > $RAW_LIST
ls "$OUT_DIR"/bams_real/* | xargs -n 1 basename > "$BAM_LIST"
ls "$OUT_DIR"/raw_data/* | xargs -n 1 basename > "$RAW_LIST"


# Output files
FULL="${OUT_DIR}/samples_full.txt"
MISSING="${OUT_DIR}/samples_missing.txt"

# Create or clear output files
> $FULL
> $MISSING

# Extract unique sample prefixes (first 13 chars like FG_LP_19T_084)
bam_prefixes=$(awk '{print substr($0,1,13)}' $BAM_LIST | sort | uniq)
raw_prefixes=$(awk '{print substr($0,1,13)}' $RAW_LIST | sort | uniq)

# Combine all sample prefixes
all_samples=$(echo -e "$bam_prefixes\n$raw_prefixes" | sort | uniq)

# Check counts for each sample
for sample in $all_samples; do
    bam_count=$(grep -c "^$sample" $BAM_LIST)
    raw_count=$(grep -c "^$sample" $RAW_LIST)

    if [[ $bam_count -eq 4 && $raw_count -eq 2 ]]; then
        echo $sample >> $FULL
    else
        echo -e "$sample\tBAM=$bam_count\tRAW=$raw_count" >> $MISSING
    fi
done

echo "Sanity check complete!"
echo "Full samples: $FULL"
echo "Missing/incomplete samples: $MISSING"
