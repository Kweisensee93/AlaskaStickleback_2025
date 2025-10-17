#!/bin/bash
#SBATCH --job-name=samplot_sv
#SBATCH --output=/storage/homefs/kw23y068/logfiles/samplot_sv_%j.log
#SBATCH --error=/storage/homefs/kw23y068/logfiles/samplot_sv_%j.err
#SBATCH --time=04:00:00
#SBATCH --partition=epyc2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G

# Load modules
module load SAMtools/1.13-GCC-10.3.0

# Activate your environment if Samplot is installed in a venv
# source /path/to/venv/bin/activate

# Directories
BAM_DIR="/path/to/bam_files"       # BAM files for your samples
VCF_DIR="/path/to/vcf_files"       # VCFs containing SVs
OUT_DIR="/path/to/samplot_output"  # Where images will be saved

mkdir -p $OUT_DIR

# Loop over VCF files
for vcf in "$VCF_DIR"/*.vcf; do
    base=$(basename "$vcf" .vcf)
    echo "Processing VCF: $base"

    # Loop over BAMs
    for bam in "$BAM_DIR"/*.bam; do
        sample=$(basename "$bam" .bam)
        sample_out="$OUT_DIR/${base}_${sample}"

        mkdir -p "$sample_out"

        # Generate Samplot images
        samplot plot \
            -n "$bam" \
            -v "$vcf" \
            -o "$sample_out" \
            -t DEL,DUP,INV,INS,TRA \
            --min-mapq 20 \
            --min-sv-len 50 \
            --maxreadlen 150 \
            --plotsize 8
    done
done

echo "Samplot generation complete. Images are in $OUT_DIR"
