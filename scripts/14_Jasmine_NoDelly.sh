#!/bin/bash
#SBATCH --job-name=jasmine_NoDelly
#SBATCH --output=/storage/homefs/kw23y068/logfiles/jasmine_NoDelly_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/jasmine_NoDelly_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=epyc2
#SBATCH --account=gratis
#SBATCH --wckey=noop

####################
# PARAMETERS TO SET
####################
# May be altered for feeding in arguments to script
# Define which batch of samples to process
BATCH_START=1    # Change to 81, 161, 241, 321, 401 for other batches
BATCH_END=80     # Change to 160, 240, 320, 400, 480 for other batches

START_PAD=$(printf "%03d" "${BATCH_START}")
END_PAD=$(printf "%03d" "${BATCH_END}")

PROJECT_DIR=/storage/scratch/iee_evol/kw23y068
OUTPUT_DIR=${PROJECT_DIR}/jasmine_NoDelly_${START_PAD}_${END_PAD}
mkdir -p ${OUTPUT_DIR}
STORAGE_DIR="/storage/research/iee_evol/Korbi"
REFERENCE="${STORAGE_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna"

# Get sample: Go to ouptut directory and copy from SV callers
cd ${OUTPUT_DIR}
# Define expected files
#DELLY="${PROJECT_DIR}/Delly_Merge_${START_PAD}_${END_PAD}/Delly_merged_filtered_${START_PAD}_${END_PAD}.vcf.gz"
MANTA="${PROJECT_DIR}/manta_${START_PAD}_${END_PAD}_merged/manta_${START_PAD}_${END_PAD}_merged.vcf.gz"
SMOOVE="${PROJECT_DIR}/smoove_${START_PAD}_${END_PAD}/smoove_${START_PAD}_${END_PAD}_filtered.vcf.gz"


# Copy files - Note: existing files will be overwritten
cp -f "$MANTA" "$SMOOVE" .

# Safety check - we all know typos and stuff
for f in "$MANTA" "$SMOOVE" "$REFERENCE"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Required file not found: $f" >&2
        echo "Check pathway and previous callers."
        exit 1
    fi
done

# Load modules for preprocessing
module load BCFtools/1.12-GCC-10.3.0

echo "======================================"
echo "Preprocessing VCF files..."
echo "======================================"

# Decompress and preprocess VCFs
for vcf_gz in *.vcf.gz; do
    vcf="${vcf_gz%.gz}"
    echo "Processing $vcf_gz"
    
    # Decompress
    gunzip -f "$vcf_gz"
    
    # Check if file exists after decompression
    if [[ ! -f "$vcf" ]]; then
        echo "ERROR: Failed to decompress $vcf_gz"
        exit 1
    fi
    
    # Fix the VCF:
    # 1. Keep header
    # 2. Filter body: only valid lines with >=8 fields
    # 3. Remove lines with comma-separated values in ID field (Manta merge artifacts)
    grep "^#" "$vcf" > "${vcf}.fixed"
    grep -v "^#" "$vcf" | \
        awk -F'\t' 'NF>=8' | \
        awk -F'\t' '$3 !~ /,/' >> "${vcf}.fixed"
    
    # Check results
    original_lines=$(grep -v "^#" "$vcf" | wc -l)
    fixed_lines=$(grep -v "^#" "${vcf}.fixed" | wc -l)
    removed_lines=$((original_lines - fixed_lines))
    
    echo "  Original variants: $original_lines"
    echo "  Fixed variants: $fixed_lines"
    echo "  Removed (malformed/multi-ID): $removed_lines"
    
    if [ "$fixed_lines" -eq 0 ]; then
        echo "  ERROR: All variants removed from $vcf"
        exit 1
    fi
    
    # Replace with fixed version
    mv "${vcf}.fixed" "$vcf"
done

# Create list of VCFs for Jasmine
ls *.vcf > Manta_Smoove.list

# Start with Jasmine

module load Anaconda3/2024.02-1
source $(conda info --base)/etc/profile.d/conda.sh
conda activate jasmine

# I leave out the module reset to avoid issues with conda
#module -q reset
# I leave the module load for BCFtools and SAMtools - may be needed by jasmine ?
#module load BCFtools/1.9-intel-2018b
#module load SAMtools/1.9-GCC-7.3.0-2.30
module load BCFtools/1.12-GCC-10.3.0
module load SAMtools/1.13-GCC-10.3.0

#input a list of vcf files (should not be zip)
BASE="Manta_Smoove"

VCFlist=${BASE}.list
OUTfile=${BASE}.vcf

# max_dist_linear may be set to 0.1 ?
# Why use min_dist=50 ? to avoid merging very close calls ? --> for smaller SVs, we may need to adapt
# original creators (the ones who got copied by the ones we copy) didnt have -- allow_intrasample
# maybe add --output_genotypes
# reference used --default_zero_genotype \
# reference used --mutual_distance \
jasmine file_list=${VCFlist} \
    out_file=${OUTfile} \
    genome_file=${REFERENCE} \
    --ignore_strand \
    --max_dist_linear=0.5 \
    --max_dist=1000 \
    --min_dist=50 \
    --use_end \
    --normalize_type \
    --allow_intrasample \
    --output_genotypes \
    threads=${SLURM_CPUS_PER_TASK}

conda deactivate

module -q reset
module load HTSlib/1.12-GCC-10.3.0
module load BCFtools/1.12-GCC-10.3.0
module load VCFtools/0.1.16-GCC-10.3.0

bgzip -f -c ${OUTfile}
tabix -f ${OUTfile}.gz
