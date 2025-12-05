#!/bin/bash
#SBATCH --job-name=jasmine
#SBATCH --output=/storage/homefs/kw23y068/logfiles/jasmine_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/jasmine_%j.err
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
BATCH_START=1
BATCH_END=80

START_PAD=$(printf "%03d" "${BATCH_START}")
END_PAD=$(printf "%03d" "${BATCH_END}")

PROJECT_DIR=/storage/scratch/iee_evol/kw23y068
OUTPUT_DIR=${PROJECT_DIR}/jasmine_${START_PAD}_${END_PAD}
mkdir -p ${OUTPUT_DIR}

STORAGE_DIR="/storage/research/iee_evol/Korbi"
REFERENCE="${STORAGE_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna"

# Get sample: Go to output directory and copy from SV callers
cd ${OUTPUT_DIR}

# Define expected files
DELLY="${PROJECT_DIR}/Delly_Merge_${START_PAD}_${END_PAD}/Delly_merged_filtered_${START_PAD}_${END_PAD}.vcf.gz"
MANTA="${PROJECT_DIR}/manta_${START_PAD}_${END_PAD}_merged/manta_${START_PAD}_${END_PAD}_merged.vcf.gz"
SMOOVE="${PROJECT_DIR}/smoove_${START_PAD}_${END_PAD}/smoove_${START_PAD}_${END_PAD}_filtered.vcf.gz"

# Clean leftover files from previous runs
rm -rf *

# Copy files
cp -f "$DELLY" "$MANTA" "$SMOOVE" .

# Safety check
for f in "$DELLY" "$MANTA" "$SMOOVE" "$REFERENCE"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Required file not found: $f" >&2
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

echo "======================================"
echo "VCF preprocessing complete"
echo "======================================"

# Create list of VCFs for Jasmine
> Manta_Smoove_Delly.list
ls *smoove*.vcf >> Manta_Smoove_Delly.list
ls *manta*.vcf >> Manta_Smoove_Delly.list
ls *Delly*.vcf >> Manta_Smoove_Delly.list

echo "VCF files to merge:"
cat Manta_Smoove_Delly.list

# Start Jasmine
#module -q reset
module load Anaconda3/2024.02-1
source $(conda info --base)/etc/profile.d/conda.sh
conda activate jasmine

module load SAMtools/1.13-GCC-10.3.0

BASE="Manta_Smoove_Delly"
VCFlist=${BASE}.list
OUTfile=${BASE}.vcf

echo "======================================"
echo "Running Jasmine..."
echo "======================================"

#excluded: --output_genotypes \
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
    threads=${SLURM_CPUS_PER_TASK}

# Check if Jasmine succeeded
if [ ! -f "${OUTfile}" ]; then
    echo "ERROR: Jasmine did not create output file ${OUTfile}"
    exit 1
fi

echo "Jasmine completed successfully"

conda deactivate

# Reload modules for post-processing
module -q reset
module load HTSlib/1.12-GCC-10.3.0
module load BCFtools/1.12-GCC-10.3.0

echo "======================================"
echo "Post-processing Jasmine output..."
echo "======================================"

# Additional cleanup: remove any remaining problematic entries
echo "Cleaning Jasmine output..."
grep "^#" "${OUTfile}" > "${OUTfile}.clean"
grep -v "^#" "${OUTfile}" | \
    awk -F'\t' 'NF>=8 && $3 !~ /,/' >> "${OUTfile}.clean"

clean_count=$(grep -v "^#" "${OUTfile}.clean" | wc -l)
echo "Variants after final cleaning: $clean_count"

if [ "$clean_count" -eq 0 ]; then
    echo "ERROR: No variants remaining after cleaning"
    exit 1
fi


mv "${OUTfile}.clean" "${OUTfile}"

echo "======================================"
echo "Sorting and compressing..."
echo "======================================"


# Sort the VCF
echo "Sorting VCF..."
bcftools sort -o ${OUTfile}.sorted ${OUTfile}

if [ ! -f "${OUTfile}.sorted" ]; then
    echo "ERROR: Sorting failed"
    exit 1
fi

# Compress
bgzip -f ${OUTfile}.sorted

if [ ! -f "${OUTfile}.sorted.gz" ]; then
    echo "ERROR: Compression failed"
    exit 1
fi

# Index
tabix -p vcf -f ${OUTfile}.sorted.gz

if [ ! -f "${OUTfile}.sorted.gz.tbi" ]; then
    echo "ERROR: Indexing failed"
    exit 1
fi

echo "======================================"
echo "Pipeline completed successfully!"
echo "Output: ${OUTPUT_DIR}/${OUTfile}.sorted.gz"
echo "======================================"
