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
# May be altered for feeding in arguments to script
# Define which batch of samples to process
BATCH_START=1    # Change to 81, 161, 241, 321, 401 for other batches
BATCH_END=80     # Change to 160, 240, 320, 400, 480 for other batches

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

# Copy files - Note: existing files will be overwritten
cp -f "$DELLY" "$MANTA" "$SMOOVE" .

# Safety check - we all know typos and stuff
for f in "$DELLY" "$MANTA" "$SMOOVE" "$REFERENCE"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Required file not found: $f" >&2
        echo "Check pathway and previous callers."
        exit 1
    fi
done

gunzip -f *.vcf.gz

echo "======================================"
echo "Checking VCF validity and filtering..."
echo "======================================"

for vcf in *.vcf; do
    echo "Processing $vcf"
    
    # Diagnostic: Check for malformed lines
    echo "  Checking for lines with too few fields..."
    malformed_count=$(grep -v "^#" "$vcf" | awk -F'\t' '{if(NF<8) print}' | wc -l)
    if [ "$malformed_count" -gt 0 ]; then
        echo "  WARNING: Found $malformed_count malformed lines in $vcf"
    fi
    
    # Diagnostic: Check for extremely long sequences
    echo "  Checking for very long variants (>5000bp)..."
    grep -v "^#" "$vcf" | awk -F'\t' '{if(length($4)>5000 || length($5)>5000) print "  Very long variant at",$1,$2,"REF length:",length($4),"ALT length:",length($5)}' | head -5
    
    # Filter: Keep only lines with at least 8 fields (valid VCF format)
    grep "^#" "$vcf" > "${vcf}.filtered"
    grep -v "^#" "$vcf" | awk -F'\t' 'NF>=8' >> "${vcf}.filtered"
    
    # Check if filtering was successful
    filtered_lines=$(grep -v "^#" "${vcf}.filtered" | wc -l)
    original_lines=$(grep -v "^#" "$vcf" | wc -l)
    
    if [ "$filtered_lines" -eq 0 ]; then
        echo "  ERROR: Filtering removed all variants from $vcf"
        exit 1
    fi
    
    echo "  Original variants: $original_lines"
    echo "  Filtered variants: $filtered_lines"
    echo "  Removed: $((original_lines - filtered_lines)) variants"
    
    # Replace original with filtered version
    mv "${vcf}.filtered" "$vcf"
done

echo "======================================"
echo "VCF validation and filtering complete"
echo "======================================"

# Create list of VCFs for Jasmine
ls *.vcf > Manta_Smoove_Delly.list

echo "VCF files to merge:"
cat Manta_Smoove_Delly.list

# Start with Jasmine
module load Anaconda3/2024.02-1
source $(conda info --base)/etc/profile.d/conda.sh
conda activate jasmine

module load BCFtools/1.12-GCC-10.3.0
module load SAMtools/1.13-GCC-10.3.0

# Input a list of vcf files (should not be zipped)
BASE="Manta_Smoove_Delly"
VCFlist=${BASE}.list
OUTfile=${BASE}.vcf

echo "======================================"
echo "Running Jasmine..."
echo "======================================"


# We get an error with --output_genotypes \  --> removed for now
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
    threads=${SLURM_CPUS_PER_TASK}

# Check if Jasmine succeeded
if [ ! -f "${OUTfile}" ]; then
    echo "ERROR: Jasmine did not create output file ${OUTfile}"
    exit 1
fi

echo "Jasmine completed successfully"

conda deactivate

module -q reset
module load HTSlib/1.12-GCC-10.3.0
module load BCFtools/1.12-GCC-10.3.0
module load VCFtools/0.1.16-GCC-10.3.0

echo "======================================"
echo "Compressing and indexing output..."
echo "======================================"


echo "Sorting VCF..."
bcftools sort -o ${OUTfile}.sorted ${OUTfile}


if [ ! -f "${OUTfile}.sorted" ]; then
    echo "ERROR: Sorting failed"
    exit 1
fi

# Compress the output
bgzip -f ${OUTfile}.sorted

# Check if compression succeeded
if [ ! -f "${OUTfile}.sorted.gz" ]; then
    echo "ERROR: bgzip failed to create ${OUTfile}.sorted.gz"
    exit 1
fi

# Index the compressed file
tabix -p vcf -f ${OUTfile}.sorted.gz

# Check if indexing succeeded
if [ ! -f "${OUTfile}.sorted.gz.tbi" ]; then
    echo "ERROR: tabix failed to create index"
    exit 1
fi

echo "======================================"
echo "Pipeline completed successfully!"
echo "Output: ${OUTPUT_DIR}/${OUTfile}.gz"
echo "======================================"
