#!/bin/bash
#SBATCH --job-name=Delly_BCFtools_merge
#SBATCH --output=/storage/homefs/kw23y068/logfiles/Delly_BCFtools_merge_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/Delly_BCFtools_merge_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --mem=16G
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


echo "Delly merge with bcftools"
echo "Date: $(date)"
echo "Host: $(hostname)"

module load BCFtools/1.12-GCC-10.3.0
module load SAMtools/1.13-GCC-10.3.0

PROJECT_DIR=/storage/research/iee_evol/Korbi
OUTPUT_DIR=/storage/scratch/iee_evol/kw23y068/

# Set up merge directory
START_PAD=$(printf "%03d" "${BATCH_START}")
END_PAD=$(printf "%03d" "${BATCH_END}")
OUT_DIR="${OUTPUT_DIR}/Delly_Merge_BCFtools_${START_PAD}_${END_PAD}"
if [ ! -d "${OUT_DIR}" ]; then
    mkdir -p "${OUT_DIR}"
fi

# Create a subdirectory for converted VCFs
VCF_DIR="${OUT_DIR}/converted_vcfs"
if [ ! -d "${VCF_DIR}" ]; then
    mkdir -p "${VCF_DIR}"
fi

# Define output merged BCF and VCF filenames
OUTPUT_VCF="${OUT_DIR}/Delly_BCFmerged_${START_PAD}_${END_PAD}.vcf"

# Construct list of input VCFs
SAMPLE_LIST=/storage/research/iee_evol/Korbi/output/bams_full.csv

VCF_FILES=()
MISSING_BCF=0
CONVERSION_FAILED=0

echo "Converting BCF files to VCF format..."

# Loop through samples in the batch
for i in $(seq ${BATCH_START} ${BATCH_END}); do
    SAMPLE_PATH=$(sed -n "${i}p" "${SAMPLE_LIST}")
    if [ -z "${SAMPLE_PATH}" ]; then
        echo "WARNING: Could not read sample at line ${i}"
        continue
    fi

    SAMPLE_FILE=$(basename "${SAMPLE_PATH}")
    SAMPLE_NAME="${SAMPLE_FILE%.fixmate.coordsorted.bam}"

    # Add warning if BCF file is missing
    BCF="${OUTPUT_DIR}/Delly/${SAMPLE_NAME}/${SAMPLE_NAME}.bcf"

        if [ ! -f "${BCF}" ]; then
        echo "WARNING: BCF file not found: ${BCF}"
        MISSING_BCF=$((MISSING_BCF + 1))
        continue
    fi

    # Define output VCF filename
    VCF="${VCF_DIR}/${SAMPLE_NAME}.vcf.gz"
    # Convert BCF to compressed VCF
    echo "  [${i}] Converting ${SAMPLE_NAME}..."
    bcftools view "${BCF}" -Oz -o "${VCF}"
    # Check if conversion was successful
    if [ $? -ne 0 ]; then
        echo "ERROR: BCF to VCF conversion failed for ${SAMPLE_NAME}"
        CONVERSION_FAILED=$((CONVERSION_FAILED + 1))
        continue
    fi
    # Index the VCF file
    bcftools index -t "${VCF}"
    if [ $? -ne 0 ]; then
        echo "ERROR: VCF indexing failed for ${SAMPLE_NAME}"
        CONVERSION_FAILED=$((CONVERSION_FAILED + 1))
        continue
    fi
    
    # Append to VCF file list
    VCF_FILES+=("${VCF}")


done

echo ""
echo "=== Conversion Summary ==="
echo "Total samples processed: $((BATCH_END - BATCH_START + 1))"
echo "Missing BCF files: ${MISSING_BCF}"
echo "Conversion failures: ${CONVERSION_FAILED}"
echo "Successfully converted VCF files: ${#VCF_FILES[@]}"
echo ""

if [ ${#VCF_FILES[@]} -eq 0 ]; then
    echo "ERROR: No VCF files available for merging"
    exit 1
fi

# Merge VCF files
echo "Merging ${#VCF_FILES[@]} VCF files..."
bcftools merge \
    --threads ${SLURM_CPUS_PER_TASK} \
    --merge none \
    --output-type z \
    --output ${OUTPUT_VCF} \
    "${VCF_FILES[@]}"

# Check if merge was successful
if [ $? -ne 0 ]; then
    echo "ERROR: VCF merging failed"
    exit 1
fi

bgzip ${OUTPUT_VCF}

# Index the merged VCF
echo "Indexing merged VCF..."
bcftools index -t ${OUTPUT_VCF}.gz

if [ $? -ne 0 ]; then
    echo "ERROR: Merged VCF indexing failed"
    exit 1
fi

echo ""
echo "=== Merge Complete ==="
echo "Output file: ${OUTPUT_VCF}"
echo "Output index: ${OUTPUT_VCF}.tbi"
echo "Conversion directory: ${VCF_DIR}"
echo "Date: $(date)"
echo ""
