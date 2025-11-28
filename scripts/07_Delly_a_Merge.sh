#!/bin/bash
#SBATCH --job-name=Delly_Merge
#SBATCH --output=/storage/homefs/kw23y068/logfiles/Delly_Merge_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/Delly_Merge_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=epyc2


####################
# PARAMETERS TO SET
####################
# May be altered for feeding in arguments to script
# Define which batch of samples to process
BATCH_START=1    # Change to 81, 161, 241, 321, 401 for other batches
BATCH_END=80     # Change to 160, 240, 320, 400, 480 for other batches


echo "Delly merge"
echo "Date: $(date)"
echo "Host: $(hostname)"

DELLY_SIF=/storage/homefs/kw23y068/software/Delly/delly_v1.5.0.sif
PROJECT_DIR=/storage/research/iee_evol/Korbi
OUTPUT_DIR=/storage/scratch/iee_evol/kw23y068/

# Set up merge directory
START_PAD=$(printf "%03d" "${BATCH_START}")
END_PAD=$(printf "%03d" "${BATCH_END}")
OUT_DIR="${OUTPUT_DIR}/Delly_Merge_${START_PAD}_${END_PAD}"
if [ ! -d "${OUT_DIR}" ]; then
    mkdir -p "${OUT_DIR}"
fi

# Define output merged BCF and VCF filenames
OUTPUT_BCF="${OUT_DIR}/Delly_merged_${START_PAD}_${END_PAD}.bcf"
OUTPUT_VCF="${OUT_DIR}/Delly_merged_${START_PAD}_${END_PAD}.vcf"

# Construct list of input VCFs
SAMPLE_LIST=/storage/research/iee_evol/Korbi/output/bams_full.csv

BCF_FILES=()
MISSING_BCF=0

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

    # Append to BCF file list
        BCF_FILES+=("${BCF}")
    echo "  [${i}] ${SAMPLE_NAME}"
done

echo ""
echo "Collected ${#BCF_FILES[@]} BCF files"

if [ ${#BCF_FILES[@]} -eq 0 ]; then
    echo "ERROR: No BCF files collected"
    exit 1
fi

apptainer exec \
    --bind ${PROJECT_DIR} \
    --bind ${OUTPUT_DIR} \
    ${DELLY_SIF} delly merge \
    -o ${OUTPUT_BCF} \
    ${BCF_FILES[@]}

# Convert BCF to VCF
module load BCFtools/1.12-GCC-10.3.0
module load SAMtools/1.13-GCC-10.3.0

bcftools view ${OUTPUT_BCF} -Ov -o ${OUTPUT_VCF}

# Check if conversion was successful
if [ $? -ne 0 ]; then
    echo "ERROR: BCF to VCF conversion failed"
    exit 1
fi

bgzip ${OUTPUT_VCF}
tabix -p vcf ${OUTPUT_VCF}.gz

echo "  BCF files merged: ${#BCF_FILES[@]}"
echo "  Missing BCF files: ${MISSING_BCF}"
echo "Done: $(date)"

# No re-genotyping step included here, which is in the Delly pipeline
# another rerun of delly call with -v on the merged VCF would be needed
# This can be added as a separate script if required
# Example command:
# delly call -g reference.fasta -v merged_sites.bcf -o sample1.geno.bcf sample1.bam
# re-merging again:
# delly merge -m samples -o final_multisample.bcf sample1.geno.bcf sample2.geno.bcf ...