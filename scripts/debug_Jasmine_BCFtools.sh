#!/bin/bash
#SBATCH --job-name=jasmine_debug
#SBATCH --output=/storage/homefs/kw23y068/logfiles/jasmine_debug_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/jasmine_debug_%j.err
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
#mkdir -p ${OUTPUT_DIR}

#STORAGE_DIR="/storage/research/iee_evol/Korbi"
#REFERENCE="${STORAGE_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna"
# Get sample: Go to output directory and copy from SV callers
cd ${OUTPUT_DIR}

BASE="Manta_Smoove_Delly"
VCFlist=${BASE}.list
OUTfile=${BASE}.vcf

# Check if Jasmine succeeded
if [ ! -f "${OUTfile}" ]; then
    echo "ERROR: Jasmine did not create output file ${OUTfile}"
    exit 1
fi


# Reload modules for post-processing
module -q reset
module load HTSlib/1.12-GCC-10.3.0
module load BCFtools/1.12-GCC-10.3.0

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
