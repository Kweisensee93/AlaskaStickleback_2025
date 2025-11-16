#!/bin/bash
#SBATCH --job-name=GRIDSS_call
#SBATCH --output=/storage/homefs/kw23y068/logfiles/GRIDSS_call_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/GRIDSS_call_%j.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --partition=epyc2

# According to GRIDSS documentation, it is optimized for 8 CPU threads and 32GB RAM - lets see

####################
# PARAMETERS TO SET
####################
# May be altered for feeding in arguments to script
# Define which batch of samples to process
BATCH_START=81    # Change to 81, 161, 241, 321, 401 for other batches
BATCH_END=160     # Change to 160, 240, 320, 400, 480 for other batches

# This script is seting up the reference and preprocessing BAM files for GRIDSS SV calling.
# This is a sub-step of 05_GRIDSS.sh

echo "========================================"
echo "GRIDSS Call"
echo "Batch range: ${BATCH_START}-${BATCH_END}"
echo "Date: $(date)"
echo "========================================"


# Load GRIDSS environment
module load Anaconda3/2024.02-1
source $(conda info --base)/etc/profile.d/conda.sh
conda activate gridss


# Paths
PROJECT_DIR=/storage/research/iee_evol/Korbi
REFERENCE=${PROJECT_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna
SAMPLE_LIST=/storage/research/iee_evol/Korbi/output/bams_full.csv

# Set up run directory based on batch; with 0 padding to 3 digits
START_PAD=$(printf "%03d" "${BATCH_START}")
END_PAD=$(printf "%03d" "${BATCH_END}")
RUN_DIR="/storage/scratch/iee_evol/kw23y068/Gridss_joint_${START_PAD}_${END_PAD}/preprocessed"

if [ ! -d "${RUN_DIR}" ]; then
    echo "ERROR: Run directory does not exist: ${RUN_DIR}"
    echo "Either Batch number may be wrong or run preprocessing first."
    exit 1
fi

# Get sample
echo "Collecting BAM files for batch ${BATCH_START}-${BATCH_END}..."
BAM_FILES=()
MISSING_PREPROCESS=0

# Loop through samples in the batch
for i in $(seq ${BATCH_START} ${BATCH_END}); do
    SAMPLE_PATH=$(sed -n "${i}p" "${SAMPLE_LIST}")
    
    if [ -z "${SAMPLE_PATH}" ]; then
        echo "WARNING: Could not read sample at line ${i}"
        continue
    fi
    
    SAMPLE_FILE=$(basename "${SAMPLE_PATH}")
    SAMPLE_NAME="${SAMPLE_FILE%.fixmate.coordsorted.bam}"
    
    BAM="${PROJECT_DIR}/bams_real/${SAMPLE_NAME}.fixmate.coordsorted.bam"
    
    # Verify BAM exists
    if [ ! -f "${BAM}" ]; then
        echo "ERROR: BAM file not found: ${BAM}"
        exit 1
    fi
    
    # Check if preprocessing was completed for this BAM
    WORKING_DIR="${RUN_DIR}/${SAMPLE_NAME}.fixmate.coordsorted.bam.gridss.working"
    if [ ! -d "${WORKING_DIR}" ]; then
        echo "WARNING: Preprocessing not found for ${SAMPLE_NAME}"
        echo "  Expected: ${WORKING_DIR}"
        MISSING_PREPROCESS=$((MISSING_PREPROCESS + 1))
        continue # So we don't fail run, if preprocessing missing for some samples
    fi
    
    BAM_FILES+=("${BAM}")
    echo "  [${i}] ${SAMPLE_NAME}"
done

echo ""
echo "Collected ${#BAM_FILES[@]} BAM files"
if [ ${#BAM_FILES[@]} -eq 0 ]; then
    echo "ERROR: No BAM files collected"
    exit 1
fi
# Short overview:
if [ ${MISSING_PREPROCESS} -gt 0 ]; then
    echo "WARNING: ${MISSING_PREPROCESS} samples missing preprocessing"
fi

# Make assembly output directory
ASSEMBLY_DIR="${RUN_DIR}/assembly"
mkdir -p "${ASSEMBLY_DIR}"

# Define output VCF
VCF_OUT="${ASSEMBLY_DIR}/joint.vcf.gz"

# STEP 4: Call variants jointly
gridss -s call \
  -r "${REFERENCE}" \
  -w "${RUN_DIR}" \
  -o "${VCF_OUT}" \
  -a "${ASSEMBLY_DIR}/batch_assembly.bam" \
  "${BAM_FILES[@]}"

CALL_EXIT_CODE=$?

echo ""
echo "========================================"
echo "Finished GRIDSS assembly"
echo "Exit code: ${CALL_EXIT_CODE}"
echo "Date: $(date)"
echo "========================================"

conda deactivate

# Get job statistics
echo ""
echo "========================================"
echo "Resource Usage Summary:"
echo "========================================"
sacct -j ${SLURM_JOB_ID} --format=JobID,MaxRSS,MaxVMSize,Elapsed,CPUTime,TotalCPU -P
