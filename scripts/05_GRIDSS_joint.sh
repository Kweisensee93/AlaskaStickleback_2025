#!/bin/bash
#SBATCH --job-name=GRIDSS_joint
#SBATCH --output=/storage/homefs/kw23y068/logfiles/GRIDSS_joint_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/GRIDSS_joint_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --accounting-stats
#SBATCH --partition=epyc2

# According to GRIDSS documentation, it is optimized for 8 CPU threads and 32GB RAM

####################
# PARAMETERS TO SET
####################
N_SAMPLES=4  # Change this to 4, 8, 16, or 80 for different test runs

echo "========================================"
echo "Starting GRIDSS Joint-Calling"
echo "Number of samples: ${N_SAMPLES}"
echo "Date: $(date)"
echo "========================================"

# Clean up any leftover gridss lock files from previous runs
rm -f /storage/research/iee_evol/Korbi/ref/*.tmp.gridsslock


# Load Anaconda module to access conda environments
module load Anaconda3/2024.02-1

# Define paths
SAMPLE_LIST=/storage/research/iee_evol/Korbi/output/bams_full.csv
GRIDSS_CONDA=/storage/homefs/kw23y068/.conda/envs/gridss/bin/gridss
PROJECT_DIR=/storage/research/iee_evol/Korbi
RUN_DIR=/storage/scratch/iee_evol/kw23y068/Gridss_joint/joint_${N_SAMPLES}samples

# Create run directory if it doesn't exist
if [ ! -d "${RUN_DIR}" ]; then
    mkdir -p "${RUN_DIR}"
fi

# Input files
REFERENCE=${PROJECT_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna

# Output files
VCF_OUT=${RUN_DIR}/joint_${N_SAMPLES}samples.vcf.gz

# Verify reference exists
[ -f "${REFERENCE}" ] || { echo "ERROR: Reference not found: ${REFERENCE}"; exit 1; }

# Read the first N_SAMPLES from the sample list
echo "Reading first ${N_SAMPLES} samples from ${SAMPLE_LIST}..."
mapfile -t SAMPLE_PATHS < <(head -n ${N_SAMPLES} "${SAMPLE_LIST}")

# Check that we got the expected number of samples
if [ ${#SAMPLE_PATHS[@]} -ne ${N_SAMPLES} ]; then
    echo "ERROR: Expected ${N_SAMPLES} samples but got ${#SAMPLE_PATHS[@]}"
    exit 1
fi

# Build arrays for BAM files and assembly BAM outputs
BAM_FILES=()
ASSEMBLY_BAMS=()
SAMPLE_NAMES=()

echo ""
echo "Processing sample list:"
for SAMPLE_PATH in "${SAMPLE_PATHS[@]}"; do
    SAMPLE_FILE=$(basename "${SAMPLE_PATH}")
    SAMPLE_NAME="${SAMPLE_FILE%.fixmate.coordsorted.bam}"
    SAMPLE_NAMES+=("${SAMPLE_NAME}")
    
    BAM_FILE=${PROJECT_DIR}/bams_real/${SAMPLE_NAME}.fixmate.coordsorted.bam
    
    # Verify BAM exists
    if [ ! -f "${BAM_FILE}" ]; then
        echo "ERROR: BAM file not found: ${BAM_FILE}"
        exit 1
    fi
    
    BAM_FILES+=("${BAM_FILE}")
    ASSEMBLY_BAMS+=("${RUN_DIR}/${SAMPLE_NAME}.assembly.bam")
    
    echo "  - ${SAMPLE_NAME}"
done

echo ""
echo "All ${N_SAMPLES} BAM files verified."
echo ""

# Activate GRIDSS environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate gridss


# Execute GRIDSS
echo ""
echo "Executing GRIDSS..."


# STEP 1: Reference setup (only once)
if [ ! -d "${REFERENCE}.gridss.working" ]; then
  echo "Running GRIDSS reference setup..."
  gridss -s setupreference -r "${REFERENCE}" -w "${RUN_DIR}"
else
  echo "Reference already set up. Skipping setupreference step."
fi


# STEP 2: Preprocess each sample
for BAM in "${BAM_FILES[@]}"; do
  gridss -s preprocess -r "${REFERENCE}" -w "${RUN_DIR}" "${BAM}"
done

# STEP 3: Assemble -Maybe outsource this to another script???
gridss -s assemble -r "${REFERENCE}" -w "${RUN_DIR}" -a "${RUN_DIR}/joint_assembly.bam" "${BAM_FILES[@]}"

# STEP 4: Call variants jointly
gridss -s call \
  -r "${REFERENCE}" \
  -w "${RUN_DIR}" \
  -o "${VCF_OUT}" \
  -a "${RUN_DIR}/joint_assembly.bam" \
  "${BAM_FILES[@]}"

GRIDSS_EXIT_CODE=$?

echo ""
echo "========================================"
echo "Finished GRIDSS joint-calling"
echo "Exit code: ${GRIDSS_EXIT_CODE}"
echo "Date: $(date)"
echo "========================================"

# Deactivate environment
conda deactivate

echo ""
echo "Output files:"
echo "  Joint VCF: ${VCF_OUT}"
echo ""
echo "Assembly BAMs:"
for i in "${!SAMPLE_NAMES[@]}"; do
    echo "  - ${SAMPLE_NAMES[$i]}: ${ASSEMBLY_BAMS[$i]}"
done
echo ""

if [ ${GRIDSS_EXIT_CODE} -eq 0 ]; then
    echo "SUCCESS: Joint-calling completed successfully for ${N_SAMPLES} samples"
else
    echo "ERROR: GRIDSS failed with exit code ${GRIDSS_EXIT_CODE}"
    exit ${GRIDSS_EXIT_CODE}
fi

# Clean up any leftover gridss lock files
rm -f /storage/research/iee_evol/Korbi/ref/*.tmp.gridsslock

# Get job statistics
echo ""
echo "========================================"
echo "Resource Usage Summary:"
echo "========================================"
sacct -j ${SLURM_JOB_ID} --format=JobID,MaxRSS,MaxVMSize,Elapsed,CPUTime,TotalCPU -P
