######
1st
########
#!/bin/bash
#SBATCH --job-name=GRIDSS_preprocess
#SBATCH --output=/storage/homefs/kw23y068/logfiles/GRIDSS_preprocess_%x_%A_%a.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/GRIDSS_preprocess_%x_%A_%a.err
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=epyc2
#SBATCH --array=1-4   # Adjust to number of samples

module load Anaconda3/2024.02-1
source $(conda info --base)/etc/profile.d/conda.sh
conda activate gridss

# Paths
PROJECT_DIR=/storage/research/iee_evol/Korbi
REFERENCE=${PROJECT_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna
SAMPLE_LIST=/storage/research/iee_evol/Korbi/output/bams_full.csv
RUN_DIR=/storage/scratch/iee_evol/kw23y068/Gridss_joint/preprocessed

mkdir -p "${RUN_DIR}"

# Get sample
SAMPLE_PATH=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_LIST})
SAMPLE_FILE=$(basename "${SAMPLE_PATH}")
SAMPLE_NAME="${SAMPLE_FILE%.fixmate.coordsorted.bam}"
BAM=${PROJECT_DIR}/bams_real/${SAMPLE_NAME}.fixmate.coordsorted.bam

echo "Processing sample ${SAMPLE_NAME}"

# Step 1: Setup reference (only once globally)
if [ ${SLURM_ARRAY_TASK_ID} -eq 1 ]; then
  gridss -s setupreference -r ${REFERENCE}
fi

# Step 2: Preprocess
gridss -s preprocess -r ${REFERENCE} "${BAM}"

conda deactivate

########
#2nd
# Assemble + Call
########

#!/bin/bash
#SBATCH --job-name=GRIDSS_joint
#SBATCH --output=/storage/homefs/kw23y068/logfiles/GRIDSS_joint_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/GRIDSS_joint_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --partition=epyc2

module load Anaconda3/2024.02-1
source $(conda info --base)/etc/profile.d/conda.sh
conda activate gridss

N_SAMPLES=4
PROJECT_DIR=/storage/research/iee_evol/Korbi
RUN_DIR=/storage/scratch/iee_evol/kw23y068/Gridss_joint/joint_${N_SAMPLES}samples
REFERENCE=${PROJECT_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna
VCF_OUT=${RUN_DIR}/joint_${N_SAMPLES}samples.vcf.gz

mkdir -p "${RUN_DIR}"

mapfile -t BAM_FILES < <(head -n ${N_SAMPLES} /storage/research/iee_evol/Korbi/output/bams_full.csv)

echo "Assembling and calling variants for ${N_SAMPLES} samples"

# Step 3: Assemble
gridss -s assemble -a ${RUN_DIR}/joint_assembly.bam -r ${REFERENCE} "${BAM_FILES[@]}"

# Step 4: Call
gridss -s call -r ${REFERENCE} -a ${RUN_DIR}/joint_assembly.bam -o ${VCF_OUT} "${BAM_FILES[@]}"

conda deactivate
