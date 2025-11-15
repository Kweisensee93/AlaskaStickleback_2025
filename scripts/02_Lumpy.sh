#!/bin/bash
#SBATCH --job-name=lumpy_run
#SBATCH --output=/storage/homefs/kw23y068/logfiles/lumpy_%j.log
#SBATCH --error=/storage/homefs/kw23y068/logfiles/lumpy_%j.err
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH --partition=epyc2

module load Anaconda3/2024.02-1
source $(conda info --base)/etc/profile.d/conda.sh
conda activate lumpy

#paths
PROJECT_DIR=/storage/research/iee_evol/Korbi
REFERENCE=${PROJECT_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna
BAM="${PROJECT_DIR}/bams_real/FG_CC_19T_031.fixmate.coordsorted.bam"
OUT_DIR="/storage/scratch/iee_evol/kw23y068/lumpy/"
#SAMPLE_LIST=/storage/research/iee_evol/Korbi/output/bams_full.csv
SAMPLE_NAME="FG_CC_19T_031"
RUN_DIR="${OUT_DIR}/${SAMPLE_NAME}"

if [ ! -d "$OUT_DIR" ]; then
    mkdir -p "$OUT_DIR"
fi
if [ ! -d "$RUN_DIR" ]; then
    mkdir -p "$RUN_DIR"
fi

# Lets follow LUMPY example workflow
# https://github.com/arq5x/lumpy-sv
# We go with the BWA-MEM and samtools route; we already hace a .bam file

# Extract discordant & split-read alignments
DISCORDANT="${RUN_DIR}/${SAMPLE_NAME}.discordants.unsorted.bam"
SPLITREAD="${RUN_DIR}/${SAMPLE_NAME}.splitters.unsorted.bam"

# Extract the discordant paired-end alignments.
samtools view -b -F 1294 $BAM > ${DISCORDANT}
# Extract the split-read alignments
samtools view -h $BAM \
    | extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > ${SPLITREAD}


# Sort both alignments
samtools sort $DISCORDANT ${SAMPLE_NAME}.discordants
samtools sort $SPLITREAD ${SAMPLE_NAME}.splitters

echo "Running LUMPY…"

lumpyexpress \
    -B ${BAM} \
    -S ${SPLITREAD} \
    -D ${DISCORDANT} \
    -o ${RUN_DIR}/${SAMPLE_NAME}_lumpy.vcf

echo "Done. Output written to: ${RUN_DIR}/${SAMPLE_NAME}_lumpy.vcf"

# Compress VCF
# the conda environment has bgzip installed
bgzip -c ${RUN_DIR}/${SAMPLE_NAME}_lumpy.vcf > ${RUN_DIR}/${SAMPLE_NAME}_lumpy.vcf.gz

conda deactivate
# No BFCtools in lumpy env, so load it separately, afte conda deactivate

module load BCFtools/1.12-GCC-10.3.0

bcftools index ${RUN_DIR}/${SAMPLE_NAME}_lumpy.vcf.gz
