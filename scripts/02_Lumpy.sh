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
#SAMPLE_LIST=/storage/research/iee_evol/Korbi/output/bams_full.csv
BAM="/storage/research/iee_evol/Korbi/bams_real/FG_CC_19T_031.fixmate.coordsorted.bam"
OUT_DIR="/storage/scratch/iee_evol/kw23y068/lumpy/"
# sample_lumpy.vcf
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
DISCORDANT="${OUT_DIR}/${SAMPLE_NAME}.discordants.unsorted.bam"
SPLITREAD="${OUT_DIR}/${SAMPLE_NAME}.splitters.unsorted.bam"

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

# lumpy \
#     -mw 4 \
#     -tt 0 \
#     -pe id:sample,bam_file:$DISCORDANT,histo_file:sample.histo,mean:300,stdev:50,read_length:150,min_non_overlap:100,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
#     -sr id:sample,bam_file:$SPLITREAD,back_distance:10,weight:1,min_mapping_threshold:20 \
#     > $OUT

lumpyexpress \
    -B ${BAM} \
    -S ${SPLITREAD}.bam \
    -D ${DISCORDANT}.bam \
    -o ${RUN_DIR}/${SAMPLE_NAME}_lumpy.vcf

echo "Done. Output written to: ${RUN_DIR}/${SAMPLE_NAME}_lumpy.vcf"
