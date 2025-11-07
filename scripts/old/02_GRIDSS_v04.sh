#!/bin/bash
#SBATCH --job-name=GRIDSS
#SBATCH --output=/storage/homefs/kw23y068/logfiles/GRIDSS_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/GRIDSS_%j.err
#SBATCH --time=02:30:00
#SBATCH --cpus-per-task=9
#SBATCH --mem=32G
#SBATCH --partition=epyc2

echo "starting GRIDSS SV calling"
echo "Date: $(date)"

# Load Anaconda module to access conda environments
module load Anaconda3/2024.02-1

SAMPLE_NAME="FG_CC_19T_031"

# Define paths
GRIDSS_CONDA=/storage/homefs/kw23y068/.conda/envs/gridss/bin/gridss
PROJECT_DIR=/storage/research/iee_evol/Korbi
RUN_DIR=/storage/scratch/iee_evol/kw23y068/Gridss/${SAMPLE_NAME}

# Create run directory if it doesn't exist
if [ ! -d "${RUN_DIR}" ]; then
    mkdir -p "${RUN_DIR}"
fi

# Input files
REFERENCE=${PROJECT_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna
BAM=${PROJECT_DIR}/bams_real/${SAMPLE_NAME}.fixmate.coordsorted.bam

# Output files
VCF_OUT=${RUN_DIR}/${SAMPLE_NAME}.vcf.gz
ASSEMBLY_BAM=${RUN_DIR}/${SAMPLE_NAME}.assembly.bam
VCF_ANNOTATED=${RUN_DIR}/${SAMPLE_NAME}.annotated.vcf
SIMPLE_BED=${RUN_DIR}/${SAMPLE_NAME}.simple.bed

# Verify inputs exist
[ -f "${REFERENCE}" ] || { echo "ERROR: Reference not found: ${REFERENCE}"; exit 1; }
[ -f "${BAM}" ] || { echo "ERROR: BAM file not found: ${BAM}"; exit 1; }

# Activate GRIDSS environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate gridss

# Run GRIDSS
echo "Running GRIDSS..."
${GRIDSS_CONDA} \
    -t $((SLURM_CPUS_PER_TASK - 1)) \
    -r ${REFERENCE} \
    -o ${VCF_OUT} \
    -a ${ASSEMBLY_BAM} \
    -w ${RUN_DIR} \
    ${BAM}

echo "Finished GRIDSS SV calling"
echo "Date: $(date)"

# Deactivate GRIDSS environment
conda deactivate
# Activate R environment for annotation
conda activate R_gridss_annotation

# Annotate with StructuralVariantAnnotation
echo "Starting SV annotation..."
echo "Date: $(date)"

# Create R script for annotation based on GRIDSS example on GitHub
# https://github.com/PapenfussLab/gridss/blob/master/example/simple-event-annotation.R
cat > ${RUN_DIR}/annotate_sv.R << 'EOF'
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
    stop("Usage: annotate_sv.R input.vcf.gz genome_name output.vcf output.bed")
}

input_vcf <- args[1]
genome_name <- args[2]
output_vcf <- args[3]
output_bed <- args[4]

# Load libraries
suppressPackageStartupMessages({
    library(VariantAnnotation)
    library(StructuralVariantAnnotation)
    library(stringr)
})

#' Simple SV type classifier
simpleEventType <- function(gr) {
  pgr = partner(gr)
  return(ifelse(seqnames(gr) != seqnames(pgr), "CTX", # inter-chromosomal
    ifelse(strand(gr) == strand(pgr), "INV",
      ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS", # TODO: improve classification of complex events
        ifelse(xor(start(gr) < start(pgr), strand(gr) == "-"), "DEL",
          "DUP")))))
}

cat("Reading VCF...\n")
vcf <- readVcf(input_vcf, genome_name)

cat("Adding SIMPLE_TYPE andto header...\n")
info(header(vcf)) = unique(as(rbind(as.data.frame(info(header(vcf))), data.frame(
    row.names=c("SIMPLE_TYPE"),
    Number=c("1"),
    Type=c("String"),
    Description=c("Simple event type annotation based purely on breakend position and orientation."))), "DataFrame"))

cat("Converting to breakpoint ranges...\n")
gr <- breakpointRanges(vcf)

cat("Classifying SV types...\n")
svtype <- simpleEventType(gr)

cat("Updating VCF with annotations...\n")
info(vcf)$SIMPLE_TYPE <- NA_character_
info(vcf[gr$sourceId])$SIMPLE_TYPE <- svtype
info(vcf[gr$sourceId])$SVLEN <- gr$svLen

cat("Writing annotated VCF...\n")
writeVcf(vcf, output_vcf)

cat("Filtering to PASS variants...\n")
# By default, GRIDSS is very sensitive but this comes at the cost of a high false discovery rate
gr <- gr[gr$FILTER == "PASS" & partner(gr)$FILTER == "PASS"] # Remove low confidence calls

cat("Creating simple BED file...\n")
simplegr <- gr[simpleEventType(gr) %in% c("INS", "INV", "DEL", "DUP")]
simplebed <- data.frame(
    chrom=seqnames(simplegr),
    # call the centre of the homology/inexact interval
    start=as.integer((start(simplegr) + end(simplegr)) / 2),
    end=as.integer((start(partner(simplegr)) + end(partner(simplegr))) / 2),
    name=simpleEventType(simplegr),
    score=simplegr$QUAL,
    strand="."
)

# Just the lower of the two breakends so we don't output everything twice
simplebed <- simplebed[simplebed$start < simplebed$end,]

write.table(simplebed, output_bed, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)

cat("Done!\n")
cat("Summary:\n")
cat("  Total breakpoints:", length(gr), "\n")
cat("  PASS filter only:", sum(gr$FILTER == "PASS" & partner(gr)$FILTER == "PASS"), "\n")
cat("  Simple SVs (DEL/DUP/INS/INV):", nrow(simplebed), "\n")
EOF

chmod +x ${RUN_DIR}/annotate_sv.R

# Run annotation
# Genome name is for VCF header
Rscript ${RUN_DIR}/annotate_sv.R \
    ${VCF_OUT} \
    "GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna" \
    ${VCF_ANNOTATED} \
    ${SIMPLE_BED}

echo "Finished SV annotation"
echo "Date: $(date)"
echo ""
echo "Output files:"
echo "  Raw GRIDSS VCF:     ${VCF_OUT}"
echo "  Annotated VCF:      ${VCF_ANNOTATED}"
echo "  Simple BED (PASS):  ${SIMPLE_BED}"
echo ""
echo "The BED file contains only PASS-filtered simple SVs (DEL/DUP/INS/INV)"