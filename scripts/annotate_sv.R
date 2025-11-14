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

# we get an error for missing SVLEN, so we adapt original script
cat("Adding SIMPLE_TYPE and SVLEN to header...\n")
info(header(vcf)) = unique(as(rbind(as.data.frame(info(header(vcf))), data.frame(
    row.names=c("SIMPLE_TYPE", "SVLEN"),
    Number=c("1", "1"),
    Type=c("String", "Integer"),
    Description=c("Simple event type annotation based purely on breakend position and orientation.",
    "Difference in length between REF and ALT alleles"))), "DataFrame"))

cat("Converting to breakpoint ranges...\n")
gr <- breakpointRanges(vcf)

cat("Classifying SV types...\n")
svtype <- simpleEventType(gr)

cat("Updating VCF with annotations...\n")
info(vcf)$SIMPLE_TYPE <- NA_character_
# We add SVLEN initialization in the same way as original script to avoid errors
info(vcf)$SVLEN <- NA_integer_
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
# add for debugging, may be removed later:
if (length(warnings()) > 0) {
  cat("Warnings:\n")
  print(warnings())
}

