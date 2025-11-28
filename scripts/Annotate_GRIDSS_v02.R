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
  pgr <- partner(gr)
  return(ifelse(seqnames(gr) != seqnames(pgr), "CTX",
    ifelse(strand(gr) == strand(pgr), "INV",
      ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS",
        ifelse(xor(start(gr) < start(pgr), strand(gr) == "-"), "DEL",
          "DUP")))))
}

cat("Reading VCF...\n")
vcf <- readVcf(input_vcf, genome_name)

# Add headers for annotation
cat("Adding SIMPLE_TYPE and SVLEN to header...\n")
info(header(vcf)) <- unique(as(rbind(as.data.frame(info(header(vcf))),
    data.frame(
        row.names=c("SIMPLE_TYPE", "SVLEN"),
        Number=c("1", "1"),
        Type=c("String", "Integer"),
        Description=c(
          "Simple event type annotation based purely on breakend position and orientation.",
          "Difference in length between REF and ALT alleles"
        )
    )), "DataFrame"))

cat("Converting to breakpoint ranges...\n")
gr <- breakpointRanges(vcf, inferMissingBreakends=TRUE)

cat("Classifying SV types...\n")
svtype <- simpleEventType(gr)

# Annotate original VCF
info(vcf)$SIMPLE_TYPE <- NA_character_
info(vcf)$SVLEN <- NA_integer_
info(vcf[gr$sourceId])$SIMPLE_TYPE <- svtype
info(vcf[gr$sourceId])$SVLEN <- gr$svLen

cat("Writing annotated VCF (with BNDs)...\n")
writeVcf(vcf, output_vcf)

# Now create unified VCF with standard SV notation
cat("Creating unified VCF (standard SVs + non-simple BNDs)...\n")

# Identify simple events (DEL/DUP/INV/INS)
simple_types <- c("DEL", "DUP", "INV", "INS")
is_simple <- svtype %in% simple_types

# For simple events, keep only one breakend per pair
simple_gr <- gr[is_simple & !duplicated(gr$sourceId)]

# Get all BND IDs that are part of simple events (both breakends)
simple_bnd_ids <- unique(c(simple_gr$sourceId, partner(simple_gr)$sourceId))

cat("  Simple events found:", length(simple_gr), "\n")
cat("  BND records to remove:", length(simple_bnd_ids), "\n")

# Create standard VCF records for simple events
if (length(simple_gr) > 0) {
    simple_vcf_records <- vcf[simple_gr$sourceId]
    
    # Update to standard SV format
    for (i in seq_along(simple_vcf_records)) {
        sv_type <- svtype[is_simple & !duplicated(gr$sourceId)][i]
        
        # Set SVTYPE in INFO
        info(simple_vcf_records[i])$SVTYPE <- sv_type
        
        # Set END position
        info(simple_vcf_records[i])$END <- start(partner(simple_gr[i]))
        
        # Change ALT to <DEL>, <DUP>, etc.
        alt(simple_vcf_records[i]) <- DNAStringSetList(DNAStringSet(paste0("<", sv_type, ">")))
    }
} else {
    simple_vcf_records <- NULL
}

# Keep BNDs that are NOT part of simple events (CTX, unpaired, etc.)
non_simple_vcf <- vcf[!names(vcf) %in% simple_bnd_ids]

cat("  Non-simple BNDs kept:", length(non_simple_vcf), "\n")

# Combine them
if (!is.null(simple_vcf_records) && length(simple_vcf_records) > 0) {
    unified_vcf <- rbind(simple_vcf_records, non_simple_vcf)
} else {
    unified_vcf <- non_simple_vcf
}

# Write unified VCF
unified_vcf_path <- sub("\\.vcf.*$", ".unified.vcf", output_vcf)
cat("Writing unified VCF:", unified_vcf_path, "\n")
writeVcf(unified_vcf, unified_vcf_path)

# Create BED file (PASS only)
cat("Creating BED file from PASS variants...\n")
gr_pass <- gr[gr$FILTER == "PASS" & partner(gr)$FILTER == "PASS"]

simplegr <- gr_pass[simpleEventType(gr_pass) %in% simple_types]
simplebed <- data.frame(
    chrom=seqnames(simplegr),
    start=as.integer((start(simplegr) + end(simplegr)) / 2),
    end=as.integer((start(partner(simplegr)) + end(partner(simplegr))) / 2),
    name=simpleEventType(simplegr),
    score=simplegr$QUAL,
    strand="."
)
simplebed <- simplebed[simplebed$start < simplebed$end,]

write.table(simplebed, output_bed, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)

cat("\n=== DONE ===\n")
cat("Annotated VCF (all BNDs):", output_vcf, "\n")
cat("Unified VCF (standard SVs + non-simple BNDs):", unified_vcf_path, "\n")
cat("BED file (PASS simple SVs):", output_bed, "\n\n")
cat("Summary:\n")
cat("  Total variants:", length(vcf), "\n")
cat("  Simple SVs (DEL/DUP/INV/INS):", length(simple_gr), "\n")
cat("  Non-simple BNDs (CTX/unpaired):", length(non_simple_vcf), "\n")
cat("  Unified VCF total:", length(unified_vcf), "\n")
cat("  BED records (PASS only):", nrow(simplebed), "\n")

if (length(warnings()) > 0) {
  cat("\nWarnings:\n")
  print(warnings())
}
