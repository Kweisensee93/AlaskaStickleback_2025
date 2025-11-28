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
gr <- breakpointRanges(vcf)

cat("Classifying SV types...\n")
svtype <- simpleEventType(gr)

cat("Updating VCF with annotations...\n")
info(vcf)$SIMPLE_TYPE <- NA_character_
info(vcf)$SVLEN <- NA_integer_
info(vcf[gr$sourceId])$SIMPLE_TYPE <- svtype
info(vcf[gr$sourceId])$SVLEN <- gr$svLen

cat("Writing annotated original VCF...\n")
writeVcf(vcf, output_vcf)

cat("Building unified VCF (standard SV notation + unpaired BNDs)...\n")

# Identify canonical SV types that can be converted
canonical_types <- c("DEL", "DUP", "INV", "INS")
is_canonical <- svtype %in% canonical_types

# Get unique event IDs for canonical SVs (keep one breakend per pair)
canonical_gr <- gr[is_canonical & !duplicated(gr$sourceId)]

cat("  Canonical SVs to convert:", length(canonical_gr), "\n")

# Build standard VCF records for canonical SVs
if (length(canonical_gr) > 0) {
    # Prepare END and SVLEN fields
    canonical_gr$END <- start(partner(canonical_gr))
    canonical_gr$SVLEN <- canonical_gr$svLen
    
    # Convert to single breakend format then to standard VCF
    vr <- breakpointRangesToSingleBreakend(canonical_gr)
    canonical_vcf <- writeVcfAsBreakpoints(vr)
    
    # Get IDs of BND records to remove (both ends of canonical events)
    bnd_ids_to_remove <- unique(c(canonical_gr$sourceId, partner(canonical_gr)$sourceId))
} else {
    canonical_vcf <- NULL
    bnd_ids_to_remove <- character(0)
}

# Keep BNDs that are NOT part of canonical events
remaining_bnd_vcf <- vcf[!names(vcf) %in% bnd_ids_to_remove]

cat("  Remaining BNDs (unpaired/CTX):", length(remaining_bnd_vcf), "\n")

# Merge canonical SVs with remaining BNDs
if (!is.null(canonical_vcf) && length(canonical_vcf) > 0) {
    cat("Merging canonical SVs with remaining BNDs...\n")
    
    combined_vcf <- remaining_bnd_vcf
    rowRanges(combined_vcf) <- c(rowRanges(remaining_bnd_vcf), rowRanges(canonical_vcf))
    geno(combined_vcf) <- c(geno(remaining_bnd_vcf), geno(canonical_vcf))
    fixed(combined_vcf) <- rbind(fixed(remaining_bnd_vcf), fixed(canonical_vcf))
    info(combined_vcf) <- rbind(info(remaining_bnd_vcf), info(canonical_vcf))
    meta(combined_vcf) <- meta(remaining_bnd_vcf)
} else {
    cat("No canonical SVs found, keeping all BNDs...\n")
    combined_vcf <- remaining_bnd_vcf
}

# Write unified VCF
unified_vcf_path <- sub("\\.vcf.*$", ".unified.vcf", output_vcf)
cat("Writing unified VCF:", unified_vcf_path, "\n")
writeVcf(combined_vcf, unified_vcf_path)

# Create BED file from PASS variants
cat("Creating BED file from PASS variants...\n")
gr_pass <- gr[gr$FILTER == "PASS" & partner(gr)$FILTER == "PASS"]

simplegr <- gr_pass[simpleEventType(gr_pass) %in% c("INS", "INV", "DEL", "DUP")]
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

cat("\n=== SUMMARY ===\n")
cat("Original VCF written to:", output_vcf, "\n")
cat("Unified VCF written to:", unified_vcf_path, "\n")
cat("BED file written to:", output_bed, "\n\n")
cat("Statistics:\n")
cat("  Total breakpoints (all):", length(gr), "\n")
cat("  Canonical SVs converted:", length(canonical_gr), "\n")
cat("  Remaining BNDs:", length(remaining_bnd_vcf), "\n")
cat("  Total records in unified VCF:", length(combined_vcf), "\n")
cat("  PASS simple SVs in BED:", nrow(simplebed), "\n")

if (length(warnings()) > 0) {
  cat("\nWarnings:\n")
  print(warnings())
}