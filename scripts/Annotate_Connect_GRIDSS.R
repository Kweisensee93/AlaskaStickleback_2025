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


# Simple SV Type Classifier
simpleEventType <- function(gr) {
  pgr <- partner(gr)
  return(ifelse(seqnames(gr) != seqnames(pgr), "CTX",
    ifelse(strand(gr) == strand(pgr), "INV",
      ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS",
        ifelse(xor(start(gr) < start(pgr), strand(gr) == "-"), "DEL",
          "DUP")))))
}

# Load VCF
cat("Reading VCF...\n")
vcf <- readVcf(input_vcf, genome_name)

# we get an error for missing SVLEN, so we adapt original script
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


# Convert to breakpoint ranges
cat("Converting to breakpoint ranges...\n")
gr <- breakpointRanges(vcf)

cat("Classifying SV types...\n")
svtype <- simpleEventType(gr)
# We leave out original assignment to gr$SVTYPE to avoid overwriting existing SVTYPE
#gr$SVTYPE <- svtype

## Add annotations back into original VCF
cat("Updating VCF with annotations...\n")
info(vcf)$SIMPLE_TYPE <- NA_character_
# We add SVLEN initialization in the same way as original script to avoid errors
info(vcf)$SVLEN <- NA_integer_
info(vcf[gr$sourceId])$SIMPLE_TYPE <- svtype
info(vcf[gr$sourceId])$SVLEN <- gr$svLen

writeVcf(vcf, output_vcf)
cat("Annotated original VCF written to:", output_vcf, "\n")

# Build canonical SV events (DEL/DUP/INV/INS)
### ADAPTED SECTION STARTS HERE ###

cat("Extracting canonical SV events...\n")
canonical_types <- c("DEL", "DUP", "INV", "INS")
is_canonical <- gr$SVTYPE %in% canonical_types

# Keep only one breakend per event
keep_one <- !duplicated(gr$sourceId) & is_canonical
gr_events <- gr[keep_one]

# Add standard fields
gr_events$END <- start(partner(gr_events))
gr_events$SVLEN <- gr_events$svLen

cat("Converting canonical events into VCF records...\n")
vr <- breakpointRangesToSingleBreakend(gr_events)
canonical_vcf <- writeVcfAsBreakpoints(vr)

# Remove BNDs for canonical events
bnd_to_remove <- unique(c(gr_events$sourceId, partner(gr_events)$sourceId))
remaining_vcf <- vcf[!names(vcf) %in% bnd_to_remove]

# Merge canonical SVs + remaining BNDs
cat("Merging canonical SVs with remaining BNDs...\n")

# rbind VCFs safely
combined_vcf <- remaining_vcf
rowRanges(combined_vcf) <- c(rowRanges(remaining_vcf), rowRanges(canonical_vcf))
geno(combined_vcf) <- c(geno(remaining_vcf), geno(canonical_vcf))
fixed(combined_vcf) <- rbind(fixed(remaining_vcf), fixed(canonical_vcf))
info(combined_vcf) <- rbind(info(remaining_vcf), info(canonical_vcf))
meta(combined_vcf) <- meta(remaining_vcf)  # keep original header

final_vcf_path <- sub("\\.vcf.*$", ".combined.vcf", output_vcf)

cat("Writing combined VCF (BND + canonical SV):", final_vcf_path, "\n")
writeVcf(combined_vcf, final_vcf_path)

# BED File Creation (same as original script)
cat("Filtering to PASS variants for BED...\n")
gr_pass <- gr[gr$FILTER == "PASS" & partner(gr)$FILTER == "PASS"]

cat("Creating simple BED file...\n")
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

############################################################

cat("Done!\n")
cat("Summary:\n")
cat("  Total breakpoints:", length(gr), "\n")
cat("  Canonical SV events:", nrow(gr_events), "\n")
cat("  Remaining BNDs:", length(remaining_vcf), "\n")
cat("  Simple SVs (DEL/DUP/INS/INV):", nrow(simplebed), "\n")

if (length(warnings()) > 0) {
  cat("Warnings:\n")
  print(warnings())
}
