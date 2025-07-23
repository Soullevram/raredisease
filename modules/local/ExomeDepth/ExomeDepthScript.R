#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ExomeDepth)
  library(optparse)
  library(GenomicRanges)
  library(rtracklayer)
})

# -----------------------------
# CLI argument parser
# -----------------------------
option_list <- list(
  make_option("--bam", type="character", help="Input BAM file"),
  make_option("--bed", type="character", help="BED file with target regions"),
  make_option("--out", type="character", help="Output filename for CNV calls")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# -----------------------------
# Load BED file as target regions
# -----------------------------
cat("Loading BED file: ", opt$bed, "\n")
target_regions <- import(opt$bed, format = "BED")

# -----------------------------
# Extract read count from BAM
# -----------------------------
cat("Reading BAM: ", opt$bam, "\n")
bam_counts <- getBamCounts(bed.frame = as.data.frame(target_regions),
                           bam.files = opt$bam,
                           include.chr = TRUE)

count_data <- as(bam_counts, 'data.frame')

# -----------------------------
# ExomeDepth CNV calling (mock ref)
# -----------------------------
# NOTE: This uses a self-reference for now because I'm still testing to see if it'll work.
# For this benchmarking, I'll replace with actual references.

exome <- new('ExomeDepth',
             test = count_data$reads,
             reference = count_data$reads,
             formula = 'cbind(test, reference) ~ 1')

cnvs <- CallCNVs(x = exome,
                 transition.probability = 10^-4,
                 chromosome = count_data$space,
                 start = count_data$start,
                 end = count_data$end,
                 name = count_data$names)

# -----------------------------
# Write CNVs to output
# -----------------------------
write.table(cnvs, file = opt$out, sep = "\t", row.names = FALSE, quote = FALSE)
cat("CNV calling complete. Results saved to", opt$out, "\n")
