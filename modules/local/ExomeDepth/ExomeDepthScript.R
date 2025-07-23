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
  make_option("--bam", type = "character", help = "Input BAM file"),
  make_option("--bed", type = "character", help = "BED file with target regions"),
  make_option("--fasta", type = "character", help = "Reference FASTA file", default = NULL),
  make_option("--out", type = "character", help = "Output file for CNV calls"),
  make_option("--readlen", type = "integer", default = 150, help = "Read length (default: 150)"),
  make_option("--transition_prob", type = "double", default = 1e-4, help = "Transition probability for CNV calling"),
  make_option("--expected_cnv_len", type = "double", default = 1e6, help = "Expected CNV length in bp")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# -----------------------------
# Load BED and add dummy gene column if missing
# -----------------------------
cat("Loading BED file:", opt$bed, "\n")
bed_df <- read.table(opt$bed, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(bed_df) <- c("chr", "start", "end", "gene")[1:ncol(bed_df)]
if (!"gene" %in% colnames(bed_df)) {
  bed_df$gene <- paste0("GENE_", seq_len(nrow(bed_df)))
}

# -----------------------------
# Count reads from BAM
# -----------------------------
cat("Reading BAM:", opt$bam, "\n")
bam_counts <- getBamCounts(
  bed.frame = bed_df,
  bam.files = opt$bam,
  include.chr = TRUE,
  read.width = opt$readlen,
  referenceFasta = opt$fasta
)

count_df <- as.data.frame(bam_counts)
colnames(count_df)[which(names(count_df) == "reads")] <- "test"

# Use the test sample as its own reference (mock)
count_df$reference <- count_df$test

# -----------------------------
# ExomeDepth CNV calling
# -----------------------------
cat("Running ExomeDepth\n")
exome <- new('ExomeDepth',
             test = count_df$test,
             reference = count_df$reference,
             formula = 'cbind(test, reference) ~ 1')

cnvs <- CallCNVs(
  x = exome,
  transition.probability = opt$transition_prob,
  expected.CNV.length = opt$expected_cnv_len,
  chromosome = count_df$space,
  start = count_df$start,
  end = count_df$end,
  name = count_df$names
)

# -----------------------------
# Add gene annotation
# -----------------------------
getGeneForPosition <- function(bedData, chr, pos) {
  gene <- bedData[bedData$chr == chr & pos >= bedData$start & pos <= bedData$end, "gene"]
  if (length(gene) == 0) return(NA)
  return(gene[1])
}

if (nrow(cnvs) > 0) {
  cnvs$gene <- mapply(getGeneForPosition, chr = cnvs$chromosome, pos = cnvs$end, MoreArgs = list(bedData = bed_df))
}

# -----------------------------
# Save results
# -----------------------------
cat("Saving CNV results to", opt$out, "\n")
write.table(cnvs, file = opt$out, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Finished.\n")
