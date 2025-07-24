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
  make_option("--bam", type = "character", help = "Test BAM file"),
  make_option("--controls", type = "character", help = "Comma-separated list of control BAM files"),
  make_option("--bed", type = "character", help = "BED file with target regions"),
  make_option("--fasta", type = "character", help = "Reference FASTA file", default = NULL),
  make_option("--out", type = "character", help = "Output file for CNV calls"),
  make_option("--readlen", type = "integer", default = 150, help = "Read length (default: 150)"),
  make_option("--transition_prob", type = "double", default = 1e-4),
  make_option("--expected_cnv_len", type = "double", default = 1e6),
  make_option("--n_bins_reduced", type = "integer", default = 10000),
  make_option("--phi_bins", type = "double", default = 1)
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# -----------------------------
# Load BED
# -----------------------------
bed_df <- read.table(opt$bed, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(bed_df) <- c("chr", "start", "end", "gene")[1:ncol(bed_df)]
if (!"gene" %in% colnames(bed_df)) {
  bed_df$gene <- paste0("GENE_", seq_len(nrow(bed_df)))
}

# -----------------------------
# Count reads from BAMs
# -----------------------------
test_bam <- opt$bam
control_bams <- strsplit(opt$controls, ",")[[1]]

bam_files <- c(test_bam, control_bams)

counts <- getBamCounts(
  bed.frame = bed_df,
  bam.files = bam_files,
  include.chr = TRUE,
  read.width = opt$readlen,
  referenceFasta = opt$fasta
)

counts_df <- as.data.frame(counts)
colnames(counts_df) <- gsub(".bam", "", colnames(counts_df))

# -----------------------------
# Split counts into test + controls
# -----------------------------
test_sample <- gsub(".bam", "", basename(test_bam))
control_samples <- gsub(".bam", "", basename(control_bams))

test_counts <- as.matrix(counts_df[, test_sample, drop = FALSE])
control_counts <- as.matrix(counts_df[, control_samples, drop = FALSE])

# -----------------------------
# Select reference set
# -----------------------------
references <- select.reference.set(
  test.counts = test_counts,
  reference.count = control_counts,
  bin.length = (counts_df$end - counts_df$start) / 1000,
  n.bins.reduced = opt$n_bins_reduced,
  phi.bins = opt$phi_bins
)

selected_controls <- apply(
  X = as.matrix(control_counts[, references$reference.choice]),
  MARGIN = 1,
  FUN = sum
)

# -----------------------------
# Run ExomeDepth
# -----------------------------
exome <- new('ExomeDepth',
             test = test_counts,
             reference = selected_controls,
             formula = 'cbind(test, reference) ~ 1')

cnvs <- CallCNVs(
  x = exome,
  transition.probability = opt$transition_prob,
  expected.CNV.length = opt$expected_cnv_len,
  chromosome = counts_df$space,
  start = counts_df$start,
  end = counts_df$end,
  name = counts_df$names
)

# -----------------------------
# Annotate gene names
# -----------------------------
getGeneForPosition <- function(bedData, chr, pos) {
  gene <- bedData[bedData$chr == chr & pos >= bedData$start & pos <= bedData$end, "gene"]
  if (length(gene) == 0) return(NA)
  return(gene[1])
}

if (nrow(cnvs) > 0) {
  cnvs$gene <- mapply(getGeneForPosition, chr = cnvs$chromosome, pos = cnvs$end, MoreArgs = list(bedData = bed_df))
  cnvs$sample <- test_sample
}

# -----------------------------
# Output
# -----------------------------
write.table(cnvs, file = opt$out, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Finished CNV calling for", test_sample, "\n")
