
# DESeq2 Analysis for GSE122897 (ruptured vs unruptured aneurysm tissue)
# Run with: Rscript run_deseq2_gse122897.R
# Raw counts: GSE122897_raw_counts_GRCh38.p13_NCBI.tsv.gz
# Gene IDs: NCBI/Entrez (integer)

if (!require("DESeq2", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("DESeq2")
}

library(DESeq2)

cat("DESeq2 analysis for GSE122897 (tissue RNA-seq)\n")
cat("============================================\n")

# Paths (run from project root)
counts_file <- "data/sc_rnaseq/GSE122897/GSE122897_raw_counts_GRCh38.p13_NCBI.tsv.gz"
meta_file <- "data/sc_rnaseq/GSE122897/sample_metadata_with_clinical.tsv"
if (!file.exists(meta_file)) meta_file <- "data/sc_rnaseq/GSE122897/sample_metadata.tsv"
out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE)

# Load counts
cat("\nLoading counts...\n")
counts <- read.delim(counts_file, row.names = 1, check.names = FALSE)
counts <- round(counts)  # ensure integers
cat("  Genes:", nrow(counts), " Samples:", ncol(counts), "\n")

# Load metadata
meta <- read.delim(meta_file, sep = "\t", stringsAsFactors = FALSE)
if (!"Group" %in% colnames(meta)) {
  # Try characteristics_ch1 for rupture status
  if ("characteristics_ch1" %in% colnames(meta)) {
    meta$Group <- ifelse(grepl("Ruptured|ruptured", meta$characteristics_ch1), "Ruptured", "Unruptured")
  } else {
    stop("No Group column and no characteristics_ch1")
  }
}
meta <- meta[meta$Group %in% c("Ruptured", "Unruptured"), ]
meta <- meta[meta$Sample %in% colnames(counts), ]

# Ensure sample order matches
counts <- counts[, meta$Sample, drop = FALSE]
meta$Group <- factor(meta$Group, levels = c("Unruptured", "Ruptured"))

cat("  Ruptured:", sum(meta$Group == "Ruptured"), " Unruptured:", sum(meta$Group == "Unruptured"), "\n")

# Filter low-count genes (keep genes with >= 10 total counts)
keep <- rowSums(counts) >= 10
counts <- counts[keep, ]
cat("  Genes after filtering (sum >= 10):", nrow(counts), "\n")

# DESeq2
cat("\nRunning DESeq2...\n")
dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ Group)
dds <- DESeq(dds)
res <- results(dds, contrast = c("Group", "Ruptured", "Unruptured"))
res <- res[order(res$padj, na.last = TRUE), ]

# Write results
out_file <- file.path(out_dir, "gse122897_deseq2_results.tsv")
res_df <- as.data.frame(res)
res_df$Gene_ID <- rownames(res_df)
res_df <- res_df[, c("Gene_ID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
write.table(res_df, out_file, sep = "\t", row.names = FALSE)
cat("\nResults saved to:", out_file, "\n")

# Summary
cat("\nSummary:\n")
cat("  Significant (padj < 0.05):", sum(res_df$padj < 0.05, na.rm = TRUE), "\n")
cat("  Significant (padj < 0.1):", sum(res_df$padj < 0.1, na.rm = TRUE), "\n")
