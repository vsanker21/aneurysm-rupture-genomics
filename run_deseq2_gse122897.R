# GSE122897 DESeq2 — PRIMARY phenotype: explicit GEO source_name, 21 vs 21.
# Ruptured intracranial aneurysm vs Unruptured intracranial aneurysm.
# Excludes cortical artery, GSM3487858 (Unknown intracranial aneurysm),
# and GSM3487892 (Unknown rupture status).
# Does not use grepl("Ruptured|ruptured", ...) substring matching.
# Frozen 22 vs 21 (as-submitted) is not this script.

if (!require("DESeq2", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("DESeq2")
}

library(DESeq2)

cat("GSE122897 DESeq2 — explicit GEO labels, 21 vs 21\n")
cat("================================================\n")

counts_file <- "data/sc_rnaseq/GSE122897/GSE122897_raw_counts_GRCh38.p13_NCBI.tsv.gz"
explicit_meta <- "metadata/gse122897_sample_metadata_explicit_21v21.tsv"
fallback_meta <- "data/sc_rnaseq/GSE122897/sample_metadata.tsv"
out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

counts <- read.delim(counts_file, row.names = 1, check.names = FALSE)
counts <- round(counts)
cat("  Genes:", nrow(counts), " Samples in matrix:", ncol(counts), "\n")

assign_explicit_group <- function(meta) {
  src <- if ("source_name_ch1" %in% colnames(meta)) {
    trimws(as.character(meta$source_name_ch1))
  } else if ("characteristics_ch1" %in% colnames(meta)) {
    trimws(as.character(meta$characteristics_ch1))
  } else {
    stop("Metadata must contain source_name_ch1 or characteristics_ch1")
  }
  group <- rep(NA_character_, length(src))
  group[src == "Ruptured intracranial aneurysm" | src == "tissue/status: Ruptured intracranial aneurysm"] <- "Ruptured"
  group[src == "Unruptured intracranial aneurysm" | src == "tissue/status: Unruptured intracranial aneurysm"] <- "Unruptured"
  meta$Group <- group
  meta
}

if (file.exists(explicit_meta)) {
  meta <- read.delim(explicit_meta, sep = "\t", stringsAsFactors = FALSE)
  cat("  Metadata:", explicit_meta, "\n")
} else if (file.exists(fallback_meta)) {
  cat("  Explicit metadata missing; assigning Group from exact GEO source_name\n")
  meta <- read.delim(fallback_meta, sep = "\t", stringsAsFactors = FALSE)
  meta <- assign_explicit_group(meta)
} else {
  stop("No 21 vs 21 sample metadata found")
}

stopifnot(all(c("Sample", "Group") %in% colnames(meta)))
meta <- meta[meta$Group %in% c("Ruptured", "Unruptured"), ]
meta <- meta[meta$Sample != "GSM3487892", ]
meta <- meta[meta$Sample %in% colnames(counts), ]
if ("GSM3487892" %in% meta$Sample) {
  stop("GSM3487892 must not enter the 21 vs 21 contrast")
}

counts <- counts[, meta$Sample, drop = FALSE]
meta$Group <- factor(meta$Group, levels = c("Unruptured", "Ruptured"))
rownames(meta) <- meta$Sample

n_r <- sum(meta$Group == "Ruptured")
n_u <- sum(meta$Group == "Unruptured")
cat("  Ruptured:", n_r, " Unruptured:", n_u, "\n")
if (n_r != 21L || n_u != 21L) {
  stop(sprintf("Expected 21 vs 21, got %d vs %d", n_r, n_u))
}

keep <- rowSums(counts) >= 10
counts <- counts[keep, ]
cat("  Genes after filtering (sum >= 10):", nrow(counts), "\n")

dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ Group)
dds <- DESeq(dds)
res <- results(dds, contrast = c("Group", "Ruptured", "Unruptured"))
res <- res[order(res$padj, na.last = TRUE), ]

res_df <- as.data.frame(res)
res_df$Gene_ID <- rownames(res_df)
res_df <- res_df[, c("Gene_ID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
out_file <- file.path(out_dir, "gse122897_deseq2_results.tsv")
write.table(res_df, out_file, sep = "\t", row.names = FALSE)

cat("\nResults saved to:", out_file, "\n")
cat("  Samples:", paste(rownames(meta), collapse = ","), "\n")
cat("  Significant (padj < 0.05):", sum(res_df$padj < 0.05, na.rm = TRUE), "\n")
cat("  Significant (padj < 0.05 & |log2FC| >= 0.5):",
    sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 0.5, na.rm = TRUE), "\n")
