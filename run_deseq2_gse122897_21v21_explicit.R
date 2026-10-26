# GSE122897 DESeq2 — primary phenotype from explicit GEO source_name labels.
# 21 ruptured vs 21 unruptured. GSM3487892 (Unknown rupture status) excluded.
# Does not overwrite the frozen 22 vs 21 as-submitted results.

if (!require("DESeq2", quietly = TRUE)) {
  stop("DESeq2 is required")
}

library(DESeq2)

base <- "."
counts_file <- file.path(base, "data/sc_rnaseq/GSE122897/GSE122897_raw_counts_GRCh38.p13_NCBI.tsv.gz")
meta_file <- file.path(base, "metadata/gse122897_sample_metadata_explicit_21v21.tsv")
if (!file.exists(meta_file)) {
  meta_file <- file.path(base, "results/deseq2_gse122897_21v21_explicit/sample_metadata_explicit_21v21.tsv")
}
out_dir <- file.path(base, "results/deseq2_gse122897_21v21_explicit")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("GSE122897 DESeq2 — explicit GEO labels, 21 vs 21\n")
cat("================================================\n")

counts <- read.delim(counts_file, row.names = 1, check.names = FALSE)
counts <- round(counts)
cat("  Genes:", nrow(counts), " Samples in matrix:", ncol(counts), "\n")

meta <- read.delim(meta_file, sep = "\t", stringsAsFactors = FALSE)
stopifnot(all(c("Sample", "Group") %in% colnames(meta)))
meta <- meta[meta$Group %in% c("Ruptured", "Unruptured"), ]
missing <- setdiff(meta$Sample, colnames(counts))
if (length(missing)) {
  stop(paste("Samples missing from counts:", paste(missing, collapse = ", ")))
}
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
out_file <- file.path(out_dir, "gse122897_deseq2_results_21v21.tsv")
write.table(res_df, out_file, sep = "\t", row.names = FALSE)

cat("\nResults saved to:", out_file, "\n")
cat("  Samples:", paste(rownames(meta), collapse = ","), "\n")
cat("  Significant (padj < 0.05):", sum(res_df$padj < 0.05, na.rm = TRUE), "\n")
cat("  Significant (padj < 0.05 & |log2FC| >= 0.5):",
    sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 0.5, na.rm = TRUE), "\n")
