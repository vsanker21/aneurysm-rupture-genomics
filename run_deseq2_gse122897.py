"""
DESeq2 analysis for GSE122897 tissue RNA-seq (ruptured vs unruptured).
Uses raw counts for proper count-based differential expression.
Creates and runs an R script; requires R and DESeq2 (Bioconductor).
"""
import os
import subprocess
import glob

def find_rscript():
    """Find Rscript executable, checking PATH and common Windows install locations."""
    import shutil
    rscript = shutil.which("Rscript")
    if rscript:
        return rscript
    # Common R install paths on Windows
    for base in [r"C:\Program Files\R", os.path.expandvars(r"%LOCALAPPDATA%\Programs\R")]:
        if not base or not os.path.exists(base):
            continue
        for p in glob.glob(os.path.join(base, "R-*", "bin", "*", "Rscript.exe")):
            if os.path.isfile(p):
                return p
        for p in glob.glob(os.path.join(base, "R-*", "bin", "Rscript.exe")):
            if os.path.isfile(p):
                return p
    return None

def check_raw_counts_and_metadata():
    """Verify raw counts and metadata exist"""
    data_dir = "data/sc_rnaseq/GSE122897"
    raw_file = os.path.join(data_dir, "GSE122897_raw_counts_GRCh38.p13_NCBI.tsv.gz")
    meta_files = [
        os.path.join(data_dir, "sample_metadata_with_clinical.tsv"),
        os.path.join(data_dir, "sample_metadata.tsv"),
    ]
    for m in meta_files:
        if os.path.exists(m):
            meta_file = m
            break
    else:
        meta_file = None

    if not os.path.exists(raw_file):
        print(f"[ERROR] Raw counts not found: {raw_file}")
        return None, None
    if meta_file is None:
        print("[ERROR] No sample metadata found")
        return None, None

    return raw_file, meta_file

def create_r_script():
    """Create R script for DESeq2 analysis of GSE122897"""
    r_script = '''
# DESeq2 Analysis for GSE122897 (ruptured vs unruptured aneurysm tissue)
# Run with: Rscript run_deseq2_gse122897.R
# Raw counts: GSE122897_raw_counts_GRCh38.p13_NCBI.tsv.gz
# Gene IDs: NCBI/Entrez (integer)

if (!require("DESeq2", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("DESeq2")
}

library(DESeq2)

cat("DESeq2 analysis for GSE122897 (tissue RNA-seq)\\n")
cat("============================================\\n")

# Paths (run from project root)
counts_file <- "data/sc_rnaseq/GSE122897/GSE122897_raw_counts_GRCh38.p13_NCBI.tsv.gz"
meta_file <- "data/sc_rnaseq/GSE122897/sample_metadata_with_clinical.tsv"
if (!file.exists(meta_file)) meta_file <- "data/sc_rnaseq/GSE122897/sample_metadata.tsv"
out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE)

# Load counts
cat("\\nLoading counts...\\n")
counts <- read.delim(counts_file, row.names = 1, check.names = FALSE)
counts <- round(counts)  # ensure integers
cat("  Genes:", nrow(counts), " Samples:", ncol(counts), "\\n")

# Load metadata
meta <- read.delim(meta_file, sep = "\\t", stringsAsFactors = FALSE)
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

cat("  Ruptured:", sum(meta$Group == "Ruptured"), " Unruptured:", sum(meta$Group == "Unruptured"), "\\n")

# Filter low-count genes (keep genes with >= 10 total counts)
keep <- rowSums(counts) >= 10
counts <- counts[keep, ]
cat("  Genes after filtering (sum >= 10):", nrow(counts), "\\n")

# DESeq2
cat("\\nRunning DESeq2...\\n")
dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ Group)
dds <- DESeq(dds)
res <- results(dds, contrast = c("Group", "Ruptured", "Unruptured"))
res <- res[order(res$padj, na.last = TRUE), ]

# Write results
out_file <- file.path(out_dir, "gse122897_deseq2_results.tsv")
res_df <- as.data.frame(res)
res_df$Gene_ID <- rownames(res_df)
res_df <- res_df[, c("Gene_ID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
write.table(res_df, out_file, sep = "\\t", row.names = FALSE)
cat("\\nResults saved to:", out_file, "\\n")

# Summary
cat("\\nSummary:\\n")
cat("  Significant (padj < 0.05):", sum(res_df$padj < 0.05, na.rm = TRUE), "\\n")
cat("  Significant (padj < 0.1):", sum(res_df$padj < 0.1, na.rm = TRUE), "\\n")
'''
    return r_script

def main():
    raw_file, meta_file = check_raw_counts_and_metadata()
    if raw_file is None:
        return

    print("GSE122897 DESeq2 Analysis")
    print("=" * 50)
    print(f"Raw counts: {raw_file}")
    print(f"Metadata: {meta_file}")

    # Write R script
    r_script = create_r_script()
    script_path = "run_deseq2_gse122897.R"
    with open(script_path, "w") as f:
        f.write(r_script)
    print(f"\nR script written to: {script_path}")

    # Try to run R
    rscript = find_rscript()
    if not rscript:
        print("\n[INFO] Rscript not found in PATH or common install locations.")
        print("  Add R to PATH, or run manually: Rscript run_deseq2_gse122897.R")
        return
    try:
        result = subprocess.run(
            [rscript, script_path],
            capture_output=True,
            text=True,
            timeout=300,
            cwd=os.getcwd()
        )
        if result.returncode == 0:
            print(result.stdout)
            print("\n[OK] DESeq2 analysis completed successfully")
            print("  Results: results/gse122897_deseq2_results.tsv")
        else:
            print("[ERROR] R script failed")
            print(result.stderr or result.stdout)
    except FileNotFoundError:
        print("\n[INFO] R not found. To run DESeq2:")
        print("  1. Install R from https://cran.r-project.org/")
        print("  2. Install DESeq2: install.packages('BiocManager'); BiocManager::install('DESeq2')")
        print("  3. Run: Rscript run_deseq2_gse122897.R")
    except subprocess.TimeoutExpired:
        print("[ERROR] DESeq2 analysis timed out (5 min)")

if __name__ == "__main__":
    main()
