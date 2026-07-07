"""
DESeq2 Analysis for Blood RNA-seq Data (GSE159610)
Proper count-based differential expression analysis
"""
import pandas as pd
import numpy as np
import os
import subprocess
import sys

def check_deseq2_availability():
    """Check if DESeq2 is available (via R or pyDESeq2)"""
    print("=" * 80)
    print("CHECKING DESEQ2 AVAILABILITY")
    print("=" * 80)
    
    # Check for pyDESeq2 (Python wrapper)
    try:
        import pydeseq2
        print("[OK] pyDESeq2 is available")
        return 'pydeseq2'
    except ImportError:
        print("[INFO] pyDESeq2 not available")
    
    # Check for R and DESeq2
    try:
        result = subprocess.run(['Rscript', '--version'], 
                              capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            print("[OK] R is available")
            # Check if DESeq2 is installed in R
            check_script = """
            if (!require("DESeq2", quietly = TRUE)) {
                cat("DESeq2_NOT_INSTALLED")
            } else {
                cat("DESeq2_AVAILABLE")
            }
            """
            result = subprocess.run(['Rscript', '-e', check_script],
                                  capture_output=True, text=True, timeout=10)
            if 'DESeq2_AVAILABLE' in result.stdout:
                print("[OK] DESeq2 R package is available")
                return 'r_deseq2'
            else:
                print("[INFO] DESeq2 R package not installed")
                return 'r_available'
        else:
            print("[INFO] R is not available")
    except (subprocess.TimeoutExpired, FileNotFoundError):
        print("[INFO] R is not available")
    
    return None

def create_r_deseq2_script():
    """Create R script for DESeq2 analysis"""
    r_script = """
# DESeq2 Analysis for GSE159610 Blood RNA-seq
# Run with: Rscript deseq2_blood_analysis.R

library(DESeq2)
library(readr)

# Load data
expr_file <- "data/blood_rnaseq/GSE159610/GSE159610_raw_gene_counts_matrix.txt.gz"
metadata_file <- "data/blood_rnaseq/GSE159610/sample_metadata_with_labels.tsv"

cat("Loading expression data...\\n")
counts <- read_delim(expr_file, delim = "\\t", escape_double = FALSE, 
                     trim_ws = TRUE)
rownames(counts) <- counts[[1]]
counts <- counts[, -1]

# Filter to DSA columns only
dsa_cols <- grep("^DSA", colnames(counts), value = TRUE)
dsa_cols <- dsa_cols[order(as.numeric(gsub("DSA", "", dsa_cols)))]
counts <- counts[, dsa_cols]

cat("Loading metadata...\\n")
metadata <- read_delim(metadata_file, delim = "\\t", escape_double = FALSE, 
                       trim_ws = TRUE)

# Create sample groups
status <- metadata$Status
names(status) <- dsa_cols[1:length(status)]

# Create DESeq2 object
cat("Creating DESeq2 object...\\n")
colData <- data.frame(
    row.names = names(status),
    condition = factor(status, levels = c("Control", "Aneurysm"))
)

# Ensure counts are integers
counts_matrix <- as.matrix(counts)
counts_matrix[is.na(counts_matrix)] <- 0
counts_matrix <- round(counts_matrix)
mode(counts_matrix) <- "integer"

# Filter low count genes
keep <- rowSums(counts_matrix >= 10) >= 3
counts_matrix <- counts_matrix[keep, ]

cat("Running DESeq2 analysis...\\n")
dds <- DESeqDataSetFromMatrix(
    countData = counts_matrix,
    colData = colData,
    design = ~ condition
)

# Run DESeq2
dds <- DESeq(dds)

# Get results
results <- results(dds, contrast = c("condition", "Aneurysm", "Control"))
results_df <- as.data.frame(results)
results_df$Gene_ID <- rownames(results_df)

# Save results
output_file <- "results/gse159610_deseq2_results.tsv"
write_delim(results_df, output_file, delim = "\\t")

cat("\\nDESeq2 analysis complete!\\n")
cat(paste("Results saved to:", output_file, "\\n"))
cat(paste("Total genes tested:", nrow(results_df), "\\n"))
cat(paste("Genes with padj < 0.1:", sum(results_df$padj < 0.1, na.rm = TRUE), "\\n"))
"""
    
    script_file = "deseq2_blood_analysis.R"
    with open(script_file, 'w') as f:
        f.write(r_script)
    
    print(f"[OK] Created R script: {script_file}")
    return script_file

def run_r_deseq2_analysis(script_file):
    """Run DESeq2 analysis using R script"""
    print("\n" + "=" * 80)
    print("RUNNING DESEQ2 ANALYSIS (R)")
    print("=" * 80)
    
    print(f"\nExecuting R script: {script_file}")
    print("This may take several minutes...")
    
    try:
        result = subprocess.run(['Rscript', script_file],
                              capture_output=True, text=True, timeout=600)
        
        if result.returncode == 0:
            print("\n[OK] DESeq2 analysis completed successfully")
            print(result.stdout)
            return True
        else:
            print("\n[ERROR] DESeq2 analysis failed")
            print("STDOUT:", result.stdout)
            print("STDERR:", result.stderr)
            return False
    except subprocess.TimeoutExpired:
        print("\n[ERROR] Analysis timed out (>10 minutes)")
        return False
    except Exception as e:
        print(f"\n[ERROR] Failed to run R script: {e}")
        return False

def analyze_deseq2_results():
    """Analyze DESeq2 results if available"""
    results_file = "results/gse159610_deseq2_results.tsv"
    
    if not os.path.exists(results_file):
        print(f"[INFO] DESeq2 results not found: {results_file}")
        return None
    
    print("\n" + "=" * 80)
    print("DESEQ2 RESULTS ANALYSIS")
    print("=" * 80)
    
    df = pd.read_csv(results_file, sep='\t')
    
    print(f"\nTotal genes tested: {len(df):,}")
    print(f"Genes with valid padj: {(~df['padj'].isna()).sum():,}")
    
    if (~df['padj'].isna()).sum() > 0:
        print(f"\nSignificance (adjusted p-value):")
        print(f"  padj < 0.05: {(df['padj'] < 0.05).sum():,}")
        print(f"  padj < 0.1: {(df['padj'] < 0.1).sum():,}")
        print(f"  padj < 0.2: {(df['padj'] < 0.2).sum():,}")
    
    print(f"\nFold change statistics:")
    abs_fc = np.abs(df['log2FoldChange'])
    print(f"  Mean |log2FC|: {abs_fc.mean():.3f}")
    print(f"  Median |log2FC|: {abs_fc.median():.3f}")
    print(f"  Genes with |log2FC| > 0.585: {(abs_fc > 0.585).sum():,}")
    print(f"  Genes with |log2FC| > 1.0: {(abs_fc > 1.0).sum():,}")
    
    # Top genes
    if (~df['padj'].isna()).sum() > 0:
        sig_genes = df[df['padj'] < 0.1].copy()
        if len(sig_genes) > 0:
            sig_genes = sig_genes.sort_values('padj')
            print(f"\nTop 20 significant genes (padj < 0.1):")
            print(sig_genes.head(20)[['Gene_ID', 'log2FoldChange', 'pvalue', 'padj']].to_string(index=False))
        else:
            # Show top by p-value
            top_genes = df.nsmallest(20, 'pvalue')
            print(f"\nTop 20 genes by p-value (no padj < 0.1):")
            print(top_genes[['Gene_ID', 'log2FoldChange', 'pvalue', 'padj']].to_string(index=False))
    
    return df

def main():
    """Main DESeq2 analysis pipeline"""
    print("=" * 80)
    print("DESEQ2 ANALYSIS FOR BLOOD RNA-SEQ")
    print("=" * 80)
    
    # Check availability
    deseq2_type = check_deseq2_availability()
    
    if deseq2_type is None:
        print("\n" + "=" * 80)
        print("DESEQ2 NOT AVAILABLE")
        print("=" * 80)
        print("\nTo use DESeq2, you need to install one of:")
        print("1. R with DESeq2 package:")
        print("   - Install R from https://www.r-project.org/")
        print("   - Install DESeq2: install.packages('BiocManager'); BiocManager::install('DESeq2')")
        print("2. pyDESeq2 (Python wrapper):")
        print("   - pip install pydeseq2")
        print("\nAlternatively, continue with t-test based analysis (already completed)")
        return
    
    if deseq2_type == 'r_deseq2':
        # Create and run R script
        script_file = create_r_deseq2_script()
        success = run_r_deseq2_analysis(script_file)
        
        if success:
            analyze_deseq2_results()
    elif deseq2_type == 'r_available':
        print("\nR is available but DESeq2 is not installed.")
        print("Installing DESeq2...")
        install_script = """
        if (!require("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager")
        }
        BiocManager::install("DESeq2")
        """
        result = subprocess.run(['Rscript', '-e', install_script],
                              capture_output=True, text=True, timeout=300)
        if result.returncode == 0:
            print("[OK] DESeq2 installed. Re-run this script.")
        else:
            print("[ERROR] Failed to install DESeq2")
            print("Please install manually: BiocManager::install('DESeq2')")
    elif deseq2_type == 'pydeseq2':
        print("\n[INFO] pyDESeq2 is available but implementation not yet complete")
        print("Please use R-based DESeq2 for now, or implement pyDESeq2 analysis")

if __name__ == "__main__":
    main()






