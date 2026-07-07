"""
Setup and Run DESeq2 Analysis
Checks for R/DESeq2 availability and provides installation instructions
"""
import subprocess
import os
import sys

def check_r_installation():
    """Check if R is installed"""
    try:
        result = subprocess.run(['R', '--version'], 
                              capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            version_line = result.stdout.split('\n')[0]
            print(f"[OK] R is installed: {version_line}")
            return True
    except (subprocess.TimeoutExpired, FileNotFoundError):
        pass
    
    print("[INFO] R is not installed")
    return False

def check_deseq2_installation():
    """Check if DESeq2 is installed in R"""
    check_script = """
    if (!require("DESeq2", quietly = TRUE)) {
        cat("NOT_INSTALLED")
    } else {
        cat("INSTALLED")
        packageVersion("DESeq2")
    }
    """
    
    try:
        result = subprocess.run(['Rscript', '-e', check_script],
                              capture_output=True, text=True, timeout=10)
        if 'INSTALLED' in result.stdout:
            print("[OK] DESeq2 is installed")
            if 'packageVersion' in result.stdout:
                print(f"  Version info in output")
            return True
        else:
            print("[INFO] DESeq2 is not installed")
            return False
    except Exception as e:
        print(f"[INFO] Could not check DESeq2: {e}")
        return False

def install_deseq2():
    """Install DESeq2 via BiocManager"""
    print("\n" + "=" * 80)
    print("INSTALLING DESEQ2")
    print("=" * 80)
    
    install_script = """
    if (!require("BiocManager", quietly = TRUE)) {
        cat("Installing BiocManager...\\n")
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    cat("Installing DESeq2...\\n")
    BiocManager::install("DESeq2", update = FALSE, ask = FALSE)
    cat("Installation complete!\\n")
    """
    
    print("\nInstalling DESeq2 (this may take 5-10 minutes)...")
    try:
        result = subprocess.run(['Rscript', '-e', install_script],
                              capture_output=True, text=True, timeout=600)
        
        if result.returncode == 0:
            print("[OK] DESeq2 installation completed")
            print(result.stdout)
            return True
        else:
            print("[ERROR] DESeq2 installation failed")
            print("STDOUT:", result.stdout)
            print("STDERR:", result.stderr)
            return False
    except subprocess.TimeoutExpired:
        print("[ERROR] Installation timed out (>10 minutes)")
        return False
    except Exception as e:
        print(f"[ERROR] Installation failed: {e}")
        return False

def create_deseq2_script():
    """Create R script for DESeq2 analysis"""
    r_script = """# DESeq2 Analysis for GSE159610 Blood RNA-seq
# Run with: Rscript deseq2_blood_analysis.R

library(DESeq2)
library(readr)

cat("=", rep("=", 79), "\\n", sep="")
cat("DESEQ2 ANALYSIS FOR BLOOD RNA-SEQ\\n")
cat("=", rep("=", 79), "\\n\\n", sep="")

# Load data
expr_file <- "data/blood_rnaseq/GSE159610/GSE159610_raw_gene_counts_matrix.txt.gz"
metadata_file <- "data/blood_rnaseq/GSE159610/sample_metadata_with_labels.tsv"

cat("Loading expression data...\\n")
counts <- read_delim(expr_file, delim = "\\t", escape_double = FALSE, 
                     trim_ws = TRUE, show_col_types = FALSE)
rownames(counts) <- counts[[1]]
counts <- counts[, -1]

# Filter to DSA columns only
dsa_cols <- grep("^DSA", colnames(counts), value = TRUE)
dsa_cols <- dsa_cols[order(as.numeric(gsub("DSA", "", dsa_cols)))]
counts <- counts[, dsa_cols]

cat("  Loaded:", nrow(counts), "genes,", ncol(counts), "samples\\n")

cat("Loading metadata...\\n")
metadata <- read_delim(metadata_file, delim = "\\t", escape_double = FALSE, 
                       trim_ws = TRUE, show_col_types = FALSE)

# Create sample groups
status <- metadata$Status
names(status) <- dsa_cols[1:length(status)]

cat("  Sample groups:\\n")
cat("    Aneurysm:", sum(status == "Aneurysm"), "\\n")
cat("    Control:", sum(status == "Control"), "\\n")

# Create DESeq2 object
cat("\\nCreating DESeq2 object...\\n")
colData <- data.frame(
    row.names = names(status),
    condition = factor(status, levels = c("Control", "Aneurysm"))
)

# Ensure counts are integers
counts_matrix <- as.matrix(counts)
counts_matrix[is.na(counts_matrix)] <- 0
counts_matrix <- round(counts_matrix)
mode(counts_matrix) <- "integer"

# Filter low count genes (keep genes with >= 10 counts in >= 3 samples)
cat("Filtering low count genes...\\n")
keep <- rowSums(counts_matrix >= 10) >= 3
counts_matrix <- counts_matrix[keep, ]
cat("  Kept", nrow(counts_matrix), "genes after filtering\\n")

# Create DESeq2 dataset
cat("\\nRunning DESeq2 analysis...\\n")
dds <- DESeqDataSetFromMatrix(
    countData = counts_matrix,
    colData = colData,
    design = ~ condition
)

# Run DESeq2
dds <- DESeq(dds)

# Get results
cat("Extracting results...\\n")
results <- results(dds, contrast = c("condition", "Aneurysm", "Control"))
results_df <- as.data.frame(results)
results_df$Gene_ID <- rownames(results_df)

# Reorder columns
results_df <- results_df[, c("Gene_ID", "baseMean", "log2FoldChange", 
                              "lfcSE", "stat", "pvalue", "padj")]

# Save results
output_file <- "results/gse159610_deseq2_results.tsv"
write_delim(results_df, output_file, delim = "\\t")

cat("\\n", "=", rep("=", 79), "\\n", sep="")
cat("DESEQ2 ANALYSIS COMPLETE!\\n")
cat("=", rep("=", 79), "\\n\\n", sep="")
cat("Results saved to:", output_file, "\\n")
cat("Total genes tested:", nrow(results_df), "\\n")
cat("Genes with padj < 0.05:", sum(results_df$padj < 0.05, na.rm = TRUE), "\\n")
cat("Genes with padj < 0.1:", sum(results_df$padj < 0.1, na.rm = TRUE), "\\n")
cat("Genes with padj < 0.2:", sum(results_df$padj < 0.2, na.rm = TRUE), "\\n")
"""
    
    script_file = "deseq2_blood_analysis.R"
    with open(script_file, 'w') as f:
        f.write(r_script)
    
    print(f"[OK] Created R script: {script_file}")
    return script_file

def run_deseq2_analysis(script_file):
    """Run DESeq2 analysis"""
    print("\n" + "=" * 80)
    print("RUNNING DESEQ2 ANALYSIS")
    print("=" * 80)
    
    print(f"\nExecuting: Rscript {script_file}")
    print("This may take 5-10 minutes...")
    
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
        print(f"\n[ERROR] Failed to run analysis: {e}")
        return False

def main():
    """Main DESeq2 setup and analysis pipeline"""
    print("=" * 80)
    print("DESEQ2 SETUP AND ANALYSIS")
    print("=" * 80)
    
    # Check R
    if not check_r_installation():
        print("\n" + "=" * 80)
        print("R IS NOT INSTALLED")
        print("=" * 80)
        print("\nTo install R:")
        print("1. Download from: https://www.r-project.org/")
        print("2. Install R for your operating system")
        print("3. Re-run this script")
        return
    
    # Check DESeq2
    if not check_deseq2_installation():
        print("\nDESeq2 is not installed. Attempting to install...")
        if not install_deseq2():
            print("\n" + "=" * 80)
            print("DESEQ2 INSTALLATION FAILED")
            print("=" * 80)
            print("\nTo install manually:")
            print("1. Open R or RStudio")
            print("2. Run:")
            print("   install.packages('BiocManager')")
            print("   BiocManager::install('DESeq2')")
            print("3. Re-run this script")
            return
    
    # Create R script
    script_file = create_deseq2_script()
    
    # Run analysis
    success = run_deseq2_analysis(script_file)
    
    if success:
        print("\n" + "=" * 80)
        print("DESEQ2 ANALYSIS SUCCESSFUL!")
        print("=" * 80)
        print("\nResults saved to: results/gse159610_deseq2_results.tsv")
        print("\nNext steps:")
        print("1. Compare DESeq2 results with t-test results")
        print("2. Use DESeq2 results for downstream analysis")

if __name__ == "__main__":
    main()






