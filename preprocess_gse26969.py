"""
Preprocess GSE26969 Microarray Data (Affymetrix CEL files)
This script normalizes CEL files and prepares expression matrix for analysis
"""
import pandas as pd
import numpy as np
import os
import gzip
from pathlib import Path
import subprocess
import sys

def check_r_availability():
    """Check if R is available for processing CEL files"""
    try:
        result = subprocess.run(['R', '--version'], capture_output=True, text=True)
        if result.returncode == 0:
            print("[OK] R is available")
            return True
    except FileNotFoundError:
        pass
    
    try:
        result = subprocess.run(['Rscript', '--version'], capture_output=True, text=True)
        if result.returncode == 0:
            print("[OK] Rscript is available")
            return True
    except FileNotFoundError:
        pass
    
    print("[WARNING] R/Rscript not found. Will use Python-based approach (limited)")
    return False

def create_r_normalization_script():
    """Create R script for CEL file normalization"""
    r_script = """
# GSE26969 CEL file normalization using R
library(affy)
library(Biobase)

# Set working directory
cel_dir <- "results/gse26969_extracted"
setwd(cel_dir)

# List CEL files
cel_files <- list.files(pattern="\\.CEL\\.gz$", full.names=TRUE)

if(length(cel_files) == 0) {
    stop("No CEL.gz files found")
}

print(paste("Found", length(cel_files), "CEL files"))

# Read CEL files
print("Reading CEL files...")
affy_data <- ReadAffy(filenames=cel_files)

# RMA normalization
print("Performing RMA normalization...")
eset <- rma(affy_data)

# Extract expression matrix
expr_matrix <- exprs(eset)

# Get sample names from CEL file names
sample_names <- gsub("\\.CEL\\.gz$", "", basename(cel_files))
colnames(expr_matrix) <- sample_names

# Convert to data frame
expr_df <- as.data.frame(expr_matrix)
expr_df$ProbeID <- rownames(expr_matrix)

# Save normalized expression matrix
output_file <- "../../results/gse26969_normalized_expression.tsv"
write.table(expr_df, file=output_file, sep="\\t", row.names=FALSE, quote=FALSE)

print(paste("Normalized expression matrix saved to:", output_file))
print(paste("Expression matrix dimensions:", nrow(expr_matrix), "probes x", ncol(expr_matrix), "samples"))
"""
    
    script_path = Path("results/gse26969_normalize.R")
    script_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(script_path, 'w') as f:
        f.write(r_script)
    
    return script_path

def normalize_cel_files_r():
    """Normalize CEL files using R and affy package"""
    print("=" * 80)
    print("NORMALIZING GSE26969 CEL FILES USING R")
    print("=" * 80)
    
    # Check if R is available
    if not check_r_availability():
        return None
    
    # Create R script
    r_script_path = create_r_normalization_script()
    print(f"\nCreated R normalization script: {r_script_path}")
    
    # Run R script
    print("\nRunning R normalization script...")
    try:
        result = subprocess.run(
            ['Rscript', str(r_script_path)],
            capture_output=True,
            text=True,
            cwd=Path("results/gse26969_extracted")
        )
        
        if result.returncode == 0:
            print("[OK] R normalization completed successfully")
            print(result.stdout)
            
            # Check if output file was created
            output_file = Path("results/gse26969_normalized_expression.tsv")
            if output_file.exists():
                print(f"\n[OK] Normalized expression matrix saved: {output_file}")
                return output_file
            else:
                print(f"\n[WARNING] Output file not found: {output_file}")
                return None
        else:
            print("[ERROR] R script failed:")
            print(result.stderr)
            return None
            
    except Exception as e:
        print(f"[ERROR] Failed to run R script: {e}")
        return None

def annotate_probe_ids(expr_df):
    """Annotate Affymetrix probe IDs with gene symbols"""
    print("\n" + "=" * 80)
    print("ANNOTATING PROBE IDS WITH GENE SYMBOLS")
    print("=" * 80)
    
    # For now, we'll keep probe IDs and note that annotation is needed
    # In practice, would use biomaRt or AnnotationDbi in R
    print("\n[NOTE] Probe ID annotation requires gene annotation database")
    print("   This can be done using:")
    print("   1. biomaRt in R")
    print("   2. AnnotationDbi with appropriate chip package")
    print("   3. Manual mapping file from GEO")
    
    return expr_df

def load_normalized_data():
    """Load normalized expression data if it exists"""
    output_file = Path("results/gse26969_normalized_expression.tsv")
    
    if output_file.exists():
        print(f"\nLoading normalized expression data: {output_file}")
        expr_df = pd.read_csv(output_file, sep='\t', index_col='ProbeID')
        print(f"[OK] Loaded expression matrix: {expr_df.shape}")
        return expr_df
    else:
        print(f"\n[WARNING] Normalized data not found: {output_file}")
        return None

def create_metadata_template():
    """Create metadata template for GSE26969 samples"""
    print("\n" + "=" * 80)
    print("CREATING METADATA TEMPLATE")
    print("=" * 80)
    
    # Sample IDs from extracted files
    extract_dir = Path("results/gse26969_extracted")
    cel_files = list(extract_dir.glob("*.CEL.gz"))
    
    sample_ids = [f.stem.replace('.CEL', '') for f in cel_files]
    
    # Create metadata template
    metadata = pd.DataFrame({
        'Sample': sample_ids,
        'Group': ['Unknown'] * len(sample_ids),  # IA or STA
        'Type': ['Unknown'] * len(sample_ids),    # Aneurysm or Control
        'Rupture_Status': ['Unknown'] * len(sample_ids)  # If applicable
    })
    
    output_file = Path("results/gse26969_metadata_template.tsv")
    metadata.to_csv(output_file, sep='\t', index=False)
    
    print(f"\n[OK] Created metadata template: {output_file}")
    print("\nPlease update metadata with:")
    print("  - Group: 'IA' (intracranial aneurysm) or 'STA' (superficial temporal artery control)")
    print("  - Type: 'Aneurysm' or 'Control'")
    print("  - Rupture_Status: 'Ruptured', 'Unruptured', or 'N/A' (for controls)")
    
    return metadata

def quality_control_microarray(expr_df):
    """Perform quality control on microarray data"""
    print("\n" + "=" * 80)
    print("QUALITY CONTROL - MICROARRAY DATA")
    print("=" * 80)
    
    # Calculate QC metrics
    mean_expr = expr_df.mean(axis=0)
    median_expr = expr_df.median(axis=0)
    detected_probes = (expr_df > 0).sum(axis=0)
    
    print(f"\nExpression statistics per sample:")
    print(f"  Mean expression:")
    print(f"    Mean: {mean_expr.mean():.2f}")
    print(f"    Range: {mean_expr.min():.2f} - {mean_expr.max():.2f}")
    
    print(f"\n  Median expression:")
    print(f"    Mean: {median_expr.mean():.2f}")
    print(f"    Range: {median_expr.min():.2f} - {median_expr.max():.2f}")
    
    print(f"\n  Detected probes per sample:")
    print(f"    Mean: {detected_probes.mean():,.0f}")
    print(f"    Range: {detected_probes.min():,.0f} - {detected_probes.max():,.0f}")
    
    # Filter lowly expressed probes
    # Keep probes with expression > 5 in at least 2 samples
    min_expr = 5
    min_samples = 2
    expressed = (expr_df > min_expr).sum(axis=1) >= min_samples
    
    print(f"\n  Probe filtering:")
    print(f"    Probes before: {len(expr_df):,}")
    print(f"    Probes after (expr > {min_expr} in >= {min_samples} samples): {expressed.sum():,}")
    print(f"    Removed: {(~expressed).sum():,} lowly expressed probes")
    
    expr_df_qc = expr_df[expressed]
    
    return expr_df_qc

def main():
    """Main preprocessing pipeline"""
    
    print("=" * 80)
    print("GSE26969 MICROARRAY PREPROCESSING PIPELINE")
    print("=" * 80)
    
    # Create results directory
    os.makedirs("results", exist_ok=True)
    
    # Check if CEL files are extracted
    extract_dir = Path("results/gse26969_extracted")
    cel_files = list(extract_dir.glob("*.CEL.gz"))
    
    if not cel_files:
        print(f"\n[ERROR] No CEL.gz files found in {extract_dir}")
        print("   Please extract GSE26969_RAW.tar first")
        return
    
    print(f"\n[OK] Found {len(cel_files)} CEL files")
    
    # Step 1: Normalize CEL files
    normalized_file = normalize_cel_files_r()
    
    if normalized_file is None:
        print("\n[WARNING] R normalization not available or failed")
        print("   Alternative: Use Python-based approach or install R with affy package")
        print("   For now, creating metadata template only")
    else:
        # Step 2: Load normalized data
        expr_df = load_normalized_data()
        
        if expr_df is not None:
            # Step 3: Quality control
            expr_df_qc = quality_control_microarray(expr_df)
            
            # Save QC'd data
            output_file = Path("results/gse26969_expression_qc.tsv")
            expr_df_qc.to_csv(output_file, sep='\t')
            print(f"\n[OK] Saved QC'd expression matrix: {output_file}")
    
    # Step 4: Create metadata template
    metadata = create_metadata_template()
    
    print("\n" + "=" * 80)
    print("PREPROCESSING COMPLETE")
    print("=" * 80)
    print("\nNext steps:")
    print("1. Update metadata template with sample annotations")
    print("2. Annotate probe IDs with gene symbols (if needed)")
    print("3. Run integrated analysis with GSE122897")

if __name__ == "__main__":
    main()






