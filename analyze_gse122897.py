"""
Analyze GSE122897 - Ruptured vs Unruptured Aneurysm RNA-seq
This script performs QC, normalization, and differential expression analysis
"""
import pandas as pd
import numpy as np
import os
import sys
from scipy import stats
from statsmodels.stats.multitest import multipletests

# Fix encoding for Windows console - use ASCII-safe characters
# Replace Unicode symbols with ASCII equivalents for Windows compatibility

def load_and_prepare_data():
    """Load GSE122897 data and prepare for analysis"""
    
    print("=" * 80)
    print("LOADING GSE122897 DATA")
    print("=" * 80)
    
    data_dir = "data/sc_rnaseq/GSE122897"
    
    # Load TPM-normalized data (already normalized, easier to work with)
    tpm_file = os.path.join(data_dir, "GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz")
    
    if not os.path.exists(tpm_file):
        print(f"[ERROR] TPM file not found: {tpm_file}")
        return None, None
    
    print(f"\nLoading TPM-normalized data: {tpm_file}")
    expr_df = pd.read_csv(tpm_file, sep='\t', compression='gzip', index_col=0)
    
    print(f"[OK] Loaded expression matrix:")
    print(f"  Genes: {expr_df.shape[0]:,}")
    print(f"  Samples: {expr_df.shape[1]}")
    
    # Get sample names
    samples = expr_df.columns.tolist()
    print(f"\nSample IDs (first 10): {samples[:10]}")
    print(f"Total samples: {len(samples)}")
    
    # Try to load metadata file
    metadata_file = os.path.join(data_dir, "sample_metadata.tsv")
    if not os.path.exists(metadata_file):
        metadata_file = "results/gse122897_sample_metadata.tsv"
    
    if os.path.exists(metadata_file):
        print(f"\nLoading sample metadata: {metadata_file}")
        metadata = pd.read_csv(metadata_file, sep='\t')
        print(f"[OK] Loaded metadata: {metadata.shape}")
        print(f"  Columns: {list(metadata.columns)}")
        
        # Check if Group column exists
        if 'Group' in metadata.columns:
            group_counts = metadata['Group'].value_counts()
            print(f"\n[OK] Sample groups:")
            print(group_counts.to_string())
        else:
            print("\n[WARNING] No 'Group' column in metadata")
            metadata['Group'] = 'Unknown'
    else:
        print("\n[WARNING] Metadata file not found. Creating placeholder...")
        metadata = pd.DataFrame({
            'Sample': samples,
            'Group': ['Unknown'] * len(samples)
        })
        print("   Please run extract_gse122897_metadata.py first to extract metadata")
    
    return expr_df, metadata

def quality_control(expr_df):
    """Perform quality control checks"""
    
    print("\n" + "=" * 80)
    print("QUALITY CONTROL")
    print("=" * 80)
    
    # Calculate QC metrics
    total_expr = expr_df.sum(axis=0)
    detected_genes = (expr_df > 0).sum(axis=0)
    pct_zeros = (expr_df == 0).sum(axis=0) / len(expr_df) * 100
    
    print(f"\nExpression statistics:")
    print(f"  Total TPM per sample:")
    print(f"    Mean: {total_expr.mean():,.0f}")
    print(f"    Median: {total_expr.median():,.0f}")
    print(f"    Range: {total_expr.min():,.0f} - {total_expr.max():,.0f}")
    
    print(f"\n  Detected genes per sample:")
    print(f"    Mean: {detected_genes.mean():,.0f}")
    print(f"    Median: {detected_genes.median():,.0f}")
    print(f"    Range: {detected_genes.min():,.0f} - {detected_genes.max():,.0f}")
    
    print(f"\n  Percentage of zeros:")
    print(f"    Mean: {pct_zeros.mean():.2f}%")
    print(f"    Median: {pct_zeros.median():.2f}%")
    
    # Filter lowly expressed genes
    # Keep genes with TPM > 1 in at least 3 samples
    min_tpm = 1
    min_samples = 3
    expressed = (expr_df > min_tpm).sum(axis=1) >= min_samples
    
    print(f"\n  Gene filtering:")
    print(f"    Genes before: {len(expr_df):,}")
    print(f"    Genes after (TPM > {min_tpm} in >= {min_samples} samples): {expressed.sum():,}")
    print(f"    Removed: {(~expressed).sum():,} lowly expressed genes")
    
    expr_df_qc = expr_df[expressed]
    
    return expr_df_qc

def differential_expression(expr_df, metadata):
    """Perform differential expression analysis"""
    
    print("\n" + "=" * 80)
    print("DIFFERENTIAL EXPRESSION ANALYSIS")
    print("=" * 80)
    
    # Get sample groups
    if 'Group' not in metadata.columns:
        print("[ERROR] No 'Group' column in metadata!")
        return None
    
    # Check if groups are defined
    unique_groups = metadata['Group'].unique()
    if len(unique_groups) < 2 or 'Unknown' in unique_groups:
        print("\n[WARNING] Sample groups not properly defined!")
        print(f"   Current groups: {unique_groups}")
        print("\n   Please update metadata with actual group labels (Ruptured/Unruptured)")
        print("   Then re-run this analysis.")
        return None
    
    # Get group 1 and group 2 samples
    group1_name = unique_groups[0]
    group2_name = unique_groups[1]
    
    group1_samples = metadata[metadata['Group'] == group1_name]['Sample'].tolist()
    group2_samples = metadata[metadata['Group'] == group2_name]['Sample'].tolist()
    
    # Align with expression matrix
    group1_samples = [s for s in group1_samples if s in expr_df.columns]
    group2_samples = [s for s in group2_samples if s in expr_df.columns]
    
    if not group1_samples or not group2_samples:
        print(f"[ERROR] No matching samples found!")
        print(f"   {group1_name}: {len(group1_samples)} samples")
        print(f"   {group2_name}: {len(group2_samples)} samples")
        return None
    
    print(f"\nSample groups:")
    print(f"  {group1_name}: {len(group1_samples)} samples")
    print(f"  {group2_name}: {len(group2_samples)} samples")
    
    print(f"\nPerforming differential expression analysis...")
    print(f"  Testing {len(expr_df):,} genes")
    print(f"  Using vectorized operations for faster computation...")
    
    # Subset expression data to only the samples we need
    expr_subset = expr_df[group1_samples + group2_samples].copy()
    
    # Separate into groups
    expr_group1 = expr_subset[group1_samples].values
    expr_group2 = expr_subset[group2_samples].values
    
    # Log2 transform
    expr_group1_log = np.log2(expr_group1 + 1)
    expr_group2_log = np.log2(expr_group2 + 1)
    
    # Vectorized t-tests
    print("  Computing t-statistics and p-values...")
    from scipy.stats import ttest_ind
    t_stats, p_values = ttest_ind(expr_group1_log, expr_group2_log, axis=1)
    
    # Calculate means and fold changes (vectorized)
    mean_group1 = np.mean(expr_group1, axis=1)
    mean_group2 = np.mean(expr_group2, axis=1)
    
    # Calculate log2 fold change with pseudocount
    log2_fc = np.log2((mean_group1 + 1e-6) / (mean_group2 + 1e-6))
    
    # Handle special cases
    log2_fc[np.isnan(log2_fc)] = 0
    log2_fc[~np.isfinite(log2_fc)] = 0
    
    # Create results DataFrame
    de_results = pd.DataFrame({
        'Gene': expr_df.index.astype(str),  # Ensure Gene is string type
        f'Mean_{group1_name}': mean_group1,
        f'Mean_{group2_name}': mean_group2,
        'Log2FC': log2_fc,
        'T_statistic': t_stats,
        'P_value': p_values
    })
    
    # Multiple testing correction
    # Handle NaN p-values by setting them to 1.0 for correction, then restore NaN
    p_values = de_results['P_value'].values.copy()
    nan_mask = np.isnan(p_values)
    
    # Replace NaN with 1.0 for correction
    p_values_for_correction = p_values.copy()
    p_values_for_correction[nan_mask] = 1.0
    
    # Apply multiple testing correction
    _, p_adj_bh, _, _ = multipletests(p_values_for_correction, method='fdr_bh')
    _, p_adj_bonf, _, _ = multipletests(p_values_for_correction, method='bonferroni')
    
    # Restore NaN values
    p_adj_bh[nan_mask] = np.nan
    p_adj_bonf[nan_mask] = np.nan
    
    de_results['P_adj_BH'] = p_adj_bh
    de_results['P_adj_Bonferroni'] = p_adj_bonf
    
    # Sort by significance
    de_results = de_results.sort_values('P_value')
    
    print(f"\n[OK] Differential expression analysis complete")
    print(f"\nResults summary:")
    print(f"  Total genes tested: {len(de_results):,}")
    print(f"  Significant (P < 0.05): {(de_results['P_value'] < 0.05).sum():,}")
    print(f"  Significant (P_adj < 0.05): {(de_results['P_adj_BH'] < 0.05).sum():,}")
    print(f"  Highly significant (P_adj < 0.01): {(de_results['P_adj_BH'] < 0.01).sum():,}")
    
    return de_results, group1_name, group2_name

def identify_top_genes(de_results, group1_name, group2_name, top_n=50):
    """Identify top differentially expressed genes"""
    
    print("\n" + "=" * 80)
    print(f"TOP {top_n} DIFFERENTIALLY EXPRESSED GENES")
    print("=" * 80)
    
    # Filter significant genes
    sig_genes = de_results[de_results['P_adj_BH'] < 0.05].copy()
    
    if len(sig_genes) == 0:
        print("\n[WARNING] No significant genes found (P_adj < 0.05)")
        print("   Showing top genes by P-value instead...")
        sig_genes = de_results.head(100).copy()
    
    # Sort by absolute log2FC and significance
    sig_genes['Abs_Log2FC'] = np.abs(sig_genes['Log2FC'])
    sig_genes = sig_genes.sort_values(['P_adj_BH', 'Abs_Log2FC'], ascending=[True, False])
    
    top_genes = sig_genes.head(top_n)
    
    # Separate upregulated and downregulated
    upregulated = top_genes[top_genes['Log2FC'] > 0].copy()
    downregulated = top_genes[top_genes['Log2FC'] < 0].copy()
    
    print(f"\nUpregulated in {group1_name} ({len(upregulated)} genes):")
    for idx, row in upregulated.head(20).iterrows():
        gene_name = str(row['Gene'])  # Ensure string type
        print(f"  {gene_name:20s} Log2FC={row['Log2FC']:7.3f}  P_adj={row['P_adj_BH']:.2e}")
    
    print(f"\nDownregulated in {group1_name} ({len(downregulated)} genes):")
    for idx, row in downregulated.head(20).iterrows():
        gene_name = str(row['Gene'])  # Ensure string type
        print(f"  {gene_name:20s} Log2FC={row['Log2FC']:7.3f}  P_adj={row['P_adj_BH']:.2e}")
    
    return top_genes

def main():
    """Main analysis pipeline"""
    
    # Create results directory
    os.makedirs("results", exist_ok=True)
    
    # Step 1: Load data
    expr_df, metadata = load_and_prepare_data()
    if expr_df is None:
        print("[ERROR] Failed to load data. Exiting.")
        return
    
    # Step 2: Quality control
    expr_df_qc = quality_control(expr_df)
    
    # Step 3: Check if metadata has group labels
    # Filter out Control and Unknown samples for rupture comparison
    metadata_aneurysm = metadata[metadata['Group'].isin(['Ruptured', 'Unruptured'])].copy()
    
    if len(metadata_aneurysm) < 2 or metadata_aneurysm['Group'].nunique() < 2:
        print("\n" + "=" * 80)
        print("[WARNING] SAMPLE GROUPS NOT PROPERLY DEFINED")
        print("=" * 80)
        print(f"\nCurrent groups: {metadata['Group'].value_counts().to_dict()}")
        print("\nTo complete the rupture comparison analysis, you need:")
        print("1. At least 2 samples in both 'Ruptured' and 'Unruptured' groups")
        print("2. Run extract_gse122897_metadata.py to extract metadata from series matrix")
        print("3. Re-run this script")
        return
    
    print(f"\n[OK] Using {len(metadata_aneurysm)} aneurysm samples for rupture comparison:")
    print(metadata_aneurysm['Group'].value_counts().to_string())
    
    # Step 4: Differential expression (ruptured vs unruptured only)
    de_results, group1, group2 = differential_expression(expr_df_qc, metadata_aneurysm)
    
    if de_results is not None:
        # Step 5: Identify top genes
        top_genes = identify_top_genes(de_results, group1, group2, top_n=50)
        
        # Save results
        output_file = "results/gse122897_differential_expression.tsv"
        de_results.to_csv(output_file, sep='\t', index=False)
        print(f"\n[OK] Saved full DE results to {output_file}")
        
        top_file = "results/gse122897_top_genes.tsv"
        top_genes.to_csv(top_file, sep='\t', index=False)
        print(f"[OK] Saved top genes to {top_file}")
        
        print("\n" + "=" * 80)
        print("ANALYSIS COMPLETE!")
        print("=" * 80)
    else:
        print("\n[WARNING] Could not complete differential expression analysis")

if __name__ == "__main__":
    main()

