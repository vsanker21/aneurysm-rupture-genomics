"""
Full Differential Expression Analysis for GSE159610 Blood RNA-seq
Runs DE analysis on all genes using proper statistical methods
"""
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import os
import sys

def load_data():
    """Load expression data and metadata"""
    print("=" * 80)
    print("LOADING GSE159610 BLOOD RNA-SEQ DATA")
    print("=" * 80)
    
    expr_file = 'data/blood_rnaseq/GSE159610/GSE159610_raw_gene_counts_matrix.txt.gz'
    metadata_file = 'data/blood_rnaseq/GSE159610/sample_metadata_with_labels.tsv'
    
    if not os.path.exists(expr_file):
        print(f"[ERROR] Expression file not found: {expr_file}")
        sys.exit(1)
    
    if not os.path.exists(metadata_file):
        print(f"[ERROR] Metadata file not found: {metadata_file}")
        sys.exit(1)
    
    print(f"\nLoading expression data...")
    expr_df = pd.read_csv(expr_file, sep='\t', compression='gzip', index_col=0)
    print(f"  [OK] Expression matrix: {expr_df.shape[0]:,} genes, {expr_df.shape[1]} samples")
    
    print(f"\nLoading metadata...")
    metadata = pd.read_csv(metadata_file, sep='\t')
    print(f"  [OK] Metadata: {metadata.shape}")
    
    # Get sample columns
    data_cols = sorted([c for c in expr_df.columns if c.startswith('DSA')], 
                       key=lambda x: int(x.replace('DSA', '')))
    
    # Identify samples
    aneurysm_samples = [data_cols[i] for i in range(len(metadata)) 
                        if metadata.iloc[i]['Status'] == 'Aneurysm']
    control_samples = [data_cols[i] for i in range(len(metadata)) 
                       if metadata.iloc[i]['Status'] == 'Control']
    
    print(f"\nSample groups:")
    print(f"  Aneurysm: {len(aneurysm_samples)} samples")
    print(f"  Control: {len(control_samples)} samples")
    
    return expr_df, aneurysm_samples, control_samples

def perform_de_analysis(expr_df, aneurysm_samples, control_samples):
    """
    Perform differential expression analysis on all genes using vectorized operations
    """
    print("\n" + "=" * 80)
    print("DIFFERENTIAL EXPRESSION ANALYSIS")
    print("=" * 80)
    
    print(f"\nTesting {len(expr_df):,} genes...")
    print("  Using vectorized operations for faster computation...")
    
    # Subset expression data to only the samples we need
    expr_subset = expr_df[aneurysm_samples + control_samples].copy()
    
    # Separate into groups
    expr_aneurysm = expr_subset[aneurysm_samples].values.astype(float)
    expr_control = expr_subset[control_samples].values.astype(float)
    
    print("  Filtering lowly expressed genes...")
    
    # Filter genes with low expression (vectorized)
    # Need at least 3 samples with expression > 0 in each group
    aneurysm_expressed = (expr_aneurysm > 0).sum(axis=1) >= 3
    control_expressed = (expr_control > 0).sum(axis=1) >= 3
    expressed_mask = aneurysm_expressed & control_expressed
    
    print(f"  Genes passing expression filter: {expressed_mask.sum():,} / {len(expr_df):,}")
    
    # Apply filter
    expr_aneurysm_filtered = expr_aneurysm[expressed_mask]
    expr_control_filtered = expr_control[expressed_mask]
    gene_ids_filtered = expr_subset.index[expressed_mask]
    
    print("  Computing statistics...")
    
    # Log2 transform (vectorized)
    expr_aneurysm_log = np.log2(expr_aneurysm_filtered + 1)
    expr_control_log = np.log2(expr_control_filtered + 1)
    
    # Vectorized t-tests
    from scipy.stats import ttest_ind
    t_stats, p_values = ttest_ind(expr_aneurysm_log, expr_control_log, axis=1)
    
    # Calculate means and fold changes (vectorized)
    mean_aneurysm = np.mean(expr_aneurysm_filtered, axis=1)
    mean_control = np.mean(expr_control_filtered, axis=1)
    median_aneurysm = np.median(expr_aneurysm_filtered, axis=1)
    median_control = np.median(expr_control_filtered, axis=1)
    
    # Calculate log2 fold change with pseudocount
    log2_fc = np.log2((mean_aneurysm + 1e-6) / (mean_control + 1e-6))
    
    # Handle special cases
    log2_fc[np.isnan(log2_fc)] = 0
    log2_fc[~np.isfinite(log2_fc)] = 0
    
    print(f"\n  [OK] Tested {len(gene_ids_filtered):,} genes (after filtering)")
    
    # Create results DataFrame
    de_df = pd.DataFrame({
        'Gene_ID': gene_ids_filtered.astype(str),
        'Mean_Aneurysm': mean_aneurysm,
        'Mean_Control': mean_control,
        'Median_Aneurysm': median_aneurysm,
        'Median_Control': median_control,
        'Log2FC': log2_fc,
        'T_statistic': t_stats,
        'P_value': p_values
    })
    
    # Apply FDR correction
    print("\nApplying multiple testing correction...")
    valid_p = de_df['P_value'].dropna()
    if len(valid_p) == 0:
        print("[ERROR] No valid p-values!")
        return None
    
    # Create adjusted p-value column
    p_adj = np.full(len(de_df), np.nan)
    valid_mask = ~de_df['P_value'].isna()
    
    _, p_adj_valid, _, _ = multipletests(de_df.loc[valid_mask, 'P_value'], method='fdr_bh')
    p_adj[valid_mask] = p_adj_valid
    
    de_df['P_adj_BH'] = p_adj
    
    # Sort by p-value
    de_df = de_df.sort_values('P_value')
    
    return de_df

def summarize_results(de_df):
    """Summarize DE analysis results"""
    print("\n" + "=" * 80)
    print("RESULTS SUMMARY")
    print("=" * 80)
    
    print(f"\nTotal genes tested: {len(de_df):,}")
    print(f"Genes with valid p-values: {(~de_df['P_value'].isna()).sum():,}")
    
    if len(de_df) == 0:
        print("[WARNING] No results to summarize!")
        return
    
    # Summary statistics
    print(f"\nP-value statistics:")
    valid_p = de_df['P_value'].dropna()
    if len(valid_p) > 0:
        print(f"  Min p-value: {valid_p.min():.6f}")
        print(f"  Median p-value: {valid_p.median():.6f}")
        print(f"  Genes with p < 0.05: {(valid_p < 0.05).sum():,}")
        print(f"  Genes with p < 0.01: {(valid_p < 0.01).sum():,}")
    
    valid_padj = de_df['P_adj_BH'].dropna()
    if len(valid_padj) > 0:
        print(f"\nAdjusted p-value (FDR) statistics:")
        print(f"  Min adjusted p-value: {valid_padj.min():.6f}")
        print(f"  Genes with FDR < 0.05: {(valid_padj < 0.05).sum():,}")
        print(f"  Genes with FDR < 0.1: {(valid_padj < 0.1).sum():,}")
        print(f"  Genes with FDR < 0.2: {(valid_padj < 0.2).sum():,}")
    
    # Fold change statistics
    print(f"\nFold change statistics:")
    abs_fc = np.abs(de_df['Log2FC'])
    print(f"  Genes with |Log2FC| > 0.585 (1.5-fold): {(abs_fc > 0.585).sum():,}")
    print(f"  Genes with |Log2FC| > 1.0 (2-fold): {(abs_fc > 1.0).sum():,}")
    print(f"  Genes with |Log2FC| > 1.585 (3-fold): {(abs_fc > 1.585).sum():,}")
    
    # Combined criteria
    print(f"\nGenes passing combined criteria:")
    for fdr_thresh in [0.05, 0.1, 0.2]:
        for fc_thresh in [0.585, 1.0]:
            count = ((de_df['P_adj_BH'] < fdr_thresh) & (abs_fc > fc_thresh)).sum()
            print(f"  FDR < {fdr_thresh} AND |Log2FC| > {fc_thresh}: {count:,}")
    
    # Top genes
    print(f"\nTop 20 genes by p-value:")
    top_20 = de_df.head(20)[['Gene_ID', 'Log2FC', 'P_value', 'P_adj_BH']]
    print(top_20.to_string(index=False))
    
    # Top upregulated and downregulated
    print(f"\nTop 10 upregulated genes (Log2FC > 0):")
    upregulated = de_df[de_df['Log2FC'] > 0].head(10)[['Gene_ID', 'Log2FC', 'P_value', 'P_adj_BH']]
    print(upregulated.to_string(index=False))
    
    print(f"\nTop 10 downregulated genes (Log2FC < 0):")
    downregulated = de_df[de_df['Log2FC'] < 0].head(10)[['Gene_ID', 'Log2FC', 'P_value', 'P_adj_BH']]
    print(downregulated.to_string(index=False))

def save_results(de_df, output_dir="results"):
    """Save results to files"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Save full results
    output_file = os.path.join(output_dir, "gse159610_full_de_results.tsv")
    de_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved full results to {output_file}")
    
    # Save significant genes at different thresholds
    for fdr_thresh in [0.05, 0.1, 0.2]:
        sig_genes = de_df[de_df['P_adj_BH'] < fdr_thresh].copy()
        if len(sig_genes) > 0:
            output_file = os.path.join(output_dir, f"gse159610_sig_genes_fdr_{fdr_thresh}.tsv")
            sig_genes.to_csv(output_file, sep='\t', index=False)
            print(f"[OK] Saved {len(sig_genes):,} significant genes (FDR < {fdr_thresh}) to {output_file}")
    
    # Save genes with substantial fold changes even if not FDR-significant
    substantial_fc = de_df[np.abs(de_df['Log2FC']) > 0.585].copy()
    if len(substantial_fc) > 0:
        output_file = os.path.join(output_dir, "gse159610_substantial_fc_genes.tsv")
        substantial_fc = substantial_fc.sort_values('P_value')
        substantial_fc.to_csv(output_file, sep='\t', index=False)
        print(f"[OK] Saved {len(substantial_fc):,} genes with |Log2FC| > 0.585 to {output_file}")

def main():
    """Main analysis pipeline"""
    # Load data
    expr_df, aneurysm_samples, control_samples = load_data()
    
    # Perform DE analysis
    de_df = perform_de_analysis(expr_df, aneurysm_samples, control_samples)
    
    if de_df is None or len(de_df) == 0:
        print("[ERROR] No results generated!")
        return
    
    # Summarize results
    summarize_results(de_df)
    
    # Save results
    save_results(de_df)
    
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE!")
    print("=" * 80)

if __name__ == "__main__":
    main()

