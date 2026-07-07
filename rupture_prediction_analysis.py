"""
Rupture Prediction Analysis using GSE122897
Identifies top rupture-associated genes from ruptured vs unruptured comparison
"""
import pandas as pd
import numpy as np
import gzip
import os
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

def load_gse122897():
    """Load GSE122897 RNA-seq data"""
    
    print("=" * 80)
    print("LOADING GSE122897 (Ruptured vs Unruptured RNA-seq)")
    print("=" * 80)
    
    data_dir = "data/sc_rnaseq/GSE122897"
    
    # Try to load the raw counts file
    raw_counts_file = os.path.join(data_dir, "GSE122897_raw_counts_GRCh38.p13_NCBI.tsv.gz")
    
    if not os.path.exists(raw_counts_file):
        print(f"⚠️  Raw counts file not found: {raw_counts_file}")
        # Try alternative file
        alt_file = os.path.join(data_dir, "GSE122897_readCounts_raw.txt.gz")
        if os.path.exists(alt_file):
            raw_counts_file = alt_file
            print(f"✓ Using alternative file: {alt_file}")
        else:
            print("❌ No raw counts file found!")
            return None, None
    
    print(f"\nLoading: {raw_counts_file}")
    
    # Load expression matrix
    print("Reading expression matrix...")
    expr_df = pd.read_csv(raw_counts_file, sep='\t', compression='gzip', index_col=0)
    
    print(f"✓ Loaded expression matrix:")
    print(f"  Genes: {expr_df.shape[0]:,}")
    print(f"  Samples: {expr_df.shape[1]}")
    print(f"  Memory usage: {expr_df.memory_usage(deep=True).sum() / 1024**2:.2f} MB")
    
    # Load metadata if available
    metadata_file = os.path.join(data_dir, "GSE122897_series_matrix.txt.gz")
    metadata = None
    
    if os.path.exists(metadata_file):
        print(f"\nLoading metadata: {metadata_file}")
        try:
            with gzip.open(metadata_file, 'rt', encoding='utf-8') as f:
                lines = f.readlines()
                # Parse series matrix format
                metadata_dict = {}
                for line in lines:
                    if line.startswith('!Sample_'):
                        parts = line.strip().split('\t')
                        if len(parts) > 1:
                            key = parts[0].replace('!Sample_', '')
                            values = parts[1:]
                            metadata_dict[key] = values
                
                if metadata_dict:
                    metadata = pd.DataFrame(metadata_dict)
                    print(f"✓ Loaded metadata: {metadata.shape}")
        except Exception as e:
            print(f"⚠️  Could not parse metadata: {e}")
    
    # If no metadata file, try to infer from sample names
    if metadata is None:
        print("\n⚠️  No metadata file found. Attempting to infer sample groups from filenames...")
        sample_names = expr_df.columns.tolist()
        
        # Try to identify ruptured vs unruptured from sample names
        # Common patterns: "Ruptured", "Unruptured", "R", "U", "Rupt", "Unrupt"
        ruptured_samples = []
        unruptured_samples = []
        
        for sample in sample_names:
            sample_lower = sample.lower()
            if any(x in sample_lower for x in ['rupt', 'r_', '_r']):
                ruptured_samples.append(sample)
            elif any(x in sample_lower for x in ['unrupt', 'u_', '_u', 'control']):
                unruptured_samples.append(sample)
        
        if ruptured_samples and unruptured_samples:
            print(f"✓ Inferred sample groups:")
            print(f"  Ruptured: {len(ruptured_samples)} samples")
            print(f"  Unruptured: {len(unruptured_samples)} samples")
            
            metadata = pd.DataFrame({
                'Sample': sample_names,
                'Group': ['Ruptured' if s in ruptured_samples else 'Unruptured' if s in unruptured_samples else 'Unknown' 
                         for s in sample_names]
            })
        else:
            print("⚠️  Could not infer sample groups. Will need manual annotation.")
            metadata = pd.DataFrame({
                'Sample': sample_names,
                'Group': ['Unknown'] * len(sample_names)
            })
    
    return expr_df, metadata

def quality_control(expr_df, metadata):
    """Perform quality control on RNA-seq data"""
    
    print("\n" + "=" * 80)
    print("QUALITY CONTROL")
    print("=" * 80)
    
    # Calculate QC metrics
    print("\nCalculating QC metrics...")
    
    # Total counts per sample
    total_counts = expr_df.sum(axis=0)
    
    # Number of detected genes (counts > 0)
    detected_genes = (expr_df > 0).sum(axis=0)
    
    # Percentage of zeros
    pct_zeros = (expr_df == 0).sum(axis=0) / len(expr_df) * 100
    
    qc_metrics = pd.DataFrame({
        'Sample': expr_df.columns,
        'Total_Counts': total_counts.values,
        'Detected_Genes': detected_genes.values,
        'Pct_Zeros': pct_zeros.values
    })
    
    # Merge with metadata
    if 'Sample' in metadata.columns:
        qc_metrics = qc_metrics.merge(metadata, on='Sample', how='left')
    
    print("\nQC Summary Statistics:")
    print(f"  Total counts per sample:")
    print(f"    Mean: {qc_metrics['Total_Counts'].mean():,.0f}")
    print(f"    Median: {qc_metrics['Total_Counts'].median():,.0f}")
    print(f"    Range: {qc_metrics['Total_Counts'].min():,.0f} - {qc_metrics['Total_Counts'].max():,.0f}")
    
    print(f"\n  Detected genes per sample:")
    print(f"    Mean: {qc_metrics['Detected_Genes'].mean():,.0f}")
    print(f"    Median: {qc_metrics['Detected_Genes'].median():,.0f}")
    print(f"    Range: {qc_metrics['Detected_Genes'].min():,.0f} - {qc_metrics['Detected_Genes'].max():,.0f}")
    
    print(f"\n  Percentage of zeros:")
    print(f"    Mean: {qc_metrics['Pct_Zeros'].mean():.2f}%")
    print(f"    Median: {qc_metrics['Pct_Zeros'].median():.2f}%")
    
    # Filter low-quality samples
    # Remove samples with very low total counts or very few detected genes
    min_total_counts = qc_metrics['Total_Counts'].quantile(0.1)  # Bottom 10%
    min_detected_genes = qc_metrics['Detected_Genes'].quantile(0.1)
    
    good_samples = qc_metrics[
        (qc_metrics['Total_Counts'] >= min_total_counts) &
        (qc_metrics['Detected_Genes'] >= min_detected_genes)
    ]['Sample'].tolist()
    
    print(f"\n✓ Quality filtering:")
    print(f"  Samples before QC: {len(expr_df.columns)}")
    print(f"  Samples after QC: {len(good_samples)}")
    print(f"  Removed: {len(expr_df.columns) - len(good_samples)} samples")
    
    # Filter expression matrix
    expr_df_qc = expr_df[good_samples]
    
    # Filter lowly expressed genes
    # Keep genes with at least 10 counts in at least 3 samples
    min_counts = 10
    min_samples = 3
    gene_counts = (expr_df_qc >= min_counts).sum(axis=1)
    expressed_genes = gene_counts >= min_samples
    
    print(f"\n  Genes before QC: {len(expr_df)}")
    print(f"  Genes after QC: {expressed_genes.sum():,}")
    print(f"  Removed: {(~expressed_genes).sum():,} lowly expressed genes")
    
    expr_df_qc = expr_df_qc[expressed_genes]
    
    return expr_df_qc, qc_metrics[qc_metrics['Sample'].isin(good_samples)]

def normalize_data(expr_df):
    """Normalize RNA-seq data using TPM or DESeq2-style normalization"""
    
    print("\n" + "=" * 80)
    print("NORMALIZATION")
    print("=" * 80)
    
    # Check if TPM file exists (pre-normalized)
    data_dir = "data/sc_rnaseq/GSE122897"
    tpm_file = os.path.join(data_dir, "GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz")
    
    if os.path.exists(tpm_file):
        print(f"\n✓ Found pre-normalized TPM file: {tpm_file}")
        print("Loading TPM-normalized data...")
        tpm_df = pd.read_csv(tpm_file, sep='\t', compression='gzip', index_col=0)
        
        # Align with QC'd expression matrix
        common_genes = expr_df.index.intersection(tpm_df.index)
        common_samples = expr_df.columns.intersection(tpm_df.columns)
        
        tpm_df_aligned = tpm_df.loc[common_genes, common_samples]
        
        print(f"✓ Loaded TPM data: {tpm_df_aligned.shape}")
        return tpm_df_aligned, 'TPM'
    
    # Otherwise, perform normalization
    print("\nPerforming normalization...")
    
    # Method 1: Log2(CPM + 1) normalization
    print("  Method: Log2(CPM + 1)")
    
    # Calculate CPM (Counts Per Million)
    total_counts = expr_df.sum(axis=0)
    cpm = expr_df.div(total_counts, axis=1) * 1e6
    
    # Log2 transform
    log2_cpm = np.log2(cpm + 1)
    
    print(f"✓ Normalized data: {log2_cpm.shape}")
    return log2_cpm, 'Log2_CPM'

def differential_expression(expr_norm, metadata):
    """Perform differential expression analysis: Ruptured vs Unruptured"""
    
    print("\n" + "=" * 80)
    print("DIFFERENTIAL EXPRESSION ANALYSIS")
    print("=" * 80)
    
    # Get sample groups
    if 'Group' in metadata.columns:
        ruptured_samples = metadata[metadata['Group'] == 'Ruptured']['Sample'].tolist()
        unruptured_samples = metadata[metadata['Group'] == 'Unruptured']['Sample'].tolist()
    else:
        print("⚠️  No group information in metadata. Using all samples as one group.")
        return None
    
    if not ruptured_samples or not unruptured_samples:
        print("⚠️  Cannot perform DE analysis: missing group information")
        print(f"  Ruptured samples: {len(ruptured_samples)}")
        print(f"  Unruptured samples: {len(unruptured_samples)}")
        return None
    
    print(f"\nSample groups:")
    print(f"  Ruptured: {len(ruptured_samples)} samples")
    print(f"  Unruptured: {len(unruptured_samples)} samples")
    
    # Align samples
    ruptured_samples = [s for s in ruptured_samples if s in expr_norm.columns]
    unruptured_samples = [s for s in unruptured_samples if s in expr_norm.columns]
    
    if not ruptured_samples or not unruptured_samples:
        print("⚠️  No matching samples found in expression matrix")
        return None
    
    print(f"\nPerforming differential expression analysis...")
    print(f"  Testing {len(expr_norm):,} genes")
    
    # Perform t-test for each gene
    results = []
    
    for gene in expr_norm.index:
        ruptured_expr = expr_norm.loc[gene, ruptured_samples].values
        unruptured_expr = expr_norm.loc[gene, unruptured_samples].values
        
        # Remove NaN values
        ruptured_expr = ruptured_expr[~np.isnan(ruptured_expr)]
        unruptured_expr = unruptured_expr[~np.isnan(unruptured_expr)]
        
        if len(ruptured_expr) < 2 or len(unruptured_expr) < 2:
            continue
        
        # T-test
        t_stat, p_value = stats.ttest_ind(ruptured_expr, unruptured_expr)
        
        # Calculate fold change
        mean_ruptured = np.mean(ruptured_expr)
        mean_unruptured = np.mean(unruptured_expr)
        
        if mean_unruptured > 0:
            log2_fc = np.log2(mean_ruptured / mean_unruptured)
        else:
            log2_fc = np.nan
        
        results.append({
            'Gene': gene,
            'Mean_Ruptured': mean_ruptured,
            'Mean_Unruptured': mean_unruptured,
            'Log2FC': log2_fc,
            'T_statistic': t_stat,
            'P_value': p_value
        })
    
    de_results = pd.DataFrame(results)
    
    # Multiple testing correction
    de_results['P_adj_BH'] = multipletests(de_results['P_value'], method='fdr_bh')[1]
    de_results['P_adj_Bonferroni'] = multipletests(de_results['P_value'], method='bonferroni')[1]
    
    # Sort by significance
    de_results = de_results.sort_values('P_value')
    
    print(f"\n✓ Differential expression analysis complete")
    print(f"\nResults summary:")
    print(f"  Total genes tested: {len(de_results):,}")
    print(f"  Significant (P < 0.05): {(de_results['P_value'] < 0.05).sum():,}")
    print(f"  Significant (P_adj < 0.05): {(de_results['P_adj_BH'] < 0.05).sum():,}")
    print(f"  Highly significant (P_adj < 0.01): {(de_results['P_adj_BH'] < 0.01).sum():,}")
    
    # Top genes
    print(f"\nTop 20 rupture-associated genes:")
    top_genes = de_results.head(20)
    for idx, row in top_genes.iterrows():
        direction = "↑" if row['Log2FC'] > 0 else "↓"
        print(f"  {row['Gene']:20s} {direction} Log2FC={row['Log2FC']:7.3f}  P={row['P_value']:.2e}  P_adj={row['P_adj_BH']:.2e}")
    
    return de_results

def identify_top_rupture_genes(de_results, top_n=50):
    """Identify top rupture-associated genes"""
    
    print("\n" + "=" * 80)
    print(f"TOP {top_n} RUPTURE-ASSOCIATED GENES")
    print("=" * 80)
    
    # Filter significant genes
    sig_genes = de_results[de_results['P_adj_BH'] < 0.05].copy()
    
    # Sort by absolute log2FC and significance
    sig_genes['Abs_Log2FC'] = np.abs(sig_genes['Log2FC'])
    sig_genes = sig_genes.sort_values(['P_adj_BH', 'Abs_Log2FC'], ascending=[True, False])
    
    top_genes = sig_genes.head(top_n)
    
    # Separate upregulated and downregulated
    upregulated = top_genes[top_genes['Log2FC'] > 0].copy()
    downregulated = top_genes[top_genes['Log2FC'] < 0].copy()
    
    print(f"\nUpregulated in ruptured aneurysms ({len(upregulated)} genes):")
    for idx, row in upregulated.iterrows():
        print(f"  {row['Gene']:20s} Log2FC={row['Log2FC']:7.3f}  P_adj={row['P_adj_BH']:.2e}")
    
    print(f"\nDownregulated in ruptured aneurysms ({len(downregulated)} genes):")
    for idx, row in downregulated.iterrows():
        print(f"  {row['Gene']:20s} Log2FC={row['Log2FC']:7.3f}  P_adj={row['P_adj_BH']:.2e}")
    
    # Save results
    output_file = "results/gse122897_top_rupture_genes.tsv"
    os.makedirs("results", exist_ok=True)
    top_genes.to_csv(output_file, sep='\t', index=False)
    print(f"\n✓ Saved top {len(top_genes)} genes to {output_file}")
    
    return top_genes

def main():
    """Main analysis pipeline"""
    
    # Create results directory
    os.makedirs("results", exist_ok=True)
    
    # Step 1: Load data
    expr_df, metadata = load_gse122897()
    if expr_df is None:
        print("❌ Failed to load data. Exiting.")
        return
    
    # Step 2: Quality control
    expr_df_qc, qc_metrics = quality_control(expr_df, metadata)
    
    # Step 3: Normalization
    expr_norm, norm_method = normalize_data(expr_df_qc)
    
    # Step 4: Differential expression
    de_results = differential_expression(expr_norm, qc_metrics)
    
    if de_results is not None:
        # Step 5: Identify top genes
        top_genes = identify_top_rupture_genes(de_results, top_n=50)
        
        # Save full DE results
        output_file = "results/gse122897_differential_expression.tsv"
        de_results.to_csv(output_file, sep='\t', index=False)
        print(f"\n✓ Saved full DE results to {output_file}")
        
        print("\n" + "=" * 80)
        print("ANALYSIS COMPLETE!")
        print("=" * 80)
    else:
        print("\n⚠️  Could not complete differential expression analysis")

if __name__ == "__main__":
    main()






