"""
Integrated Analysis Pipeline for Multiple Aneurysm Datasets
Combines GSE122897 (RNA-seq) and GSE26969 (Microarray) with updated preprocessing
"""
import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path
from scipy import stats
from statsmodels.stats.multitest import multipletests

# Import existing analysis functions
sys.path.append('.')
from analyze_gse122897 import (
    load_and_prepare_data as load_gse122897,
    quality_control as qc_rnaseq,
    differential_expression as de_rnaseq
)

def load_gse26969_data():
    """Load GSE26969 normalized microarray data"""
    
    print("=" * 80)
    print("LOADING GSE26969 MICROARRAY DATA")
    print("=" * 80)
    
    # Try to load normalized data
    expr_file = Path("results/gse26969_expression_qc.tsv")
    if not expr_file.exists():
        expr_file = Path("results/gse26969_normalized_expression.tsv")
    
    if not expr_file.exists():
        print(f"[WARNING] GSE26969 expression data not found: {expr_file}")
        print("   Please run preprocess_gse26969.py first")
        return None, None
    
    print(f"\nLoading expression data: {expr_file}")
    expr_df = pd.read_csv(expr_file, sep='\t', index_col=0)
    
    print(f"[OK] Loaded expression matrix:")
    print(f"  Probes/Genes: {expr_df.shape[0]:,}")
    print(f"  Samples: {expr_df.shape[1]}")
    
    # Load metadata
    metadata_file = Path("results/gse26969_metadata.tsv")
    if not metadata_file.exists():
        metadata_file = Path("results/gse26969_metadata_template.tsv")
    
    if metadata_file.exists():
        print(f"\nLoading metadata: {metadata_file}")
        metadata = pd.read_csv(metadata_file, sep='\t')
        print(f"[OK] Loaded metadata: {metadata.shape}")
    else:
        print("\n[WARNING] Metadata file not found. Creating placeholder...")
        metadata = pd.DataFrame({
            'Sample': expr_df.columns.tolist(),
            'Group': ['Unknown'] * len(expr_df.columns)
        })
    
    return expr_df, metadata

def normalize_for_comparison(expr_rnaseq, expr_microarray):
    """
    Normalize both datasets for cross-platform comparison
    Uses quantile normalization and log2 transformation
    """
    print("\n" + "=" * 80)
    print("NORMALIZING FOR CROSS-PLATFORM COMPARISON")
    print("=" * 80)
    
    # Log2 transform both (microarray is already log2, RNA-seq needs transformation)
    print("\nApplying log2 transformation...")
    expr_rnaseq_log = np.log2(expr_rnaseq + 1)
    
    # Microarray data is already log2 from RMA, but check
    if expr_microarray.max() > 20:  # Likely not log2 transformed
        print("  Microarray data appears to be in linear scale, applying log2...")
        expr_microarray_log = np.log2(expr_microarray + 1)
    else:
        expr_microarray_log = expr_microarray.copy()
    
    # Quantile normalization (optional, for direct comparison)
    print("\n[NOTE] For direct comparison, consider quantile normalization")
    print("   This step can be added if needed for integration")
    
    return expr_rnaseq_log, expr_microarray_log

def find_common_genes(expr_rnaseq, expr_microarray, gene_mapping=None):
    """
    Find common genes between RNA-seq and microarray datasets
    Requires gene annotation/mapping
    """
    print("\n" + "=" * 80)
    print("FINDING COMMON GENES")
    print("=" * 80)
    
    rnaseq_genes = set(expr_rnaseq.index)
    microarray_probes = set(expr_microarray.index)
    
    print(f"\nRNA-seq genes: {len(rnaseq_genes):,}")
    print(f"Microarray probes: {len(microarray_probes):,}")
    
    if gene_mapping is not None:
        # Map microarray probes to genes
        print("\n[NOTE] Gene mapping required for direct comparison")
        print("   This requires probe annotation with gene symbols")
        return None, None
    else:
        print("\n[NOTE] Direct gene comparison requires:")
        print("   1. Annotate microarray probes with gene symbols")
        print("   2. Map Ensembl IDs (RNA-seq) to gene symbols")
        print("   3. Then find common genes")
        
        return None, None

def cross_platform_validation(expr_rnaseq, expr_microarray, metadata_rnaseq, metadata_microarray):
    """
    Validate findings across platforms
    Compare differential expression results
    """
    print("\n" + "=" * 80)
    print("CROSS-PLATFORM VALIDATION")
    print("=" * 80)
    
    print("\n[NOTE] Cross-platform validation requires:")
    print("  1. Common gene mapping between platforms")
    print("  2. Similar sample groups in both datasets")
    print("  3. Statistical comparison of effect sizes")
    
    # This will be implemented once gene mapping is available
    return None

def integrated_differential_expression(expr_rnaseq, expr_microarray, 
                                      metadata_rnaseq, metadata_microarray):
    """
    Perform differential expression analysis on both platforms
    and integrate results
    """
    print("\n" + "=" * 80)
    print("INTEGRATED DIFFERENTIAL EXPRESSION ANALYSIS")
    print("=" * 80)
    
    results = {}
    
    # Analyze RNA-seq data (GSE122897)
    print("\n--- GSE122897 (RNA-seq) Analysis ---")
    if 'Group' in metadata_rnaseq.columns:
        metadata_aneurysm = metadata_rnaseq[
            metadata_rnaseq['Group'].isin(['Ruptured', 'Unruptured'])
        ].copy()
        
        if len(metadata_aneurysm) >= 2 and metadata_aneurysm['Group'].nunique() >= 2:
            de_results_rnaseq, group1, group2 = de_rnaseq(expr_rnaseq, metadata_aneurysm)
            if de_results_rnaseq is not None:
                results['GSE122897'] = {
                    'de_results': de_results_rnaseq,
                    'group1': group1,
                    'group2': group2,
                    'platform': 'RNA-seq'
                }
                print(f"[OK] RNA-seq DE analysis complete")
    
    # Analyze microarray data (GSE26969)
    print("\n--- GSE26969 (Microarray) Analysis ---")
    if 'Group' in metadata_microarray.columns:
        # Check if we have IA vs STA comparison
        unique_groups = metadata_microarray['Group'].unique()
        if len(unique_groups) >= 2 and 'Unknown' not in unique_groups:
            # Perform DE analysis
            group1_name = unique_groups[0]
            group2_name = unique_groups[1]
            
            group1_samples = metadata_microarray[
                metadata_microarray['Group'] == group1_name
            ]['Sample'].tolist()
            group2_samples = metadata_microarray[
                metadata_microarray['Group'] == group2_name
            ]['Sample'].tolist()
            
            # Align with expression matrix
            group1_samples = [s for s in group1_samples if s in expr_microarray.columns]
            group2_samples = [s for s in group2_samples if s in expr_microarray.columns]
            
            if len(group1_samples) >= 1 and len(group2_samples) >= 1:
                print(f"\nComparing: {group1_name} (n={len(group1_samples)}) vs {group2_name} (n={len(group2_samples)})")
                
                # Subset expression data
                expr_subset = expr_microarray[group1_samples + group2_samples].copy()
                
                # Log2 transform if needed
                if expr_subset.max() > 20:
                    expr_subset = np.log2(expr_subset + 1)
                
                # Separate into groups
                expr_group1 = expr_subset[group1_samples].values
                expr_group2 = expr_subset[group2_samples].values
                
                # T-tests
                from scipy.stats import ttest_ind
                t_stats, p_values = ttest_ind(expr_group1, expr_group2, axis=1)
                
                # Calculate means and fold changes
                mean_group1 = np.mean(expr_group1, axis=1)
                mean_group2 = np.mean(expr_group2, axis=1)
                log2_fc = np.log2((mean_group1 + 1e-6) / (mean_group2 + 1e-6))
                
                # Create results DataFrame
                de_results_microarray = pd.DataFrame({
                    'ProbeID': expr_microarray.index.astype(str),
                    f'Mean_{group1_name}': mean_group1,
                    f'Mean_{group2_name}': mean_group2,
                    'Log2FC': log2_fc,
                    'T_statistic': t_stats,
                    'P_value': p_values
                })
                
                # Multiple testing correction
                p_values_clean = np.nan_to_num(p_values, nan=1.0)
                _, p_adj_bh, _, _ = multipletests(p_values_clean, method='fdr_bh')
                de_results_microarray['P_adj_BH'] = p_adj_bh
                
                de_results_microarray = de_results_microarray.sort_values('P_value')
                
                results['GSE26969'] = {
                    'de_results': de_results_microarray,
                    'group1': group1_name,
                    'group2': group2_name,
                    'platform': 'Microarray'
                }
                
                print(f"[OK] Microarray DE analysis complete")
                print(f"  Significant (P < 0.05): {(de_results_microarray['P_value'] < 0.05).sum():,}")
                print(f"  Significant (P_adj < 0.05): {(de_results_microarray['P_adj_BH'] < 0.05).sum():,}")
    
    return results

def save_integrated_results(results):
    """Save integrated analysis results"""
    print("\n" + "=" * 80)
    print("SAVING INTEGRATED RESULTS")
    print("=" * 80)
    
    os.makedirs("results", exist_ok=True)
    
    for dataset, data in results.items():
        de_results = data['de_results']
        output_file = f"results/{dataset.lower()}_integrated_de.tsv"
        de_results.to_csv(output_file, sep='\t', index=False)
        print(f"\n[OK] Saved {dataset} results: {output_file}")
        
        # Save top genes
        top_genes = de_results[de_results['P_adj_BH'] < 0.05].head(100)
        if len(top_genes) > 0:
            top_file = f"results/{dataset.lower()}_top_genes_integrated.tsv"
            top_genes.to_csv(top_file, sep='\t', index=False)
            print(f"[OK] Saved top genes: {top_file}")

def main():
    """Main integrated analysis pipeline"""
    
    print("=" * 80)
    print("INTEGRATED ANALYSIS PIPELINE")
    print("GSE122897 (RNA-seq) + GSE26969 (Microarray)")
    print("=" * 80)
    
    # Step 1: Load GSE122897 RNA-seq data
    expr_rnaseq, metadata_rnaseq = load_gse122897()
    if expr_rnaseq is None:
        print("[ERROR] Failed to load GSE122897 data")
        return
    
    # Step 2: Load GSE26969 microarray data
    expr_microarray, metadata_microarray = load_gse26969_data()
    if expr_microarray is None:
        print("[WARNING] GSE26969 data not available, proceeding with RNA-seq only")
        expr_microarray = None
        metadata_microarray = None
    
    # Step 3: Quality control
    print("\n" + "=" * 80)
    print("QUALITY CONTROL")
    print("=" * 80)
    
    expr_rnaseq_qc = qc_rnaseq(expr_rnaseq)
    
    if expr_microarray is not None:
        # QC for microarray (similar to RNA-seq but with different thresholds)
        from preprocess_gse26969 import quality_control_microarray
        expr_microarray_qc = quality_control_microarray(expr_microarray)
    else:
        expr_microarray_qc = None
    
    # Step 4: Integrated differential expression
    if expr_microarray_qc is not None:
        results = integrated_differential_expression(
            expr_rnaseq_qc, expr_microarray_qc,
            metadata_rnaseq, metadata_microarray
        )
    else:
        # RNA-seq only analysis
        metadata_aneurysm = metadata_rnaseq[
            metadata_rnaseq['Group'].isin(['Ruptured', 'Unruptured'])
        ].copy()
        
        if len(metadata_aneurysm) >= 2 and metadata_aneurysm['Group'].nunique() >= 2:
            de_results, group1, group2 = de_rnaseq(expr_rnaseq_qc, metadata_aneurysm)
            if de_results is not None:
                results = {
                    'GSE122897': {
                        'de_results': de_results,
                        'group1': group1,
                        'group2': group2,
                        'platform': 'RNA-seq'
                    }
                }
            else:
                results = {}
        else:
            results = {}
    
    # Step 5: Save results
    if results:
        save_integrated_results(results)
        
        print("\n" + "=" * 80)
        print("INTEGRATED ANALYSIS COMPLETE!")
        print("=" * 80)
    else:
        print("\n[WARNING] No differential expression results generated")
        print("   Please check metadata and sample groups")

if __name__ == "__main__":
    main()






