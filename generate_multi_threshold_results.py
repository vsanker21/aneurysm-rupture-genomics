"""
Generate Results Tables at Multiple FDR Thresholds
Creates separate result files for different significance levels
"""
import pandas as pd
import numpy as np
import os

def load_de_results():
    """Load differential expression results"""
    de_file = "results/gse122897_differential_expression.tsv"
    if not os.path.exists(de_file):
        print(f"[ERROR] DE results file not found: {de_file}")
        return None
    
    df = pd.read_csv(de_file, sep='\t')
    return df

def generate_threshold_results(df, thresholds=[0.05, 0.1, 0.2]):
    """
    Generate separate result files for different FDR thresholds
    """
    print("=" * 80)
    print("GENERATING MULTI-THRESHOLD RESULTS")
    print("=" * 80)
    
    results_summary = []
    
    for threshold in thresholds:
        print(f"\nProcessing FDR < {threshold}...")
        
        # Get significant genes
        sig_genes = df[df['P_adj_BH'] < threshold].copy()
        
        if len(sig_genes) == 0:
            print(f"  No genes pass FDR < {threshold}")
            results_summary.append({
                'FDR_threshold': threshold,
                'Total_genes': 0,
                'Upregulated': 0,
                'Downregulated': 0,
                'Mean_Log2FC': np.nan,
                'Max_Log2FC': np.nan,
                'Min_Log2FC': np.nan
            })
            continue
        
        # Separate upregulated and downregulated
        upregulated = sig_genes[sig_genes['Log2FC'] > 0].copy()
        downregulated = sig_genes[sig_genes['Log2FC'] < 0].copy()
        
        # Sort by significance
        sig_genes = sig_genes.sort_values('P_adj_BH')
        upregulated = upregulated.sort_values('P_adj_BH')
        downregulated = downregulated.sort_values('P_adj_BH')
        
        # Save all significant genes
        output_file = f"results/gse122897_sig_genes_fdr_{threshold}.tsv"
        sig_genes.to_csv(output_file, sep='\t', index=False)
        print(f"  [OK] Saved {len(sig_genes)} significant genes to {output_file}")
        
        # Save upregulated genes
        up_file = f"results/gse122897_upregulated_fdr_{threshold}.tsv"
        upregulated.to_csv(up_file, sep='\t', index=False)
        print(f"  [OK] Saved {len(upregulated)} upregulated genes to {up_file}")
        
        # Save downregulated genes
        down_file = f"results/gse122897_downregulated_fdr_{threshold}.tsv"
        downregulated.to_csv(down_file, sep='\t', index=False)
        print(f"  [OK] Saved {len(downregulated)} downregulated genes to {down_file}")
        
        # Statistics
        results_summary.append({
            'FDR_threshold': threshold,
            'Total_genes': len(sig_genes),
            'Upregulated': len(upregulated),
            'Downregulated': len(downregulated),
            'Mean_Log2FC': sig_genes['Log2FC'].mean(),
            'Mean_|Log2FC|': sig_genes['Log2FC'].abs().mean(),
            'Max_Log2FC': sig_genes['Log2FC'].max(),
            'Min_Log2FC': sig_genes['Log2FC'].min(),
            'Median_P_adj': sig_genes['P_adj_BH'].median()
        })
        
        # Print summary
        print(f"  Summary:")
        print(f"    Total: {len(sig_genes)} genes")
        print(f"    Upregulated: {len(upregulated)}")
        print(f"    Downregulated: {len(downregulated)}")
        print(f"    Mean |Log2FC|: {sig_genes['Log2FC'].abs().mean():.3f}")
    
    # Save summary
    summary_df = pd.DataFrame(results_summary)
    summary_file = "results/multi_threshold_summary.tsv"
    summary_df.to_csv(summary_file, sep='\t', index=False)
    print(f"\n[OK] Saved summary to {summary_file}")
    
    print("\n" + "=" * 80)
    print("MULTI-THRESHOLD RESULTS SUMMARY")
    print("=" * 80)
    print(summary_df.to_string(index=False))
    
    return summary_df

def main():
    """Main pipeline"""
    # Load DE results
    df = load_de_results()
    if df is None:
        return
    
    # Generate results at multiple thresholds
    summary = generate_threshold_results(df, thresholds=[0.05, 0.1, 0.2])
    
    print("\n" + "=" * 80)
    print("MULTI-THRESHOLD RESULTS GENERATION COMPLETE")
    print("=" * 80)

if __name__ == "__main__":
    main()






