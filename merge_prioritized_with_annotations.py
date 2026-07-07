"""
Merge Prioritized Genes with Annotations
Combines prioritized genes list with annotated gene information
"""
import pandas as pd
import numpy as np
import os

def merge_prioritized_with_annotations():
    """Merge prioritized genes with annotations"""
    
    print("=" * 80)
    print("MERGING PRIORITIZED GENES WITH ANNOTATIONS")
    print("=" * 80)
    
    # Load prioritized genes
    prioritized_file = "results/prioritized_genes_summary.tsv"
    if not os.path.exists(prioritized_file):
        print(f"[ERROR] Prioritized genes file not found: {prioritized_file}")
        return None
    
    prioritized_df = pd.read_csv(prioritized_file, sep='\t')
    print(f"\nLoaded prioritized genes: {len(prioritized_df)} genes")
    
    # Load annotated genes
    annotated_file = "results/gse122897_annotated_fdr_0.1.tsv"
    if not os.path.exists(annotated_file):
        print(f"[ERROR] Annotated genes file not found: {annotated_file}")
        return None
    
    annotated_df = pd.read_csv(annotated_file, sep='\t')
    print(f"Loaded annotated genes: {len(annotated_df)} genes")
    
    # Prepare for merging
    prioritized_df['Gene_ID_Str'] = prioritized_df['Gene_ID'].astype(str).str.rstrip('.0')
    annotated_df['Gene_ID_Str'] = annotated_df['Gene'].astype(str).str.rstrip('.0')
    
    # Merge on Gene_ID
    merged_df = prioritized_df.merge(
        annotated_df,
        left_on='Gene_ID_Str',
        right_on='Gene_ID_Str',
        how='left',
        suffixes=('_prioritized', '_annotated')
    )
    
    print(f"\nMerged results: {len(merged_df)} genes")
    print(f"  Genes with annotations: {merged_df['Gene'].notna().sum()}")
    print(f"  Genes without annotations: {merged_df['Gene'].isna().sum()}")
    
    # Reorder columns - prioritize key information
    priority_cols = [
        'Gene_ID', 'Gene_Symbol', 'Gene_Name', 'Entrez_ID', 'Ensembl_ID',
        'Gene_Type', 'Aliases', 'Log2FC_prioritized', 'P_adj_BH_prioritized',
        'Priority_Score', 'In_GWAS_Locus', 'GWAS_Association',
        'GO_Terms', 'Pathways', 'OMIM_Disease', 'Summary'
    ]
    
    # Select available columns
    available_cols = [col for col in priority_cols if col in merged_df.columns]
    
    # Add remaining columns
    remaining_cols = [col for col in merged_df.columns if col not in available_cols]
    final_cols = available_cols + remaining_cols
    
    merged_df = merged_df[final_cols]
    
    # Rename columns for clarity
    rename_dict = {
        'Log2FC_prioritized': 'Log2FC',
        'P_adj_BH_prioritized': 'P_adj_BH',
        'Log2FC_annotated': 'Log2FC_annotated',
        'P_adj_BH_annotated': 'P_adj_BH_annotated'
    }
    merged_df = merged_df.rename(columns=rename_dict)
    
    # Save merged results
    output_file = "results/prioritized_genes_annotated_merged.tsv"
    merged_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved merged results to {output_file}")
    
    # Print top annotated prioritized genes
    annotated_merged = merged_df[merged_df['Gene'].notna()].copy()
    annotated_merged = annotated_merged.sort_values('Priority_Score', ascending=False)
    
    print(f"\n\nTop 20 Prioritized Genes with Annotations:")
    print(f"{'Rank':<6} {'Symbol':<15} {'Name':<40} {'Log2FC':<10} {'P_adj_BH':<12} {'Priority':<12}")
    print("-" * 110)
    
    for idx, (i, row) in enumerate(annotated_merged.head(20).iterrows(), 1):
        symbol = row.get('Gene_Symbol', str(row['Gene_ID']).rstrip('.0'))
        name = row.get('Gene_Name', '')
        if pd.notna(name) and isinstance(name, str):
            if len(name) > 38:
                name = name[:35] + "..."
        else:
            name = ''
        log2fc = row.get('Log2FC', np.nan)
        p_adj = row.get('P_adj_BH', np.nan)
        priority = row.get('Priority_Score', np.nan)
        
        symbol_str = symbol[:13] if isinstance(symbol, str) and len(symbol) > 13 else symbol
        name_str = name[:38] if len(name) > 38 else name
        
        log2fc_str = f"{log2fc:.3f}" if pd.notna(log2fc) else "N/A"
        p_adj_str = f"{p_adj:.4f}" if pd.notna(p_adj) else "N/A"
        priority_str = f"{priority:.2f}" if pd.notna(priority) else "N/A"
        
        print(f"{idx:<6} {str(symbol_str):<15} {name_str:<40} {log2fc_str:>9}  {p_adj_str:>11}  {priority_str:>11}")
    
    # Summary statistics
    print(f"\n\nMerged Results Summary:")
    print(f"  Total prioritized genes: {len(prioritized_df)}")
    print(f"  Genes with annotations: {annotated_merged['Gene'].notna().sum()}")
    print(f"  Genes with symbols: {annotated_merged['Gene_Symbol'].notna().sum() if 'Gene_Symbol' in annotated_merged.columns else 0}")
    if 'GO_Terms' in annotated_merged.columns:
        print(f"  Genes with GO terms: {annotated_merged['GO_Terms'].notna().sum()}")
    if 'Pathways' in annotated_merged.columns:
        print(f"  Genes with pathways: {annotated_merged['Pathways'].notna().sum()}")
    
    return merged_df

def main():
    """Main function"""
    merged_df = merge_prioritized_with_annotations()
    
    print("\n" + "=" * 80)
    print("MERGE COMPLETE")
    print("=" * 80)
    
    if merged_df is not None:
        print(f"\n[NOTE] Merged results saved to results/prioritized_genes_annotated_merged.tsv")
        print("  Next steps:")
        print("  1. Review top annotated prioritized genes")
        print("  2. Use annotations for literature searches")
        print("  3. Identify functional studies for validation")

if __name__ == "__main__":
    main()

