"""
Re-run Blood Integration Analysis with Refined Metadata
Identifies biomarker candidates by comparing tissue DE genes with blood transcriptomics
"""
import pandas as pd
import numpy as np
import os
from scipy import stats
from scipy.stats import ttest_ind

def load_tissue_de_genes(fdr_threshold=0.1):
    """Load significant tissue DE genes"""
    de_file = "results/gse122897_differential_expression.tsv"
    if not os.path.exists(de_file):
        print(f"[ERROR] DE results file not found: {de_file}")
        return None
    
    df = pd.read_csv(de_file, sep='\t')
    sig_df = df[df['P_adj_BH'] < fdr_threshold].copy()
    
    print(f"Loaded {len(sig_df)} significant tissue DE genes (FDR < {fdr_threshold})")
    return sig_df

def load_blood_data():
    """Load blood RNA-seq expression data"""
    blood_file = "data/blood_rnaseq/GSE159610/GSE159610_raw_gene_counts_matrix.txt.gz"
    
    if not os.path.exists(blood_file):
        print(f"[ERROR] Blood transcriptomics file not found: {blood_file}")
        return None, None
    
    print(f"\nLoading blood transcriptomics data...")
    expr_df = pd.read_csv(blood_file, sep='\t', compression='gzip', index_col=0)
    print(f"  [OK] Expression matrix: {expr_df.shape[0]} genes, {expr_df.shape[1]} samples")
    
    # Load refined metadata
    metadata_file = "data/blood_rnaseq/GSE159610/sample_metadata_with_labels.tsv"
    if not os.path.exists(metadata_file):
        print(f"[ERROR] Refined metadata not found: {metadata_file}")
        return None, None
    
    metadata = pd.read_csv(metadata_file, sep='\t')
    print(f"  [OK] Metadata: {metadata.shape}")
    
    # Check sample groups
    if 'Status' in metadata.columns:
        status_counts = metadata['Status'].value_counts()
        print(f"  Sample groups: {status_counts.to_dict()}")
    
    return expr_df, metadata

def map_gene_ids_to_symbols(gene_ids):
    """Map gene IDs to symbols using annotated top genes"""
    # Load annotated genes if available
    annotated_file = "results/gse122897_top_genes_annotated.tsv"
    mapping = {}
    
    if os.path.exists(annotated_file):
        annotated_df = pd.read_csv(annotated_file, sep='\t')
        for idx, row in annotated_df.iterrows():
            gene_id = str(row['Gene']).rstrip('.0')
            symbol = row.get('Gene_Symbol', gene_id)
            mapping[gene_id] = symbol
    
    # For unmapped genes, use ID as symbol
    for gene_id in gene_ids:
        gene_id_str = str(gene_id).rstrip('.0')
        if gene_id_str not in mapping:
            mapping[gene_id_str] = gene_id_str
    
    return mapping

def identify_biomarker_candidates(tissue_de_df, blood_expr_df, blood_metadata, fdr_threshold=0.1):
    """
    Identify biomarker candidates by finding tissue DE genes that are also
    differentially expressed in blood between aneurysm and control samples
    """
    print(f"\n" + "=" * 80)
    print("IDENTIFYING BIOMARKER CANDIDATES")
    print("=" * 80)
    
    # Get gene ID mapping
    tissue_gene_ids = tissue_de_df['Gene'].astype(str).str.rstrip('.0').tolist()
    gene_mapping = map_gene_ids_to_symbols(tissue_gene_ids)
    
    # Get blood sample groups
    if 'Status' not in blood_metadata.columns:
        print("[ERROR] No 'Status' column in blood metadata")
        return None
    
    aneurysm_samples = blood_metadata[blood_metadata['Status'] == 'Aneurysm']['Sample'].tolist()
    control_samples = blood_metadata[blood_metadata['Status'] == 'Control']['Sample'].tolist()
    
    if not aneurysm_samples or not control_samples:
        print("[ERROR] Could not identify aneurysm and control samples")
        return None
    
    print(f"\nBlood sample groups:")
    print(f"  Aneurysm: {len(aneurysm_samples)} samples")
    print(f"  Control: {len(control_samples)} samples")
    
    # Find overlapping genes between tissue and blood
    blood_gene_ids = [str(gid).rstrip('.0') for gid in blood_expr_df.index]
    overlapping_genes = [gid for gid in tissue_gene_ids if gid in blood_gene_ids]
    
    print(f"\nOverlapping genes: {len(overlapping_genes)} out of {len(tissue_gene_ids)} tissue DE genes")
    
    # Test differential expression in blood for overlapping genes
    biomarker_candidates = []
    
    print(f"\nTesting differential expression in blood for {len(overlapping_genes)} genes...")
    
    for i, gene_id in enumerate(overlapping_genes[:100], 1):  # Limit to top 100 for speed
        if i % 20 == 0:
            print(f"  Progress: {i}/{min(len(overlapping_genes), 100)}")
        
        gene_id_str = str(gene_id).rstrip('.0')
        
        # Get expression values
        if gene_id_str not in blood_expr_df.index:
            continue
        
        aneurysm_expr = blood_expr_df.loc[gene_id_str, aneurysm_samples].values
        control_expr = blood_expr_df.loc[gene_id_str, control_samples].values
        
        # Filter out genes with too many zeros
        if np.sum(aneurysm_expr > 0) < 3 or np.sum(control_expr > 0) < 3:
            continue
        
        # Log transform
        aneurysm_expr_log = np.log2(aneurysm_expr + 1)
        control_expr_log = np.log2(control_expr + 1)
        
        # T-test
        try:
            t_stat, p_value = ttest_ind(aneurysm_expr_log, control_expr_log)
            
            # Calculate fold change
            mean_aneurysm = np.mean(aneurysm_expr)
            mean_control = np.mean(control_expr)
            log2fc_blood = np.log2((mean_aneurysm + 1e-6) / (mean_control + 1e-6))
            
            # Get tissue DE info
            tissue_row = tissue_de_df[tissue_de_df['Gene'].astype(str).str.rstrip('.0') == gene_id_str].iloc[0]
            log2fc_tissue = tissue_row['Log2FC']
            p_adj_tissue = tissue_row['P_adj_BH']
            
            # Calculate biomarker score (higher is better)
            # Consider both tissue and blood significance, and consistency of direction
            direction_consistent = (log2fc_tissue > 0 and log2fc_blood > 0) or (log2fc_tissue < 0 and log2fc_blood < 0)
            biomarker_score = abs(log2fc_tissue) * abs(log2fc_blood) * (-np.log10(max(p_value, 1e-10))) * (1 if direction_consistent else 0.5)
            
            biomarker_candidates.append({
                'Gene_ID': gene_id_str,
                'Gene_Symbol': gene_mapping.get(gene_id_str, gene_id_str),
                'Log2FC_Tissue': log2fc_tissue,
                'P_adj_Tissue': p_adj_tissue,
                'Log2FC_Blood': log2fc_blood,
                'P_value_Blood': p_value,
                'Mean_Aneurysm_Blood': mean_aneurysm,
                'Mean_Control_Blood': mean_control,
                'Direction_Consistent': direction_consistent,
                'Biomarker_Score': biomarker_score
            })
        except Exception as e:
            continue
    
    if not biomarker_candidates:
        print("[WARNING] No biomarker candidates found")
        return None
    
    candidates_df = pd.DataFrame(biomarker_candidates)
    candidates_df = candidates_df.sort_values('Biomarker_Score', ascending=False)
    
    print(f"\n[OK] Identified {len(candidates_df)} biomarker candidates")
    
    # Filter by significance
    significant_candidates = candidates_df[candidates_df['P_value_Blood'] < 0.05]
    print(f"  Significant in blood (P < 0.05): {len(significant_candidates)}")
    
    # Filter by direction consistency
    consistent_candidates = candidates_df[candidates_df['Direction_Consistent'] == True]
    print(f"  Direction consistent: {len(consistent_candidates)}")
    
    return candidates_df

def main():
    """Main analysis pipeline"""
    
    print("=" * 80)
    print("BLOOD INTEGRATION ANALYSIS WITH REFINED METADATA")
    print("=" * 80)
    
    # Create results directory
    os.makedirs("results", exist_ok=True)
    
    # Load tissue DE genes
    tissue_de_df = load_tissue_de_genes(fdr_threshold=0.1)
    if tissue_de_df is None:
        return
    
    # Load blood data
    blood_expr_df, blood_metadata = load_blood_data()
    if blood_expr_df is None or blood_metadata is None:
        return
    
    # Identify biomarker candidates
    candidates_df = identify_biomarker_candidates(tissue_de_df, blood_expr_df, blood_metadata, fdr_threshold=0.1)
    
    if candidates_df is not None:
        # Save results
        output_file = "results/blood_tissue_integration_refined.tsv"
        candidates_df.to_csv(output_file, sep='\t', index=False)
        print(f"\n[OK] Saved biomarker candidates to {output_file}")
        
        # Print top candidates
        print(f"\n\nTop 10 Biomarker Candidates:")
        print(f"{'Rank':<6} {'Symbol':<15} {'Tissue_Log2FC':<15} {'Blood_Log2FC':<15} {'Blood_P':<12} {'Score':<12}")
        print("-" * 80)
        
        for rank, (idx, row) in enumerate(candidates_df.head(10).iterrows(), 1):
            print(f"{rank:<6} {row['Gene_Symbol']:<15} {row['Log2FC_Tissue']:>14.3f} {row['Log2FC_Blood']:>14.3f} {row['P_value_Blood']:>11.4e} {row['Biomarker_Score']:>11.2f}")
        
        print("\n" + "=" * 80)
        print("BLOOD INTEGRATION ANALYSIS COMPLETE")
        print("=" * 80)
    
    return candidates_df

if __name__ == "__main__":
    main()






