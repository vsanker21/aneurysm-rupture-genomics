"""
Blood Transcriptomics Integration Script
Compares tissue DE genes with blood transcriptomics signatures
"""
import pandas as pd
import numpy as np
import gzip
import os
from typing import Dict, List, Tuple, Optional

def load_de_results(fdr_threshold=0.1):
    """Load differential expression results from tissue analysis"""
    de_file = "results/gse122897_differential_expression.tsv"
    if not os.path.exists(de_file):
        print(f"[ERROR] DE results file not found: {de_file}")
        return None
    
    df = pd.read_csv(de_file, sep='\t')
    sig_df = df[df['P_adj_BH'] < fdr_threshold].copy()
    
    print(f"Loaded {len(sig_df)} significant tissue DE genes (FDR < {fdr_threshold})")
    return sig_df

def load_blood_transcriptomics_data():
    """
    Load blood RNA-seq data from GSE159610
    """
    blood_file = "data/blood_rnaseq/GSE159610/GSE159610_raw_gene_counts_matrix.txt.gz"
    
    if not os.path.exists(blood_file):
        print(f"[WARNING] Blood transcriptomics file not found: {blood_file}")
        return None, None
    
    print(f"\nLoading blood transcriptomics data: {blood_file}")
    
    try:
        # Load expression matrix
        expr_df = pd.read_csv(blood_file, sep='\t', compression='gzip', index_col=0)
        print(f"  [OK] Loaded expression matrix: {expr_df.shape[0]} genes, {expr_df.shape[1]} samples")
        
        # Try to load parsed metadata first, otherwise parse series matrix
        # Try refined metadata first, then original parsed, then parse from series matrix
        refined_metadata_file = "data/blood_rnaseq/GSE159610/sample_metadata_with_labels.tsv"
        parsed_metadata_file = "data/blood_rnaseq/GSE159610/sample_metadata.tsv"
        metadata_file = "data/blood_rnaseq/GSE159610/GSE159610_series_matrix.txt.gz"
        metadata = None
        
        # Try refined metadata first (from refine_blood_metadata.py)
        if os.path.exists(refined_metadata_file):
            try:
                metadata = pd.read_csv(refined_metadata_file, sep='\t')
                print(f"  [OK] Loaded refined metadata: {metadata.shape}")
            except Exception as e:
                print(f"  [WARNING] Could not load refined metadata: {e}")
        
        # Try parsed metadata if refined not available
        if metadata is None and os.path.exists(parsed_metadata_file):
            try:
                metadata = pd.read_csv(parsed_metadata_file, sep='\t')
                print(f"  [OK] Loaded parsed metadata: {metadata.shape}")
            except Exception as e:
                print(f"  [WARNING] Could not load parsed metadata: {e}")
        
        # If no parsed metadata, try parsing series matrix
        if metadata is None and os.path.exists(metadata_file):
            try:
                # Parse GEO series matrix format
                with gzip.open(metadata_file, 'rt', encoding='utf-8', errors='ignore') as f:
                    lines = f.readlines()
                    
                    metadata_dict = {}
                    sample_names = []
                    for line in lines:
                        line = line.strip()
                        if not line or not line.startswith('!'):
                            continue
                        
                        if line.startswith('!Sample_'):
                            parts = line.split('\t', 1)
                            if len(parts) >= 2:
                                key = parts[0].replace('!Sample_', '').strip()
                                values = parts[1].split('\t')
                                # Clean values
                                values = [v.strip('"').strip("'").strip() for v in values]
                                metadata_dict[key] = values
                                
                                # Get sample names
                                if key == 'geo_accession':
                                    sample_names = values
                                elif key == 'title' and not sample_names:
                                    sample_names = values
                    
                    # Convert to DataFrame
                    if metadata_dict and sample_names:
                        metadata = pd.DataFrame({'Sample': sample_names})
                        for key, values in metadata_dict.items():
                            if len(values) == len(sample_names):
                                metadata[key] = values
                        print(f"  [OK] Parsed metadata from series matrix: {metadata.shape}")
            except Exception as e:
                print(f"  [WARNING] Could not parse series matrix: {e}")
        
        # Create placeholder metadata if needed
        if metadata is None:
            samples = expr_df.columns.tolist()
            metadata = pd.DataFrame({
                'Sample': samples,
                'Status': ['Unknown'] * len(samples)
            })
            print("  [NOTE] Using placeholder metadata")
        
        return expr_df, metadata
        
    except Exception as e:
        print(f"[ERROR] Failed to load blood data: {e}")
        return None, None

def analyze_blood_signature(expr_df: pd.DataFrame, metadata: pd.DataFrame):
    """
    Analyze blood transcriptomics to identify aneurysm-associated genes
    """
    print("\nAnalyzing blood transcriptomics signature...")
    
    # Identify samples by status if available
    # Check multiple possible column names
    status_col = None
    for col in ['Status', 'status', 'Group', 'group', 'Sample_characteristics_ch1', 
                'characteristics_ch1', 'Sample_title', 'title']:
        if col in metadata.columns:
            status_col = col
            break
    
    if status_col:
        # Clean up status values (may have quotes or other formatting)
        unique_statuses = metadata[status_col].unique()
        # Remove quotes and clean
        unique_statuses = [str(s).strip('"').strip("'").strip() for s in unique_statuses]
        unique_statuses = [s for s in unique_statuses if s and s != 'nan']
        unique_statuses = list(set(unique_statuses))
        
        print(f"  Sample groups: {unique_statuses}")
        
        if len(unique_statuses) >= 2:
            # Perform differential expression in blood
            # Match cleaned values
            cleaned_metadata = metadata.copy()
            cleaned_metadata[status_col] = cleaned_metadata[status_col].astype(str).str.strip('"').str.strip("'").str.strip()
            
            group1 = cleaned_metadata[cleaned_metadata[status_col] == unique_statuses[0]]['Sample'].tolist()
            group2 = cleaned_metadata[cleaned_metadata[status_col] == unique_statuses[1]]['Sample'].tolist()
            
            # Filter to samples present in expression data
            group1 = [s for s in group1 if s in expr_df.columns]
            group2 = [s for s in group2 if s in expr_df.columns]
            
            if len(group1) > 0 and len(group2) > 0:
                print(f"  Comparing {len(group1)} vs {len(group2)} samples")
                
                # Simple DE analysis
                from scipy.stats import ttest_ind
                
                blood_de_results = []
                
                for gene in expr_df.index[:5000]:  # Limit for performance
                    expr_group1 = expr_df.loc[gene, group1].values
                    expr_group2 = expr_df.loc[gene, group2].values
                    
                    # Log transform
                    expr_group1_log = np.log2(expr_group1 + 1)
                    expr_group2_log = np.log2(expr_group2 + 1)
                    
                    # T-test
                    t_stat, p_val = ttest_ind(expr_group1_log, expr_group2_log)
                    
                    # Fold change
                    mean1 = np.mean(expr_group1)
                    mean2 = np.mean(expr_group2)
                    log2fc = np.log2((mean1 + 1e-6) / (mean2 + 1e-6))
                    
                    blood_de_results.append({
                        'Gene': str(gene),
                        'Log2FC_Blood': log2fc,
                        'P_value_Blood': p_val,
                        'Mean_Group1': mean1,
                        'Mean_Group2': mean2
                    })
                
                blood_de_df = pd.DataFrame(blood_de_results)
                
                # Multiple testing correction
                from statsmodels.stats.multitest import multipletests
                _, p_adj, _, _ = multipletests(blood_de_df['P_value_Blood'], method='fdr_bh')
                blood_de_df['P_adj_BH_Blood'] = p_adj
                
                blood_de_df = blood_de_df.sort_values('P_value_Blood')
                
                print(f"  [OK] Analyzed {len(blood_de_df)} genes")
                return blood_de_df
    
    print("  [NOTE] Could not perform blood DE analysis (metadata or samples insufficient)")
    return None

def integrate_tissue_blood_signatures(tissue_de_df: pd.DataFrame, blood_de_df: pd.DataFrame = None):
    """
    Integrate tissue DE genes with blood transcriptomics
    Identifies genes that are dysregulated in both tissue and blood
    """
    print("\nIntegrating tissue and blood signatures...")
    
    # Convert gene IDs to common format (simplified - would use proper mapping)
    tissue_genes = tissue_de_df['Gene'].astype(str).tolist()
    
    integration_results = []
    
    if blood_de_df is not None:
        # Map blood genes to tissue genes
        blood_genes = blood_de_df['Gene'].astype(str).tolist()
        
        # Find overlapping genes
        common_genes = set(tissue_genes) & set(blood_genes)
        print(f"  Common genes in tissue and blood: {len(common_genes)}")
        
        # Create integration table
        for gene_id in tissue_genes:
            tissue_row = tissue_de_df[tissue_de_df['Gene'] == gene_id]
            blood_row = blood_de_df[blood_de_df['Gene'] == gene_id] if str(gene_id) in blood_genes else None
            
            if len(tissue_row) > 0:
                tissue_info = tissue_row.iloc[0]
                
                integration_results.append({
                    'Gene_ID': gene_id,
                    'Log2FC_Tissue': tissue_info['Log2FC'],
                    'P_adj_BH_Tissue': tissue_info['P_adj_BH'],
                    'Log2FC_Blood': blood_row['Log2FC_Blood'].values[0] if blood_row is not None and len(blood_row) > 0 else np.nan,
                    'P_adj_BH_Blood': blood_row['P_adj_BH_Blood'].values[0] if blood_row is not None and len(blood_row) > 0 else np.nan,
                    'In_Blood_Dataset': str(gene_id) in blood_genes,
                    'Dysregulated_in_Both': (blood_row is not None and len(blood_row) > 0 and 
                                           blood_row['P_adj_BH_Blood'].values[0] < 0.1) if blood_row is not None else False,
                    'Concordant_Direction': (np.sign(tissue_info['Log2FC']) == np.sign(blood_row['Log2FC_Blood'].values[0])) 
                                          if blood_row is not None and len(blood_row) > 0 and pd.notna(blood_row['Log2FC_Blood'].values[0]) else False,
                    'Biomarker_Score': abs(tissue_info['Log2FC']) * (-np.log10(tissue_info['P_adj_BH'] + 1e-300))
                })
    else:
        # No blood data available - create tissue-only structure
        print("  [NOTE] Blood data not available - creating tissue-only structure")
        for gene_id in tissue_genes:
            tissue_row = tissue_de_df[tissue_de_df['Gene'] == gene_id]
            if len(tissue_row) > 0:
                tissue_info = tissue_row.iloc[0]
                integration_results.append({
                    'Gene_ID': gene_id,
                    'Log2FC_Tissue': tissue_info['Log2FC'],
                    'P_adj_BH_Tissue': tissue_info['P_adj_BH'],
                    'Log2FC_Blood': np.nan,
                    'P_adj_BH_Blood': np.nan,
                    'In_Blood_Dataset': False,
                    'Dysregulated_in_Both': False,
                    'Concordant_Direction': False,
                    'Biomarker_Score': abs(tissue_info['Log2FC']) * (-np.log10(tissue_info['P_adj_BH'] + 1e-300))
                })
    
    integration_df = pd.DataFrame(integration_results)
    
    # Only sort if Biomarker_Score column exists
    if len(integration_df) > 0 and 'Biomarker_Score' in integration_df.columns:
        integration_df = integration_df.sort_values('Biomarker_Score', ascending=False)
    
    return integration_df

def identify_biomarker_candidates(integration_df: pd.DataFrame):
    """
    Identify top biomarker candidates based on integration criteria
    """
    print("\nIdentifying biomarker candidates...")
    
    # Criteria for biomarker candidacy:
    # 1. Significant in tissue
    # 2. Also dysregulated in blood (if available)
    # 3. Concordant direction of change
    # 4. High biomarker score
    
    # Check if required columns exist
    required_cols = ['P_adj_BH_Tissue', 'Dysregulated_in_Both', 'In_Blood_Dataset']
    if all(col in integration_df.columns for col in required_cols):
        biomarker_candidates = integration_df[
            (integration_df['P_adj_BH_Tissue'] < 0.1) &
            (
                (integration_df['Dysregulated_in_Both'] == True) |
                (integration_df['In_Blood_Dataset'] == False)  # Include if blood data unavailable
            )
        ].copy()
    else:
        # If columns don't exist, use all genes (fallback)
        print("  [WARNING] Required columns not found, using all genes as candidates")
        biomarker_candidates = integration_df.copy()
    
    # Sort by biomarker score if available
    if 'Biomarker_Score' in biomarker_candidates.columns:
        biomarker_candidates = biomarker_candidates.sort_values('Biomarker_Score', ascending=False)
    elif 'Log2FC_Tissue' in biomarker_candidates.columns:
        biomarker_candidates = biomarker_candidates.sort_values('Log2FC_Tissue', key=abs, ascending=False)
    
    print(f"  Identified {len(biomarker_candidates)} biomarker candidates")
    
    if len(biomarker_candidates) > 0:
        print(f"  Top candidates:")
        top_candidates = biomarker_candidates.head(20)
        for idx, row in top_candidates.iterrows():
            both = "Tissue+Blood" if row['Dysregulated_in_Both'] else "Tissue only"
            print(f"    {row['Gene_ID']:15s} Log2FC_T={row['Log2FC_Tissue']:7.3f}  Log2FC_B={row['Log2FC_Blood']:7.3f if pd.notna(row['Log2FC_Blood']) else 'N/A':>7}  {both}")
    
    return biomarker_candidates

def create_integration_report(fdr_threshold=0.1):
    """
    Create comprehensive blood transcriptomics integration report
    """
    print("\n" + "=" * 80)
    print("BLOOD TRANSCRIPTOMICS INTEGRATION")
    print("=" * 80)
    
    # Load tissue DE results
    tissue_de_df = load_de_results(fdr_threshold=fdr_threshold)
    if tissue_de_df is None:
        return None
    
    # Load blood transcriptomics data
    blood_expr_df, blood_metadata = load_blood_transcriptomics_data()
    
    # Analyze blood signature
    blood_de_df = None
    if blood_expr_df is not None:
        blood_de_df = analyze_blood_signature(blood_expr_df, blood_metadata)
        
        if blood_de_df is not None:
            # Save blood DE results
            blood_output = f"results/blood_de_analysis_fdr_{fdr_threshold}.tsv"
            blood_de_df.to_csv(blood_output, sep='\t', index=False)
            print(f"\n[OK] Saved blood DE results to {blood_output}")
    
    # Integrate signatures
    integration_df = integrate_tissue_blood_signatures(tissue_de_df, blood_de_df)
    
    # Identify biomarker candidates
    biomarker_candidates = identify_biomarker_candidates(integration_df)
    
    # Save integration results
    output_file = f"results/blood_tissue_integration_fdr_{fdr_threshold}.tsv"
    integration_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved integration results to {output_file}")
    
    if len(biomarker_candidates) > 0:
        biomarker_file = f"results/biomarker_candidates_fdr_{fdr_threshold}.tsv"
        biomarker_candidates.to_csv(biomarker_file, sep='\t', index=False)
        print(f"[OK] Saved biomarker candidates to {biomarker_file}")
    
    # Print summary
    print("\nIntegration Summary:")
    print(f"  Tissue DE genes: {len(tissue_de_df)}")
    if blood_de_df is not None:
        print(f"  Blood DE genes: {(blood_de_df['P_adj_BH_Blood'] < 0.1).sum()}")
        print(f"  Common genes: {integration_df['In_Blood_Dataset'].sum()}")
        print(f"  Dysregulated in both: {integration_df['Dysregulated_in_Both'].sum()}")
        print(f"  Concordant direction: {integration_df['Concordant_Direction'].sum()}")
    print(f"  Biomarker candidates: {len(biomarker_candidates)}")
    
    return integration_df, biomarker_candidates

def main():
    """Main blood transcriptomics integration pipeline"""
    
    os.makedirs("results", exist_ok=True)
    
    # Analyze at different FDR thresholds
    for fdr_threshold in [0.05, 0.1, 0.2]:
        print(f"\n{'='*80}")
        print(f"BLOOD TRANSCRIPTOMICS INTEGRATION - FDR < {fdr_threshold}")
        print(f"{'='*80}")
        
        # Create integration report
        integration_df, biomarker_candidates = create_integration_report(fdr_threshold=fdr_threshold)
    
    print("\n" + "=" * 80)
    print("BLOOD TRANSCRIPTOMICS INTEGRATION COMPLETE")
    print("=" * 80)

if __name__ == "__main__":
    main()

