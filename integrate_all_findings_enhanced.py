"""
Enhanced Comprehensive Integration of All Analysis Results
Combines blood (GSE159610), tissue (GSE122897), GWAS, pathway enrichment, and other analyses
to create a multi-omics prioritized gene list
"""
import pandas as pd
import numpy as np
import os
from typing import Dict, List, Optional

def load_blood_results():
    """Load blood DE results (GSE159610) - prioritize annotated files"""
    files_to_try = [
        "results/gse159610_top_candidate_genes_annotated.tsv",
        "results/gse159610_substantial_fc_genes_annotated.tsv",
        "results/gse159610_full_de_results_annotated.tsv",
        "results/gse159610_top_candidate_genes.tsv",
        "results/gse159610_substantial_fc_genes.tsv",
        "results/gse159610_full_de_results.tsv"
    ]
    
    for file in files_to_try:
        if os.path.exists(file):
            df = pd.read_csv(file, sep='\t')
            # Add dataset identifier
            df['Dataset'] = 'Blood_GSE159610'
            print(f"Loaded blood DE results from: {file}")
            print(f"  Genes: {len(df):,}")
            if 'Gene_Symbol' in df.columns:
                annotated_count = df['Gene_Symbol'].notna().sum() + (df['Gene_Symbol'] != '').sum()
                print(f"  Annotated genes: {annotated_count:,}")
            return df
    
    print("[WARNING] No blood DE results found")
    return None

def load_tissue_results():
    """Load tissue DE results (GSE122897 - ruptured vs unruptured) - prioritize annotated files"""
    files_to_try = [
        "results/gse122897_top_genes_annotated.tsv",
        "results/gse122897_annotated_fdr_0.05.tsv",
        "results/gse122897_annotated_fdr_0.1.tsv",
        "results/gse122897_annotated_fdr_0.2.tsv",
        "results/gse122897_top_genes.tsv",
        "results/gse122897_differential_expression.tsv"
    ]
    
    for file in files_to_try:
        if os.path.exists(file):
            df = pd.read_csv(file, sep='\t')
            # Add dataset identifier
            df['Dataset'] = 'Tissue_GSE122897'
            # Ensure Gene_ID column exists (may be 'Gene' in GSE122897)
            if 'Gene' in df.columns and 'Gene_ID' not in df.columns:
                df['Gene_ID'] = df['Gene'].astype(str)
            print(f"Loaded tissue DE results from: {file}")
            print(f"  Genes: {len(df):,}")
            if 'Gene_Symbol' in df.columns:
                annotated_count = df['Gene_Symbol'].notna().sum() + (df['Gene_Symbol'] != '').sum()
                print(f"  Annotated genes: {annotated_count:,}")
            return df
    
    print("[WARNING] No tissue DE results found")
    return None

def load_gwas_integration():
    """Load GWAS integration results"""
    gwas_file = "results/gwas_integration_report.tsv"
    
    if os.path.exists(gwas_file):
        df = pd.read_csv(gwas_file, sep='\t')
        print(f"Loaded GWAS integration: {len(df):,} genes")
        return df
    
    print("[WARNING] No GWAS integration results found")
    return None

def load_blood_tissue_integration():
    """Load blood-tissue cross-dataset integration"""
    integration_file = "results/blood_tissue_integration_fdr_0.1.tsv"
    
    if os.path.exists(integration_file):
        df = pd.read_csv(integration_file, sep='\t')
        print(f"Loaded blood-tissue integration: {len(df):,} genes")
        return df
    
    print("[WARNING] No blood-tissue integration results found")
    return None

def load_pathway_results():
    """Load pathway enrichment results"""
    pathway_files = [
        "results/pathway_enrichment_fdr_0.05.tsv",
        "results/pathway_enrichment_fdr_0.1.tsv"
    ]
    
    for file in pathway_files:
        if os.path.exists(file):
            df = pd.read_csv(file, sep='\t')
            print(f"Loaded pathway enrichment from: {file}")
            return df
    
    print("[WARNING] No pathway enrichment results found")
    return None

def load_annotated_genes():
    """Load genes with annotations"""
    files_to_try = [
        "results/gse159610_top_candidate_genes_annotated.tsv",
        "results/gse122897_top_genes_annotated.tsv"
    ]
    
    for file in files_to_try:
        if os.path.exists(file):
            df = pd.read_csv(file, sep='\t')
            print(f"Loaded annotated genes from: {file}")
            return df
    
    return None

def create_unified_gene_list(blood_df, tissue_df):
    """Create a unified gene list combining blood and tissue results using gene symbols"""
    print("\n" + "=" * 80)
    print("CREATING UNIFIED GENE LIST (GENE SYMBOL-BASED MATCHING)")
    print("=" * 80)
    
    unified_genes = {}
    
    # Process blood results - use Gene_Symbol as key if available, otherwise Gene_ID
    if blood_df is not None:
        for idx, row in blood_df.iterrows():
            # Try to get gene symbol first, fall back to ID
            gene_symbol = str(row.get('Gene_Symbol', '')).strip()
            gene_id = str(row.get('Gene_ID', '')).strip()
            
            # Use symbol as primary key if available, otherwise use ID
            if gene_symbol and gene_symbol != 'nan' and gene_symbol != '':
                key = gene_symbol.upper()  # Normalize to uppercase for matching
                key_type = 'Symbol'
            elif gene_id and gene_id != 'nan' and gene_id != '':
                key = gene_id
                key_type = 'ID'
            else:
                continue
            
            if key not in unified_genes:
                unified_genes[key] = {
                    'Gene_Symbol': gene_symbol if gene_symbol != 'nan' else np.nan,
                    'Blood_Gene_ID': gene_id if gene_id != 'nan' else np.nan,
                    'Entrez_ID': row.get('Entrez_ID', np.nan),
                    'Blood_Log2FC': row.get('Log2FC', np.nan),
                    'Blood_P_value': row.get('P_value', np.nan),
                    'Blood_P_adj': row.get('P_adj_BH', np.nan),
                    'Blood_Mean_Aneurysm': row.get('Mean_Aneurysm', np.nan),
                    'Blood_Mean_Control': row.get('Mean_Control', np.nan),
                    'In_Blood': True,
                    'In_Tissue': False,
                    'Gene_Name': row.get('Gene_Name', np.nan)
                }
            else:
                unified_genes[key]['Blood_Log2FC'] = row.get('Log2FC', np.nan)
                unified_genes[key]['Blood_P_value'] = row.get('P_value', np.nan)
                unified_genes[key]['Blood_P_adj'] = row.get('P_adj_BH', np.nan)
                unified_genes[key]['In_Blood'] = True
                if pd.isna(unified_genes[key]['Gene_Symbol']) or unified_genes[key]['Gene_Symbol'] == '':
                    unified_genes[key]['Gene_Symbol'] = gene_symbol if gene_symbol != 'nan' else np.nan
                if pd.isna(unified_genes[key]['Gene_Name']) or unified_genes[key]['Gene_Name'] == '':
                    unified_genes[key]['Gene_Name'] = row.get('Gene_Name', np.nan)
                # Add Entrez_ID if not present
                if pd.isna(unified_genes[key].get('Entrez_ID', np.nan)) or unified_genes[key].get('Entrez_ID', '') == '':
                    unified_genes[key]['Entrez_ID'] = row.get('Entrez_ID', np.nan)
    
    # Process tissue results - use Gene_Symbol as key if available, otherwise Gene_ID
    if tissue_df is not None:
        for idx, row in tissue_df.iterrows():
            # Try to get gene symbol first, fall back to ID
            gene_symbol = str(row.get('Gene_Symbol', '')).strip()
            gene_id = str(row.get('Gene_ID', row.get('Gene', ''))).strip()
            
            # Use symbol as primary key if available, otherwise use ID
            if gene_symbol and gene_symbol != 'nan' and gene_symbol != '':
                key = gene_symbol.upper()  # Normalize to uppercase for matching
                key_type = 'Symbol'
            elif gene_id and gene_id != 'nan' and gene_id != '':
                key = gene_id
                key_type = 'ID'
            else:
                continue
            
            if key not in unified_genes:
                unified_genes[key] = {
                    'Gene_Symbol': gene_symbol if gene_symbol != 'nan' else np.nan,
                    'Tissue_Gene_ID': gene_id if gene_id != 'nan' else np.nan,
                    'Entrez_ID': row.get('Entrez_ID', np.nan),
                    'Tissue_Log2FC': row.get('Log2FC', np.nan),
                    'Tissue_P_value': row.get('P_value', np.nan),
                    'Tissue_P_adj': row.get('P_adj_BH', np.nan),
                    'Tissue_Mean_Ruptured': row.get('Mean_Ruptured', np.nan),
                    'Tissue_Mean_Unruptured': row.get('Mean_Unruptured', np.nan),
                    'In_Blood': False,
                    'In_Tissue': True,
                    'Gene_Name': row.get('Gene_Name', np.nan)
                }
            else:
                unified_genes[key]['Tissue_Log2FC'] = row.get('Log2FC', np.nan)
                unified_genes[key]['Tissue_P_value'] = row.get('P_value', np.nan)
                unified_genes[key]['Tissue_P_adj'] = row.get('P_adj_BH', np.nan)
                unified_genes[key]['In_Tissue'] = True
                if pd.isna(unified_genes[key]['Gene_Symbol']) or unified_genes[key]['Gene_Symbol'] == '':
                    unified_genes[key]['Gene_Symbol'] = gene_symbol if gene_symbol != 'nan' else np.nan
                if pd.isna(unified_genes[key]['Gene_Name']) or unified_genes[key]['Gene_Name'] == '':
                    unified_genes[key]['Gene_Name'] = row.get('Gene_Name', np.nan)
                # Add tissue gene ID if not present
                if 'Tissue_Gene_ID' not in unified_genes[key] or pd.isna(unified_genes[key].get('Tissue_Gene_ID', np.nan)):
                    unified_genes[key]['Tissue_Gene_ID'] = gene_id if gene_id != 'nan' else np.nan
                # Add Entrez_ID if not present
                if pd.isna(unified_genes[key].get('Entrez_ID', np.nan)) or unified_genes[key].get('Entrez_ID', '') == '':
                    unified_genes[key]['Entrez_ID'] = row.get('Entrez_ID', np.nan)
    
    # Convert to DataFrame
    unified_df = pd.DataFrame(list(unified_genes.values()))
    
    # Create a unified Gene_ID column (prefer symbol, fall back to IDs)
    if 'Gene_Symbol' in unified_df.columns:
        unified_df['Gene_ID'] = unified_df['Gene_Symbol'].fillna(
            unified_df.get('Blood_Gene_ID', pd.Series())).fillna(
            unified_df.get('Tissue_Gene_ID', pd.Series()))
    else:
        unified_df['Gene_ID'] = unified_df.get('Blood_Gene_ID', pd.Series()).fillna(
            unified_df.get('Tissue_Gene_ID', pd.Series()))
    
    print(f"\nUnified gene list: {len(unified_df):,} unique genes")
    print(f"  In blood only: {(unified_df['In_Blood'] & ~unified_df['In_Tissue']).sum():,}")
    print(f"  In tissue only: {(unified_df['In_Tissue'] & ~unified_df['In_Blood']).sum():,}")
    print(f"  In both: {(unified_df['In_Blood'] & unified_df['In_Tissue']).sum():,}")
    
    # Show some examples of cross-validated genes
    cross_validated = unified_df[(unified_df['In_Blood'] & unified_df['In_Tissue'])]
    if len(cross_validated) > 0:
        print(f"\n  Examples of cross-validated genes (in both datasets):")
        for idx, row in cross_validated.head(5).iterrows():
            symbol = row.get('Gene_Symbol', 'N/A')
            blood_fc = row.get('Blood_Log2FC', np.nan)
            tissue_fc = row.get('Tissue_Log2FC', np.nan)
            print(f"    {symbol}: Blood Log2FC={blood_fc:.2f}, Tissue Log2FC={tissue_fc:.2f}")
    
    return unified_df

def add_gwas_integration(unified_df, gwas_df):
    """Add GWAS integration information using gene symbols or IDs"""
    if gwas_df is None:
        unified_df['In_GWAS'] = False
        unified_df['GWAS_Priority_Score'] = 0.0
        unified_df['In_GWAS_Locus'] = False
        return unified_df
    
    print("\nAdding GWAS integration data...")
    
    # Create mapping from Gene_Symbol and Gene_ID to GWAS data
    gwas_map_by_symbol = {}
    gwas_map_by_id = {}
    
    for idx, row in gwas_df.iterrows():
        gene_id = str(row.get('Gene_ID', '')).strip()
        gene_symbol = str(row.get('Gene_Symbol', '')).strip()
        
        gwas_data = {
            'In_GWAS': True,
            'GWAS_Priority_Score': row.get('Priority_Score', 0.0),
            'In_GWAS_Locus': row.get('In_GWAS_Locus', False)
        }
        
        # Map by symbol if available
        if gene_symbol and gene_symbol != 'nan' and gene_symbol != '':
            gwas_map_by_symbol[gene_symbol.upper()] = gwas_data
        
        # Map by ID as backup (normalize by removing decimals)
        if gene_id and gene_id != 'nan' and gene_id != '':
            # Remove decimal if present (e.g., 406979.0 -> 406979)
            gene_id_clean = gene_id.split('.')[0]
            gwas_map_by_id[gene_id_clean] = gwas_data
    
    # Add GWAS columns - try symbol first, then Entrez_ID, then Gene_ID
    def get_gwas_data(row):
        symbol = str(row.get('Gene_Symbol', '')).strip().upper() if pd.notna(row.get('Gene_Symbol', '')) else ''
        entrez_id = str(row.get('Entrez_ID', '')).strip() if pd.notna(row.get('Entrez_ID', '')) else ''
        gene_id = str(row.get('Gene_ID', '')).strip() if pd.notna(row.get('Gene_ID', '')) else ''
        
        # Try symbol first
        if symbol and symbol in gwas_map_by_symbol:
            return gwas_map_by_symbol[symbol]
        
        # Try Entrez_ID (GWAS file uses Entrez IDs)
        if entrez_id and entrez_id != 'nan' and entrez_id != '':
            # Remove decimal if present (e.g., 157627.0 -> 157627)
            entrez_clean = entrez_id.split('.')[0]
            if entrez_clean in gwas_map_by_id:
                return gwas_map_by_id[entrez_clean]
        
        # Try Gene_ID as backup
        if gene_id and gene_id != 'nan' and gene_id != '':
            gene_id_clean = gene_id.split('.')[0]
            if gene_id_clean in gwas_map_by_id:
                return gwas_map_by_id[gene_id_clean]
        
        # Default: no match
        return {'In_GWAS': False, 'GWAS_Priority_Score': 0.0, 'In_GWAS_Locus': False}
    
    gwas_results = unified_df.apply(get_gwas_data, axis=1)
    unified_df['In_GWAS'] = [r['In_GWAS'] for r in gwas_results]
    unified_df['GWAS_Priority_Score'] = [r['GWAS_Priority_Score'] for r in gwas_results]
    unified_df['In_GWAS_Locus'] = [r['In_GWAS_Locus'] for r in gwas_results]
    
    in_gwas = unified_df['In_GWAS'].sum()
    print(f"  Genes with GWAS data: {in_gwas:,}")
    
    return unified_df

def add_pathway_integration(unified_df, pathway_df):
    """Add pathway enrichment information"""
    if pathway_df is None:
        unified_df['In_Pathway'] = False
        unified_df['Pathway_Term'] = np.nan
        return unified_df
    
    print("\nAdding pathway enrichment data...")
    
    # Extract genes from pathway results
    pathway_genes = {}
    if 'Genes' in pathway_df.columns:
        for idx, row in pathway_df.iterrows():
            term = row.get('Term', 'Unknown')
            genes_str = row.get('Genes', '')
            if pd.notna(genes_str):
                # Genes might be separated by comma, semicolon, or space
                genes = str(genes_str).replace(', ', ',').replace('; ', ',').split(',')
                for gene in genes:
                    gene = gene.strip()
                    if gene and gene != 'None':
                        if gene not in pathway_genes:
                            pathway_genes[gene] = []
                        pathway_genes[gene].append(term)
    
    # Try to match by Gene_Symbol first, then Gene_ID
    unified_df['In_Pathway'] = False
    unified_df['Pathway_Term'] = np.nan
    
    for idx, row in unified_df.iterrows():
        gene_symbol = str(row.get('Gene_Symbol', ''))
        gene_id = str(row.get('Gene_ID', ''))
        
        # Try symbol match
        if gene_symbol in pathway_genes:
            unified_df.at[idx, 'In_Pathway'] = True
            unified_df.at[idx, 'Pathway_Term'] = '; '.join(pathway_genes[gene_symbol])
        # Try ID match
        elif gene_id in pathway_genes:
            unified_df.at[idx, 'In_Pathway'] = True
            unified_df.at[idx, 'Pathway_Term'] = '; '.join(pathway_genes[gene_id])
    
    in_pathway = unified_df['In_Pathway'].sum()
    print(f"  Genes in enriched pathways: {in_pathway:,}")
    
    return unified_df

def calculate_multi_omics_priority_score(row):
    """Calculate comprehensive multi-omics priority score"""
    score = 0.0
    
    # Blood evidence (GSE159610)
    if row.get('In_Blood', False):
        blood_fc = abs(row.get('Blood_Log2FC', 0))
        blood_p = row.get('Blood_P_value', 1.0)
        blood_padj = row.get('Blood_P_adj', 1.0)
        
        # Effect size in blood
        if blood_fc > 2.0:
            score += 8
        elif blood_fc > 1.585:
            score += 6
        elif blood_fc > 1.0:
            score += 4
        elif blood_fc > 0.585:
            score += 2
        
        # Statistical significance in blood
        if blood_padj < 0.05:
            score += 5
        elif blood_padj < 0.1:
            score += 3
        elif blood_p < 0.05:
            score += 2
    
    # Tissue evidence (GSE122897)
    if row.get('In_Tissue', False):
        tissue_fc = abs(row.get('Tissue_Log2FC', 0))
        tissue_p = row.get('Tissue_P_value', 1.0)
        tissue_padj = row.get('Tissue_P_adj', 1.0)
        
        # Effect size in tissue
        if tissue_fc > 2.0:
            score += 8
        elif tissue_fc > 1.585:
            score += 6
        elif tissue_fc > 1.0:
            score += 4
        elif tissue_fc > 0.585:
            score += 2
        
        # Statistical significance in tissue
        if tissue_padj < 0.05:
            score += 5
        elif tissue_padj < 0.1:
            score += 3
        elif tissue_p < 0.05:
            score += 2
    
    # Cross-validation bonus (present in both datasets)
    if row.get('In_Blood', False) and row.get('In_Tissue', False):
        # Check if same direction
        blood_fc = row.get('Blood_Log2FC', 0)
        tissue_fc = row.get('Tissue_Log2FC', 0)
        if (blood_fc > 0 and tissue_fc > 0) or (blood_fc < 0 and tissue_fc < 0):
            score += 10  # Same direction in both datasets
        else:
            score += 5   # Present in both but different direction
    
    # GWAS evidence
    if row.get('In_GWAS', False):
        gwas_score = row.get('GWAS_Priority_Score', 0.0)
        score += min(gwas_score * 0.5, 5)  # Cap at 5 points
        
        if row.get('In_GWAS_Locus', False):
            score += 3
    
    # Pathway enrichment
    if row.get('In_Pathway', False):
        score += 3
    
    return score

def integrate_all_results():
    """Integrate all multi-omics analysis results"""
    print("=" * 80)
    print("ENHANCED MULTI-OMICS GENE PRIORITIZATION")
    print("=" * 80)
    
    # Load all data sources
    print("\nLoading data sources...")
    blood_df = load_blood_results()
    tissue_df = load_tissue_results()
    gwas_df = load_gwas_integration()
    pathway_df = load_pathway_results()
    annotated_df = load_annotated_genes()
    
    # Create unified gene list
    if blood_df is None and tissue_df is None:
        print("[ERROR] No DE results found!")
        return None
    
    unified_df = create_unified_gene_list(blood_df, tissue_df)
    
    # Add annotations if available
    if annotated_df is not None:
        if 'Gene_Symbol' in annotated_df.columns:
            symbol_map = {}
            if 'Gene_ID' in annotated_df.columns:
                symbol_map = dict(zip(annotated_df['Gene_ID'].astype(str), annotated_df['Gene_Symbol']))
            # Update missing symbols
            mask = unified_df['Gene_Symbol'].isna() | (unified_df['Gene_Symbol'] == '')
            unified_df.loc[mask, 'Gene_Symbol'] = unified_df.loc[mask, 'Gene_ID'].astype(str).map(
                lambda x: symbol_map.get(x, np.nan)
            )
    
    # Add GWAS integration
    unified_df = add_gwas_integration(unified_df, gwas_df)
    
    # Add pathway integration
    unified_df = add_pathway_integration(unified_df, pathway_df)
    
    # Calculate multi-omics priority scores
    print("\nCalculating multi-omics priority scores...")
    unified_df['Multi_Omics_Priority_Score'] = unified_df.apply(
        calculate_multi_omics_priority_score, axis=1
    )
    
    # Sort by priority score
    unified_df = unified_df.sort_values('Multi_Omics_Priority_Score', ascending=False)
    
    return unified_df

def create_prioritized_lists(integrated_df):
    """Create different prioritized gene lists"""
    print("\n" + "=" * 80)
    print("CREATING PRIORITIZED GENE LISTS")
    print("=" * 80)
    
    os.makedirs("results", exist_ok=True)
    
    # Top 50 overall
    top_50 = integrated_df.head(50).copy()
    output_file = "results/multi_omics_prioritized_top50.tsv"
    top_50.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved top 50 multi-omics prioritized genes to {output_file}")
    
    # Top 20 cross-validated (in both blood and tissue)
    cross_validated = integrated_df[
        (integrated_df['In_Blood'] == True) & 
        (integrated_df['In_Tissue'] == True)
    ].head(20).copy()
    if len(cross_validated) > 0:
        output_file = "results/multi_omics_cross_validated_genes.tsv"
        cross_validated.to_csv(output_file, sep='\t', index=False)
        print(f"[OK] Saved {len(cross_validated)} cross-validated genes to {output_file}")
    
    # Top 20 with GWAS support
    with_gwas = integrated_df[integrated_df['In_GWAS'] == True].head(20).copy()
    if len(with_gwas) > 0:
        output_file = "results/multi_omics_gwas_supported_genes.tsv"
        with_gwas.to_csv(output_file, sep='\t', index=False)
        print(f"[OK] Saved {len(with_gwas)} GWAS-supported genes to {output_file}")
    
    # Top 20 with pathway enrichment
    in_pathway = integrated_df[integrated_df['In_Pathway'] == True].head(20).copy()
    if len(in_pathway) > 0:
        output_file = "results/multi_omics_pathway_genes.tsv"
        in_pathway.to_csv(output_file, sep='\t', index=False)
        print(f"[OK] Saved {len(in_pathway)} pathway-enriched genes to {output_file}")
    
    # Top 20 blood-specific
    blood_only = integrated_df[
        (integrated_df['In_Blood'] == True) & 
        (integrated_df['In_Tissue'] == False)
    ].head(20).copy()
    if len(blood_only) > 0:
        output_file = "results/multi_omics_blood_specific_genes.tsv"
        blood_only.to_csv(output_file, sep='\t', index=False)
        print(f"[OK] Saved {len(blood_only)} blood-specific genes to {output_file}")
    
    # Top 20 tissue-specific
    tissue_only = integrated_df[
        (integrated_df['In_Tissue'] == True) & 
        (integrated_df['In_Blood'] == False)
    ].head(20).copy()
    if len(tissue_only) > 0:
        output_file = "results/multi_omics_tissue_specific_genes.tsv"
        tissue_only.to_csv(output_file, sep='\t', index=False)
        print(f"[OK] Saved {len(tissue_only)} tissue-specific genes to {output_file}")

def generate_summary_report(integrated_df):
    """Generate comprehensive summary report"""
    print("\n" + "=" * 80)
    print("GENERATING SUMMARY REPORT")
    print("=" * 80)
    
    report = []
    report.append("=" * 80)
    report.append("ENHANCED MULTI-OMICS GENE PRIORITIZATION SUMMARY")
    report.append("=" * 80)
    report.append("")
    
    report.append(f"Total unique genes analyzed: {len(integrated_df):,}")
    report.append("")
    
    # Dataset coverage
    in_blood = integrated_df['In_Blood'].sum()
    in_tissue = integrated_df['In_Tissue'].sum()
    in_both = ((integrated_df['In_Blood'] == True) & (integrated_df['In_Tissue'] == True)).sum()
    in_gwas = integrated_df['In_GWAS'].sum()
    in_pathway = integrated_df['In_Pathway'].sum()
    
    report.append("Dataset Coverage:")
    report.append(f"  In blood (GSE159610): {in_blood:,} ({in_blood/len(integrated_df)*100:.1f}%)")
    report.append(f"  In tissue (GSE122897): {in_tissue:,} ({in_tissue/len(integrated_df)*100:.1f}%)")
    report.append(f"  In both datasets: {in_both:,} ({in_both/len(integrated_df)*100:.1f}%)")
    report.append(f"  With GWAS support: {in_gwas:,} ({in_gwas/len(integrated_df)*100:.1f}%)")
    report.append(f"  In enriched pathways: {in_pathway:,} ({in_pathway/len(integrated_df)*100:.1f}%)")
    report.append("")
    
    # Priority score statistics
    report.append("Multi-Omics Priority Score Statistics:")
    report.append(f"  Mean: {integrated_df['Multi_Omics_Priority_Score'].mean():.2f}")
    report.append(f"  Median: {integrated_df['Multi_Omics_Priority_Score'].median():.2f}")
    report.append(f"  Max: {integrated_df['Multi_Omics_Priority_Score'].max():.2f}")
    report.append(f"  75th percentile: {integrated_df['Multi_Omics_Priority_Score'].quantile(0.75):.2f}")
    report.append("")
    
    # Top 20 prioritized genes
    report.append("Top 20 Multi-Omics Prioritized Genes:")
    report.append("-" * 80)
    
    top_20 = integrated_df.head(20)
    for i, (idx, row) in enumerate(top_20.iterrows(), 1):
        symbol = row.get('Gene_Symbol', '')
        if pd.isna(symbol) or symbol == '':
            symbol = str(row.get('Gene_ID', 'Unknown'))
        else:
            symbol = str(symbol)
        
        score = row.get('Multi_Omics_Priority_Score', 0)
        
        # Build evidence string
        evidence = []
        if row.get('In_Blood', False):
            blood_fc = row.get('Blood_Log2FC', 0)
            evidence.append(f"Blood:Log2FC={blood_fc:.2f}")
        if row.get('In_Tissue', False):
            tissue_fc = row.get('Tissue_Log2FC', 0)
            evidence.append(f"Tissue:Log2FC={tissue_fc:.2f}")
        if row.get('In_GWAS', False):
            evidence.append("GWAS")
        if row.get('In_Pathway', False):
            evidence.append("Pathway")
        
        evidence_str = ' | '.join(evidence) if evidence else 'No evidence'
        
        report.append(f"{i:2d}. {symbol:20s} | Score: {score:5.1f} | {evidence_str}")
    
    report.append("")
    report.append("=" * 80)
    
    # Save report
    report_text = "\n".join(report)
    output_file = "results/multi_omics_prioritization_report.txt"
    with open(output_file, 'w') as f:
        f.write(report_text)
    
    print(report_text)
    print(f"\n[OK] Saved summary report to {output_file}")

def main():
    """Main integration pipeline"""
    # Integrate all results
    integrated_df = integrate_all_results()
    
    if integrated_df is None:
        print("[ERROR] Could not integrate results")
        return
    
    # Save full integrated results
    output_file = "results/multi_omics_comprehensive_prioritized_genes.tsv"
    integrated_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved full multi-omics integrated results to {output_file}")
    
    # Create prioritized lists
    create_prioritized_lists(integrated_df)
    
    # Generate summary report
    generate_summary_report(integrated_df)
    
    print("\n" + "=" * 80)
    print("ENHANCED MULTI-OMICS INTEGRATION COMPLETE!")
    print("=" * 80)

if __name__ == "__main__":
    main()





