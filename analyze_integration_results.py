"""
Comprehensive Analysis of Integration Results
Identifies prioritized genes, biomarker candidates, and enriched pathways
"""
import pandas as pd
import numpy as np
import os

def load_results():
    """Load all integration result files"""
    results = {}
    
    # Load GWAS integration
    gwas_file = "results/gwas_integration_report.tsv"
    if os.path.exists(gwas_file):
        results['gwas'] = pd.read_csv(gwas_file, sep='\t')
        print(f"[OK] Loaded GWAS integration: {len(results['gwas'])} genes")
    else:
        print(f"[WARNING] GWAS integration file not found")
        results['gwas'] = None
    
    # Load blood-tissue integration
    blood_files = {
        '0.05': "results/blood_tissue_integration_fdr_0.05.tsv",
        '0.1': "results/blood_tissue_integration_fdr_0.1.tsv",
        '0.2': "results/blood_tissue_integration_fdr_0.2.tsv"
    }
    results['blood'] = {}
    for threshold, file_path in blood_files.items():
        if os.path.exists(file_path):
            try:
                # Check if file is empty or has content
                file_size = os.path.getsize(file_path)
                if file_size > 10:  # At least 10 bytes (not just headers)
                    df = pd.read_csv(file_path, sep='\t')
                    if len(df) > 0:
                        results['blood'][threshold] = df
                        print(f"[OK] Loaded blood integration (FDR {threshold}): {len(df)} genes")
                    else:
                        print(f"[NOTE] Blood integration file (FDR {threshold}) is empty")
                else:
                    print(f"[NOTE] Blood integration file (FDR {threshold}) appears empty")
            except Exception as e:
                print(f"[WARNING] Could not load blood integration (FDR {threshold}): {e}")
    
    # Load pathway enrichment
    pathway_file = "results/pathway_enrichment_fdr_0.05.tsv"
    if os.path.exists(pathway_file):
        results['pathways'] = pd.read_csv(pathway_file, sep='\t')
        print(f"[OK] Loaded pathway enrichment: {len(results['pathways'])} pathways")
    else:
        results['pathways'] = None
    
    # Load annotated genes for context
    annotated_file = "results/gse122897_annotated_fdr_0.1.tsv"
    if os.path.exists(annotated_file):
        results['annotated'] = pd.read_csv(annotated_file, sep='\t')
        print(f"[OK] Loaded annotated genes: {len(results['annotated'])} genes")
    else:
        results['annotated'] = None
    
    # Load DE results for full context
    de_file = "results/gse122897_differential_expression.tsv"
    if os.path.exists(de_file):
        results['de'] = pd.read_csv(de_file, sep='\t')
        print(f"[OK] Loaded DE results: {len(results['de'])} genes")
    else:
        results['de'] = None
    
    return results

def identify_prioritized_genes(results):
    """Identify top prioritized genes from GWAS integration"""
    
    print("\n" + "=" * 80)
    print("PRIORITIZED GENES (GWAS INTEGRATION)")
    print("=" * 80)
    
    if results['gwas'] is None or len(results['gwas']) == 0:
        print("[WARNING] No GWAS integration results available")
        return None
    
    gwas_df = results['gwas'].copy()
    
    # Sort by priority score
    gwas_df = gwas_df.sort_values('Priority_Score', ascending=False)
    
    # Get top 50 prioritized genes
    top_genes = gwas_df.head(50)
    
    print(f"\nTop 50 Prioritized Genes (by Priority Score):")
    print(f"{'Rank':<6} {'Gene_ID':<12} {'Log2FC':<10} {'P_adj_BH':<12} {'Priority':<12} {'GWAS':<8}")
    print("-" * 80)
    
    for idx, (i, row) in enumerate(top_genes.iterrows(), 1):
        gene_id = str(row['Gene_ID']).rstrip('.0')
        log2fc = row['Log2FC']
        p_adj = row['P_adj_BH']
        priority = row['Priority_Score']
        gwas_status = "Yes" if row['In_GWAS_Locus'] else "No"
        
        print(f"{idx:<6} {gene_id:<12} {log2fc:>9.3f}  {p_adj:>11.4f}  {priority:>11.2f}  {gwas_status:<8}")
    
    # Get genes with highest effect sizes
    print(f"\n\nTop 20 Genes by Absolute Log2FC:")
    top_fc = gwas_df.reindex(gwas_df['Log2FC'].abs().nlargest(20).index)
    print(f"{'Rank':<6} {'Gene_ID':<12} {'Log2FC':<10} {'P_adj_BH':<12} {'Priority':<12}")
    print("-" * 80)
    
    for idx, (i, row) in enumerate(top_fc.iterrows(), 1):
        gene_id = str(row['Gene_ID']).rstrip('.0')
        log2fc = row['Log2FC']
        p_adj = row['P_adj_BH']
        priority = row['Priority_Score']
        
        print(f"{idx:<6} {gene_id:<12} {log2fc:>9.3f}  {p_adj:>11.4f}  {priority:>11.2f}")
    
    return top_genes

def analyze_biomarker_candidates(results):
    """Analyze biomarker candidates from blood-tissue integration"""
    
    print("\n" + "=" * 80)
    print("BIOMARKER CANDIDATES (BLOOD-TISSUE INTEGRATION)")
    print("=" * 80)
    
    if not results['blood'] or len(results['blood']) == 0:
        print("[WARNING] No blood-tissue integration results available")
        return None
    
    # Use FDR 0.1 threshold for analysis (good balance)
    if '0.1' in results['blood']:
        blood_df = results['blood']['0.1'].copy()
        
        print(f"\nBlood-Tissue Integration Results (FDR < 0.1):")
        print(f"  Total genes: {len(blood_df)}")
        
        # Check for biomarker-related columns
        if 'Biomarker_Score' in blood_df.columns:
            blood_df = blood_df.sort_values('Biomarker_Score', ascending=False)
            print(f"  Genes with biomarker scores: {blood_df['Biomarker_Score'].notna().sum()}")
        
        if 'Dysregulated_in_Both' in blood_df.columns:
            both = blood_df['Dysregulated_in_Both'].sum()
            print(f"  Dysregulated in both tissue and blood: {both}")
        
        if 'Concordant_Direction' in blood_df.columns:
            concordant = blood_df['Concordant_Direction'].sum()
            print(f"  Concordant direction of change: {concordant}")
        
        # Get top biomarker candidates
        print(f"\n\nTop 30 Biomarker Candidates:")
        
        if 'Biomarker_Score' in blood_df.columns:
            top_candidates = blood_df.head(30)
            
            print(f"{'Rank':<6} {'Gene_ID':<12} {'Log2FC_T':<12} {'Log2FC_B':<12} {'Score':<12} {'Both':<8}")
            print("-" * 80)
            
            for idx, (i, row) in enumerate(top_candidates.iterrows(), 1):
                gene_id = str(row.get('Gene_ID', 'N/A'))
                log2fc_t = row.get('Log2FC_Tissue', np.nan)
                log2fc_b = row.get('Log2FC_Blood', np.nan)
                score = row.get('Biomarker_Score', np.nan)
                both = "Yes" if row.get('Dysregulated_in_Both', False) else "No"
                
                fc_t_str = f"{log2fc_t:.3f}" if pd.notna(log2fc_t) else "N/A"
                fc_b_str = f"{log2fc_b:.3f}" if pd.notna(log2fc_b) else "N/A"
                score_str = f"{score:.2f}" if pd.notna(score) else "N/A"
                
                print(f"{idx:<6} {gene_id:<12} {fc_t_str:>11}  {fc_b_str:>11}  {score_str:>11}  {both:<8}")
        else:
            print("  [NOTE] Biomarker score column not found in results")
            # Use Log2FC_Tissue instead
            if 'Log2FC_Tissue' in blood_df.columns:
                top_candidates = blood_df.reindex(
                    blood_df['Log2FC_Tissue'].abs().nlargest(30).index
                )
                print(f"{'Rank':<6} {'Gene_ID':<12} {'Log2FC_Tissue':<15} {'P_adj_BH':<12}")
                print("-" * 80)
                
                for idx, (i, row) in enumerate(top_candidates.iterrows(), 1):
                    gene_id = str(row.get('Gene_ID', 'N/A'))
                    log2fc_t = row.get('Log2FC_Tissue', np.nan)
                    p_adj = row.get('P_adj_BH_Tissue', np.nan)
                    
                    fc_t_str = f"{log2fc_t:.3f}" if pd.notna(log2fc_t) else "N/A"
                    p_adj_str = f"{p_adj:.4f}" if pd.notna(p_adj) else "N/A"
                    
                    print(f"{idx:<6} {gene_id:<12} {fc_t_str:>14}  {p_adj_str:>11}")
        
        return blood_df
    
    return None

def analyze_enriched_pathways(results):
    """Analyze enriched pathways for biological insights"""
    
    print("\n" + "=" * 80)
    print("ENRICHED PATHWAYS (BIOLOGICAL INSIGHTS)")
    print("=" * 80)
    
    if results['pathways'] is None or len(results['pathways']) == 0:
        print("[WARNING] No pathway enrichment results available")
        return None
    
    pathways_df = results['pathways'].copy()
    
    # Sort by significance
    pathways_df = pathways_df.sort_values('P_value')
    
    print(f"\nEnriched Pathways (sorted by P-value):")
    print(f"  Total pathways found: {len(pathways_df)}")
    
    print(f"\n\nAll Enriched Pathways:")
    print(f"{'Rank':<6} {'Pathway':<50} {'Source':<10} {'P_value':<12} {'Size':<8}")
    print("-" * 100)
    
    for idx, (i, row) in enumerate(pathways_df.iterrows(), 1):
        pathway = row.get('Term', 'N/A')[:48]
        source = row.get('Source', 'N/A')
        p_val = row.get('P_value', np.nan)
        term_size = row.get('Term_size', 0)
        
        p_val_str = f"{p_val:.6f}" if pd.notna(p_val) else "N/A"
        
        print(f"{idx:<6} {pathway:<50} {source:<10} {p_val_str:>11}  {term_size:>7}")
    
    # Group by source
    print(f"\n\nPathways by Database:")
    if 'Source' in pathways_df.columns:
        source_counts = pathways_df['Source'].value_counts()
        for source, count in source_counts.items():
            print(f"  {source}: {count} pathways")
    
    # Identify key biological themes
    print(f"\n\nKey Biological Themes:")
    print("  [NOTE] Manual review recommended for biological interpretation")
    print("  - Review pathway terms for aneurysm-relevant processes:")
    print("    * Vascular remodeling")
    print("    * Extracellular matrix")
    print("    * Inflammation")
    print("    * Endothelial function")
    print("    * Apoptosis")
    print("    * Cell adhesion")
    
    return pathways_df

def create_integrated_summary(results):
    """Create integrated summary of all findings"""
    
    print("\n" + "=" * 80)
    print("INTEGRATED SUMMARY")
    print("=" * 80)
    
    summary = {
        'prioritized_genes': [],
        'biomarker_candidates': [],
        'enriched_pathways': []
    }
    
    # Get top 10 prioritized genes
    if results['gwas'] is not None and len(results['gwas']) > 0:
        top_gwas = results['gwas'].nlargest(10, 'Priority_Score')
        summary['prioritized_genes'] = [
            {
                'Gene_ID': str(row['Gene_ID']).rstrip('.0'),
                'Log2FC': row['Log2FC'],
                'P_adj_BH': row['P_adj_BH'],
                'Priority_Score': row['Priority_Score']
            }
            for _, row in top_gwas.iterrows()
        ]
    
    # Get top biomarker candidates
    if '0.1' in results['blood'] and len(results['blood']['0.1']) > 0:
        blood_df = results['blood']['0.1'].copy()
        if 'Biomarker_Score' in blood_df.columns:
            top_blood = blood_df.nlargest(10, 'Biomarker_Score')
            summary['biomarker_candidates'] = [
                {
                    'Gene_ID': str(row.get('Gene_ID', 'N/A')),
                    'Log2FC_Tissue': row.get('Log2FC_Tissue', np.nan),
                    'Log2FC_Blood': row.get('Log2FC_Blood', np.nan),
                    'Biomarker_Score': row.get('Biomarker_Score', np.nan)
                }
                for _, row in top_blood.iterrows()
            ]
    
    # Get enriched pathways
    if results['pathways'] is not None and len(results['pathways']) > 0:
        summary['enriched_pathways'] = [
            {
                'Pathway': row.get('Term', 'N/A'),
                'Source': row.get('Source', 'N/A'),
                'P_value': row.get('P_value', np.nan)
            }
            for _, row in results['pathways'].iterrows()
        ]
    
    return summary

def main():
    """Main analysis pipeline"""
    
    print("=" * 80)
    print("COMPREHENSIVE INTEGRATION RESULTS ANALYSIS")
    print("=" * 80)
    
    # Load all results
    results = load_results()
    
    # Analyze prioritized genes
    prioritized_genes = identify_prioritized_genes(results)
    
    # Analyze biomarker candidates
    biomarker_candidates = analyze_biomarker_candidates(results)
    
    # Analyze enriched pathways
    enriched_pathways = analyze_enriched_pathways(results)
    
    # Create integrated summary
    summary = create_integrated_summary(results)
    
    # Save summary report
    os.makedirs("results", exist_ok=True)
    
    # Save prioritized genes
    if prioritized_genes is not None:
        output_file = "results/prioritized_genes_summary.tsv"
        prioritized_genes.to_csv(output_file, sep='\t', index=False)
        print(f"\n[OK] Saved prioritized genes to {output_file}")
    
    # Save biomarker candidates
    if biomarker_candidates is not None:
        output_file = "results/biomarker_candidates_summary.tsv"
        biomarker_candidates.to_csv(output_file, sep='\t', index=False)
        print(f"[OK] Saved biomarker candidates to {output_file}")
    
    # Save enriched pathways
    if enriched_pathways is not None:
        output_file = "results/enriched_pathways_summary.tsv"
        enriched_pathways.to_csv(output_file, sep='\t', index=False)
        print(f"[OK] Saved enriched pathways to {output_file}")
    
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    
    print("\n[RECOMMENDATIONS]")
    print("1. Review top prioritized genes for validation experiments")
    print("2. Validate biomarker candidates in independent cohorts")
    print("3. Perform literature review on enriched pathways")
    print("4. Cross-reference findings with known IA biology")
    print("5. Consider functional experiments for top candidates")

if __name__ == "__main__":
    main()

