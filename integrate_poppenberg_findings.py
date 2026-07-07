"""
Integrate Poppenberg et al. 2023 Findings with Existing Analyses
Combines the 39-gene growth trajectory signature with existing rupture prediction analyses
"""
import pandas as pd
import numpy as np
from pathlib import Path
import json

def load_poppenberg_genes():
    """Load the extracted Poppenberg et al. 2023 gene list"""
    gene_file = Path("results/poppenberg_2023_integration/poppenberg_2023_39_genes.tsv")
    if gene_file.exists():
        return pd.read_csv(gene_file, sep='\t')
    else:
        print(f"[ERROR] Gene file not found: {gene_file}")
        return None

def load_existing_analyses():
    """Load existing differential expression analyses"""
    results = {}
    
    # GSE122897 rupture analysis
    de_file = Path("results/gse122897_differential_expression.tsv")
    if de_file.exists():
        results['gse122897_rupture'] = pd.read_csv(de_file, sep='\t')
        print(f"[OK] Loaded GSE122897 DE results: {len(results['gse122897_rupture'])} genes")
    
    # Top genes from GSE122897
    top_file = Path("results/gse122897_top_genes.tsv")
    if top_file.exists():
        results['gse122897_top'] = pd.read_csv(top_file, sep='\t')
        print(f"[OK] Loaded GSE122897 top genes: {len(results['gse122897_top'])} genes")
    
    return results

def compare_gene_lists(poppenberg_df, existing_results):
    """Compare Poppenberg genes with existing analyses"""
    
    print("\n" + "=" * 80)
    print("GENE LIST COMPARISON")
    print("=" * 80)
    
    poppenberg_genes = set(poppenberg_df['Gene_Symbol'].str.upper())
    
    comparisons = {}
    
    for analysis_name, df in existing_results.items():
        print(f"\n[INFO] Comparing with {analysis_name}...")
        
        # Find gene column
        gene_col = None
        for col in ['Gene', 'gene', 'Gene_Symbol', 'Symbol']:
            if col in df.columns:
                gene_col = col
                break
        
        if gene_col is None:
            print(f"  [WARNING] Could not find gene column in {analysis_name}")
            continue
        
        existing_genes = set(df[gene_col].astype(str).str.upper())
        overlap = poppenberg_genes & existing_genes
        
        comparisons[analysis_name] = {
            'total_existing': len(existing_genes),
            'overlapping_genes': sorted(list(overlap)),
            'overlap_count': len(overlap),
            'overlap_pct': len(overlap) / len(poppenberg_genes) * 100 if poppenberg_genes else 0
        }
        
        if overlap:
            print(f"  [OK] Found {len(overlap)} overlapping genes ({len(overlap)/len(poppenberg_genes)*100:.1f}% of Poppenberg genes)")
            print(f"       Genes: {', '.join(sorted(overlap))}")
        else:
            print(f"  [INFO] No overlapping genes found")
    
    return comparisons

def create_integrated_gene_list(poppenberg_df, existing_results, comparisons):
    """Create an integrated gene list combining all findings"""
    
    print("\n" + "=" * 80)
    print("CREATING INTEGRATED GENE LIST")
    print("=" * 80)
    
    # Start with Poppenberg genes
    integrated = poppenberg_df.copy()
    
    # Add columns for other analyses
    integrated['Found_in_GSE122897'] = False
    integrated['GSE122897_Log2FC'] = np.nan
    integrated['GSE122897_P_adj'] = np.nan
    integrated['Consistent_Direction'] = np.nan  # True if same direction in both analyses
    
    # Check GSE122897 overlap
    if 'gse122897_rupture' in existing_results:
        gse_df = existing_results['gse122897_rupture']
        
        # Find gene column
        gene_col = None
        for col in ['Gene', 'gene', 'Gene_Symbol', 'Symbol']:
            if col in gse_df.columns:
                gene_col = col
                break
        
        if gene_col:
            gse_genes_upper = gse_df[gene_col].astype(str).str.upper()
            
            for idx, row in integrated.iterrows():
                gene = row['Gene_Symbol'].upper()
                matches = gse_df[gse_genes_upper == gene]
                
                if not matches.empty:
                    integrated.loc[idx, 'Found_in_GSE122897'] = True
                    
                    # Get Log2FC and P_adj
                    if 'Log2FC' in matches.columns:
                        integrated.loc[idx, 'GSE122897_Log2FC'] = matches['Log2FC'].iloc[0]
                    elif 'log2FoldChange' in matches.columns:
                        integrated.loc[idx, 'GSE122897_Log2FC'] = matches['log2FoldChange'].iloc[0]
                    
                    if 'P_adj_BH' in matches.columns:
                        integrated.loc[idx, 'GSE122897_P_adj'] = matches['P_adj_BH'].iloc[0]
                    elif 'padj' in matches.columns:
                        integrated.loc[idx, 'GSE122897_P_adj'] = matches['padj'].iloc[0]
                    
                    # Check if direction is consistent
                    poppenberg_fc = row['Log2FC']
                    gse_fc = integrated.loc[idx, 'GSE122897_Log2FC']
                    
                    if not pd.isna(gse_fc):
                        # Both positive or both negative
                        consistent = (poppenberg_fc > 0 and gse_fc > 0) or (poppenberg_fc < 0 and gse_fc < 0)
                        integrated.loc[idx, 'Consistent_Direction'] = consistent
    
    # Add priority score based on multiple factors
    integrated['Priority_Score'] = 0
    
    # Points for significance
    integrated.loc[integrated['Significant_0001'], 'Priority_Score'] += 5
    integrated.loc[integrated['Significant_001'], 'Priority_Score'] += 3
    integrated.loc[integrated['Significant_01'], 'Priority_Score'] += 1
    
    # Points for high fold-change
    integrated.loc[integrated['FC_Category'] == 'High', 'Priority_Score'] += 3
    
    # Points for being in other analyses
    integrated.loc[integrated['Found_in_GSE122897'], 'Priority_Score'] += 2
    
    # Points for consistent direction
    integrated.loc[integrated['Consistent_Direction'] == True, 'Priority_Score'] += 3
    
    # Points for priority status
    integrated.loc[integrated['Priority_Status'] != 'Standard', 'Priority_Score'] += 5
    
    # Sort by priority score
    integrated = integrated.sort_values('Priority_Score', ascending=False)
    
    print(f"\n[OK] Created integrated gene list with {len(integrated)} genes")
    print(f"  - Found in GSE122897: {integrated['Found_in_GSE122897'].sum()}")
    print(f"  - Consistent direction: {integrated['Consistent_Direction'].sum()}")
    print(f"\nTop 10 priority genes:")
    for idx, row in integrated.head(10).iterrows():
        status = row['Priority_Status']
        score = row['Priority_Score']
        in_gse = "YES" if row['Found_in_GSE122897'] else "NO"
        print(f"  {row['Gene_Symbol']:12s} Priority={score:2d} | Status={status:30s} | GSE122897={in_gse}")
    
    return integrated

def generate_integration_summary(poppenberg_df, integrated_df, comparisons, pat_info):
    """Generate comprehensive integration summary"""
    
    summary = {
        'poppenberg_study': {
            'paper': 'RNA Expression Signatures of Intracranial Aneurysm Growth Trajectory Identified in Circulating Whole Blood',
            'authors': 'Poppenberg, K.E. et al.',
            'year': 2023,
            'journal': 'Journal of Personalized Medicine',
            'doi': '10.3390/jpm13020266',
            'n_samples': 66,
            'n_genes': len(poppenberg_df),
            'model_auc': 0.86
        },
        'gene_statistics': {
            'total_genes': len(integrated_df),
            'decreased_in_growing': (integrated_df['Direction'] == 'Decreased').sum(),
            'increased_in_growing': (integrated_df['Direction'] == 'Increased').sum(),
            'high_fold_change': (integrated_df['FC_Category'] == 'High').sum(),
            'found_in_gse122897': integrated_df['Found_in_GSE122897'].sum(),
            'consistent_direction': integrated_df['Consistent_Direction'].sum(),
            'high_priority_genes': (integrated_df['Priority_Score'] >= 10).sum()
        },
        'priority_genes': {
            'top_10': integrated_df.head(10)[['Gene_Symbol', 'Priority_Score', 'Priority_Status', 'Log2FC', 'Found_in_GSE122897']].to_dict('records')
        },
        'comparisons': comparisons,
        'pat_score': pat_info,
        'recommendations': {
            'top_biomarker_candidates': integrated_df[integrated_df['Priority_Score'] >= 10]['Gene_Symbol'].tolist(),
            'validation_priorities': integrated_df[integrated_df['Found_in_GSE122897'] & (integrated_df['Consistent_Direction'] == True)]['Gene_Symbol'].tolist(),
            'novel_findings': integrated_df[~integrated_df['Found_in_GSE122897']]['Gene_Symbol'].tolist()
        }
    }
    
    return summary

def main():
    """Main integration function"""
    
    print("=" * 80)
    print("INTEGRATING POPPENBERG ET AL. 2023 FINDINGS")
    print("=" * 80)
    
    # Load Poppenberg genes
    print("\n[1] Loading Poppenberg et al. 2023 gene list...")
    poppenberg_df = load_poppenberg_genes()
    if poppenberg_df is None:
        return
    
    print(f"  [OK] Loaded {len(poppenberg_df)} genes from Poppenberg et al. 2023")
    
    # Load PAT score info
    pat_file = Path("results/poppenberg_2023_integration/pat_score_formula_info.json")
    pat_info = {}
    if pat_file.exists():
        with open(pat_file, 'r') as f:
            pat_info = json.load(f)
        print(f"  [OK] Loaded PAT score information")
    
    # Load existing analyses
    print("\n[2] Loading existing analyses...")
    existing_results = load_existing_analyses()
    
    # Compare gene lists
    print("\n[3] Comparing gene lists...")
    comparisons = compare_gene_lists(poppenberg_df, existing_results)
    
    # Create integrated gene list
    print("\n[4] Creating integrated gene list...")
    integrated_df = create_integrated_gene_list(poppenberg_df, existing_results, comparisons)
    
    # Generate summary
    print("\n[5] Generating integration summary...")
    summary = generate_integration_summary(poppenberg_df, integrated_df, comparisons, pat_info)
    
    # Save results
    output_dir = Path("results/poppenberg_2023_integration")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save integrated gene list
    integrated_file = output_dir / "integrated_gene_list.tsv"
    integrated_df.to_csv(integrated_file, sep='\t', index=False)
    print(f"\n[OK] Saved integrated gene list to: {integrated_file}")
    
    # Save summary
    summary_file = output_dir / "integration_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"[OK] Saved integration summary to: {summary_file}")
    
    # Print final summary
    print("\n" + "=" * 80)
    print("INTEGRATION SUMMARY")
    print("=" * 80)
    print(f"\nPoppenberg et al. 2023 Study:")
    print(f"  - {summary['poppenberg_study']['n_samples']} samples, {summary['poppenberg_study']['n_genes']} genes")
    print(f"  - Model AUC: {summary['poppenberg_study']['model_auc']}")
    
    print(f"\nGene Statistics:")
    print(f"  - Total genes: {summary['gene_statistics']['total_genes']}")
    print(f"  - Decreased in growing: {summary['gene_statistics']['decreased_in_growing']}")
    print(f"  - Increased in growing: {summary['gene_statistics']['increased_in_growing']}")
    print(f"  - High fold-change (>1.0): {summary['gene_statistics']['high_fold_change']}")
    print(f"  - Found in GSE122897: {summary['gene_statistics']['found_in_gse122897']}")
    print(f"  - Consistent direction: {summary['gene_statistics']['consistent_direction']}")
    print(f"  - High priority (score >= 10): {summary['gene_statistics']['high_priority_genes']}")
    
    print(f"\nTop Biomarker Candidates (Priority Score >= 10):")
    for gene in summary['recommendations']['top_biomarker_candidates'][:10]:
        print(f"  - {gene}")
    
    print(f"\nGenes for Validation (found in GSE122897 with consistent direction):")
    for gene in summary['recommendations']['validation_priorities'][:10]:
        print(f"  - {gene}")
    
    print("\n" + "=" * 80)
    print("INTEGRATION COMPLETE!")
    print("=" * 80)
    print(f"\nAll results saved to: {output_dir}/")

if __name__ == "__main__":
    main()






