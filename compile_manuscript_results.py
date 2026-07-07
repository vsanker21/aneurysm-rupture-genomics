"""
Compile Results for Manuscript
Creates a comprehensive summary of all findings
"""
import pandas as pd
import numpy as np
import os
from datetime import datetime

def load_gene_annotations():
    """Load annotated top genes"""
    annotated_file = "results/gse122897_top_genes_annotated.tsv"
    if os.path.exists(annotated_file):
        return pd.read_csv(annotated_file, sep='\t')
    return None

def load_literature_search():
    """Load literature search results if available"""
    lit_file = "results/pubmed_literature_search_results.tsv"
    if os.path.exists(lit_file):
        return pd.read_csv(lit_file, sep='\t')
    return None

def load_blood_integration():
    """Load blood integration results"""
    blood_file = "results/blood_tissue_integration_refined.tsv"
    if os.path.exists(blood_file):
        try:
            df = pd.read_csv(blood_file, sep='\t')
            if len(df) > 0:
                return df
        except:
            pass
    
    # Try alternative files
    for fdr in [0.05, 0.1, 0.2]:
        alt_file = f"results/blood_tissue_integration_fdr_{fdr}.tsv"
        if os.path.exists(alt_file):
            try:
                df = pd.read_csv(alt_file, sep='\t')
                if len(df) > 0:
                    return df
            except:
                continue
    return None

def load_pathway_enrichment():
    """Load pathway enrichment results"""
    pathway_file = "results/enriched_pathways_summary.tsv"
    if os.path.exists(pathway_file):
        try:
            df = pd.read_csv(pathway_file, sep='\t')
            if len(df) > 0:
                return df
        except:
            pass
    
    # Try alternative file
    alt_file = "results/pathway_enrichment_fdr_0.05.tsv"
    if os.path.exists(alt_file):
        try:
            df = pd.read_csv(alt_file, sep='\t')
            if len(df) > 0:
                return df
        except:
            pass
    return None

def create_summary_report():
    """Create comprehensive summary report"""
    
    print("=" * 80)
    print("COMPILING MANUSCRIPT RESULTS")
    print("=" * 80)
    
    # Load all data
    print("\nLoading analysis results...")
    
    annotated_genes = load_gene_annotations()
    literature_results = load_literature_search()
    blood_integration = load_blood_integration()
    pathway_enrichment = load_pathway_enrichment()
    
    # Create summary sections
    summary_sections = []
    
    # 1. Top Genes Summary
    summary_sections.append("=" * 80)
    summary_sections.append("1. TOP DIFFERENTIALLY EXPRESSED GENES")
    summary_sections.append("=" * 80)
    
    if annotated_genes is not None:
        top_10 = annotated_genes.head(10)
        summary_sections.append(f"\nTop 10 genes identified from GSE122897 analysis:")
        summary_sections.append(f"{'Rank':<6} {'Symbol':<15} {'Name':<50} {'Log2FC':<12} {'P_adj':<12}")
        summary_sections.append("-" * 100)
        
        for rank, (idx, row) in enumerate(top_10.iterrows(), 1):
            # Use direct column access for pandas Series
            symbol = str(row['Gene_Symbol']) if 'Gene_Symbol' in row and pd.notna(row['Gene_Symbol']) else str(row.get('Gene', ''))
            name = str(row['Gene_Name']) if 'Gene_Name' in row and pd.notna(row['Gene_Name']) else ''
            if len(name) > 48:
                name = name[:45] + "..."
            
            # Get Log2FC
            log2fc = 0.0
            if 'Log2FC' in row:
                log2fc = float(row['Log2FC']) if pd.notna(row['Log2FC']) else 0.0
            
            # Get P_adj
            p_adj = 1.0
            if 'P_adj_BH' in row:
                p_adj = float(row['P_adj_BH']) if pd.notna(row['P_adj_BH']) else 1.0
            
            summary_sections.append(f"{rank:<6} {symbol:<15} {name:<50} {log2fc:>11.3f}  {p_adj:>11.4e}")
    else:
        summary_sections.append("\n[WARNING] Gene annotations not available")
    
    # 2. Literature Review Summary
    summary_sections.append("\n\n" + "=" * 80)
    summary_sections.append("2. LITERATURE REVIEW")
    summary_sections.append("=" * 80)
    
    if literature_results is not None:
        summary_sections.append(f"\nLiterature search completed for {len(literature_results)} genes")
        
        # Genes with aneurysm-related publications
        aneurysm_genes = literature_results[literature_results['Query_2_Count'] > 0]
        summary_sections.append(f"\nGenes with aneurysm-related publications: {len(aneurysm_genes)}")
        
        for idx, row in aneurysm_genes.iterrows():
            symbol = row['Gene_Symbol']
            count = row['Query_2_Count']
            summary_sections.append(f"  {symbol}: {count} publications")
    else:
        summary_sections.append("\n[INFO] Literature search results not available (may need to run search_pubmed_literature.py)")
    
    # 3. Blood Integration Summary
    summary_sections.append("\n\n" + "=" * 80)
    summary_sections.append("3. BLOOD INTEGRATION - BIOMARKER CANDIDATES")
    summary_sections.append("=" * 80)
    
    if blood_integration is not None:
        top_biomarkers = blood_integration.head(10)
        summary_sections.append(f"\nTop 10 biomarker candidates (tissue + blood concordance):")
        summary_sections.append(f"{'Rank':<6} {'Symbol':<15} {'Tissue_Log2FC':<15} {'Blood_Log2FC':<15} {'Blood_P':<12} {'Score':<12}")
        summary_sections.append("-" * 80)
        
        for rank, (idx, row) in enumerate(top_biomarkers.iterrows(), 1):
            symbol = row.get('Gene_Symbol', str(row.get('Gene_ID', '')))
            tissue_fc = row.get('Log2FC_Tissue', 0)
            blood_fc = row.get('Log2FC_Blood', 0)
            blood_p = row.get('P_value_Blood', 1.0)
            score = row.get('Biomarker_Score', 0)
            summary_sections.append(f"{rank:<6} {symbol:<15} {tissue_fc:>14.3f} {blood_fc:>14.3f} {blood_p:>11.4e} {score:>11.2f}")
    else:
        summary_sections.append("\n[WARNING] Blood integration results not available")
    
    # 4. Pathway Enrichment Summary
    summary_sections.append("\n\n" + "=" * 80)
    summary_sections.append("4. PATHWAY ENRICHMENT ANALYSIS")
    summary_sections.append("=" * 80)
    
    if pathway_enrichment is not None:
        top_pathways = pathway_enrichment.head(20)
        summary_sections.append(f"\nTop enriched pathways:")
        summary_sections.append(f"{'Rank':<6} {'Pathway':<60} {'P_value':<12} {'Gene_Count':<12}")
        summary_sections.append("-" * 100)
        
        for rank, (idx, row) in enumerate(top_pathways.iterrows(), 1):
            pathway = row.get('Pathway_Name', row.get('term_name', 'N/A'))
            if len(pathway) > 58:
                pathway = pathway[:55] + "..."
            p_value = row.get('p_value', row.get('pvalue', 1.0))
            gene_count = row.get('Gene_Count', row.get('intersection_size', 0))
            summary_sections.append(f"{rank:<6} {pathway:<60} {p_value:>11.4e} {gene_count:>11}")
    else:
        summary_sections.append("\n[WARNING] Pathway enrichment results not available")
    
    # 5. Key Findings
    summary_sections.append("\n\n" + "=" * 80)
    summary_sections.append("5. KEY FINDINGS SUMMARY")
    summary_sections.append("=" * 80)
    
    key_findings = []
    
    if annotated_genes is not None:
        key_findings.append(f"- Identified {len(annotated_genes)} top differentially expressed genes")
        upregulated = annotated_genes[annotated_genes['Log2FC'] > 0]
        downregulated = annotated_genes[annotated_genes['Log2FC'] < 0]
        key_findings.append(f"  - Upregulated in ruptured: {len(upregulated)} genes")
        key_findings.append(f"  - Downregulated in ruptured: {len(downregulated)} genes")
    
    if blood_integration is not None:
        significant_blood = blood_integration[blood_integration['P_value_Blood'] < 0.05]
        consistent = blood_integration[blood_integration.get('Direction_Consistent', False) == True]
        key_findings.append(f"- Identified {len(blood_integration)} biomarker candidates from blood integration")
        key_findings.append(f"  - Significant in blood (P < 0.05): {len(significant_blood)} genes")
        key_findings.append(f"  - Direction consistent: {len(consistent)} genes")
    
    if pathway_enrichment is not None:
        key_findings.append(f"- Identified {len(pathway_enrichment)} enriched pathways")
    
    if key_findings:
        summary_sections.append("\n" + "\n".join(key_findings))
    else:
        summary_sections.append("\n[INFO] Key findings summary not available")
    
    # Combine all sections
    full_summary = "\n".join(summary_sections)
    
    # Add header
    header = f"""
================================================================================
CEREBRAL ANEURYSM GENOMICS PROJECT - ANALYSIS SUMMARY
Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
================================================================================

This summary compiles findings from:
1. Differential expression analysis (GSE122897)
2. Gene annotation and mapping
3. Literature review (PubMed)
4. Blood transcriptomics integration
5. Pathway enrichment analysis

"""
    
    full_summary = header + full_summary
    
    # Save summary
    output_file = "results/MANUSCRIPT_SUMMARY.txt"
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(full_summary)
    
    print(f"\n[OK] Saved manuscript summary to {output_file}")
    
    # Print summary to console
    print("\n" + full_summary)
    
    return full_summary

def main():
    """Main compilation pipeline"""
    summary = create_summary_report()
    return summary

if __name__ == "__main__":
    main()

