"""
Comprehensive Integration of All Analysis Results
Combines DE analysis, gene annotations, pathway enrichment, ML, and literature
to create a prioritized gene list
"""
import pandas as pd
import numpy as np
import os
from typing import Dict, List, Optional

def load_de_results():
    """Load differential expression results"""
    files_to_try = [
        "results/gse159610_top_candidate_genes.tsv",
        "results/gse159610_substantial_fc_genes.tsv",
        "results/gse159610_full_de_results.tsv"
    ]
    
    for file in files_to_try:
        if os.path.exists(file):
            df = pd.read_csv(file, sep='\t')
            print(f"Loaded DE results from: {file}")
            print(f"  Genes: {len(df):,}")
            return df
    
    print("[ERROR] No DE results found")
    return None

def load_annotated_genes():
    """Load genes with annotations (symbols, etc.)"""
    files_to_try = [
        "results/gse159610_top_candidate_genes_annotated.tsv",
        "results/gse159610_substantial_fc_genes_annotated.tsv",
        "results/gse159610_top_candidate_genes.tsv"
    ]
    
    for file in files_to_try:
        if os.path.exists(file):
            df = pd.read_csv(file, sep='\t')
            print(f"Loaded annotated genes from: {file}")
            return df
    
    return None

def load_literature_annotations():
    """Load literature review annotations"""
    lit_file = "results/gse159610_top_genes_literature_annotated.tsv"
    
    if os.path.exists(lit_file):
        df = pd.read_csv(lit_file, sep='\t')
        print(f"Loaded literature annotations: {len(df):,} genes")
        return df
    
    return None

def load_ml_results():
    """Load machine learning feature importance"""
    ml_file = "results/blood_ml_feature_importance.tsv"
    
    if os.path.exists(ml_file):
        df = pd.read_csv(ml_file, sep='\t')
        print(f"Loaded ML feature importance: {len(df):,} features")
        return df
    
    return None

def load_pathway_results():
    """Load pathway enrichment results"""
    pathway_files = [
        "results/blood_pathway_enrichment_symbols_fdr_0.05.tsv",
        "results/blood_pathway_enrichment_symbols_fdr_0.1.tsv",
        "results/blood_pathway_enrichment_symbols_all.tsv"
    ]
    
    for file in pathway_files:
        if os.path.exists(file):
            df = pd.read_csv(file, sep='\t')
            print(f"Loaded pathway enrichment from: {file}")
            return df
    
    return None

def load_deseq2_results():
    """Load DESeq2 results if available"""
    deseq2_file = "results/gse159610_deseq2_results.tsv"
    
    if os.path.exists(deseq2_file):
        df = pd.read_csv(deseq2_file, sep='\t')
        print(f"Loaded DESeq2 results: {len(df):,} genes")
        return df
    
    return None

def calculate_priority_score(gene_row, ml_importance=None, has_literature=False, 
                            in_pathway=False, deseq2_padj=None):
    """Calculate comprehensive priority score for a gene"""
    score = 0.0
    
    # Effect size (Log2FC)
    abs_fc = abs(gene_row.get('Log2FC', 0))
    if abs_fc > 2.0:
        score += 10
    elif abs_fc > 1.585:
        score += 7
    elif abs_fc > 1.0:
        score += 5
    elif abs_fc > 0.585:
        score += 3
    
    # Statistical significance (p-value)
    p_val = gene_row.get('P_value', 1.0)
    if p_val < 0.001:
        score += 5
    elif p_val < 0.01:
        score += 3
    elif p_val < 0.05:
        score += 1
    
    # ML feature importance
    if ml_importance is not None:
        score += ml_importance * 5
    
    # Literature support
    if has_literature:
        score += 3
    
    # Pathway enrichment
    if in_pathway:
        score += 2
    
    # DESeq2 significance
    if deseq2_padj is not None and deseq2_padj < 0.1:
        score += 5
    elif deseq2_padj is not None and deseq2_padj < 0.2:
        score += 2
    
    # Expression level (higher is better for reliability)
    mean_expr = (gene_row.get('Mean_Aneurysm', 0) + gene_row.get('Mean_Control', 0)) / 2
    if mean_expr > 100:
        score += 2
    elif mean_expr > 50:
        score += 1
    
    return score

def integrate_all_results():
    """Integrate all analysis results"""
    print("=" * 80)
    print("COMPREHENSIVE GENE PRIORITIZATION")
    print("=" * 80)
    
    # Load all data sources
    print("\nLoading data sources...")
    de_df = load_de_results()
    if de_df is None:
        return None
    
    annotated_df = load_annotated_genes()
    lit_df = load_literature_annotations()
    ml_df = load_ml_results()
    pathway_df = load_pathway_results()
    deseq2_df = load_deseq2_results()
    
    # Start with DE results
    integrated_df = de_df.copy()
    
    # Add gene symbols if available
    if annotated_df is not None:
        if 'Gene_Symbol' in annotated_df.columns:
            symbol_map = dict(zip(annotated_df['Gene_ID'], annotated_df['Gene_Symbol']))
            integrated_df['Gene_Symbol'] = integrated_df['Gene_ID'].map(symbol_map)
        
        if 'Gene_Name' in annotated_df.columns:
            name_map = dict(zip(annotated_df['Gene_ID'], annotated_df['Gene_Name']))
            integrated_df['Gene_Name'] = integrated_df['Gene_ID'].map(name_map)
    
    # Add literature annotations
    if lit_df is not None:
        lit_genes = set(lit_df['Gene_ID'].astype(str))
        integrated_df['Has_Literature'] = integrated_df['Gene_ID'].astype(str).isin(lit_genes)
        
        if 'PubMed_Count' in lit_df.columns:
            pubmed_map = dict(zip(lit_df['Gene_ID'].astype(str), lit_df['PubMed_Count']))
            integrated_df['PubMed_Count'] = integrated_df['Gene_ID'].astype(str).map(pubmed_map).fillna(0)
    else:
        integrated_df['Has_Literature'] = False
        integrated_df['PubMed_Count'] = 0
    
    # Add ML importance
    if ml_df is not None:
        ml_map = dict(zip(ml_df['Gene_ID'].astype(str), ml_df['Importance']))
        integrated_df['ML_Importance'] = integrated_df['Gene_ID'].astype(str).map(ml_map).fillna(0)
    else:
        integrated_df['ML_Importance'] = 0
    
    # Add pathway information
    if pathway_df is not None:
        # Extract genes from pathway results
        pathway_genes = set()
        if 'Genes' in pathway_df.columns:
            for genes_str in pathway_df['Genes'].dropna():
                pathway_genes.update(genes_str.split(', '))
        
        integrated_df['In_Pathway'] = integrated_df.get('Gene_Symbol', '').astype(str).isin(pathway_genes)
    else:
        integrated_df['In_Pathway'] = False
    
    # Add DESeq2 results
    if deseq2_df is not None:
        deseq2_map = dict(zip(deseq2_df['Gene_ID'].astype(str), deseq2_df['padj']))
        integrated_df['DESeq2_padj'] = integrated_df['Gene_ID'].astype(str).map(deseq2_map)
        
        deseq2_fc_map = dict(zip(deseq2_df['Gene_ID'].astype(str), deseq2_df['log2FoldChange']))
        integrated_df['DESeq2_Log2FC'] = integrated_df['Gene_ID'].astype(str).map(deseq2_fc_map)
    else:
        integrated_df['DESeq2_padj'] = np.nan
        integrated_df['DESeq2_Log2FC'] = np.nan
    
    # Calculate priority scores
    print("\nCalculating priority scores...")
    integrated_df['Priority_Score'] = integrated_df.apply(
        lambda row: calculate_priority_score(
            row,
            ml_importance=row.get('ML_Importance', 0),
            has_literature=row.get('Has_Literature', False),
            in_pathway=row.get('In_Pathway', False),
            deseq2_padj=row.get('DESeq2_padj')
        ), axis=1
    )
    
    # Sort by priority score
    integrated_df = integrated_df.sort_values('Priority_Score', ascending=False)
    
    return integrated_df

def create_prioritized_lists(integrated_df):
    """Create different prioritized gene lists"""
    print("\n" + "=" * 80)
    print("CREATING PRIORITIZED GENE LISTS")
    print("=" * 80)
    
    os.makedirs("results", exist_ok=True)
    
    # Top 50 overall
    top_50 = integrated_df.head(50).copy()
    output_file = "results/prioritized_genes_top50.tsv"
    top_50.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved top 50 prioritized genes to {output_file}")
    
    # Top 20 with literature support
    with_lit = integrated_df[integrated_df['Has_Literature'] == True].head(20).copy()
    if len(with_lit) > 0:
        output_file = "results/prioritized_genes_with_literature.tsv"
        with_lit.to_csv(output_file, sep='\t', index=False)
        print(f"[OK] Saved {len(with_lit)} genes with literature support to {output_file}")
    
    # Top 20 with ML importance
    high_ml = integrated_df[integrated_df['ML_Importance'] > 0].nlargest(20, 'ML_Importance').copy()
    if len(high_ml) > 0:
        output_file = "results/prioritized_genes_ml_important.tsv"
        high_ml.to_csv(output_file, sep='\t', index=False)
        print(f"[OK] Saved {len(high_ml)} ML-important genes to {output_file}")
    
    # Top 20 with pathway enrichment
    in_pathway = integrated_df[integrated_df['In_Pathway'] == True].head(20).copy()
    if len(in_pathway) > 0:
        output_file = "results/prioritized_genes_in_pathways.tsv"
        in_pathway.to_csv(output_file, sep='\t', index=False)
        print(f"[OK] Saved {len(in_pathway)} pathway-enriched genes to {output_file}")
    
    # Top 20 upregulated
    upregulated = integrated_df[integrated_df['Log2FC'] > 0].head(20).copy()
    output_file = "results/prioritized_genes_upregulated.tsv"
    upregulated.to_csv(output_file, sep='\t', index=False)
    print(f"[OK] Saved top 20 upregulated genes to {output_file}")
    
    # Top 20 downregulated
    downregulated = integrated_df[integrated_df['Log2FC'] < 0].head(20).copy()
    output_file = "results/prioritized_genes_downregulated.tsv"
    downregulated.to_csv(output_file, sep='\t', index=False)
    print(f"[OK] Saved top 20 downregulated genes to {output_file}")

def generate_summary_report(integrated_df):
    """Generate comprehensive summary report"""
    print("\n" + "=" * 80)
    print("GENERATING SUMMARY REPORT")
    print("=" * 80)
    
    report = []
    report.append("=" * 80)
    report.append("COMPREHENSIVE GENE PRIORITIZATION SUMMARY")
    report.append("=" * 80)
    report.append("")
    
    report.append(f"Total genes analyzed: {len(integrated_df):,}")
    report.append("")
    
    # Statistics
    report.append("Priority Score Statistics:")
    report.append(f"  Mean: {integrated_df['Priority_Score'].mean():.2f}")
    report.append(f"  Median: {integrated_df['Priority_Score'].median():.2f}")
    report.append(f"  Max: {integrated_df['Priority_Score'].max():.2f}")
    report.append("")
    
    # Gene symbols
    if 'Gene_Symbol' in integrated_df.columns:
        with_symbols = integrated_df['Gene_Symbol'].notna().sum()
        report.append(f"Genes with symbols: {with_symbols:,} ({with_symbols/len(integrated_df)*100:.1f}%)")
    
    # Literature
    with_lit = integrated_df['Has_Literature'].sum()
    report.append(f"Genes with literature support: {with_lit:,} ({with_lit/len(integrated_df)*100:.1f}%)")
    
    # ML importance
    with_ml = (integrated_df['ML_Importance'] > 0).sum()
    report.append(f"Genes with ML importance: {with_ml:,} ({with_ml/len(integrated_df)*100:.1f}%)")
    
    # Pathways
    in_pathway = integrated_df['In_Pathway'].sum()
    report.append(f"Genes in enriched pathways: {in_pathway:,} ({in_pathway/len(integrated_df)*100:.1f}%)")
    
    # DESeq2
    if 'DESeq2_padj' in integrated_df.columns:
        with_deseq2 = integrated_df['DESeq2_padj'].notna().sum()
        deseq2_sig = (integrated_df['DESeq2_padj'] < 0.1).sum()
        report.append(f"Genes with DESeq2 results: {with_deseq2:,}")
        report.append(f"Genes with DESeq2 padj < 0.1: {deseq2_sig:,}")
    
    report.append("")
    report.append("Top 20 Prioritized Genes:")
    report.append("-" * 80)
    
    top_20 = integrated_df.head(20)
    for idx, row in top_20.iterrows():
        symbol = row.get('Gene_Symbol', row.get('Gene_ID', 'Unknown'))
        # Ensure symbol is a string (handle NaN values)
        if pd.isna(symbol):
            symbol = str(row.get('Gene_ID', 'Unknown'))
        else:
            symbol = str(symbol)
        fc = row.get('Log2FC', 0)
        score = row.get('Priority_Score', 0)
        p_val = row.get('P_value', 1.0)
        
        report.append(f"{idx+1:2d}. {symbol:20s} | Log2FC: {fc:7.3f} | P: {p_val:.4f} | Score: {score:.1f}")
    
    report.append("")
    report.append("=" * 80)
    
    # Save report
    report_text = "\n".join(report)
    output_file = "results/comprehensive_gene_prioritization_report.txt"
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
    output_file = "results/comprehensive_prioritized_genes.tsv"
    integrated_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved full integrated results to {output_file}")
    
    # Create prioritized lists
    create_prioritized_lists(integrated_df)
    
    # Generate summary report
    generate_summary_report(integrated_df)
    
    print("\n" + "=" * 80)
    print("INTEGRATION COMPLETE!")
    print("=" * 80)

if __name__ == "__main__":
    main()

