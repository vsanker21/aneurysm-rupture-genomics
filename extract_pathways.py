"""
Extract Pathways from Annotated Genes (FDR < 0.1)
Summarizes GO terms and pathways from annotated DE genes
"""
import pandas as pd
import numpy as np
import os
from collections import Counter

def extract_pathways_from_annotated():
    """Extract pathways from annotated genes file"""
    
    print("=" * 80)
    print("PATHWAY EXTRACTION FROM ANNOTATED GENES")
    print("=" * 80)
    
    # Load annotated genes file
    annotated_file = "results/gse122897_annotated_fdr_0.1.tsv"
    
    if not os.path.exists(annotated_file):
        print(f"[ERROR] Annotated genes file not found: {annotated_file}")
        return None
    
    df = pd.read_csv(annotated_file, sep='\t')
    print(f"\nLoaded annotated genes: {len(df)} genes")
    print(f"Columns: {list(df.columns)}")
    
    # Extract pathways from GO_Terms and Pathways columns
    all_go_terms = []
    all_pathways = []
    gene_pathways = []
    
    go_col = None
    pathway_col = None
    
    for col in df.columns:
        if 'GO' in col.upper() or 'GO_TERMS' in col.upper():
            go_col = col
        if 'PATHWAY' in col.upper() and 'COUNT' not in col.upper():
            pathway_col = col
    
    print(f"\nUsing columns:")
    print(f"  GO Terms: {go_col if go_col else 'Not found'}")
    print(f"  Pathways: {pathway_col if pathway_col else 'Not found'}")
    
    if not go_col and not pathway_col:
        print("[WARNING] No GO terms or pathways columns found")
        print("  Available columns with 'GO' or 'PATHWAY':")
        for col in df.columns:
            if 'go' in col.lower() or 'pathway' in col.lower():
                print(f"    {col}")
        return None
    
    # Extract GO terms
    if go_col:
        for idx, row in df.iterrows():
            go_terms = str(row.get(go_col, ''))
            if pd.notna(go_terms) and go_terms and go_terms != 'nan':
                # Split by semicolon or comma
                terms = [t.strip() for t in go_terms.replace(';', ',').split(',') if t.strip()]
                all_go_terms.extend(terms)
                
                gene_id = str(row.get('Gene', row.get('Gene_ID', idx)))
                gene_symbol = str(row.get('Gene_Symbol', gene_id))
                
                for term in terms:
                    gene_pathways.append({
                        'Gene_ID': gene_id,
                        'Gene_Symbol': gene_symbol,
                        'Pathway/Term': term,
                        'Type': 'GO_Term',
                        'Log2FC': row.get('Log2FC', np.nan),
                        'P_adj_BH': row.get('P_adj_BH', np.nan)
                    })
    
    # Extract pathways
    if pathway_col:
        for idx, row in df.iterrows():
            pathways = str(row.get(pathway_col, ''))
            if pd.notna(pathways) and pathways and pathways != 'nan':
                # Split by semicolon or comma
                paths = [p.strip() for p in pathways.replace(';', ',').split(',') if p.strip()]
                all_pathways.extend(paths)
                
                gene_id = str(row.get('Gene', row.get('Gene_ID', idx)))
                gene_symbol = str(row.get('Gene_Symbol', gene_id))
                
                for path in paths:
                    gene_pathways.append({
                        'Gene_ID': gene_id,
                        'Gene_Symbol': gene_symbol,
                        'Pathway/Term': path,
                        'Type': 'Pathway',
                        'Log2FC': row.get('Log2FC', np.nan),
                        'P_adj_BH': row.get('P_adj_BH', np.nan)
                    })
    
    # Create pathway summary
    print(f"\n\nPathway Extraction Summary:")
    print(f"  Total GO terms found: {len(all_go_terms)}")
    print(f"  Total pathways found: {len(all_pathways)}")
    print(f"  Unique GO terms: {len(set(all_go_terms))}")
    print(f"  Unique pathways: {len(set(all_pathways))}")
    
    # Count GO terms
    if all_go_terms:
        go_counter = Counter(all_go_terms)
        top_go_terms = go_counter.most_common(30)
        
        print(f"\n\nTop 30 GO Terms (by frequency):")
        print(f"{'Rank':<6} {'GO Term':<60} {'Count':<8} {'Genes'}")
        print("-" * 100)
        
        for idx, (term, count) in enumerate(top_go_terms, 1):
            # Find genes with this term
            genes_with_term = [gp['Gene_Symbol'] for gp in gene_pathways 
                             if gp['Pathway/Term'] == term and gp['Type'] == 'GO_Term']
            genes_str = ', '.join(genes_with_term[:5])  # First 5 genes
            if len(genes_with_term) > 5:
                genes_str += f" ... (+{len(genes_with_term)-5} more)"
            
            term_display = term[:58] if len(term) > 58 else term
            print(f"{idx:<6} {term_display:<60} {count:<8} {genes_str}")
    
    # Count pathways
    if all_pathways:
        pathway_counter = Counter(all_pathways)
        top_pathways = pathway_counter.most_common(30)
        
        print(f"\n\nTop 30 Pathways (by frequency):")
        print(f"{'Rank':<6} {'Pathway':<60} {'Count':<8} {'Genes'}")
        print("-" * 100)
        
        for idx, (pathway, count) in enumerate(top_pathways, 1):
            # Find genes with this pathway
            genes_with_pathway = [gp['Gene_Symbol'] for gp in gene_pathways 
                                 if gp['Pathway/Term'] == pathway and gp['Type'] == 'Pathway']
            genes_str = ', '.join(genes_with_pathway[:5])  # First 5 genes
            if len(genes_with_pathway) > 5:
                genes_str += f" ... (+{len(genes_with_pathway)-5} more)"
            
            pathway_display = pathway[:58] if len(pathway) > 58 else pathway
            print(f"{idx:<6} {pathway_display:<60} {count:<8} {genes_str}")
    
    # Save pathway summary
    if gene_pathways:
        pathways_df = pd.DataFrame(gene_pathways)
        output_file = "results/pathways_extracted_fdr_0.1.tsv"
        pathways_df.to_csv(output_file, sep='\t', index=False)
        print(f"\n[OK] Saved pathway details to {output_file}")
    
    # Save summary tables
    if all_go_terms:
        go_summary = pd.DataFrame([
            {'GO_Term': term, 'Count': count, 'Percentage': count/len(df)*100}
            for term, count in Counter(all_go_terms).most_common(50)
        ])
        go_output = "results/go_terms_summary_fdr_0.1.tsv"
        go_summary.to_csv(go_output, sep='\t', index=False)
        print(f"[OK] Saved GO terms summary to {go_output}")
    
    if all_pathways:
        pathway_summary = pd.DataFrame([
            {'Pathway': pathway, 'Count': count, 'Percentage': count/len(df)*100}
            for pathway, count in Counter(all_pathways).most_common(50)
        ])
        pathway_output = "results/pathways_summary_fdr_0.1.tsv"
        pathway_summary.to_csv(pathway_output, sep='\t', index=False)
        print(f"[OK] Saved pathways summary to {pathway_output}")
    
    return gene_pathways

def main():
    """Main function"""
    pathways = extract_pathways_from_annotated()
    
    print("\n" + "=" * 80)
    print("PATHWAY EXTRACTION COMPLETE")
    print("=" * 80)
    
    if pathways:
        annotated_count = len(pd.read_csv('results/gse122897_annotated_fdr_0.1.tsv', sep='\t'))
        print(f"\n[NOTE] Pathways extracted from {annotated_count} annotated genes")
        print("  Next steps:")
        print("  1. Review top pathways for biological relevance")
        print("  2. Map pathways to known IA biology")
        print("  3. Identify key biological themes")

if __name__ == "__main__":
    main()

