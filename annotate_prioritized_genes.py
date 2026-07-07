"""
Annotate Top Prioritized Genes with Symbols, Names, and Functions
Uses MyGene.info API and cross-references with DE results
"""
import pandas as pd
import numpy as np
import requests
import time
import os

def annotate_genes_api(gene_ids):
    """
    Annotate genes using MyGene.info API
    """
    print(f"\nAnnotating {len(gene_ids)} genes using MyGene.info API...")
    
    # Convert to strings and limit for API
    query_ids = [str(gid).rstrip('.0') for gid in gene_ids[:100]]
    
    url = "https://mygene.info/v3/query"
    
    annotations = {}
    
    # Try multiple ID types
    id_types = ['entrezgene', 'ensembl.gene', 'symbol']
    
    for id_type in id_types:
        try:
            params = {
                'q': ','.join(query_ids),
                'scopes': id_type,
                'fields': 'symbol,name,entrezgene,ensembl.gene,alias,summary,type_of_gene,go,pathway,omim',
                'species': 'human'
            }
            
            response = requests.get(url, params=params, timeout=30)
            if response.status_code == 200:
                results = response.json()
                for item in results:
                    if 'query' in item:
                        query_id = str(item['query']).rstrip('.0')
                        symbol = item.get('symbol', query_id)
                        
                        # Extract GO terms
                        go_terms = []
                        if 'go' in item and isinstance(item['go'], dict):
                            for go_type in ['BP', 'CC', 'MF']:
                                if go_type in item['go']:
                                    terms = item['go'][go_type]
                                    if isinstance(terms, list):
                                        go_terms.extend([t.get('term', '') for t in terms[:5]])
                        
                        # Extract pathways
                        pathways = []
                        if 'pathway' in item:
                            if isinstance(item['pathway'], dict):
                                if 'kegg' in item['pathway']:
                                    pathways.extend([p.get('name', '') for p in item['pathway']['kegg'][:5]])
                                if 'reactome' in item['pathway']:
                                    pathways.extend([p.get('name', '') for p in item['pathway']['reactome'][:5]])
                        
                        annotations[query_id] = {
                            'Gene_ID': query_id,
                            'Symbol': symbol,
                            'Name': item.get('name', ''),
                            'Entrez_ID': item.get('entrezgene', query_id),
                            'Ensembl_ID': item.get('ensembl', {}).get('gene', '') if isinstance(item.get('ensembl'), dict) else '',
                            'Aliases': ', '.join(item.get('alias', [])) if item.get('alias') else '',
                            'Gene_Type': item.get('type_of_gene', ''),
                            'Summary': item.get('summary', '')[:500] if item.get('summary') else '',
                            'GO_Terms': '; '.join(go_terms[:10]),
                            'Pathways': '; '.join(pathways[:10])
                        }
                
                if annotations:
                    print(f"  [OK] Annotated {len(annotations)} genes using {id_type}")
                    break
        except Exception as e:
            print(f"  [WARNING] Annotation failed for {id_type}: {e}")
            continue
        
        time.sleep(1)  # Rate limiting
    
    return annotations

def main():
    """Main annotation pipeline"""
    
    print("=" * 80)
    print("ANNOTATING PRIORITIZED GENES")
    print("=" * 80)
    
    # Load prioritized genes
    prioritized_file = "results/prioritized_genes_summary.tsv"
    if not os.path.exists(prioritized_file):
        print(f"[ERROR] Prioritized genes file not found: {prioritized_file}")
        return
    
    prioritized_df = pd.read_csv(prioritized_file, sep='\t')
    print(f"\nLoaded {len(prioritized_df)} prioritized genes")
    
    # Get top 50 genes for annotation
    top_50 = prioritized_df.head(50).copy()
    gene_ids = top_50['Gene_ID'].astype(str).str.rstrip('.0').tolist()
    
    print(f"\nAnnotating top {len(gene_ids)} prioritized genes...")
    
    # Annotate genes
    annotations = annotate_genes_api(gene_ids)
    
    # Merge annotations with prioritized genes
    annotated_df = top_50.copy()
    annotated_df['Gene_ID_Str'] = annotated_df['Gene_ID'].astype(str).str.rstrip('.0')
    
    # Add annotation columns
    for idx, row in annotated_df.iterrows():
        gene_id_str = str(row['Gene_ID']).rstrip('.0')
        annotation = annotations.get(gene_id_str, {})
        
        annotated_df.at[idx, 'Gene_Symbol'] = annotation.get('Symbol', gene_id_str)
        annotated_df.at[idx, 'Gene_Name'] = annotation.get('Name', '')
        annotated_df.at[idx, 'Entrez_ID'] = annotation.get('Entrez_ID', gene_id_str)
        annotated_df.at[idx, 'Ensembl_ID'] = annotation.get('Ensembl_ID', '')
        annotated_df.at[idx, 'Gene_Type'] = annotation.get('Gene_Type', '')
        annotated_df.at[idx, 'Aliases'] = annotation.get('Aliases', '')
        annotated_df.at[idx, 'Summary'] = annotation.get('Summary', '')
        annotated_df.at[idx, 'GO_Terms'] = annotation.get('GO_Terms', '')
        annotated_df.at[idx, 'Pathways'] = annotation.get('Pathways', '')
    
    # Reorder columns
    cols = ['Gene_ID', 'Gene_Symbol', 'Gene_Name', 'Entrez_ID', 'Ensembl_ID',
            'Gene_Type', 'Aliases', 'Log2FC', 'P_adj_BH', 'Priority_Score',
            'In_GWAS_Locus', 'GWAS_Association', 'GO_Terms', 'Pathways', 'Summary']
    cols = [c for c in cols if c in annotated_df.columns]
    annotated_df = annotated_df[cols]
    
    # Save annotated results
    output_file = "results/prioritized_genes_annotated.tsv"
    annotated_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved annotated prioritized genes to {output_file}")
    
    # Print top annotated genes
    print(f"\nTop 20 Annotated Prioritized Genes:")
    print(f"{'Rank':<6} {'Symbol':<15} {'Name':<40} {'Log2FC':<10} {'P_adj_BH':<12}")
    print("-" * 100)
    
    for idx, (i, row) in enumerate(annotated_df.head(20).iterrows(), 1):
        symbol = row.get('Gene_Symbol', str(row['Gene_ID']).rstrip('.0'))
        name = row.get('Gene_Name', '')
        if len(name) > 38:
            name = name[:35] + "..."
        log2fc = row['Log2FC']
        p_adj = row['P_adj_BH']
        
        print(f"{idx:<6} {symbol:<15} {name:<40} {log2fc:>9.3f}  {p_adj:>11.4f}")
    
    print("\n" + "=" * 80)
    print("ANNOTATION COMPLETE")
    print("=" * 80)

if __name__ == "__main__":
    main()






