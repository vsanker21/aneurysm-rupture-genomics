"""
Gene Function Annotation for Top DE Genes
Attempts to retrieve gene names and functions from public databases
"""
import pandas as pd
import numpy as np
import requests
import time
import os

def load_de_results():
    """Load differential expression results"""
    de_file = "results/gse122897_differential_expression.tsv"
    if not os.path.exists(de_file):
        print(f"[ERROR] DE results file not found: {de_file}")
        return None
    
    df = pd.read_csv(de_file, sep='\t')
    return df

def annotate_genes_mygene(gene_ids):
    """
    Annotate genes using MyGene.info API with comprehensive annotations
    """
    print(f"\nAnnotating {len(gene_ids)} genes using MyGene.info...")
    
    # MyGene.info API endpoint
    url = "https://mygene.info/v3/query"
    
    # Convert to list of strings
    query_ids = [str(gid) for gid in gene_ids[:100]]  # Limit to 100
    
    # Try different ID types
    id_types = ['entrezgene', 'ensembl.gene', 'symbol']
    
    annotations = {}
    
    # Comprehensive fields to retrieve
    fields = 'symbol,name,entrezgene,ensembl.gene,summary,alias,type_of_gene,go,pathway,kegg,uniprot,omim,pharmgkb'
    
    for id_type in id_types:
        try:
            params = {
                'q': ','.join(query_ids),
                'scopes': id_type,
                'fields': fields,
                'species': 'human'
            }
            
            response = requests.get(url, params=params, timeout=30)
            if response.status_code == 200:
                results = response.json()
                for item in results:
                    if 'query' in item and 'symbol' in item:
                        query_id = str(item['query'])
                        
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
                        
                        # Extract OMIM (disease associations)
                        omim_diseases = []
                        if 'omim' in item and isinstance(item['omim'], dict):
                            omim_diseases = [item['omim'].get('name', '')] if 'name' in item['omim'] else []
                        
                        annotations[query_id] = {
                            'Symbol': item.get('symbol', query_id),
                            'Name': item.get('name', ''),
                            'Entrez_ID': item.get('entrezgene', ''),
                            'Ensembl_ID': item.get('ensembl', {}).get('gene', '') if isinstance(item.get('ensembl'), dict) else '',
                            'Aliases': ', '.join(item.get('alias', [])) if item.get('alias') else '',
                            'Gene_Type': item.get('type_of_gene', ''),
                            'Summary': item.get('summary', '')[:300] if item.get('summary') else '',
                            'GO_Terms': '; '.join(go_terms[:10]),
                            'Pathways': '; '.join(pathways[:10]),
                            'OMIM_Disease': '; '.join(omim_diseases) if omim_diseases else '',
                            'UniProt_ID': item.get('uniprot', {}).get('Swiss-Prot', '') if isinstance(item.get('uniprot'), dict) else ''
                        }
                if annotations:
                    print(f"[OK] Annotated {len(annotations)} genes using {id_type}")
                    break
        except Exception as e:
            print(f"[WARNING] Annotation failed for {id_type}: {e}")
            continue
        
        time.sleep(1)  # Rate limiting
    
    return annotations

def create_annotated_results_table(df, fdr_thresholds=[0.05, 0.1, 0.2]):
    """
    Create annotated results tables for different FDR thresholds
    """
    print("=" * 80)
    print("GENE FUNCTION ANNOTATION")
    print("=" * 80)
    
    # Get top genes to annotate (limit to top 100 per threshold)
    all_genes_to_annotate = set()
    for threshold in fdr_thresholds:
        sig_genes = df[df['P_adj_BH'] < threshold].copy()
        if len(sig_genes) > 0:
            top_genes = sig_genes.nsmallest(100, 'P_adj_BH')['Gene'].tolist()
            all_genes_to_annotate.update(top_genes)
    
    print(f"\nAnnotating {len(all_genes_to_annotate)} unique genes...")
    
    # Annotate genes
    annotations = annotate_genes_mygene(list(all_genes_to_annotate))
    
    # Create annotation mapping
    annotation_map = {}
    for gene_id, info in annotations.items():
        annotation_map[str(gene_id)] = info
    
    # Create annotated results for each threshold
    for threshold in fdr_thresholds:
        sig_genes = df[df['P_adj_BH'] < threshold].copy()
        
        if len(sig_genes) == 0:
            continue
        
        print(f"\nCreating annotated table for FDR < {threshold} ({len(sig_genes)} genes)...")
        
        # Add annotation columns
        annotated_df = sig_genes.copy()
        annotated_df['Gene_Symbol'] = annotated_df['Gene'].astype(str).map(
            lambda x: annotation_map.get(x, {}).get('Symbol', x)
        )
        annotated_df['Gene_Name'] = annotated_df['Gene'].astype(str).map(
            lambda x: annotation_map.get(x, {}).get('Name', '')
        )
        annotated_df['Entrez_ID'] = annotated_df['Gene'].astype(str).map(
            lambda x: annotation_map.get(x, {}).get('Entrez_ID', '')
        )
        annotated_df['Ensembl_ID'] = annotated_df['Gene'].astype(str).map(
            lambda x: annotation_map.get(x, {}).get('Ensembl_ID', '')
        )
        annotated_df['Gene_Summary'] = annotated_df['Gene'].astype(str).map(
            lambda x: annotation_map.get(x, {}).get('Summary', '')
        )
        annotated_df['Gene_Type'] = annotated_df['Gene'].astype(str).map(
            lambda x: annotation_map.get(x, {}).get('Gene_Type', '')
        )
        annotated_df['Aliases'] = annotated_df['Gene'].astype(str).map(
            lambda x: annotation_map.get(x, {}).get('Aliases', '')
        )
        annotated_df['GO_Terms'] = annotated_df['Gene'].astype(str).map(
            lambda x: annotation_map.get(x, {}).get('GO_Terms', '')
        )
        annotated_df['Pathways'] = annotated_df['Gene'].astype(str).map(
            lambda x: annotation_map.get(x, {}).get('Pathways', '')
        )
        annotated_df['OMIM_Disease'] = annotated_df['Gene'].astype(str).map(
            lambda x: annotation_map.get(x, {}).get('OMIM_Disease', '')
        )
        annotated_df['UniProt_ID'] = annotated_df['Gene'].astype(str).map(
            lambda x: annotation_map.get(x, {}).get('UniProt_ID', '')
        )
        
        # Reorder columns
        cols = ['Gene', 'Gene_Symbol', 'Gene_Name', 'Entrez_ID', 'Ensembl_ID', 'UniProt_ID',
                'Gene_Type', 'Aliases', 'GO_Terms', 'Pathways', 'OMIM_Disease', 'Gene_Summary',
                'Log2FC', 'P_value', 'P_adj_BH', 'Mean_Ruptured', 'Mean_Unruptured']
        cols = [c for c in cols if c in annotated_df.columns]
        cols.extend([c for c in annotated_df.columns if c not in cols])
        annotated_df = annotated_df[cols]
        
        # Save annotated results
        output_file = f"results/gse122897_annotated_fdr_{threshold}.tsv"
        annotated_df.to_csv(output_file, sep='\t', index=False)
        print(f"[OK] Saved annotated results to {output_file}")
        
        # Print top annotated genes
        print(f"\nTop 10 annotated genes:")
        top_10 = annotated_df.head(10)
        for idx, row in top_10.iterrows():
            symbol = row['Gene_Symbol'] if pd.notna(row['Gene_Symbol']) else row['Gene']
            name = row['Gene_Name'] if pd.notna(row['Gene_Name']) else ''
            print(f"  {symbol:15s} {name[:50]:50s} Log2FC={row['Log2FC']:7.3f}  P_adj={row['P_adj_BH']:.3f}")

def main():
    """Main annotation pipeline"""
    import os
    
    # Load DE results
    df = load_de_results()
    if df is None:
        return
    
    # Create annotated results tables
    create_annotated_results_table(df, fdr_thresholds=[0.05, 0.1, 0.2])
    
    print("\n" + "=" * 80)
    print("GENE ANNOTATION COMPLETE")
    print("=" * 80)

if __name__ == "__main__":
    main()

