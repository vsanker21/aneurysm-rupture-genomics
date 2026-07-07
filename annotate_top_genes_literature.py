"""
Annotate Top Candidate Genes with Literature Review
Uses Entrez IDs and gene symbols to search for relevant information
"""
import pandas as pd
import numpy as np
import requests
import json
import time
import os
from typing import Dict, List

def load_annotated_genes():
    """Load genes with annotations (symbols, Entrez IDs)"""
    # Try annotated top candidates first
    annotated_file = "results/gse159610_top_candidate_genes_annotated.tsv"
    
    if os.path.exists(annotated_file):
        df = pd.read_csv(annotated_file, sep='\t')
        print(f"Loaded {len(df):,} annotated top candidate genes")
        return df
    
    # Try annotated substantial FC genes
    fc_annotated = "results/gse159610_substantial_fc_genes_annotated.tsv"
    if os.path.exists(fc_annotated):
        df = pd.read_csv(fc_annotated, sep='\t')
        # Take top 100 by effect score
        df['Effect_Score'] = np.abs(df['Log2FC']) * (-np.log10(df['P_value'] + 1e-300))
        df = df.nlargest(100, 'Effect_Score')
        print(f"Loaded top 100 from annotated substantial FC genes")
        return df
    
    print("[ERROR] No annotated genes file found")
    print("Please run convert_ensembl_to_symbols.py first")
    return None

def get_gene_info_from_mygene(gene_symbol: str) -> Dict:
    """Get gene information from mygene API"""
    if not gene_symbol or pd.isna(gene_symbol):
        return {}
    
    url = f"https://mygene.info/v3/query"
    params = {
        'q': gene_symbol,
        'scopes': 'symbol',
        'fields': 'symbol,name,summary,alias,entrezgene,ensembl',
        'species': 'human',
        'limit': 1
    }
    
    try:
        response = requests.get(url, params=params, timeout=10)
        if response.status_code == 200:
            results = response.json()
            if 'hits' in results and len(results['hits']) > 0:
                hit = results['hits'][0]
                return {
                    'Gene_Symbol': hit.get('symbol', ''),
                    'Gene_Name': hit.get('name', ''),
                    'Summary': hit.get('summary', ''),
                    'Aliases': ', '.join(hit.get('alias', [])),
                    'Entrez_ID': hit.get('entrezgene', '')
                }
    except Exception as e:
        pass
    
    return {}

def search_pubmed_aneurysm(gene_symbol: str, max_results: int = 5) -> List[Dict]:
    """Search PubMed for aneurysm-related publications"""
    if not gene_symbol or pd.isna(gene_symbol):
        return []
    
    # Search query: gene symbol AND (aneurysm OR cerebral aneurysm OR intracranial aneurysm)
    query = f'({gene_symbol}[Title/Abstract]) AND (aneurysm[Title/Abstract] OR "cerebral aneurysm"[Title/Abstract] OR "intracranial aneurysm"[Title/Abstract])'
    
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        'db': 'pubmed',
        'term': query,
        'retmax': max_results,
        'retmode': 'json'
    }
    
    try:
        response = requests.get(url, params=params, timeout=10)
        if response.status_code == 200:
            data = response.json()
            if 'esearchresult' in data and 'idlist' in data['esearchresult']:
                pmids = data['esearchresult']['idlist']
                return [{'PMID': pmid} for pmid in pmids[:max_results]]
    except Exception as e:
        pass
    
    return []

def annotate_genes(df: pd.DataFrame, top_n: int = 50):
    """Annotate top genes with literature and gene information"""
    print("=" * 80)
    print("ANNOTATING TOP GENES WITH LITERATURE REVIEW")
    print("=" * 80)
    
    # Sort by effect score
    if 'Effect_Score' not in df.columns:
        df['Effect_Score'] = np.abs(df['Log2FC']) * (-np.log10(df['P_value'] + 1e-300))
    
    top_genes = df.nlargest(top_n, 'Effect_Score').copy()
    
    print(f"\nAnnotating top {len(top_genes)} genes...")
    print("This may take several minutes...")
    
    annotations = []
    
    for idx, row in top_genes.iterrows():
        gene_symbol = row.get('Gene_Symbol', '')
        gene_id = row.get('Gene_ID', '')
        
        if not gene_symbol or pd.isna(gene_symbol):
            print(f"  Skipping {gene_id} (no symbol)")
            continue
        
        print(f"  Annotating {gene_symbol} ({idx+1}/{len(top_genes)})...")
        
        # Get gene information
        gene_info = get_gene_info_from_mygene(gene_symbol)
        
        # Search PubMed
        pubmed_results = search_pubmed_aneurysm(gene_symbol, max_results=3)
        
        annotation = {
            'Gene_ID': gene_id,
            'Gene_Symbol': gene_symbol,
            'Log2FC': row.get('Log2FC', np.nan),
            'P_value': row.get('P_value', np.nan),
            'Effect_Score': row.get('Effect_Score', np.nan),
            'Gene_Name': gene_info.get('Gene_Name', row.get('Gene_Name', '')),
            'Gene_Summary': gene_info.get('Summary', ''),
            'Gene_Aliases': gene_info.get('Aliases', ''),
            'PubMed_Count': len(pubmed_results),
            'PubMed_PMIDs': ', '.join([p['PMID'] for p in pubmed_results]) if pubmed_results else ''
        }
        
        annotations.append(annotation)
        
        # Rate limiting
        time.sleep(0.5)
    
    if len(annotations) == 0:
        print(f"\n[WARNING] No genes were successfully annotated")
        # Return empty DataFrame with expected columns
        return pd.DataFrame(columns=['Gene_ID', 'Gene_Symbol', 'Log2FC', 'P_value', 'Effect_Score', 
                                    'Gene_Name', 'Gene_Summary', 'Gene_Aliases', 'PubMed_Count', 'PubMed_PMIDs'])
    
    annotations_df = pd.DataFrame(annotations)
    
    print(f"\nAnnotation complete:")
    print(f"  Genes annotated: {len(annotations_df)}")
    if len(annotations_df) > 0 and 'PubMed_Count' in annotations_df.columns:
        print(f"  Genes with PubMed results: {(annotations_df['PubMed_Count'] > 0).sum()}")
    
    return annotations_df

def save_annotations(annotations_df: pd.DataFrame):
    """Save annotated genes"""
    os.makedirs("results", exist_ok=True)
    
    output_file = "results/gse159610_top_genes_literature_annotated.tsv"
    annotations_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved annotated genes to {output_file}")
    
    # Save genes with PubMed results
    genes_with_pubmed = annotations_df[annotations_df['PubMed_Count'] > 0].copy()
    if len(genes_with_pubmed) > 0:
        pubmed_file = "results/gse159610_genes_with_pubmed_results.tsv"
        genes_with_pubmed = genes_with_pubmed.sort_values('PubMed_Count', ascending=False)
        genes_with_pubmed.to_csv(pubmed_file, sep='\t', index=False)
        print(f"[OK] Saved {len(genes_with_pubmed)} genes with PubMed results to {pubmed_file}")
        
        print(f"\nTop genes with aneurysm-related publications:")
        print(genes_with_pubmed.head(10)[['Gene_Symbol', 'Gene_Name', 'Log2FC', 'PubMed_Count']].to_string(index=False))

def main():
    """Main annotation pipeline"""
    # Load annotated genes
    df = load_annotated_genes()
    if df is None:
        return
    
    # Annotate top genes
    annotations_df = annotate_genes(df, top_n=50)
    
    # Save results
    save_annotations(annotations_df)
    
    print("\n" + "=" * 80)
    print("LITERATURE ANNOTATION COMPLETE!")
    print("=" * 80)
    print("\nNext steps:")
    print("1. Review annotated genes in results/gse159610_top_genes_literature_annotated.tsv")
    print("2. Check PubMed results for genes with aneurysm-related publications")
    print("3. Use gene summaries and aliases for further investigation")

if __name__ == "__main__":
    main()

