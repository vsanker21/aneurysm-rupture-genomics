"""
Enhanced Literature Annotation for Multi-Omics Prioritized Genes
Annotates top prioritized genes with PubMed references for aneurysm research
"""
import pandas as pd
import numpy as np
import requests
import time
import os
from typing import Dict, List

def get_gene_info_mygene(gene_symbol: str) -> Dict:
    """Get gene information from MyGene.info API"""
    if not gene_symbol or pd.isna(gene_symbol) or gene_symbol == '':
        return {}
    
    url = "https://mygene.info/v3/query"
    params = {
        'q': str(gene_symbol).strip(),
        'scopes': 'symbol',
        'fields': 'symbol,name,summary,alias,entrezgene,ensembl',
        'species': 'human',
        'limit': 1
    }
    
    try:
        response = requests.get(url, params=params, timeout=10)
        if response.status_code == 200:
            results = response.json()
            if results and 'hits' in results and len(results['hits']) > 0:
                hit = results['hits'][0]
                return {
                    'Gene_Symbol': hit.get('symbol', gene_symbol),
                    'Gene_Name': hit.get('name', ''),
                    'Gene_Summary': hit.get('summary', ''),
                    'Gene_Aliases': ', '.join(hit.get('alias', [])),
                    'Entrez_ID': hit.get('entrezgene', ''),
                    'Ensembl_ID': hit.get('ensembl', {}).get('gene', '') if isinstance(hit.get('ensembl'), dict) else ''
                }
    except Exception as e:
        print(f"    Error fetching gene info for {gene_symbol}: {e}")
    
    return {'Gene_Symbol': gene_symbol}

def search_pubmed_aneurysm(gene_symbol: str, max_results: int = 5) -> Dict:
    """Search PubMed for aneurysm-related publications"""
    if not gene_symbol or pd.isna(gene_symbol) or gene_symbol == '':
        return {'count': 0, 'pmids': [], 'query': ''}
    
    # Single most relevant search strategy (reduced from 3 to 1 for speed)
    queries = [
        f'({gene_symbol}[Title/Abstract]) AND (aneurysm[Title/Abstract] OR "cerebral aneurysm"[Title/Abstract] OR "intracranial aneurysm"[Title/Abstract])'
    ]
    
    all_pmids = []
    
    for query in queries:
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
                    all_pmids.extend(pmids)
        except Exception as e:
            continue
        
        time.sleep(0.2)  # Reduced rate limiting
    
    # Get unique PMIDs
    unique_pmids = list(set(all_pmids))[:max_results]
    
    return {
        'count': len(unique_pmids),
        'pmids': unique_pmids,
        'query': queries[0] if queries else ''
    }

def get_pubmed_details(pmids: List[str], max_details: int = 5) -> List[Dict]:
    """Get details for top PubMed articles"""
    if not pmids or len(pmids) == 0:
        return []
    
    pmids_to_fetch = pmids[:max_details]
    
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        'db': 'pubmed',
        'id': ','.join(pmids_to_fetch),
        'retmode': 'xml',
        'rettype': 'abstract'
    }
    
    try:
        response = requests.get(url, params=params, timeout=15)
        if response.status_code == 200:
            # Simple parsing - in production would use proper XML parser
            content = response.text
            details = []
            for pmid in pmids_to_fetch:
                details.append({'PMID': pmid, 'Available': True})
            return details
    except Exception as e:
        pass
    
    return [{'PMID': pmid, 'Available': False} for pmid in pmids_to_fetch]

def annotate_multi_omics_genes(df, top_n=20):
    """Annotate top multi-omics prioritized genes"""
    print("=" * 80)
    print("LITERATURE ANNOTATION FOR MULTI-OMICS PRIORITIZED GENES")
    print("=" * 80)
    
    # Determine priority column
    priority_col = None
    for col in ['Integrated_Priority_Score', 'Multi_Omics_Priority_Score', 'Priority_Score']:
        if col in df.columns:
            priority_col = col
            break
    
    if priority_col is None:
        print("[WARNING] No priority score column found, using first {top_n} genes")
        top_genes = df.head(top_n).copy()
    else:
        top_genes = df.nlargest(top_n, priority_col).copy()
    
    print(f"\nAnnotating top {len(top_genes)} genes...")
    print("This may take 2-3 minutes due to API rate limiting...")
    print("(Processing approximately 0.5 seconds per gene)")
    
    annotations = []
    
    for idx, (row_idx, row) in enumerate(top_genes.iterrows(), 1):
        # Get gene symbol
        gene_symbol = row.get('Gene_Symbol', '')
        if pd.isna(gene_symbol) or gene_symbol == '':
            # Try alternative columns
            gene_symbol = row.get('Gene', str(row.get('Gene_ID', '')))
        
        if not gene_symbol or gene_symbol == 'nan':
            print(f"  [{idx}/{len(top_genes)}] Skipping row {row_idx} (no gene identifier)")
            continue
        
        print(f"  [{idx}/{len(top_genes)}] Annotating {gene_symbol}...")
        
        # Get gene information
        try:
            gene_info = get_gene_info_mygene(gene_symbol)
            time.sleep(0.3)  # Reduced wait time
        except Exception as e:
            print(f"      Warning: Could not fetch gene info: {e}")
            gene_info = {}
        
        # Search PubMed
        try:
            pubmed_results = search_pubmed_aneurysm(gene_symbol, max_results=5)  # Reduced from 10 to 5
            time.sleep(0.3)  # Reduced wait time
        except Exception as e:
            print(f"      Warning: Could not search PubMed: {e}")
            pubmed_results = {'count': 0, 'pmids': []}
        
        # Combine with original row data
        annotation = row.to_dict()
        annotation.update({
            'Annotated_Gene_Symbol': gene_info.get('Gene_Symbol', gene_symbol),
            'Annotated_Gene_Name': gene_info.get('Gene_Name', ''),
            'Gene_Summary': gene_info.get('Gene_Summary', ''),
            'Gene_Aliases': gene_info.get('Gene_Aliases', ''),
            'PubMed_Count': pubmed_results['count'],
            'PubMed_PMIDs': ', '.join(pubmed_results['pmids']),
            'Has_Literature_Support': pubmed_results['count'] > 0
        })
        
        # Add gene info IDs if missing
        if 'Entrez_ID' not in annotation or pd.isna(annotation.get('Entrez_ID')):
            annotation['Entrez_ID'] = gene_info.get('Entrez_ID', '')
        if 'Ensembl_ID' not in annotation or pd.isna(annotation.get('Ensembl_ID')):
            annotation['Ensembl_ID'] = gene_info.get('Ensembl_ID', '')
        
        annotations.append(annotation)
    
    annotated_df = pd.DataFrame(annotations)
    
    # Summary
    print("\n" + "=" * 80)
    print("ANNOTATION SUMMARY")
    print("=" * 80)
    print(f"Total genes annotated: {len(annotated_df)}")
    print(f"Genes with literature support: {annotated_df['Has_Literature_Support'].sum()}")
    print(f"Total PubMed references: {annotated_df['PubMed_Count'].sum()}")
    print(f"Average references per gene: {annotated_df['PubMed_Count'].mean():.1f}")
    
    # Top genes with literature
    with_lit = annotated_df[annotated_df['Has_Literature_Support'] == True].nlargest(10, 'PubMed_Count')
    if len(with_lit) > 0:
        print("\nTop genes with literature support:")
        for idx, row in with_lit.iterrows():
            symbol = row.get('Annotated_Gene_Symbol', 'Unknown')
            count = row.get('PubMed_Count', 0)
            print(f"  {symbol}: {count} references")
    
    return annotated_df

def main():
    """Main annotation pipeline"""
    # Try to load multi-omics integrated data
    files_to_try = [
        "results/integrated_cross_dataset_genes.tsv",
        "results/multi_omics_comprehensive_prioritized_genes.tsv",
        "results/comprehensive_prioritized_genes.tsv",
        "results/gse159610_top_candidate_genes_annotated.tsv"
    ]
    
    df = None
    for file in files_to_try:
        if os.path.exists(file):
            df = pd.read_csv(file, sep='\t')
            print(f"Loaded data from: {file}")
            print(f"  Genes: {len(df):,}")
            break
    
    if df is None:
        print("[ERROR] No prioritized genes file found!")
        return
    
    # Annotate top genes (reduced to 20 for faster execution)
    annotated_df = annotate_multi_omics_genes(df, top_n=20)
    
    # Save annotated results
    os.makedirs("results", exist_ok=True)
    output_file = "results/multi_omics_prioritized_genes_literature_annotated.tsv"
    annotated_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved annotated results to {output_file}")
    
    # Save top genes with literature
    with_lit = annotated_df[annotated_df['Has_Literature_Support'] == True].copy()
    if len(with_lit) > 0:
        lit_file = "results/top_genes_with_literature_support.tsv"
        with_lit.to_csv(lit_file, sep='\t', index=False)
        print(f"[OK] Saved {len(with_lit)} genes with literature support to {lit_file}")

if __name__ == "__main__":
    main()

