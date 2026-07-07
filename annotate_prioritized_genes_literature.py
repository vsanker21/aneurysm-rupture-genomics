"""
Literature Annotation for Top Prioritized Genes
Searches PubMed for aneurysm-related literature for each gene
"""
import pandas as pd
import numpy as np
import requests
import time
import os
from typing import Dict, List

def search_pubmed_aneurysm(gene_symbol: str, max_results: int = 5):
    """
    Search PubMed for literature linking gene to aneurysm
    Uses NCBI E-utilities API
    """
    if not gene_symbol or gene_symbol == 'nan' or gene_symbol.startswith('ENSG'):
        return []
    
    # Clean gene symbol
    gene_symbol = str(gene_symbol).strip()
    
    # Search query: gene AND (aneurysm OR intracranial aneurysm OR cerebral aneurysm)
    query = f'({gene_symbol}[Gene Name] OR {gene_symbol}[Title/Abstract]) AND (aneurysm OR "intracranial aneurysm" OR "cerebral aneurysm" OR "brain aneurysm")'
    
    # ESearch API
    esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        'db': 'pubmed',
        'term': query,
        'retmax': max_results,
        'retmode': 'json'
    }
    
    try:
        response = requests.get(esearch_url, params=params, timeout=10)
        if response.status_code == 200:
            data = response.json()
            if 'esearchresult' in data and 'idlist' in data['esearchresult']:
                pmids = data['esearchresult']['idlist']
                return pmids[:max_results]
    except Exception as e:
        print(f"    [WARNING] PubMed search failed for {gene_symbol}: {e}")
    
    return []

def get_pubmed_summaries(pmids: List[str], max_results: int = 3):
    """
    Get summaries for PubMed IDs
    """
    if not pmids:
        return []
    
    # Limit to max_results
    pmids = pmids[:max_results]
    
    # EFetch API
    efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        'db': 'pubmed',
        'id': ','.join(pmids),
        'retmode': 'xml',
        'rettype': 'abstract'
    }
    
    try:
        response = requests.get(efetch_url, params=params, timeout=15)
        if response.status_code == 200:
            # Parse XML to extract titles and years
            import xml.etree.ElementTree as ET
            root = ET.fromstring(response.text)
            
            summaries = []
            for article in root.findall('.//PubmedArticle'):
                title_elem = article.find('.//ArticleTitle')
                year_elem = article.find('.//PubDate/Year')
                pmid_elem = article.find('.//PMID')
                
                title = title_elem.text if title_elem is not None else 'N/A'
                year = year_elem.text if year_elem is not None else 'N/A'
                pmid = pmid_elem.text if pmid_elem is not None else 'N/A'
                
                summaries.append({
                    'PMID': pmid,
                    'Title': title[:100] + '...' if len(title) > 100 else title,
                    'Year': year
                })
            
            return summaries
    except Exception as e:
        print(f"    [WARNING] Failed to fetch summaries: {e}")
    
    return []

def annotate_gene_literature(gene_symbol: str, gene_name: str = ''):
    """
    Annotate a single gene with literature information
    """
    print(f"  Annotating {gene_symbol}...", end='', flush=True)
    
    # Search PubMed
    pmids = search_pubmed_aneurysm(gene_symbol, max_results=5)
    
    if pmids:
        summaries = get_pubmed_summaries(pmids, max_results=3)
        
        # Format literature info
        literature_count = len(pmids)
        literature_pmids = ', '.join(pmids[:3])
        literature_titles = '; '.join([f"{s['Title']} ({s['Year']})" for s in summaries[:2]])
        
        print(f" [OK] {literature_count} publications")
        
        return {
            'Literature_Count': literature_count,
            'Literature_PMIDs': literature_pmids,
            'Literature_Titles': literature_titles,
            'Has_Literature': True
        }
    else:
        print(f" [None]")
        return {
            'Literature_Count': 0,
            'Literature_PMIDs': '',
            'Literature_Titles': '',
            'Has_Literature': False
        }

def annotate_prioritized_genes(input_file: str, top_n: int = 20):
    """
    Annotate top prioritized genes with literature
    """
    print("=" * 80)
    print("LITERATURE ANNOTATION FOR PRIORITIZED GENES")
    print("=" * 80)
    
    if not os.path.exists(input_file):
        print(f"[ERROR] Input file not found: {input_file}")
        return None
    
    # Load prioritized genes
    df = pd.read_csv(input_file, sep='\t')
    print(f"\nLoaded {len(df):,} prioritized genes")
    
    # Get top N genes with valid symbols
    top_genes = df.head(top_n).copy()
    
    # Filter to genes with valid symbols
    valid_genes = top_genes[
        (top_genes['Gene_Symbol'].notna()) & 
        (~top_genes['Gene_Symbol'].astype(str).str.startswith('ENSG')) &
        (top_genes['Gene_Symbol'].astype(str) != 'nan')
    ].copy()
    
    print(f"Annotating {len(valid_genes)} genes with valid symbols...")
    
    # Annotate each gene
    literature_data = []
    for idx, row in valid_genes.iterrows():
        gene_symbol = str(row['Gene_Symbol']).strip()
        gene_name = str(row.get('Gene_Name', '')).strip()
        
        lit_info = annotate_gene_literature(gene_symbol, gene_name)
        lit_info['Gene_Symbol'] = gene_symbol
        lit_info['Priority_Score'] = row.get('Multi_Omics_Priority_Score', 0)
        literature_data.append(lit_info)
        
        time.sleep(0.5)  # Rate limiting for API
    
    # Merge with original dataframe
    lit_df = pd.DataFrame(literature_data)
    
    # Merge back to original
    df = df.merge(lit_df[['Gene_Symbol', 'Literature_Count', 'Literature_PMIDs', 
                          'Literature_Titles', 'Has_Literature']], 
                  on='Gene_Symbol', how='left')
    
    # Fill NaN values
    df['Literature_Count'] = df['Literature_Count'].fillna(0).astype(int)
    df['Has_Literature'] = df['Has_Literature'].fillna(False)
    df['Literature_PMIDs'] = df['Literature_PMIDs'].fillna('')
    df['Literature_Titles'] = df['Literature_Titles'].fillna('')
    
    # Save annotated results
    output_file = input_file.replace('.tsv', '_literature_annotated.tsv')
    df.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved literature-annotated genes to: {output_file}")
    
    # Summary
    genes_with_lit = df['Has_Literature'].sum()
    print(f"\nLiterature Annotation Summary:")
    print(f"  Genes with literature: {genes_with_lit}/{len(valid_genes)} ({genes_with_lit/len(valid_genes)*100:.1f}%)")
    print(f"  Total publications found: {df['Literature_Count'].sum()}")
    
    # Show top genes with literature
    print(f"\nTop genes with literature support:")
    with_lit = df[df['Has_Literature'] == True].head(10)
    for idx, row in with_lit.iterrows():
        symbol = row['Gene_Symbol']
        count = row['Literature_Count']
        score = row.get('Multi_Omics_Priority_Score', 0)
        print(f"  {symbol:15s} | Score: {score:5.1f} | {count} publications")
    
    return df

def main():
    """Main annotation pipeline"""
    
    # Annotate top 50 prioritized genes
    input_file = "results/multi_omics_prioritized_top50.tsv"
    
    if not os.path.exists(input_file):
        print(f"[ERROR] File not found: {input_file}")
        return
    
    annotated_df = annotate_prioritized_genes(input_file, top_n=20)
    
    if annotated_df is not None:
        print("\n" + "=" * 80)
        print("LITERATURE ANNOTATION COMPLETE")
        print("=" * 80)

if __name__ == "__main__":
    main()


