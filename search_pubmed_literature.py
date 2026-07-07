"""
Search PubMed for top 10 genes from literature search template
Uses the NCBI E-utilities API to search PubMed
"""
import pandas as pd
import requests
import time
import os
from urllib.parse import quote

def search_pubmed(query, max_results=100):
    """
    Search PubMed using E-utilities API
    Returns: list of PMIDs and count
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    
    params = {
        'db': 'pubmed',
        'term': query,
        'retmax': max_results,
        'retmode': 'json',
        'sort': 'relevance'
    }
    
    try:
        response = requests.get(base_url, params=params, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if 'esearchresult' in data:
                pmids = data['esearchresult'].get('idlist', [])
                count = int(data['esearchresult'].get('count', 0))
                return pmids, count
        return [], 0
    except Exception as e:
        print(f"    [ERROR] PubMed search failed: {e}")
        return [], 0

def get_pubmed_summary(pmids, max_results=10):
    """
    Get summary information for PubMed IDs
    """
    if not pmids:
        return []
    
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    
    # Limit to top results
    pmids_to_fetch = pmids[:max_results]
    
    params = {
        'db': 'pubmed',
        'id': ','.join(pmids_to_fetch),
        'retmode': 'json'
    }
    
    try:
        response = requests.get(base_url, params=params, timeout=30)
        if response.status_code == 200:
            data = response.json()
            if 'result' in data:
                summaries = []
                for pmid in pmids_to_fetch:
                    if pmid in data['result']:
                        article = data['result'][pmid]
                        summaries.append({
                            'PMID': pmid,
                            'Title': article.get('title', 'N/A'),
                            'Authors': article.get('authors', [{}])[0].get('name', 'N/A') if article.get('authors') else 'N/A',
                            'Journal': article.get('source', 'N/A'),
                            'Year': article.get('pubdate', 'N/A').split()[0] if article.get('pubdate') else 'N/A',
                            'DOI': article.get('elocationid', 'N/A')
                        })
                return summaries
        return []
    except Exception as e:
        print(f"    [ERROR] Failed to fetch summaries: {e}")
        return []

def main():
    """Main literature search pipeline"""
    
    print("=" * 80)
    print("PUBMED LITERATURE SEARCH FOR TOP 10 GENES")
    print("=" * 80)
    
    # Load literature search template
    template_file = "results/literature_search_template_top10.tsv"
    if not os.path.exists(template_file):
        print(f"[ERROR] Literature search template not found: {template_file}")
        print("   Please run annotate_top_genes_for_literature.py first")
        return
    
    template_df = pd.read_csv(template_file, sep='\t')
    print(f"\nLoaded {len(template_df)} genes from template")
    
    # Search PubMed for each gene
    results = []
    
    for idx, row in template_df.iterrows():
        rank = row['Rank']
        gene_symbol = row['Gene_Symbol']
        gene_name = row['Gene_Name']
        
        print(f"\n[{rank}/10] Searching PubMed for {gene_symbol} ({gene_name})...")
        
        # Search with each query type
        query_results = {}
        
        query_mapping = {
            1: ('Search_Query_1_Direct', 'Direct'),
            2: ('Search_Query_2_Aneurysm', 'Aneurysm'),
            3: ('Search_Query_3_Vascular', 'Vascular'),
            4: ('Search_Query_4_Stroke', 'Stroke'),
            5: ('Search_Query_5_Biomarker', 'Biomarker')
        }
        
        for query_num in range(1, 6):
            query_col, query_name = query_mapping[query_num]
            
            query = row[query_col]
            print(f"  Query {query_num} ({query_name}): {query[:60]}...", end=' ')
            
            pmids, count = search_pubmed(query, max_results=50)
            query_results[f'Query_{query_num}_Count'] = count
            query_results[f'Query_{query_num}_PMIDs'] = ';'.join(pmids[:10])  # Store top 10 PMIDs
            
            print(f"{count} results")
            
            time.sleep(0.5)  # Rate limiting
        
        # Get summaries for the most relevant query (aneurysm query)
        aneurysm_query = row['Search_Query_2_Aneurysm']
        pmids, _ = search_pubmed(aneurysm_query, max_results=10)
        
        summaries = get_pubmed_summary(pmids, max_results=5)
        query_results['Top_5_PMIDs'] = ';'.join([s['PMID'] for s in summaries])
        query_results['Top_5_Titles'] = ' | '.join([s['Title'] for s in summaries])
        
        # Combine with gene info
        result_row = {
            'Rank': rank,
            'Gene_ID': row['Gene_ID'],
            'Gene_Symbol': gene_symbol,
            'Gene_Name': gene_name,
            'Log2FC': row['Log2FC'],
            'P_adj_BH': row['P_adj_BH'],
            **query_results
        }
        results.append(result_row)
        
        time.sleep(1)  # Rate limiting between genes
    
    # Create results DataFrame
    results_df = pd.DataFrame(results)
    
    # Save results
    output_file = "results/pubmed_literature_search_results.tsv"
    results_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved PubMed search results to {output_file}")
    
    # Print summary
    print(f"\n\nLiterature Search Summary:")
    print(f"{'Rank':<6} {'Symbol':<12} {'Direct':<8} {'Aneurysm':<10} {'Vascular':<10} {'Stroke':<8} {'Biomarker':<10}")
    print("-" * 80)
    
    for idx, row in results_df.iterrows():
        print(f"{row['Rank']:<6} {row['Gene_Symbol']:<12} {row['Query_1_Count']:<8} {row['Query_2_Count']:<10} {row['Query_3_Count']:<10} {row['Query_4_Count']:<8} {row['Query_5_Count']:<10}")
    
    # Genes with aneurysm-related publications
    aneurysm_genes = results_df[results_df['Query_2_Count'] > 0]
    print(f"\n\nGenes with aneurysm-related publications: {len(aneurysm_genes)}")
    if len(aneurysm_genes) > 0:
        for idx, row in aneurysm_genes.iterrows():
            print(f"  {row['Gene_Symbol']}: {row['Query_2_Count']} publications")
            if pd.notna(row['Top_5_Titles']) and row['Top_5_Titles']:
                titles = row['Top_5_Titles'].split(' | ')[:2]  # Show first 2
                for title in titles:
                    if title and title != 'N/A':
                        print(f"    - {title[:70]}...")
    
    print("\n" + "=" * 80)
    print("LITERATURE SEARCH COMPLETE")
    print("=" * 80)
    
    return results_df

if __name__ == "__main__":
    main()

