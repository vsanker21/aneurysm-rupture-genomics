"""
Pathway Enrichment Analysis using Gene Symbols
Uses g:Profiler API with gene symbols for better pathway matching
"""
import pandas as pd
import numpy as np
import requests
import json
import time
import os

def load_genes_with_symbols():
    """Load genes with gene symbols"""
    symbol_file = "results/gse159610_genes_with_symbols.tsv"
    
    if not os.path.exists(symbol_file):
        print(f"[ERROR] File not found: {symbol_file}")
        print("Please run convert_ensembl_to_symbols.py first")
        return None
    
    df = pd.read_csv(symbol_file, sep='\t')
    print(f"Loaded {len(df):,} genes with symbols")
    return df

def prepare_gene_list(df, max_genes=500):
    """Prepare gene symbol list for pathway enrichment"""
    # Sort by effect score
    df['Effect_Score'] = np.abs(df['Log2FC']) * (-np.log10(df['P_value'] + 1e-300))
    df_sorted = df.sort_values('Effect_Score', ascending=False)
    
    # Take top genes with symbols
    top_genes = df_sorted.head(max_genes)
    gene_symbols = top_genes['Gene_Symbol'].dropna().unique().tolist()
    
    print(f"\nPrepared {len(gene_symbols)} unique gene symbols for pathway enrichment")
    print(f"  Top 10 symbols: {gene_symbols[:10]}")
    
    return gene_symbols, top_genes

def perform_pathway_enrichment(gene_symbols, threshold=0.05):
    """Perform pathway enrichment using g:Profiler API with gene symbols"""
    print("\n" + "=" * 80)
    print("PATHWAY ENRICHMENT ANALYSIS (GENE SYMBOLS)")
    print("=" * 80)
    
    print(f"\nQuerying g:Profiler API for {len(gene_symbols)} gene symbols...")
    print(f"  FDR threshold: {threshold}")
    
    url = "https://biit.cs.ut.ee/gprofiler/api/gost/profile"
    
    # Limit to 500 genes for API
    query_genes = gene_symbols[:500]
    
    payload = {
        "organism": "hsapiens",
        "query": query_genes,
        "sources": ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"],
        "user_threshold": threshold,
        "all_results": False,
        "ordered": True,
        "no_iea": False,
        "domain_scope": "annotated"
    }
    
    try:
        print("  Sending request to g:Profiler API...")
        response = requests.post(url, json=payload, timeout=120)
        print(f"  API response status: {response.status_code}")
        
        if response.status_code == 200:
            result = response.json()
            if 'result' in result:
                pathways = result['result']
                print(f"  [OK] Found {len(pathways)} enriched pathways")
                return pathways
            else:
                print(f"  [WARNING] No 'result' key in API response")
                return None
        else:
            print(f"[ERROR] API returned status {response.status_code}")
            print(f"  Response: {response.text[:500]}")
            return None
    except Exception as e:
        print(f"[ERROR] API request failed: {e}")
        return None

def process_pathway_results(pathways, gene_df):
    """Process and format pathway enrichment results"""
    if pathways is None or len(pathways) == 0:
        print("[WARNING] No pathway results to process")
        return None
    
    print("\nProcessing pathway results...")
    
    pathway_data = []
    for pathway in pathways:
        pathway_data.append({
            'Pathway_ID': pathway.get('native', ''),
            'Pathway_Name': pathway.get('name', ''),
            'Source': pathway.get('source', ''),
            'P_value': pathway.get('p_value', np.nan),
            'FDR': pathway.get('adjusted_p_value', np.nan),
            'Gene_Count': pathway.get('intersection_size', 0),
            'Genes': ', '.join(pathway.get('intersection', [])),
            'Term_Size': pathway.get('term_size', 0),
            'Query_Size': pathway.get('query_size', 0)
        })
    
    pathway_df = pd.DataFrame(pathway_data)
    pathway_df = pathway_df.sort_values('FDR')
    
    print(f"\nPathway enrichment summary:")
    print(f"  Total pathways: {len(pathway_df)}")
    print(f"  FDR < 0.05: {(pathway_df['FDR'] < 0.05).sum()}")
    print(f"  FDR < 0.1: {(pathway_df['FDR'] < 0.1).sum()}")
    print(f"  FDR < 0.2: {(pathway_df['FDR'] < 0.2).sum()}")
    
    # Group by source
    print(f"\nBy database:")
    print(pathway_df['Source'].value_counts().to_string())
    
    # Top pathways
    print(f"\nTop 20 enriched pathways:")
    top_pathways = pathway_df.head(20)
    print(top_pathways[['Pathway_Name', 'Source', 'FDR', 'Gene_Count']].to_string(index=False))
    
    return pathway_df

def save_pathway_results(pathway_df):
    """Save pathway enrichment results"""
    os.makedirs("results", exist_ok=True)
    
    if pathway_df is None:
        print("[WARNING] No pathway results to save")
        return
    
    # Save all pathways
    output_file = "results/blood_pathway_enrichment_symbols_all.tsv"
    pathway_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved all pathways to {output_file}")
    
    # Save significant pathways
    for fdr_thresh in [0.05, 0.1, 0.2]:
        sig_pathways = pathway_df[pathway_df['FDR'] < fdr_thresh].copy()
        if len(sig_pathways) > 0:
            output_file = f"results/blood_pathway_enrichment_symbols_fdr_{fdr_thresh}.tsv"
            sig_pathways.to_csv(output_file, sep='\t', index=False)
            print(f"[OK] Saved {len(sig_pathways)} pathways (FDR < {fdr_thresh}) to {output_file}")

def main():
    """Main pathway enrichment pipeline"""
    print("=" * 80)
    print("PATHWAY ENRICHMENT WITH GENE SYMBOLS")
    print("=" * 80)
    
    # Load genes with symbols
    gene_df = load_genes_with_symbols()
    if gene_df is None:
        return
    
    # Prepare gene list
    gene_symbols, top_genes = prepare_gene_list(gene_df, max_genes=500)
    
    # Perform enrichment
    pathways = perform_pathway_enrichment(gene_symbols, threshold=0.05)
    
    # Process results
    pathway_df = process_pathway_results(pathways, top_genes)
    
    # Save results
    save_pathway_results(pathway_df)
    
    print("\n" + "=" * 80)
    print("PATHWAY ENRICHMENT COMPLETE!")
    print("=" * 80)

if __name__ == "__main__":
    main()






