"""
Pathway Enrichment Analysis for Blood DE Genes
Uses g:Profiler API to identify enriched pathways
"""
import pandas as pd
import numpy as np
import requests
import json
import time
import os

def load_blood_genes():
    """Load blood DE genes for pathway enrichment"""
    # Try substantial FC genes first
    fc_file = "results/gse159610_substantial_fc_genes.tsv"
    
    if os.path.exists(fc_file):
        df = pd.read_csv(fc_file, sep='\t')
        print(f"Loaded {len(df):,} genes with substantial FC")
        return df
    else:
        # Fall back to full results
        full_file = "results/gse159610_full_de_results.tsv"
        if os.path.exists(full_file):
            df = pd.read_csv(full_file, sep='\t')
            # Filter to top genes
            df = df.nsmallest(500, 'P_value')
            print(f"Loaded top 500 genes from full results")
            return df
        else:
            print("[ERROR] No blood DE results found")
            return None

def prepare_gene_list(df, max_genes=500):
    """Prepare gene list for pathway enrichment"""
    # Sort by effect score (combining FC and p-value)
    df['Effect_Score'] = np.abs(df['Log2FC']) * (-np.log10(df['P_value'] + 1e-300))
    df_sorted = df.sort_values('Effect_Score', ascending=False)
    
    # Take top genes
    top_genes = df_sorted.head(max_genes)
    
    # Extract gene IDs (Ensembl format)
    gene_list = top_genes['Gene_ID'].astype(str).tolist()
    
    print(f"\nPrepared {len(gene_list)} genes for pathway enrichment")
    print(f"  Top 10 genes: {gene_list[:10]}")
    
    return gene_list, top_genes

def perform_pathway_enrichment(gene_list, threshold=0.05):
    """
    Perform pathway enrichment using g:Profiler API
    """
    print("\n" + "=" * 80)
    print("PATHWAY ENRICHMENT ANALYSIS")
    print("=" * 80)
    
    print(f"\nQuerying g:Profiler API for {len(gene_list)} genes...")
    print(f"  FDR threshold: {threshold}")
    
    url = "https://biit.cs.ut.ee/gprofiler/api/gost/profile"
    
    # Limit to 500 genes for API
    query_genes = gene_list[:500]
    
    payload = {
        "organism": "hsapiens",
        "query": query_genes,
        "sources": ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"],  # Multiple databases
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
    except requests.exceptions.Timeout:
        print("[ERROR] Request timed out")
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
    
    # Convert to DataFrame
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

def save_pathway_results(pathway_df, gene_df):
    """Save pathway enrichment results"""
    os.makedirs("results", exist_ok=True)
    
    if pathway_df is None:
        print("[WARNING] No pathway results to save")
        return
    
    # Save all pathways
    output_file = "results/blood_pathway_enrichment_all.tsv"
    pathway_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved all pathways to {output_file}")
    
    # Save significant pathways at different thresholds
    for fdr_thresh in [0.05, 0.1, 0.2]:
        sig_pathways = pathway_df[pathway_df['FDR'] < fdr_thresh].copy()
        if len(sig_pathways) > 0:
            output_file = f"results/blood_pathway_enrichment_fdr_{fdr_thresh}.tsv"
            sig_pathways.to_csv(output_file, sep='\t', index=False)
            print(f"[OK] Saved {len(sig_pathways)} pathways (FDR < {fdr_thresh}) to {output_file}")
    
    # Save by database
    for source in pathway_df['Source'].unique():
        source_pathways = pathway_df[pathway_df['Source'] == source].copy()
        if len(source_pathways) > 0:
            output_file = f"results/blood_pathway_enrichment_{source}.tsv"
            source_pathways.to_csv(output_file, sep='\t', index=False)
            print(f"[OK] Saved {len(source_pathways)} {source} pathways to {output_file}")

def main():
    """Main pathway enrichment pipeline"""
    print("=" * 80)
    print("PATHWAY ENRICHMENT FOR BLOOD DE GENES")
    print("=" * 80)
    
    # Load genes
    gene_df = load_blood_genes()
    if gene_df is None:
        return
    
    # Prepare gene list
    gene_list, top_genes = prepare_gene_list(gene_df, max_genes=500)
    
    # Perform enrichment
    pathways = perform_pathway_enrichment(gene_list, threshold=0.05)
    
    # Process results
    pathway_df = process_pathway_results(pathways, top_genes)
    
    # Save results
    save_pathway_results(pathway_df, top_genes)
    
    print("\n" + "=" * 80)
    print("PATHWAY ENRICHMENT COMPLETE!")
    print("=" * 80)

if __name__ == "__main__":
    main()






