"""
Pathway Enrichment Analysis for Multi-Omics Prioritized Genes
Uses g:Profiler API for comprehensive pathway enrichment
"""
import pandas as pd
import numpy as np
import requests
import json
import time
import os

def load_prioritized_genes():
    """Load multi-omics prioritized genes"""
    files_to_try = [
        "results/multi_omics_prioritized_top50.tsv",
        "results/multi_omics_comprehensive_prioritized_genes.tsv"
    ]
    
    for file in files_to_try:
        if os.path.exists(file):
            df = pd.read_csv(file, sep='\t')
            print(f"Loaded prioritized genes from: {file}")
            print(f"  Total genes: {len(df):,}")
            return df
    
    print("[ERROR] No prioritized genes file found")
    return None

def extract_gene_symbols(df):
    """Extract valid gene symbols from prioritized genes"""
    gene_symbols = []
    
    for idx, row in df.iterrows():
        symbol = str(row.get('Gene_Symbol', '')).strip()
        
        # Skip if empty, nan, or Ensembl ID
        if symbol and symbol != 'nan' and symbol != '':
            # Skip if it looks like an Ensembl ID
            if not symbol.startswith('ENSG'):
                gene_symbols.append(symbol)
    
    # Remove duplicates while preserving order
    unique_symbols = []
    seen = set()
    for symbol in gene_symbols:
        if symbol.upper() not in seen:
            unique_symbols.append(symbol)
            seen.add(symbol.upper())
    
    print(f"\nExtracted {len(unique_symbols)} unique gene symbols")
    print(f"  Examples: {', '.join(unique_symbols[:10])}")
    
    return unique_symbols

def perform_pathway_enrichment(gene_symbols, threshold=0.05):
    """
    Perform pathway enrichment using g:Profiler API
    """
    print(f"\nPerforming pathway enrichment for {len(gene_symbols)} genes...")
    print(f"  Using FDR threshold: {threshold}")
    
    url = "https://biit.cs.ut.ee/gprofiler/api/gost/profile"
    
    # Limit to 500 genes for API
    query_genes = gene_symbols[:500]
    print(f"  Querying {len(query_genes)} genes (from {len(gene_symbols)} total)...")
    
    payload = {
        "organism": "hsapiens",
        "query": query_genes,
        "sources": ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP", "HPA", "CORUM"],  # Multiple databases
        "user_threshold": threshold,
        "all_results": False,
        "ordered": True,
        "no_iea": False,
        "domain_scope": "annotated"
    }
    
    try:
        print("  Sending request to g:Profiler API...")
        response = requests.post(url, json=payload, timeout=90)
        print(f"  API response status: {response.status_code}")
        
        if response.status_code == 200:
            result = response.json()
            if 'result' in result:
                num_pathways = len(result['result'])
                print(f"  [OK] Found {num_pathways} enriched pathways")
                return result['result']
            else:
                print(f"  [WARNING] No 'result' key in API response")
                if 'meta' in result:
                    print(f"  Meta info: {result['meta']}")
                return None
        else:
            print(f"[ERROR] API returned status {response.status_code}")
            print(f"  Response text: {response.text[:500]}")
            return None
    except requests.exceptions.Timeout:
        print(f"[ERROR] Pathway enrichment timed out after 90 seconds")
        return None
    except Exception as e:
        print(f"[ERROR] Pathway enrichment failed: {e}")
        import traceback
        traceback.print_exc()
        return None

def process_enrichment_results(enrichment_results, gene_symbols):
    """Process and format enrichment results"""
    pathways = []
    
    for item in enrichment_results:
        # Get intersection genes
        intersections = item.get('intersections', [])
        gene_list = []
        
        if intersections:
            # g:Profiler returns intersections as list of lists
            for entry in intersections:
                if isinstance(entry, list) and len(entry) > 0:
                    gene_list.append(str(entry[0]))
                elif isinstance(entry, str):
                    gene_list.append(entry)
        
        # Limit to 20 genes for display
        genes_str = ', '.join(gene_list[:20]) if gene_list else ''
        if len(gene_list) > 20:
            genes_str += f" ... (+{len(gene_list) - 20} more)"
        
        pathways.append({
            'Term': item.get('name', ''),
            'Source': item.get('source', ''),
            'Term_ID': item.get('native', ''),
            'P_value': item.get('p_value', np.nan),
            'P_adj': item.get('p_value', np.nan),  # g:Profiler already adjusts
            'Term_size': item.get('term_size', 0),
            'Query_size': item.get('query_size', 0),
            'Intersection_size': item.get('intersection_size', 0),
            'Genes': genes_str,
            'Gene_Count': len(gene_list)
        })
    
    pathways_df = pd.DataFrame(pathways)
    
    # Calculate enrichment ratio
    pathways_df['Enrichment_Ratio'] = (
        pathways_df['Intersection_size'] / pathways_df['Query_size']
    ) / (pathways_df['Term_size'] / 20000)  # Approximate human genome size
    
    # Sort by p-value
    pathways_df = pathways_df.sort_values('P_value')
    
    return pathways_df

def analyze_pathways_by_priority(df, top_n=50):
    """Analyze pathways for top N prioritized genes"""
    
    print("\n" + "=" * 80)
    print("PATHWAY ENRICHMENT FOR MULTI-OMICS PRIORITIZED GENES")
    print("=" * 80)
    
    # Get top N genes
    top_genes = df.head(top_n).copy()
    print(f"\nAnalyzing top {len(top_genes)} prioritized genes...")
    
    # Extract gene symbols
    gene_symbols = extract_gene_symbols(top_genes)
    
    if len(gene_symbols) == 0:
        print("[ERROR] No valid gene symbols found")
        return None
    
    # Perform enrichment at different thresholds
    thresholds = [0.05, 0.1, 0.2]
    all_results = {}
    
    for threshold in thresholds:
        print(f"\n{'='*80}")
        print(f"PATHWAY ENRICHMENT - FDR < {threshold}")
        print(f"{'='*80}")
        
        enrichment_results = perform_pathway_enrichment(gene_symbols, threshold=threshold)
        
        if enrichment_results:
            pathways_df = process_enrichment_results(enrichment_results, gene_symbols)
            
            # Save results
            output_file = f"results/pathway_enrichment_prioritized_fdr_{threshold}.tsv"
            pathways_df.to_csv(output_file, sep='\t', index=False)
            print(f"\n[OK] Saved pathway enrichment results to {output_file}")
            
            # Print top pathways
            print(f"\nTop 20 enriched pathways (FDR < {threshold}):")
            print("-" * 80)
            top_pathways = pathways_df.head(20)
            for idx, row in top_pathways.iterrows():
                print(f"{row['Term'][:60]:60s} | {row['Source']:8s} | p={row['P_value']:.2e} | {row['Intersection_size']:3d} genes")
            
            all_results[threshold] = pathways_df
        else:
            print("[WARNING] No pathway enrichment results obtained")
            all_results[threshold] = None
        
        time.sleep(2)  # Rate limiting
    
    return all_results

def create_pathway_summary(all_results):
    """Create summary of pathway enrichment across thresholds"""
    print("\n" + "=" * 80)
    print("PATHWAY ENRICHMENT SUMMARY")
    print("=" * 80)
    
    summary_data = []
    
    for threshold, pathways_df in all_results.items():
        if pathways_df is not None and len(pathways_df) > 0:
            summary_data.append({
                'FDR_Threshold': threshold,
                'Total_Pathways': len(pathways_df),
                'Top_Pathway': pathways_df.iloc[0]['Term'],
                'Top_Source': pathways_df.iloc[0]['Source'],
                'Top_P_value': pathways_df.iloc[0]['P_value'],
                'Top_Intersection': pathways_df.iloc[0]['Intersection_size']
            })
        else:
            summary_data.append({
                'FDR_Threshold': threshold,
                'Total_Pathways': 0,
                'Top_Pathway': 'None',
                'Top_Source': 'None',
                'Top_P_value': np.nan,
                'Top_Intersection': 0
            })
    
    summary_df = pd.DataFrame(summary_data)
    print("\n" + summary_df.to_string(index=False))
    
    # Save summary
    summary_df.to_csv("results/pathway_enrichment_prioritized_summary.tsv", sep='\t', index=False)
    print(f"\n[OK] Saved summary to results/pathway_enrichment_prioritized_summary.tsv")
    
    return summary_df

def main():
    """Main pathway enrichment pipeline for prioritized genes"""
    
    print("=" * 80)
    print("PATHWAY ENRICHMENT FOR MULTI-OMICS PRIORITIZED GENES")
    print("=" * 80)
    
    # Load prioritized genes
    df = load_prioritized_genes()
    if df is None:
        return
    
    # Analyze pathways
    all_results = analyze_pathways_by_priority(df, top_n=50)
    
    if all_results:
        # Create summary
        summary_df = create_pathway_summary(all_results)
        
        print("\n" + "=" * 80)
        print("PATHWAY ENRICHMENT COMPLETE")
        print("=" * 80)
    else:
        print("\n[ERROR] Pathway enrichment failed")

if __name__ == "__main__":
    main()


