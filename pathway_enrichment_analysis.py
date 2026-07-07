"""
Pathway Enrichment Analysis for GSE122897 DE Genes
Uses g:Profiler API for pathway enrichment
"""
import pandas as pd
import numpy as np
import requests
import json
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

def get_gene_symbols(gene_ids):
    """
    Convert gene IDs to gene symbols using g:Profiler API
    Gene IDs appear to be Entrez IDs or similar numeric IDs
    """
    print(f"\nConverting {len(gene_ids)} gene IDs to symbols...")
    
    # Try to map using g:Profiler
    # First, try as Entrez IDs
    gene_list = [str(gid) for gid in gene_ids[:100]]  # Limit to 100 for API
    
    url = "https://biit.cs.ut.ee/gprofiler/api/convert/convert"
    
    payload = {
        "organism": "hsapiens",
        "target_namespace": "ENTREZGENE_ACC",
        "query": gene_list
    }
    
    try:
        response = requests.post(url, json=payload, timeout=30)
        if response.status_code == 200:
            result = response.json()
            if 'result' in result:
                mapping = {}
                for item in result['result']:
                    if 'converted' in item and len(item['converted']) > 0:
                        mapping[item['input']] = item['converted'][0]['name']
                return mapping
    except Exception as e:
        print(f"[WARNING] Could not convert via API: {e}")
    
    # Fallback: return IDs as strings (they might already be symbols)
    return {str(gid): str(gid) for gid in gene_ids}

def perform_pathway_enrichment(gene_list, threshold=0.1):
    """
    Perform pathway enrichment using g:Profiler API
    """
    print(f"\nPerforming pathway enrichment for {len(gene_list)} genes...")
    print(f"  Using FDR threshold: {threshold}")
    
    url = "https://biit.cs.ut.ee/gprofiler/api/gost/profile"
    
    # Convert gene list to strings
    query_genes = [str(g) for g in gene_list[:500]]  # Limit to 500 genes
    print(f"  Querying {len(query_genes)} genes (limited from {len(gene_list)})...")
    
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
        response = requests.post(url, json=payload, timeout=60)
        print(f"  API response status: {response.status_code}")
        
        if response.status_code == 200:
            result = response.json()
            if 'result' in result:
                num_pathways = len(result['result'])
                print(f"  [OK] Found {num_pathways} enriched pathways")
                return result['result']
            else:
                print(f"  [WARNING] No 'result' key in API response")
                print(f"  Response keys: {list(result.keys()) if isinstance(result, dict) else 'Not a dict'}")
                return None
        else:
            print(f"[ERROR] API returned status {response.status_code}")
            print(f"  Response text: {response.text[:500]}")
            return None
    except requests.exceptions.Timeout:
        print(f"[ERROR] Pathway enrichment timed out after 60 seconds")
        return None
    except Exception as e:
        print(f"[ERROR] Pathway enrichment failed: {e}")
        import traceback
        traceback.print_exc()
        return None

def analyze_pathways_by_fdr_threshold(df, thresholds=None):
    """Analyze pathways at different FDR thresholds"""
    
    if thresholds is None:
        thresholds = [0.05, 0.1, 0.2]
    
    results_summary = {}
    
    for fdr_threshold in thresholds:
        print(f"\n{'='*80}")
        print(f"PATHWAY ENRICHMENT ANALYSIS - FDR < {fdr_threshold}")
        print(f"{'='*80}")
        
        # Get significant genes
        sig_genes = df[df['P_adj_BH'] < fdr_threshold].copy()
        
        if len(sig_genes) == 0:
            print(f"No genes pass FDR < {fdr_threshold}")
            continue
        
        print(f"\nSignificant genes: {len(sig_genes)}")
        
        # Separate upregulated and downregulated
        upregulated = sig_genes[sig_genes['Log2FC'] > 0]['Gene'].tolist()
        downregulated = sig_genes[sig_genes['Log2FC'] < 0]['Gene'].tolist()
        
        print(f"  Upregulated: {len(upregulated)}")
        print(f"  Downregulated: {len(downregulated)}")
        
        # Perform enrichment for all significant genes
        all_sig_genes = sig_genes['Gene'].tolist()
        
        print(f"\nEnriching all significant genes...")
        print(f"  Total significant genes: {len(all_sig_genes)}")
        enrichment_results = perform_pathway_enrichment(all_sig_genes, threshold=fdr_threshold)
        
        if enrichment_results:
            # Convert to DataFrame
            pathways = []
            for item in enrichment_results:
                # Handle intersections - g:Profiler returns dict with gene IDs as keys
                intersections = item.get('intersections', [])
                gene_list = []
                
                # g:Profiler returns intersections as a list
                # Each item can be a string (gene ID) or a list
                if intersections:
                    for entry in intersections[:10]:  # First 10 genes
                        if isinstance(entry, list):
                            # If it's a list, take the first element (usually the gene ID)
                            if entry:
                                gene_list.append(str(entry[0]))
                        elif isinstance(entry, str):
                            # If it's a string, use it directly
                            gene_list.append(entry)
                        else:
                            # Convert to string
                            gene_list.append(str(entry))
                    genes_str = ', '.join(gene_list[:10])
                else:
                    # Fallback: try to get intersection_size and use query genes
                    intersection_size = item.get('intersection_size', 0)
                    if intersection_size > 0:
                        genes_str = f"{intersection_size} genes"
                    else:
                        genes_str = ''
                
                pathways.append({
                    'Term': item.get('name', ''),
                    'Source': item.get('source', ''),
                    'P_value': item.get('p_value', np.nan),
                    'P_adj': item.get('p_value', np.nan),  # g:Profiler already adjusts
                    'Term_size': item.get('term_size', 0),
                    'Query_size': item.get('query_size', 0),
                    'Intersection_size': item.get('intersection_size', 0),
                    'Genes': genes_str
                })
            
            pathways_df = pd.DataFrame(pathways)
            pathways_df = pathways_df.sort_values('P_value')
            
            # Save results
            output_file = f"results/pathway_enrichment_fdr_{fdr_threshold}.tsv"
            pathways_df.to_csv(output_file, sep='\t', index=False)
            print(f"\n[OK] Saved pathway enrichment results to {output_file}")
            
            # Print top pathways
            print(f"\nTop 20 enriched pathways:")
            print(pathways_df.head(20)[['Term', 'Source', 'P_value', 'Intersection_size']].to_string(index=False))
            
            results_summary[fdr_threshold] = {
                'total_genes': len(sig_genes),
                'pathways_found': len(pathways_df),
                'top_pathway': pathways_df.iloc[0]['Term'] if len(pathways_df) > 0 else None
            }
        else:
            print("[WARNING] No pathway enrichment results obtained")
            results_summary[fdr_threshold] = {
                'total_genes': len(sig_genes),
                'pathways_found': 0,
                'top_pathway': None
            }
        
        time.sleep(2)  # Rate limiting
    
    return results_summary

def main():
    """Main pathway enrichment pipeline"""
    
    print("=" * 80)
    print("PATHWAY ENRICHMENT ANALYSIS")
    print("=" * 80)
    
    # Load DE results
    df = load_de_results()
    if df is None:
        return
    
    # Analyze pathways at different FDR thresholds
    summary = analyze_pathways_by_fdr_threshold(df)
    
    # Save summary
    summary_df = pd.DataFrame(summary).T
    summary_df.index.name = 'FDR_threshold'
    summary_df.to_csv("results/pathway_enrichment_summary.tsv", sep='\t')
    
    print("\n" + "=" * 80)
    print("PATHWAY ENRICHMENT COMPLETE")
    print("=" * 80)
    print("\nSummary:")
    print(summary_df.to_string())

if __name__ == "__main__":
    main()

