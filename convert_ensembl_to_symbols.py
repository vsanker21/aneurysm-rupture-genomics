"""
Convert Ensembl IDs to Gene Symbols using mygene API
This enables better pathway enrichment and gene annotation
"""
import pandas as pd
import numpy as np
import requests
import json
import time
import os
from typing import Dict, List

def convert_ensembl_to_symbols_batch(gene_ids: List[str], batch_size: int = 1000):
    """
    Convert Ensembl IDs to gene symbols using mygene API
    """
    print(f"Converting {len(gene_ids)} Ensembl IDs to gene symbols...")
    print(f"Using batch size: {batch_size}")
    
    url = "https://mygene.info/v3/query"
    
    all_results = {}
    failed_ids = []
    
    # Process in batches
    for i in range(0, len(gene_ids), batch_size):
        batch = gene_ids[i:i+batch_size]
        print(f"  Processing batch {i//batch_size + 1}/{(len(gene_ids)-1)//batch_size + 1} ({len(batch)} genes)...")
        
        # Query mygene API
        params = {
            'q': ','.join(batch),
            'scopes': 'ensembl.gene',
            'fields': 'symbol,entrezgene,name,ensembl',
            'species': 'human'
        }
        
        try:
            response = requests.get(url, params=params, timeout=60)
            if response.status_code == 200:
                results = response.json()
                
                # Process results
                if isinstance(results, list):
                    for item in results:
                        if 'query' in item and 'symbol' in item:
                            query_id = item['query']
                            symbol = item.get('symbol', '')
                            entrez = item.get('entrezgene', '')
                            name = item.get('name', '')
                            
                            all_results[query_id] = {
                                'Gene_Symbol': symbol,
                                'Entrez_ID': entrez,
                                'Gene_Name': name
                            }
                        elif 'query' in item:
                            # No symbol found
                            failed_ids.append(item['query'])
                else:
                    print(f"    [WARNING] Unexpected response format")
                
                # Rate limiting
                time.sleep(0.5)
            else:
                print(f"    [ERROR] API returned status {response.status_code}")
                failed_ids.extend(batch)
        except Exception as e:
            print(f"    [ERROR] Batch failed: {e}")
            failed_ids.extend(batch)
    
    print(f"\nConversion complete:")
    print(f"  Successfully converted: {len(all_results)}")
    print(f"  Failed: {len(failed_ids)}")
    
    if len(failed_ids) > 0 and len(failed_ids) <= 20:
        print(f"  Failed IDs: {failed_ids[:20]}")
    
    return all_results, failed_ids

def convert_blood_genes():
    """Convert blood DE gene IDs to symbols"""
    print("=" * 80)
    print("CONVERTING BLOOD DE GENE IDs TO SYMBOLS")
    print("=" * 80)
    
    # Load substantial FC genes
    fc_file = "results/gse159610_substantial_fc_genes.tsv"
    if not os.path.exists(fc_file):
        print(f"[ERROR] File not found: {fc_file}")
        return None
    
    df = pd.read_csv(fc_file, sep='\t')
    print(f"\nLoaded {len(df):,} genes")
    
    # Extract Ensembl IDs (remove version numbers for better matching)
    gene_ids = df['Gene_ID'].astype(str).tolist()
    gene_ids_clean = [gid.split('.')[0] for gid in gene_ids]
    
    # Convert to symbols
    mapping, failed = convert_ensembl_to_symbols_batch(gene_ids_clean, batch_size=1000)
    
    # Add symbols to dataframe
    df['Gene_ID_Clean'] = gene_ids_clean
    df['Gene_Symbol'] = df['Gene_ID_Clean'].map(lambda x: mapping.get(x, {}).get('Gene_Symbol', ''))
    df['Entrez_ID'] = df['Gene_ID_Clean'].map(lambda x: mapping.get(x, {}).get('Entrez_ID', ''))
    df['Gene_Name'] = df['Gene_ID_Clean'].map(lambda x: mapping.get(x, {}).get('Gene_Name', ''))
    
    # Save annotated results
    output_file = "results/gse159610_substantial_fc_genes_annotated.tsv"
    df.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved annotated genes to {output_file}")
    
    # Statistics
    print(f"\nAnnotation statistics:")
    print(f"  Genes with symbols: {(df['Gene_Symbol'] != '').sum():,}")
    print(f"  Genes with Entrez IDs: {(df['Entrez_ID'] != '').sum():,}")
    print(f"  Genes with names: {(df['Gene_Name'] != '').sum():,}")
    
    # Save genes with symbols for pathway analysis
    genes_with_symbols = df[df['Gene_Symbol'] != ''].copy()
    if len(genes_with_symbols) > 0:
        symbol_file = "results/gse159610_genes_with_symbols.tsv"
        genes_with_symbols.to_csv(symbol_file, sep='\t', index=False)
        print(f"[OK] Saved {len(genes_with_symbols):,} genes with symbols to {symbol_file}")
        
        # Create gene symbol list for pathway analysis
        symbol_list = genes_with_symbols['Gene_Symbol'].dropna().unique().tolist()
        symbol_list_file = "results/gse159610_gene_symbols_list.txt"
        with open(symbol_list_file, 'w') as f:
            f.write('\n'.join(symbol_list))
        print(f"[OK] Saved {len(symbol_list)} unique gene symbols to {symbol_list_file}")
    
    return df

def convert_top_candidates():
    """Convert top candidate gene IDs to symbols"""
    print("\n" + "=" * 80)
    print("CONVERTING TOP CANDIDATE GENE IDs")
    print("=" * 80)
    
    top_file = "results/gse159610_top_candidate_genes.tsv"
    if not os.path.exists(top_file):
        print(f"[INFO] Top candidates file not found: {top_file}")
        return None
    
    df = pd.read_csv(top_file, sep='\t')
    print(f"\nLoaded {len(df):,} top candidate genes")
    
    # Extract Ensembl IDs
    gene_ids = df['Gene_ID'].astype(str).tolist()
    gene_ids_clean = [gid.split('.')[0] for gid in gene_ids]
    
    # Convert to symbols
    mapping, failed = convert_ensembl_to_symbols_batch(gene_ids_clean, batch_size=100)
    
    # Add symbols
    df['Gene_ID_Clean'] = gene_ids_clean
    df['Gene_Symbol'] = df['Gene_ID_Clean'].map(lambda x: mapping.get(x, {}).get('Gene_Symbol', ''))
    df['Entrez_ID'] = df['Gene_ID_Clean'].map(lambda x: mapping.get(x, {}).get('Entrez_ID', ''))
    df['Gene_Name'] = df['Gene_ID_Clean'].map(lambda x: mapping.get(x, {}).get('Gene_Name', ''))
    
    # Save
    output_file = "results/gse159610_top_candidate_genes_annotated.tsv"
    df.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved annotated top candidates to {output_file}")
    
    # Show top genes with symbols
    genes_with_symbols = df[df['Gene_Symbol'] != ''].copy()
    if len(genes_with_symbols) > 0:
        print(f"\nTop 20 genes with symbols:")
        top_20 = genes_with_symbols.head(20)
        print(top_20[['Gene_Symbol', 'Gene_ID', 'Log2FC', 'P_value', 'Gene_Name']].to_string(index=False))
    
    return df

def main():
    """Main conversion pipeline"""
    # Convert substantial FC genes
    df_all = convert_blood_genes()
    
    # Convert top candidates
    df_top = convert_top_candidates()
    
    print("\n" + "=" * 80)
    print("GENE ID CONVERSION COMPLETE!")
    print("=" * 80)
    print("\nNext steps:")
    print("1. Run pathway enrichment with gene symbols:")
    print("   python pathway_enrichment_blood_genes_symbols.py")
    print("2. Use annotated genes for literature review")

if __name__ == "__main__":
    main()






