"""
Annotate Blood Dataset (GSE159610) with Gene Symbols using MyGene.info API
Converts Ensembl IDs to Gene Symbols, Entrez IDs, and Gene Names
"""
import pandas as pd
import numpy as np
import requests
import time
import os
from typing import Dict, List

def annotate_ensembl_ids_mygene(ensembl_ids: List[str], batch_size: int = 100):
    """
    Annotate Ensembl IDs with gene symbols, Entrez IDs, and gene names using MyGene.info API
    
    Args:
        ensembl_ids: List of Ensembl gene IDs (with or without version numbers)
        batch_size: Number of IDs to process per API call (max 1000)
    
    Returns:
        Dictionary mapping Ensembl ID to annotation data
    """
    print(f"\nAnnotating {len(ensembl_ids)} Ensembl IDs using MyGene.info API...")
    
    # Clean Ensembl IDs (remove version numbers for matching)
    clean_ids = []
    id_mapping = {}  # Map clean ID back to original
    
    for eid in ensembl_ids:
        eid_str = str(eid).strip()
        if eid_str and eid_str != 'nan' and eid_str != '':
            # Remove version number if present (e.g., ENSG00000123456.1 -> ENSG00000123456)
            clean_id = eid_str.split('.')[0] if '.' in eid_str else eid_str
            clean_ids.append(clean_id)
            id_mapping[clean_id] = eid_str
    
    unique_ids = list(set(clean_ids))
    print(f"  Unique Ensembl IDs: {len(unique_ids)}")
    
    url = "https://mygene.info/v3/query"
    annotation_map = {}
    
    # Process in batches
    total_batches = (len(unique_ids) + batch_size - 1) // batch_size
    for batch_idx in range(0, len(unique_ids), batch_size):
        batch_end = min(batch_idx + batch_size, len(unique_ids))
        batch_ids = unique_ids[batch_idx:batch_end]
        batch_num = batch_idx // batch_size + 1
        
        print(f"  Processing batch {batch_num}/{total_batches}: {len(batch_ids)} IDs...", end='', flush=True)
        
        try:
            # Use POST for batch queries (more reliable)
            data = {
                'q': batch_ids,
                'scopes': 'ensembl.gene',
                'fields': 'symbol,entrezgene,name,alias,ensembl',
                'species': 'human'
            }
            response = requests.post(url, json=data, timeout=60)
            
            if response.status_code == 200:
                results = response.json()
                
                batch_annotated = 0
                # POST with list returns a list where each item corresponds to query order
                if isinstance(results, list):
                    # Results are in same order as batch_ids
                    for i, item in enumerate(results):
                        if item and 'notfound' not in item:
                            query_id = batch_ids[i]  # Match by position
                            original_id = id_mapping.get(query_id, query_id)
                            
                            # Extract Ensembl ID from response to verify match
                            ensembl_id_from_response = None
                            if 'ensembl' in item:
                                ens = item['ensembl']
                                if isinstance(ens, dict) and 'gene' in ens:
                                    ensembl_id_from_response = str(ens['gene'])
                                elif isinstance(ens, list) and len(ens) > 0:
                                    if isinstance(ens[0], dict) and 'gene' in ens[0]:
                                        ensembl_id_from_response = str(ens[0]['gene'])
                            
                            annotation = {
                                'Ensembl_ID': original_id,
                                'Ensembl_ID_Clean': query_id.split('.')[0] if '.' in query_id else query_id,
                                'Gene_Symbol': item.get('symbol', ''),
                                'Entrez_ID': item.get('entrezgene', item.get('_id', '')),
                                'Gene_Name': item.get('name', ''),
                                'Aliases': ','.join(item.get('alias', [])) if isinstance(item.get('alias'), list) else ''
                            }
                            
                            # Handle Entrez ID if it's a list
                            if isinstance(annotation['Entrez_ID'], list):
                                annotation['Entrez_ID'] = annotation['Entrez_ID'][0] if annotation['Entrez_ID'] else ''
                            
                            annotation_map[original_id] = annotation
                            batch_annotated += 1
                
                print(f" [OK] {batch_annotated} annotated")
            else:
                print(f" [ERROR] HTTP {response.status_code}")
                if response.status_code == 414:
                    print(f"    URL too long - trying smaller batch size...")
                try:
                    error_text = response.text[:200]
                    print(f"    Response: {error_text}")
                except:
                    pass
            
            # Rate limiting
            time.sleep(0.5)
            
        except Exception as e:
            print(f" [ERROR] {e}")
            continue
    
    print(f"\n  Successfully annotated {len(annotation_map)} IDs")
    return annotation_map

def annotate_blood_dataset(input_file: str, output_file: str = None):
    """
    Annotate blood dataset with gene symbols and other identifiers
    
    Args:
        input_file: Path to blood dataset TSV file
        output_file: Path to save annotated dataset (default: adds '_annotated' to input filename)
    """
    print("=" * 80)
    print("ANNOTATING BLOOD DATASET (GSE159610) WITH GENE SYMBOLS")
    print("=" * 80)
    
    if not os.path.exists(input_file):
        print(f"[ERROR] Input file not found: {input_file}")
        return None
    
    # Load blood dataset
    print(f"\nLoading blood dataset: {input_file}")
    df = pd.read_csv(input_file, sep='\t')
    print(f"  Loaded {len(df):,} genes")
    
    # Check if already annotated
    if 'Gene_Symbol' in df.columns and df['Gene_Symbol'].notna().sum() > 0:
        print(f"\n[NOTE] Dataset already has {df['Gene_Symbol'].notna().sum()} genes with symbols")
        response = input("  Re-annotate anyway? (y/n): ").strip().lower()
        if response != 'y':
            print("  Skipping annotation")
            return df
    
    # Get Ensembl IDs
    if 'Gene_ID' not in df.columns:
        print("[ERROR] 'Gene_ID' column not found in dataset")
        return None
    
    ensembl_ids = df['Gene_ID'].astype(str).unique().tolist()
    print(f"  Found {len(ensembl_ids)} unique Ensembl IDs")
    
    # Annotate using MyGene.info
    annotation_map = annotate_ensembl_ids_mygene(ensembl_ids)
    
    # Add annotations to dataframe
    print("\nAdding annotations to dataset...")
    
    # Initialize annotation columns if they don't exist
    for col in ['Gene_Symbol', 'Entrez_ID', 'Gene_Name', 'Aliases']:
        if col not in df.columns:
            df[col] = ''
    
    annotated_count = 0
    for idx, row in df.iterrows():
        gene_id = str(row['Gene_ID'])
        if gene_id in annotation_map:
            annotation = annotation_map[gene_id]
            df.at[idx, 'Gene_Symbol'] = annotation['Gene_Symbol']
            df.at[idx, 'Entrez_ID'] = annotation['Entrez_ID']
            df.at[idx, 'Gene_Name'] = annotation['Gene_Name']
            df.at[idx, 'Aliases'] = annotation['Aliases']
            annotated_count += 1
    
    print(f"  Annotated {annotated_count:,} genes ({annotated_count/len(df)*100:.1f}%)")
    
    # Statistics
    print("\nAnnotation Statistics:")
    print(f"  Genes with symbols: {df['Gene_Symbol'].notna().sum() + (df['Gene_Symbol'] != '').sum():,}")
    print(f"  Genes with Entrez IDs: {df['Entrez_ID'].notna().sum() + (df['Entrez_ID'] != '').sum():,}")
    print(f"  Genes with names: {df['Gene_Name'].notna().sum() + (df['Gene_Name'] != '').sum():,}")
    
    # Save annotated dataset
    if output_file is None:
        base_name = os.path.splitext(input_file)[0]
        output_file = f"{base_name}_annotated.tsv"
    
    os.makedirs(os.path.dirname(output_file) if os.path.dirname(output_file) else 'results', exist_ok=True)
    df.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved annotated dataset to: {output_file}")
    
    # Show some examples
    print("\nExample annotations (first 10 with symbols):")
    annotated = df[df['Gene_Symbol'].notna() & (df['Gene_Symbol'] != '')].head(10)
    for idx, row in annotated.iterrows():
        print(f"  {row['Gene_ID']:25s} → {row['Gene_Symbol']:15s} ({row.get('Entrez_ID', 'N/A')}) - {row.get('Gene_Name', 'N/A')[:50]}")
    
    return df

def main():
    """Main annotation pipeline"""
    
    # Files to annotate (in priority order)
    files_to_annotate = [
        "results/gse159610_substantial_fc_genes.tsv",
        "results/gse159610_top_candidate_genes.tsv",
        "results/gse159610_full_de_results.tsv"
    ]
    
    annotated_files = []
    
    for input_file in files_to_annotate:
        if os.path.exists(input_file):
            print(f"\n{'='*80}")
            print(f"Processing: {input_file}")
            print(f"{'='*80}")
            
            annotated_df = annotate_blood_dataset(input_file)
            if annotated_df is not None:
                annotated_files.append(input_file)
        else:
            print(f"\n[SKIP] File not found: {input_file}")
    
    if annotated_files:
        print(f"\n{'='*80}")
        print("ANNOTATION COMPLETE")
        print(f"{'='*80}")
        print(f"\nAnnotated {len(annotated_files)} file(s):")
        for f in annotated_files:
            base_name = os.path.splitext(f)[0]
            print(f"  [OK] {base_name}_annotated.tsv")
        print("\nYou can now re-run integrate_all_findings_enhanced.py to see improved cross-dataset matching!")
    else:
        print("\n[WARNING] No files were annotated")

if __name__ == "__main__":
    main()

