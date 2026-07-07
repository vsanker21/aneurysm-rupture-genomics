"""
Gene ID Conversion and Enhanced Cross-Dataset Integration
Converts between Ensembl and Entrez IDs for better integration
"""
import pandas as pd
import numpy as np
import requests
import time
import os
from typing import Dict, List, Set

def convert_gene_ids_mygene(gene_ids: List[str], from_type: str = 'auto', to_types: List[str] = ['symbol', 'entrezgene', 'ensembl.gene']):
    """
    Convert gene IDs using MyGene.info API
    from_type: 'ensembl', 'entrezgene', 'symbol', or 'auto' (detect automatically)
    """
    print(f"\nConverting {len(gene_ids)} gene IDs using MyGene.info API...")
    
    # Clean gene IDs
    clean_ids = []
    for gid in gene_ids:
        gid_str = str(gid).strip()
        if gid_str and gid_str != 'nan' and gid_str != '':
            # Remove version numbers from Ensembl IDs
            if '.' in gid_str and gid_str.startswith('ENSG'):
                gid_str = gid_str.split('.')[0]
            clean_ids.append(gid_str)
    
    if not clean_ids:
        print("  No valid gene IDs to convert")
        return {}
    
    unique_ids = list(set(clean_ids))
    print(f"  Unique IDs: {len(unique_ids)}")
    
    # Determine scope based on ID format - check ALL IDs, not just first 10
    if from_type == 'auto':
        has_ensembl = any(id.startswith('ENSG') for id in unique_ids)
        all_numeric = all(id.replace('.', '').replace('-', '').isdigit() for id in unique_ids[:50])
        
        if has_ensembl:
            from_type = 'ensembl.gene'
        elif all_numeric:
            from_type = 'entrezgene'
        else:
            from_type = 'symbol'
    
    print(f"  Detected ID type: {from_type}")
    
    url = "https://mygene.info/v3/query"
    conversion_map = {}
    batch_size = 500  # MyGene.info allows up to 1000
    
    # Process in batches
    for batch_start in range(0, len(unique_ids), batch_size):
        batch_end = min(batch_start + batch_size, len(unique_ids))
        batch_ids = unique_ids[batch_start:batch_end]
        
        print(f"  Processing batch {batch_start//batch_size + 1}: {len(batch_ids)} IDs")
        
        try:
            params = {
                'q': ','.join(batch_ids),
                'scopes': from_type,
                'fields': ','.join(['symbol'] + to_types),
                'species': 'human'
            }
            
            response = requests.get(url, params=params, timeout=60)
            if response.status_code == 200:
                results = response.json()
                
                for item in results:
                    if 'query' in item and 'notfound' not in item:
                        query_id = str(item['query'])
                        converted = {'original_id': query_id}
                        
                        # Extract symbol
                        if 'symbol' in item:
                            converted['symbol'] = item['symbol']
                        
                        # Extract Entrez ID
                        if 'entrezgene' in item:
                            entrez = item['entrezgene']
                            if isinstance(entrez, list):
                                entrez = entrez[0] if entrez else None
                            converted['entrezgene'] = entrez
                        
                        # Extract Ensembl ID
                        if 'ensembl' in item:
                            ensembl = item['ensembl']
                            if isinstance(ensembl, dict) and 'gene' in ensembl:
                                converted['ensembl_gene'] = ensembl['gene']
                            elif isinstance(ensembl, list) and ensembl:
                                converted['ensembl_gene'] = ensembl[0] if isinstance(ensembl[0], str) else ensembl[0].get('gene', None)
                        
                        conversion_map[query_id] = converted
            
            time.sleep(0.5)  # Rate limiting
            
        except Exception as e:
            print(f"    Error processing batch: {e}")
            continue
    
    print(f"  Successfully converted {len(conversion_map)} IDs")
    return conversion_map

def create_unified_gene_mapping(blood_df, tissue_df):
    """Create unified gene mapping across datasets"""
    print("\n" + "=" * 80)
    print("CREATING UNIFIED GENE MAPPING")
    print("=" * 80)
    
    # Collect all gene IDs
    blood_ids = []
    if blood_df is not None and 'Gene_ID' in blood_df.columns:
        blood_ids = blood_df['Gene_ID'].astype(str).unique().tolist()
        print(f"Blood dataset (GSE159610): {len(blood_ids)} unique Ensembl IDs")
    
    tissue_ids = []
    if tissue_df is not None:
        if 'Gene_ID' in tissue_df.columns:
            tissue_ids = tissue_df['Gene_ID'].astype(str).unique().tolist()
        elif 'Gene' in tissue_df.columns:
            tissue_ids = tissue_df['Gene'].astype(str).unique().tolist()
        print(f"Tissue dataset (GSE122897): {len(tissue_ids)} unique Entrez IDs")
    
    all_ids = list(set(blood_ids + tissue_ids))
    print(f"Total unique IDs to convert: {len(all_ids)}")
    
    # Convert IDs
    conversion_map = convert_gene_ids_mygene(all_ids, from_type='auto')
    
    # Create mapping DataFrames
    mapping_data = []
    for original_id, converted in conversion_map.items():
        mapping_data.append({
            'Original_ID': original_id,
            'Gene_Symbol': converted.get('symbol', ''),
            'Entrez_ID': converted.get('entrezgene', ''),
            'Ensembl_ID': converted.get('ensembl_gene', '')
        })
    
    if len(mapping_data) > 0:
        mapping_df = pd.DataFrame(mapping_data)
    else:
        # Create empty DataFrame with correct columns
        mapping_df = pd.DataFrame(columns=['Original_ID', 'Gene_Symbol', 'Entrez_ID', 'Ensembl_ID'])
        print("  [WARNING] No conversions successful, creating empty mapping")
    
    # Save mapping
    os.makedirs("results", exist_ok=True)
    mapping_file = "results/gene_id_mapping_unified.tsv"
    mapping_df.to_csv(mapping_file, sep='\t', index=False)
    print(f"\n[OK] Saved gene ID mapping to {mapping_file}")
    
    return mapping_df, conversion_map

def integrate_with_gene_symbols(blood_df, tissue_df, mapping_df):
    """Integrate blood and tissue results using gene symbols as common identifier"""
    print("\n" + "=" * 80)
    print("INTEGRATING WITH GENE SYMBOLS")
    print("=" * 80)
    
    # If mapping is empty, try to use existing annotations
    if len(mapping_df) == 0:
        print("  [WARNING] No gene ID mapping available, using existing annotations")
        mapping_df = None
    
    # Create symbol mappings
    blood_symbol_map = {}
    if blood_df is not None:
        # Check if Gene_Symbol already exists
        if 'Gene_Symbol' in blood_df.columns:
            for idx, row in blood_df.iterrows():
                symbol = str(row.get('Gene_Symbol', ''))
                if symbol and symbol != 'nan' and pd.notna(symbol):
                    blood_symbol_map[symbol] = row.to_dict()
        elif mapping_df is not None and len(mapping_df) > 0:
            for idx, row in blood_df.iterrows():
                gene_id = str(row.get('Gene_ID', ''))
                mapping_row = mapping_df[mapping_df['Original_ID'] == gene_id]
                if len(mapping_row) > 0:
                    symbol = mapping_row.iloc[0]['Gene_Symbol']
                    if pd.notna(symbol) and symbol != '':
                        blood_symbol_map[symbol] = row.to_dict()
    
    tissue_symbol_map = {}
    if tissue_df is not None:
        gene_col = 'Gene_ID' if 'Gene_ID' in tissue_df.columns else 'Gene'
        
        # Check if Gene_Symbol already exists
        if 'Gene_Symbol' in tissue_df.columns:
            for idx, row in tissue_df.iterrows():
                symbol = str(row.get('Gene_Symbol', ''))
                if symbol and symbol != 'nan' and pd.notna(symbol):
                    tissue_symbol_map[symbol] = row.to_dict()
        elif mapping_df is not None and len(mapping_df) > 0:
            for idx, row in tissue_df.iterrows():
                gene_id = str(row.get(gene_col, ''))
                mapping_row = mapping_df[mapping_df['Original_ID'] == gene_id]
                if len(mapping_row) > 0:
                    symbol = mapping_row.iloc[0]['Gene_Symbol']
                    if pd.notna(symbol) and symbol != '':
                        tissue_symbol_map[symbol] = row.to_dict()
    
    print(f"Blood genes with symbols: {len(blood_symbol_map)}")
    print(f"Tissue genes with symbols: {len(tissue_symbol_map)}")
    
    # Find common symbols
    common_symbols = set(blood_symbol_map.keys()) & set(tissue_symbol_map.keys())
    print(f"Genes found in both datasets: {len(common_symbols)}")
    
    # Create integrated DataFrame
    all_symbols = set(blood_symbol_map.keys()) | set(tissue_symbol_map.keys())
    integrated_data = []
    
    for symbol in all_symbols:
        gene_data = {
            'Gene_Symbol': symbol,
            'In_Blood': symbol in blood_symbol_map,
            'In_Tissue': symbol in tissue_symbol_map
        }
        
        # Add blood data
        if symbol in blood_symbol_map:
            blood_data = blood_symbol_map[symbol]
            gene_data.update({
                'Blood_Log2FC': blood_data.get('Log2FC', np.nan),
                'Blood_P_value': blood_data.get('P_value', np.nan),
                'Blood_P_adj': blood_data.get('P_adj_BH', np.nan),
                'Blood_Mean_Aneurysm': blood_data.get('Mean_Aneurysm', np.nan),
                'Blood_Mean_Control': blood_data.get('Mean_Control', np.nan),
                'Blood_Ensembl_ID': blood_data.get('Gene_ID', '')
            })
        
        # Add tissue data
        if symbol in tissue_symbol_map:
            tissue_data = tissue_symbol_map[symbol]
            gene_data.update({
                'Tissue_Log2FC': tissue_data.get('Log2FC', np.nan),
                'Tissue_P_value': tissue_data.get('P_value', np.nan),
                'Tissue_P_adj': tissue_data.get('P_adj_BH', np.nan),
                'Tissue_Mean_Ruptured': tissue_data.get('Mean_Ruptured', np.nan),
                'Tissue_Mean_Unruptured': tissue_data.get('Mean_Unruptured', np.nan),
                'Tissue_Entrez_ID': tissue_data.get('Gene' if 'Gene' in tissue_data else 'Gene_ID', '')
            })
        
        integrated_data.append(gene_data)
    
    integrated_df = pd.DataFrame(integrated_data)
    
    # Calculate cross-validation score
    def calc_integrated_score(row):
        score = 0.0
        
        if row['In_Blood']:
            blood_fc = abs(row.get('Blood_Log2FC', 0))
            if blood_fc > 2.0:
                score += 8
            elif blood_fc > 1.585:
                score += 6
            elif blood_fc > 1.0:
                score += 4
            
            if row.get('Blood_P_adj', 1.0) < 0.05:
                score += 5
            elif row.get('Blood_P_value', 1.0) < 0.05:
                score += 2
        
        if row['In_Tissue']:
            tissue_fc = abs(row.get('Tissue_Log2FC', 0))
            if tissue_fc > 2.0:
                score += 8
            elif tissue_fc > 1.585:
                score += 6
            elif tissue_fc > 1.0:
                score += 4
            
            if row.get('Tissue_P_adj', 1.0) < 0.05:
                score += 5
            elif row.get('Tissue_P_value', 1.0) < 0.05:
                score += 2
        
        # Cross-validation bonus
        if row['In_Blood'] and row['In_Tissue']:
            score += 15  # Significant bonus for being in both
            # Check same direction
            blood_fc = row.get('Blood_Log2FC', 0)
            tissue_fc = row.get('Tissue_Log2FC', 0)
            if (blood_fc > 0 and tissue_fc > 0) or (blood_fc < 0 and tissue_fc < 0):
                score += 10  # Same direction
        
        return score
    
    integrated_df['Integrated_Priority_Score'] = integrated_df.apply(calc_integrated_score, axis=1)
    integrated_df = integrated_df.sort_values('Integrated_Priority_Score', ascending=False)
    
    return integrated_df

def main():
    """Main pipeline"""
    print("=" * 80)
    print("GENE ID CONVERSION AND CROSS-DATASET INTEGRATION")
    print("=" * 80)
    
    # Load datasets - prefer annotated versions
    print("\nLoading datasets...")
    
    # Try annotated blood file first
    blood_file = "results/gse159610_top_candidate_genes_annotated.tsv"
    if not os.path.exists(blood_file):
        blood_file = "results/gse159610_top_candidate_genes.tsv"
    
    blood_df = None
    if os.path.exists(blood_file):
        blood_df = pd.read_csv(blood_file, sep='\t')
        print(f"Loaded blood data: {len(blood_df)} genes from {blood_file}")
    
    # Try annotated tissue file first
    tissue_file = "results/gse122897_top_genes_annotated.tsv"
    if not os.path.exists(tissue_file):
        tissue_file = "results/gse122897_differential_expression.tsv"
    
    tissue_df = None
    if os.path.exists(tissue_file):
        tissue_df = pd.read_csv(tissue_file, sep='\t')
        # Limit to top significant genes for efficiency if not already filtered
        if 'P_adj_BH' in tissue_df.columns:
            tissue_df = tissue_df[tissue_df['P_adj_BH'] < 0.1].copy()
        # Take top 500 if too many
        if len(tissue_df) > 500:
            tissue_df = tissue_df.head(500).copy()
        print(f"Loaded tissue data: {len(tissue_df)} significant genes from {tissue_file}")
    
    if blood_df is None and tissue_df is None:
        print("[ERROR] No datasets found!")
        return
    
    # Create unified mapping
    mapping_df, conversion_map = create_unified_gene_mapping(blood_df, tissue_df)
    
    # Integrate using gene symbols
    integrated_df = integrate_with_gene_symbols(blood_df, tissue_df, mapping_df)
    
    # Save integrated results
    output_file = "results/integrated_cross_dataset_genes.tsv"
    integrated_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved integrated results to {output_file}")
    
    # Summary
    print("\n" + "=" * 80)
    print("INTEGRATION SUMMARY")
    print("=" * 80)
    print(f"Total unique genes: {len(integrated_df)}")
    print(f"  In blood only: {(integrated_df['In_Blood'] & ~integrated_df['In_Tissue']).sum()}")
    print(f"  In tissue only: {(integrated_df['In_Tissue'] & ~integrated_df['In_Blood']).sum()}")
    print(f"  In both datasets: {(integrated_df['In_Blood'] & integrated_df['In_Tissue']).sum()}")
    
    cross_validated = integrated_df[(integrated_df['In_Blood'] & integrated_df['In_Tissue'])].head(20)
    if len(cross_validated) > 0:
        cross_file = "results/cross_validated_genes_symbol_based.tsv"
        cross_validated.to_csv(cross_file, sep='\t', index=False)
        print(f"\n[OK] Saved top 20 cross-validated genes to {cross_file}")
        print("\nTop cross-validated genes:")
        for idx, row in cross_validated.head(10).iterrows():
            symbol = row['Gene_Symbol']
            blood_fc = row.get('Blood_Log2FC', 0)
            tissue_fc = row.get('Tissue_Log2FC', 0)
            score = row['Integrated_Priority_Score']
            print(f"  {symbol:15s} | Blood FC: {blood_fc:6.2f} | Tissue FC: {tissue_fc:6.2f} | Score: {score:.1f}")

if __name__ == "__main__":
    main()

