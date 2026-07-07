"""
Helper script to create gene symbol to Entrez ID mapping
This will help you create the mapping file needed for feature calculation
"""
import pandas as pd
import requests
import time

# Genes that need to be mapped
INFLAMMATION_GENES = ['IL6', 'TNF', 'IL1B', 'CRP', 'CXCL8', 'IL8', 'CXCL1', 'CXCL2']
ECM_GENES = ['MMP9', 'MMP2', 'TIMP1', 'COL3A1', 'COL1A1', 'ELN', 'FBN1', 'LAMA3']
EC_DYSFUNCTION_GENES = ['ANGPT2', 'ICAM1', 'VCAM1', 'THBS1', 'PERP', 'CXCL8', 'CDH5', 'ENG', 'KDR', 'TEK', 'VWF', 'PECAM1', 'LYVE1', 'SELE', 'SELP', 'NOS3', 'NOS2']

# Known mappings (from common knowledge)
KNOWN_MAPPINGS = {
    'IL6': 3569,
    'TNF': 7124,
    'IL1B': 3553,
    'CRP': 1401,
    'CXCL8': 3576,
    'IL8': 3576,  # Same as CXCL8
    'CXCL1': 2919,
    'CXCL2': 2920,
    'MMP9': 4318,
    'MMP2': 4313,
    'TIMP1': 7076,
    'COL3A1': 1281,
    'COL1A1': 1277,
    'ELN': 2006,
    'FBN1': 2200,
    'LAMA3': 3909,
    'ANGPT2': 285,
    'ICAM1': 3383,
    'VCAM1': 7412,
    'THBS1': 7057,
    'PERP': 64065,
    'CDH5': 1003,
    'ENG': 2022,
    'KDR': 3791,
    'TEK': 7010,
    'VWF': 7450,
    'PECAM1': 5175,
    'LYVE1': 10894,
    'SELE': 6401,
    'SELP': 6403,
    'NOS3': 4846,
    'NOS2': 4843,
}

def get_gene_info_from_api(gene_symbol):
    """Get Entrez ID from MyGene.info API"""
    try:
        url = "https://mygene.info/v3/query"
        params = {
            'q': gene_symbol,
            'scopes': 'symbol',
            'fields': 'entrezgene,symbol,name',
            'species': 'human'
        }
        response = requests.get(url, params=params, timeout=10)
        if response.status_code == 200:
            data = response.json()
            if 'hits' in data and len(data['hits']) > 0:
                hit = data['hits'][0]
                return {
                    'Entrez_ID': hit.get('entrezgene', ''),
                    'Gene_Name': hit.get('name', ''),
                    'Symbol': hit.get('symbol', gene_symbol)
                }
        time.sleep(0.2)  # Rate limiting
    except Exception as e:
        print(f"  [WARNING] API call failed for {gene_symbol}: {e}")
    return None

def create_mapping_file():
    """Create gene symbol to Entrez ID mapping file"""
    print("=" * 80)
    print("CREATING GENE SYMBOL TO ENTREZ ID MAPPING")
    print("=" * 80)
    
    all_genes = list(set(INFLAMMATION_GENES + ECM_GENES + EC_DYSFUNCTION_GENES))
    
    mappings = []
    
    print(f"\nProcessing {len(all_genes)} genes...")
    print("Using known mappings first, then API for missing ones...\n")
    
    for gene_symbol in sorted(all_genes):
        if gene_symbol in KNOWN_MAPPINGS:
            entrez_id = KNOWN_MAPPINGS[gene_symbol]
            print(f"  ✓ {gene_symbol} → {entrez_id} (known mapping)")
            mappings.append({
                'Gene_Symbol': gene_symbol,
                'Entrez_ID': entrez_id,
                'Gene_Name': '',  # Can be filled manually
                'Source': 'Known'
            })
        else:
            print(f"  ? {gene_symbol} → Looking up...")
            info = get_gene_info_from_api(gene_symbol)
            if info and info.get('Entrez_ID'):
                print(f"    → {info['Entrez_ID']} ({info.get('Gene_Name', '')})")
                mappings.append({
                    'Gene_Symbol': gene_symbol,
                    'Entrez_ID': info['Entrez_ID'],
                    'Gene_Name': info.get('Gene_Name', ''),
                    'Source': 'API'
                })
            else:
                print(f"    → NOT FOUND - needs manual lookup")
                mappings.append({
                    'Gene_Symbol': gene_symbol,
                    'Entrez_ID': '',
                    'Gene_Name': '',
                    'Source': 'Manual_Needed'
                })
    
    # Create DataFrame
    df = pd.DataFrame(mappings)
    
    # Save to file
    output_file = "results/gene_symbol_to_entrez_mapping.tsv"
    df.to_csv(output_file, sep='\t', index=False)
    
    print("\n" + "=" * 80)
    print("MAPPING FILE CREATED")
    print("=" * 80)
    print(f"\nSaved to: {output_file}")
    print(f"\nTotal genes: {len(df)}")
    print(f"  Mapped: {(df['Entrez_ID'] != '').sum()}")
    print(f"  Need manual: {(df['Entrez_ID'] == '').sum()}")
    
    if (df['Entrez_ID'] == '').sum() > 0:
        print("\n⚠️  Genes needing manual mapping:")
        manual = df[df['Entrez_ID'] == '']
        for _, row in manual.iterrows():
            print(f"  - {row['Gene_Symbol']}: Look up at https://www.ncbi.nlm.nih.gov/gene")
    
    print("\n" + "=" * 80)
    print("NEXT STEPS")
    print("=" * 80)
    print("\n1. Review the mapping file: results/gene_symbol_to_entrez_mapping.tsv")
    print("2. Fill in any missing Entrez IDs manually")
    print("3. Verify Entrez IDs exist in your expression data")
    print("4. Let me know when ready - I'll update the scripts to use this mapping")
    
    return df

if __name__ == "__main__":
    create_mapping_file()


