"""Verify gene mapping and update scripts to use it"""
import pandas as pd
import os

def verify_mapping():
    """Verify that mapped Entrez IDs exist in expression data"""
    print("=" * 80)
    print("VERIFYING GENE MAPPING")
    print("=" * 80)
    
    # Load mapping
    mapping_file = "results/gene_symbol_to_entrez_mapping.tsv"
    if not os.path.exists(mapping_file):
        print(f"[ERROR] Mapping file not found: {mapping_file}")
        return None
    
    mapping = pd.read_csv(mapping_file, sep='\t')
    print(f"\nLoaded mapping: {len(mapping)} genes")
    
    # Load expression data (just to check IDs)
    expr_file = "data/sc_rnaseq/GSE122897/GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
    if not os.path.exists(expr_file):
        print(f"[ERROR] Expression file not found: {expr_file}")
        return None
    
    # Read just the index (gene IDs)
    expr_df = pd.read_csv(expr_file, sep='\t', compression='gzip', index_col=0, nrows=0)
    expr_ids = set(expr_df.index.astype(str))
    
    print(f"Expression data has {len(expr_ids):,} genes")
    
    # Check which mapped IDs exist
    mapped_ids = set(mapping['Entrez_ID'].astype(str).astype(int).astype(str))
    found = mapped_ids.intersection(expr_ids)
    missing = mapped_ids - expr_ids
    
    print(f"\nMapping verification:")
    print(f"  Total genes to map: {len(mapped_ids)}")
    print(f"  Found in expression data: {len(found)}")
    print(f"  Missing: {len(missing)}")
    
    if len(missing) > 0:
        print(f"\n⚠️  Missing Entrez IDs (may need alternative IDs):")
        missing_genes = mapping[mapping['Entrez_ID'].astype(str).isin(missing)]
        for _, row in missing_genes.iterrows():
            print(f"  - {row['Gene_Symbol']} (Entrez ID: {row['Entrez_ID']})")
    
    if len(found) > 0:
        print(f"\n✅ Found genes (will work):")
        found_genes = mapping[mapping['Entrez_ID'].astype(str).isin(found)]
        print(f"  {len(found_genes)} genes ready to use")
        print(f"\n  Examples: {', '.join(found_genes['Gene_Symbol'].head(10).tolist())}")
    
    return mapping, found

if __name__ == "__main__":
    verify_mapping()


