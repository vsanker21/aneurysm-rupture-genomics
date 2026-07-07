"""Properly verify gene mapping"""
import pandas as pd

# Load mapping
mapping = pd.read_csv('results/gene_symbol_to_entrez_mapping.tsv', sep='\t')
print(f"Loaded mapping: {len(mapping)} genes")

# Load expression data (just index)
expr_df = pd.read_csv('data/sc_rnaseq/GSE122897/GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz', 
                      sep='\t', compression='gzip', index_col=0)
expr_ids = set(expr_df.index.astype(str))
print(f"Expression data has {len(expr_ids):,} genes")

# Check mapping
mapped_ids = set(mapping['Entrez_ID'].astype(str))
found = mapped_ids.intersection(expr_ids)
missing = mapped_ids - expr_ids

print(f"\nMapping verification:")
print(f"  Total genes to map: {len(mapped_ids)}")
print(f"  Found in expression data: {len(found)}")
print(f"  Missing: {len(missing)}")

if len(found) > 0:
    print(f"\n✅ Found genes ({len(found)}):")
    found_genes = mapping[mapping['Entrez_ID'].astype(str).isin(found)]
    for _, row in found_genes.iterrows():
        print(f"  {row['Gene_Symbol']} → {row['Entrez_ID']}")

if len(missing) > 0:
    print(f"\n⚠️  Missing genes ({len(missing)}):")
    missing_genes = mapping[mapping['Entrez_ID'].astype(str).isin(missing)]
    for _, row in missing_genes.iterrows():
        print(f"  {row['Gene_Symbol']} → {row['Entrez_ID']} (not in expression data)")

# Check specific examples
print(f"\n\nChecking specific examples:")
test_genes = ['IL6', 'TNF', 'MMP9', 'ANGPT2', 'LAMA3']
for gene in test_genes:
    gene_row = mapping[mapping['Gene_Symbol'] == gene]
    if len(gene_row) > 0:
        entrez_id = str(gene_row.iloc[0]['Entrez_ID'])
        found_status = "✅ FOUND" if entrez_id in expr_ids else "❌ MISSING"
        print(f"  {gene} (Entrez {entrez_id}): {found_status}")


