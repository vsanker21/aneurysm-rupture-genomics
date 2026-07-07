"""Quick test of integration script"""
import os
import pandas as pd

print("=" * 80)
print("TESTING DATA AVAILABILITY")
print("=" * 80)

# Check blood data
blood_file = "data/blood_rnaseq/GSE159610/GSE159610_raw_gene_counts_matrix.txt.gz"
print(f"\nBlood file exists: {os.path.exists(blood_file)}")
if os.path.exists(blood_file):
    df = pd.read_csv(blood_file, sep='\t', compression='gzip', nrows=5, index_col=0)
    print(f"  Shape: {df.shape}")
    print(f"  Sample columns: {list(df.columns[:5])}")

# Check tissue data
tissue_file = "data/sc_rnaseq/GSE122897/GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
print(f"\nTissue file exists: {os.path.exists(tissue_file)}")
if os.path.exists(tissue_file):
    df = pd.read_csv(tissue_file, sep='\t', compression='gzip', nrows=5, index_col=0)
    print(f"  Shape: {df.shape}")
    print(f"  Sample columns: {list(df.columns[:5])}")

# Check DE results
blood_de = "results/gse159610_full_de_results.tsv"
tissue_de = "results/gse122897_differential_expression.tsv"

print(f"\nBlood DE exists: {os.path.exists(blood_de)}")
print(f"Tissue DE exists: {os.path.exists(tissue_de)}")

print("\n" + "=" * 80)
print("READY TO RUN INTEGRATION")
print("=" * 80)


