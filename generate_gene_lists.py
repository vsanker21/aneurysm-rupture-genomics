"""Generate Gene Lists - Individual Script"""
import sys
import os
import traceback

if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

print("=" * 80)
print("GENERATING GENE LISTS")
print("=" * 80)
sys.stdout.flush()

try:
    import pandas as pd
    
    print("✓ Imports successful")
    sys.stdout.flush()
    
    os.makedirs("results/gene_lists", exist_ok=True)
    
    gene_lists = []
    
    # Load tissue signature genes
    print("\nLoading tissue signature genes...")
    sys.stdout.flush()
    if os.path.exists('results/multi_omics_tissue_signature_genes.tsv'):
        tissue_sig = pd.read_csv('results/multi_omics_tissue_signature_genes.tsv', sep='\t')
        print(f"✓ Loaded {len(tissue_sig)} tissue genes")
        print(f"  Columns: {list(tissue_sig.columns)}")
        sys.stdout.flush()
        
        for _, row in tissue_sig.iterrows():
            gene_lists.append({
                'Gene_ID': str(row.get('Gene', '')),
                'Gene_Symbol': str(row.get('Gene_Symbol', '')),
                'Dataset': 'Tissue (GSE122897)',
                'Log2FC': float(row.get('Log2FC', 0)),
                'P_adj': float(row.get('P_adj_BH', row.get('P_adj', row.get('P_value', 1)))),
                'Direction': 'Up' if float(row.get('Log2FC', 0)) > 0 else 'Down'
            })
    else:
        print("✗ File not found: results/multi_omics_tissue_signature_genes.tsv")
    
    # Load blood signature genes
    print("\nLoading blood signature genes...")
    sys.stdout.flush()
    if os.path.exists('results/multi_omics_blood_signature_genes.tsv'):
        blood_sig = pd.read_csv('results/multi_omics_blood_signature_genes.tsv', sep='\t')
        print(f"✓ Loaded {len(blood_sig)} blood genes")
        print(f"  Columns: {list(blood_sig.columns)}")
        sys.stdout.flush()
        
        for _, row in blood_sig.iterrows():
            gene_lists.append({
                'Gene_ID': str(row.get('Gene', row.get('Gene_ID', ''))),
                'Gene_Symbol': str(row.get('Gene_Symbol', '')),
                'Dataset': 'Blood (GSE159610)',
                'Log2FC': float(row.get('Log2FC', 0)),
                'P_adj': float(row.get('P_adj_BH', row.get('P_adj', row.get('P_value', 1)))),
                'Direction': 'Up' if float(row.get('Log2FC', 0)) > 0 else 'Down'
            })
    else:
        print("✗ File not found: results/multi_omics_blood_signature_genes.tsv")
    
    # Create comprehensive gene list
    if len(gene_lists) > 0:
        print(f"\nCreating comprehensive gene list ({len(gene_lists)} genes)...")
        sys.stdout.flush()
        
        all_genes_df = pd.DataFrame(gene_lists)
        output_file = 'results/gene_lists/comprehensive_signature_genes_annotated.tsv'
        all_genes_df.to_csv(output_file, sep='\t', index=False)
        
        if os.path.exists(output_file):
            file_size = os.path.getsize(output_file)
            print(f"✓ SUCCESS: Created {output_file} ({file_size/1024:.2f} KB)")
        else:
            print(f"✗ ERROR: File not created: {output_file}")
        
        # Create top 50 lists
        tissue_genes_df = all_genes_df[all_genes_df['Dataset'] == 'Tissue (GSE122897)'].head(50)
        blood_genes_df = all_genes_df[all_genes_df['Dataset'] == 'Blood (GSE159610)'].head(50)
        
        if len(tissue_genes_df) > 0:
            tissue_file = 'results/gene_lists/tissue_top50_signature_genes.tsv'
            tissue_genes_df.to_csv(tissue_file, sep='\t', index=False)
            if os.path.exists(tissue_file):
                print(f"✓ SUCCESS: Created {tissue_file} ({len(tissue_genes_df)} genes)")
        
        if len(blood_genes_df) > 0:
            blood_file = 'results/gene_lists/blood_top50_signature_genes.tsv'
            blood_genes_df.to_csv(blood_file, sep='\t', index=False)
            if os.path.exists(blood_file):
                print(f"✓ SUCCESS: Created {blood_file} ({len(blood_genes_df)} genes)")
        else:
            # Create empty file if no blood genes
            blood_file = 'results/gene_lists/blood_top50_signature_genes.tsv'
            empty_df = pd.DataFrame(columns=['Gene_ID', 'Gene_Symbol', 'Dataset', 'Log2FC', 'P_adj', 'Direction'])
            empty_df.to_csv(blood_file, sep='\t', index=False)
            print(f"⚠ Created empty {blood_file} (no blood genes found)")
    else:
        print("\n✗ ERROR: No genes found to export")
    
    sys.stdout.flush()
    
except Exception as e:
    print(f"\n✗ ERROR: {e}")
    print("\nFull traceback:")
    traceback.print_exc()
    sys.stdout.flush()
    sys.exit(1)

print("\n" + "=" * 80)
print("COMPLETE")
print("=" * 80)

