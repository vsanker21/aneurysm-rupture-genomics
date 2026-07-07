"""
Verify and organize GSE158558 files (Human Intracranial Aneurysm RNA-seq)
"""

from pathlib import Path
import shutil

downloads_dir = Path("C:/Users/vsanker/Downloads")
project_dir = Path(".")
data_dir = project_dir / "data" / "sc_rnaseq" / "GSE158558"
data_dir.mkdir(parents=True, exist_ok=True)

# GSE158558 files to look for
gse158558_files = [
    "GSE158558_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz",
    "GSE158558_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz",
    "GSE158558_raw_counts_GRCh38.p13_NCBI.tsv.gz",
    "GSE158558_mRNA_lncRNA_Expression_Transcript.xlsx",
    "Human.GRCh38.p13.annot.tsv.gz"
]

print("="*80)
print("VERIFYING GSE158558 FILES (Human Intracranial Aneurysm RNA-seq)")
print("="*80)

found_files = []
missing_files = []

for filename in gse158558_files:
    source = downloads_dir / filename
    if source.exists():
        size_mb = source.stat().st_size / (1024 * 1024)
        found_files.append((filename, size_mb, source))
        print(f"✅ {filename} ({size_mb:.2f} MB)")
    else:
        missing_files.append(filename)
        print(f"❌ {filename} - Not found")

print("\n" + "="*80)
print("DATASET INFORMATION")
print("="*80)

print("\nGSE158558: Intracranial aneurysm sequencing data")
print("  - Organism: Homo sapiens (Human)")
print("  - Type: Bulk RNA-seq (not single-cell)")
print("  - Expression profiling by high throughput sequencing")
print("  - This is an ALTERNATIVE to single-cell data")
print("  - Contains: mRNA and lncRNA expression data")

print("\n" + "="*80)
print("FILE TYPES")
print("="*80)

if found_files:
    print("\nFiles found:")
    for filename, size_mb, path in found_files:
        file_type = "Unknown"
        if "raw_counts" in filename:
            file_type = "Raw gene counts"
        elif "TPM" in filename:
            file_type = "TPM normalized expression"
        elif "FPKM" in filename:
            file_type = "FPKM normalized expression"
        elif "mRNA_lncRNA" in filename:
            file_type = "mRNA and lncRNA expression matrix"
        elif "annot" in filename:
            file_type = "Gene annotation file"
        
        print(f"  - {filename}")
        print(f"    Type: {file_type}")
        print(f"    Size: {size_mb:.2f} MB")
        print()

print("="*80)
print("RECOMMENDATION")
print("="*80)

print("\n✅ This dataset CAN be included in the project!")
print("\nWhy it's valuable:")
print("  - Human intracranial aneurysm tissue (correct organism)")
print("  - Bulk RNA-seq data (alternative to single-cell)")
print("  - Contains both mRNA and lncRNA expression")
print("  - Multiple normalization formats (TPM, FPKM, raw counts)")
print("  - Gene annotation file included")

print("\nHow it fits the project:")
print("  - Can be used for endothelial cell gene expression analysis")
print("  - Alternative to single-cell data for identifying dysregulated genes")
print("  - Can identify aneurysm-associated gene signatures")
print("  - Complements blood RNA-seq (GSE159610) and GWAS data")

print("\n" + "="*80)
print("NEXT STEPS")
print("="*80)

print("\n1. Move files to: data/sc_rnaseq/GSE158558/")
print("2. Update documentation to include GSE158558")
print("3. Note: This is bulk RNA-seq (not single-cell), but still valuable")
print("4. Can be used for endothelial dysfunction gene analysis")

if found_files:
    response = input("\nMove files to project directory? (y/n): ").strip().lower()
    if response == 'y':
        print("\nMoving files...")
        moved_count = 0
        for filename, size_mb, source in found_files:
            dest = data_dir / filename
            try:
                shutil.move(str(source), str(dest))
                print(f"  ✅ Moved: {filename}")
                moved_count += 1
            except Exception as e:
                print(f"  ❌ Error moving {filename}: {e}")
        
        print(f"\n✅ Moved {moved_count}/{len(found_files)} files")
        print(f"   Location: {data_dir}")
    else:
        print("\nFiles remain in Downloads folder")






