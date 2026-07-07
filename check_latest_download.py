"""
Check the most recently downloaded file for relevance to aneurysm project
"""

from pathlib import Path
from datetime import datetime, timedelta

downloads_dir = Path("C:/Users/vsanker/Downloads")

# Get all files with modification times
files_with_time = []
for file_path in downloads_dir.iterdir():
    if file_path.is_file():
        mtime = datetime.fromtimestamp(file_path.stat().st_mtime)
        size = file_path.stat().st_size
        files_with_time.append((file_path, mtime, size))

# Sort by modification time (most recent first)
files_with_time.sort(key=lambda x: x[1], reverse=True)

print("="*80)
print("CHECKING MOST RECENTLY DOWNLOADED FILES")
print("="*80)

# Show top 5 most recent
print("\nTop 5 most recently modified files:\n")

for i, (file_path, mtime, size) in enumerate(files_with_time[:5], 1):
    size_mb = size / (1024 * 1024)
    time_ago = datetime.now() - mtime
    
    # Check relevance
    name_lower = file_path.name.lower()
    relevant_keywords = ['gse', 'geo', 'aneurysm', 'cerebral', 'intracranial', 
                        'rna', 'seq', 'expression', 'gene', 'counts', 'gwas',
                        'genome', 'genomic', 'tissue', 'vascular', 'bulk', 'single', 'cell']
    
    is_relevant = any(kw in name_lower for kw in relevant_keywords)
    
    relevance = "✅ RELEVANT" if is_relevant else "⚠️  Check"
    
    print(f"{i}. {relevance}: {file_path.name}")
    print(f"   Size: {size_mb:.2f} MB")
    print(f"   Modified: {mtime.strftime('%Y-%m-%d %H:%M:%S')} ({time_ago.total_seconds()/60:.1f} minutes ago)")
    print()

# Check the most recent file in detail
if files_with_time:
    most_recent = files_with_time[0]
    file_path, mtime, size = most_recent
    
    print("="*80)
    print(f"DETAILED ANALYSIS: {file_path.name}")
    print("="*80)
    
    size_mb = size / (1024 * 1024)
    name_lower = file_path.name.lower()
    
    print(f"\nFile: {file_path.name}")
    print(f"Size: {size_mb:.2f} MB ({size:,} bytes)")
    print(f"Modified: {mtime.strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Check for GEO/GSE pattern
    if 'gse' in name_lower:
        gse_match = [s for s in name_lower.split() if 'gse' in s]
        if gse_match:
            print(f"\n✅ Contains GEO accession: {gse_match[0].upper()}")
    
    # Check file type
    if file_path.suffix == '.gz':
        print("✅ Compressed file (.gz)")
    elif file_path.suffix in ['.tsv', '.txt', '.csv']:
        print(f"✅ Text data file ({file_path.suffix})")
    elif file_path.suffix in ['.xlsx', '.xls']:
        print(f"✅ Excel file ({file_path.suffix})")
    elif file_path.suffix in ['.h5', '.h5ad']:
        print(f"✅ HDF5/AnnData file ({file_path.suffix}) - Single-cell format")
    elif file_path.suffix == '.zip':
        print("⚠️  ZIP archive - May contain multiple files")
    
    # Relevance check
    print("\n" + "="*80)
    print("RELEVANCE ASSESSMENT")
    print("="*80)
    
    relevance_factors = []
    
    if 'gse' in name_lower or 'geo' in name_lower:
        relevance_factors.append("✅ GEO dataset file")
    
    if 'aneurysm' in name_lower or 'cerebral' in name_lower or 'intracranial' in name_lower:
        relevance_factors.append("✅ Aneurysm-related")
    
    if 'rna' in name_lower or 'seq' in name_lower or 'expression' in name_lower:
        relevance_factors.append("✅ RNA-seq/expression data")
    
    if 'counts' in name_lower or 'tpm' in name_lower or 'fpkm' in name_lower:
        relevance_factors.append("✅ Expression counts/normalized data")
    
    if 'gene' in name_lower or 'transcript' in name_lower:
        relevance_factors.append("✅ Gene/transcript data")
    
    if 'gwas' in name_lower or 'genome' in name_lower:
        relevance_factors.append("✅ GWAS/genomic data")
    
    if 'tissue' in name_lower or 'vascular' in name_lower:
        relevance_factors.append("✅ Tissue/vascular data")
    
    if relevance_factors:
        print("\n✅ RELEVANT FACTORS FOUND:")
        for factor in relevance_factors:
            print(f"   {factor}")
        print("\n✅ This file CAN be included in the project!")
    else:
        print("\n⚠️  No obvious relevance factors found")
        print("   File may still be relevant - manual review recommended")
    
    print("\n" + "="*80)
    print("RECOMMENDATION")
    print("="*80)
    
    if relevance_factors:
        print("\n✅ RECOMMENDED: Include this file in the project")
        print("\nNext steps:")
        print("1. Verify file contents match expected format")
        print("2. Move to appropriate project directory")
        print("3. Update documentation")
    else:
        print("\n⚠️  Review file manually to determine relevance")






