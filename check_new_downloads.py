"""
Check for newly downloaded files that might be relevant to the aneurysm project
"""

from pathlib import Path
import os
from datetime import datetime, timedelta

downloads_dir = Path("C:/Users/vsanker/Downloads")

# Keywords that might indicate relevant datasets
relevant_keywords = [
    'gse', 'geo', 'aneurysm', 'cerebral', 'intracranial', 'vascular',
    'rna', 'seq', 'expression', 'transcript', 'gene', 'counts',
    'gwas', 'genome', 'genomic', 'biobank', 'gtex', 'tissue',
    'single', 'cell', 'scRNA', 'bulk', 'tsv', 'csv', 'gz', 'zip'
]

print("="*80)
print("CHECKING FOR NEWLY DOWNLOADED FILES")
print("="*80)

# Get files modified in last 24 hours
recent_time = datetime.now() - timedelta(hours=24)
recent_files = []

for file_path in downloads_dir.iterdir():
    if file_path.is_file():
        mtime = datetime.fromtimestamp(file_path.stat().st_mtime)
        if mtime > recent_time:
            recent_files.append((file_path, mtime, file_path.stat().st_size))

# Sort by modification time
recent_files.sort(key=lambda x: x[1], reverse=True)

if recent_files:
    print(f"\nFound {len(recent_files)} file(s) modified in last 24 hours:\n")
    
    for file_path, mtime, size in recent_files:
        size_mb = size / (1024 * 1024)
        name_lower = file_path.name.lower()
        
        # Check if file might be relevant
        is_relevant = any(kw in name_lower for kw in relevant_keywords)
        
        relevance = "✅ RELEVANT" if is_relevant else "⚠️  Check"
        
        print(f"{relevance}: {file_path.name}")
        print(f"  Size: {size_mb:.2f} MB")
        print(f"  Modified: {mtime.strftime('%Y-%m-%d %H:%M:%S')}")
        print()
else:
    print("\nNo files modified in last 24 hours.")
    print("\nChecking most recently modified files (last 10):")
    
    all_files = []
    for file_path in downloads_dir.iterdir():
        if file_path.is_file():
            mtime = datetime.fromtimestamp(file_path.stat().st_mtime)
            all_files.append((file_path, mtime, file_path.stat().st_size))
    
    all_files.sort(key=lambda x: x[1], reverse=True)
    
    for file_path, mtime, size in all_files[:10]:
        size_mb = size / (1024 * 1024)
        name_lower = file_path.name.lower()
        is_relevant = any(kw in name_lower for kw in relevant_keywords)
        relevance = "✅ RELEVANT" if is_relevant else "⚠️  Check"
        
        print(f"{relevance}: {file_path.name}")
        print(f"  Size: {size_mb:.2f} MB")
        print(f"  Modified: {mtime.strftime('%Y-%m-%d %H:%M:%S')}")
        print()

# Check for specific patterns
print("="*80)
print("CHECKING FOR SPECIFIC DATASET PATTERNS")
print("="*80)

patterns_to_check = [
    ('GSE', 'GEO dataset'),
    ('GSE158558', 'Intracranial aneurysm dataset'),
    ('GSE159610', 'Blood RNA-seq dataset'),
    ('GTEx', 'GTEx expression data'),
    ('gwas', 'GWAS data'),
    ('biobank', 'UK Biobank data'),
    ('aneurysm', 'Aneurysm-related data'),
    ('cerebral', 'Cerebral/vascular data'),
    ('.tsv', 'Tab-separated data'),
    ('.csv', 'CSV data'),
    ('.gz', 'Compressed data'),
    ('.h5', 'HDF5 format (single-cell)'),
    ('.h5ad', 'AnnData format (single-cell)'),
    ('.mtx', 'Matrix format (single-cell)'),
]

found_patterns = []
for pattern, description in patterns_to_check:
    matching_files = [f for f in downloads_dir.iterdir() 
                     if f.is_file() and pattern.lower() in f.name.lower()]
    if matching_files:
        found_patterns.append((pattern, description, matching_files))

if found_patterns:
    print("\nFiles matching relevant patterns:\n")
    for pattern, description, files in found_patterns:
        print(f"  {pattern} ({description}):")
        for f in files[:5]:  # Show first 5 matches
            size_mb = f.stat().st_size / (1024 * 1024)
            mtime = datetime.fromtimestamp(f.stat().st_mtime)
            print(f"    - {f.name} ({size_mb:.2f} MB, {mtime.strftime('%Y-%m-%d')})")
        if len(files) > 5:
            print(f"    ... and {len(files) - 5} more")
        print()

print("="*80)






