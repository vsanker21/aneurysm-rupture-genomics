"""
Verify GTEx downloaded files and check if they match requirements
"""

from pathlib import Path
import os
import gzip
import json

# Required vessel tissues
REQUIRED_TISSUES = [
    "Artery_Aorta",
    "Artery_Coronary", 
    "Artery_Tibial"
]

# Common GTEx file patterns
GTEX_FILE_PATTERNS = [
    "*gtex*",
    "*GTEx*",
    "*Artery*",
    "*Aorta*",
    "*Coronary*",
    "*Tibial*",
    "*.gct",
    "*.gct.gz",
    "*.txt",
    "*.tsv",
    "*.csv"
]

def find_gtex_files(search_dirs=None):
    """Find all potential GTEx files"""
    if search_dirs is None:
        search_dirs = [
            Path("data/gtex"),
            Path("."),
            Path.home() / "Downloads",
        ]
    
    found_files = []
    for search_dir in search_dirs:
        if not search_dir.exists():
            continue
        
        print(f"\nSearching in: {search_dir}")
        for pattern in GTEX_FILE_PATTERNS:
            for file_path in search_dir.rglob(pattern):
                if file_path.is_file() and file_path.suffix in ['.gz', '.gct', '.txt', '.tsv', '.csv', '.json']:
                    if file_path not in found_files:
                        found_files.append(file_path)
    
    return found_files

def analyze_file(file_path):
    """Analyze a GTEx file to determine its type and content"""
    info = {
        'path': file_path,
        'name': file_path.name,
        'size': file_path.stat().st_size,
        'size_mb': file_path.stat().st_size / (1024 * 1024),
        'type': 'unknown',
        'tissue': None,
        'format': None,
        'sample_count': None,
        'gene_count': None
    }
    
    # Check file name for tissue type
    name_lower = file_path.name.lower()
    if 'aorta' in name_lower:
        info['tissue'] = 'Artery_Aorta'
    elif 'coronary' in name_lower:
        info['tissue'] = 'Artery_Coronary'
    elif 'tibial' in name_lower:
        info['tissue'] = 'Artery_Tibial'
    
    # Check file extension
    if file_path.suffix == '.gz':
        info['format'] = 'compressed'
        if '.gct' in file_path.name:
            info['type'] = 'GCT (Gene Cluster Text)'
        elif '.txt' in file_path.name or '.tsv' in file_path.name:
            info['type'] = 'Text/TSV'
    elif file_path.suffix == '.gct':
        info['type'] = 'GCT (Gene Cluster Text)'
        info['format'] = 'uncompressed'
    elif file_path.suffix in ['.txt', '.tsv', '.csv']:
        info['type'] = 'Text/TSV/CSV'
        info['format'] = 'uncompressed'
    elif file_path.suffix == '.json':
        info['type'] = 'JSON'
        info['format'] = 'uncompressed'
    
    # Try to peek at file contents
    try:
        if file_path.suffix == '.gz':
            with gzip.open(file_path, 'rt', encoding='utf-8', errors='ignore') as f:
                first_lines = [f.readline() for _ in range(5)]
        else:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                first_lines = [f.readline() for _ in range(5)]
        
        # Check if it's a GCT file
        if first_lines and '#1.2' in first_lines[0]:
            info['type'] = 'GCT (Gene Cluster Text)'
            if len(first_lines) > 1:
                try:
                    parts = first_lines[1].strip().split('\t')
                    if len(parts) >= 2:
                        info['gene_count'] = int(parts[0])
                        info['sample_count'] = int(parts[1])
                except:
                    pass
        
        # Check for JSON
        elif first_lines and first_lines[0].strip().startswith('['):
            info['type'] = 'JSON'
            try:
                with open(file_path, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                    if isinstance(data, list) and len(data) > 0:
                        info['gene_count'] = len(data)
            except:
                pass
        
        # Check for tab/CSV format
        elif first_lines and '\t' in first_lines[0]:
            info['type'] = 'Tab-separated'
            parts = first_lines[0].strip().split('\t')
            info['sample_count'] = len(parts) - 1  # Subtract gene ID column
        elif first_lines and ',' in first_lines[0]:
            info['type'] = 'CSV'
            parts = first_lines[0].strip().split(',')
            info['sample_count'] = len(parts) - 1
        
    except Exception as e:
        info['error'] = str(e)
    
    return info

def main():
    print("="*80)
    print("GTEx File Verification")
    print("="*80)
    
    # Find files
    print("\n[1] Searching for GTEx files...")
    files = find_gtex_files()
    
    if not files:
        print("\n❌ No GTEx files found in common locations")
        print("\nPlease specify the location of your downloaded files, or")
        print("move them to: data/gtex/")
        return
    
    print(f"\n✅ Found {len(files)} potential GTEx file(s)")
    
    # Analyze files
    print("\n[2] Analyzing files...")
    file_info = []
    for file_path in files:
        print(f"\n  Analyzing: {file_path.name}")
        info = analyze_file(file_path)
        file_info.append(info)
    
    # Display results
    print("\n" + "="*80)
    print("FILE ANALYSIS RESULTS")
    print("="*80)
    
    vessel_files = []
    other_files = []
    
    for info in file_info:
        if info.get('tissue') in REQUIRED_TISSUES:
            vessel_files.append(info)
        else:
            other_files.append(info)
    
    print(f"\n✅ VESSEL TISSUE FILES ({len(vessel_files)}):")
    print("-" * 80)
    for info in vessel_files:
        print(f"\n  File: {info['name']}")
        print(f"    Location: {info['path']}")
        print(f"    Size: {info['size_mb']:.2f} MB ({info['size']:,} bytes)")
        print(f"    Tissue: {info.get('tissue', 'Unknown')}")
        print(f"    Type: {info.get('type', 'Unknown')}")
        if info.get('gene_count'):
            print(f"    Genes: {info['gene_count']:,}")
        if info.get('sample_count'):
            print(f"    Samples: {info['sample_count']:,}")
    
    if other_files:
        print(f"\n📋 OTHER FILES ({len(other_files)}):")
        print("-" * 80)
        for info in other_files:
            print(f"\n  File: {info['name']}")
            print(f"    Location: {info['path']}")
            print(f"    Size: {info['size_mb']:.2f} MB")
            print(f"    Type: {info.get('type', 'Unknown')}")
    
    # Check requirements
    print("\n" + "="*80)
    print("REQUIREMENTS CHECK")
    print("="*80)
    
    found_tissues = [info.get('tissue') for info in vessel_files if info.get('tissue')]
    
    print("\nRequired vessel tissues:")
    for tissue in REQUIRED_TISSUES:
        status = "✅" if tissue in found_tissues else "❌"
        print(f"  {status} {tissue}")
    
    missing = [t for t in REQUIRED_TISSUES if t not in found_tissues]
    
    if not missing:
        print("\n✅ All required vessel tissues found!")
    else:
        print(f"\n⚠️  Missing {len(missing)} tissue(s):")
        for tissue in missing:
            print(f"   - {tissue}")
        print("\nYou may need to download additional files from GTEx Portal.")
    
    # Recommendations
    print("\n" + "="*80)
    print("RECOMMENDATIONS")
    print("="*80)
    
    if len(files) < 3:
        print("\n📋 You have fewer than 3 vessel tissue files.")
        print("   Expected: 3 vessel tissues (Aorta, Coronary, Tibial)")
        print("   You may also want:")
        print("   - Sample metadata files")
        print("   - Gene annotation files")
    
    if vessel_files:
        print("\n✅ Move vessel tissue files to: data/gtex/")
        print("   This will organize them for analysis.")
    
    print("\n" + "="*80)

if __name__ == "__main__":
    main()






