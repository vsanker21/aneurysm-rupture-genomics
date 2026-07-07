"""
Quick check for GTEx downloaded files
Run this and specify the directory where your 6 GTEx files are located
"""

from pathlib import Path
import sys

# Required vessel tissues
REQUIRED_TISSUES = {
    "Artery_Aorta": ["aorta"],
    "Artery_Coronary": ["coronary"],
    "Artery_Tibial": ["tibial"]
}

def check_files_in_directory(directory):
    """Check files in a specific directory"""
    dir_path = Path(directory)
    if not dir_path.exists():
        print(f"❌ Directory not found: {directory}")
        return
    
    print(f"\n📁 Checking directory: {directory}")
    print("="*80)
    
    # Find all data files (exclude README, .md, etc.)
    files = [f for f in dir_path.iterdir() 
             if f.is_file() 
             and f.suffix in ['.gz', '.gct', '.txt', '.tsv', '.csv', '.json']
             and not f.name.endswith('.md')]
    
    if not files:
        print("❌ No GTEx data files found in this directory")
        return
    
    print(f"\n✅ Found {len(files)} file(s):\n")
    
    found_tissues = {}
    other_files = []
    
    for file_path in sorted(files):
        name_lower = file_path.name.lower()
        size_mb = file_path.stat().st_size / (1024 * 1024)
        
        # Check which tissue this file belongs to
        tissue_found = None
        for tissue, keywords in REQUIRED_TISSUES.items():
            if any(kw in name_lower for kw in keywords):
                tissue_found = tissue
                break
        
        file_info = {
            'name': file_path.name,
            'size_mb': size_mb,
            'tissue': tissue_found
        }
        
        if tissue_found:
            if tissue_found not in found_tissues:
                found_tissues[tissue_found] = []
            found_tissues[tissue_found].append(file_info)
        else:
            other_files.append(file_info)
    
    # Display results
    print("VESSEL TISSUE FILES:")
    print("-" * 80)
    for tissue in REQUIRED_TISSUES.keys():
        if tissue in found_tissues:
            print(f"\n  ✅ {tissue}:")
            for f in found_tissues[tissue]:
                print(f"     - {f['name']} ({f['size_mb']:.2f} MB)")
        else:
            print(f"\n  ❌ {tissue}: NOT FOUND")
    
    if other_files:
        print(f"\n\nOTHER FILES ({len(other_files)}):")
        print("-" * 80)
        for f in other_files:
            print(f"  - {f['name']} ({f['size_mb']:.2f} MB)")
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    
    missing = [t for t in REQUIRED_TISSUES.keys() if t not in found_tissues]
    
    if not missing:
        print("\n✅ EXCELLENT! All 3 required vessel tissues found!")
        print(f"   You have {len(files)} file(s) total")
    else:
        print(f"\n⚠️  Missing {len(missing)} required tissue(s):")
        for tissue in missing:
            print(f"   - {tissue}")
        print(f"\n   You have {len(files)} file(s), but need files for all 3 vessel tissues")
    
    if len(files) == 6:
        print("\n📋 You have exactly 6 files - this matches expected count!")
        print("   Likely includes: 3 vessel tissues + metadata/annotation files")
    
    print("\n" + "="*80)
    print("RECOMMENDATIONS")
    print("="*80)
    
    if not missing:
        print("\n✅ All required files present!")
        print("   Next step: Move files to data/gtex/ directory")
        print("   Command: Move-Item <files> data/gtex/")
    else:
        print(f"\n⚠️  You need to download {len(missing)} more file(s)")
        print("   Visit: https://gtexportal.org/home/datasets")
        print("   Download expression data for:", ", ".join(missing))

if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Directory provided as argument
        check_files_in_directory(sys.argv[1])
    else:
        # Check common locations
        print("GTEx File Checker")
        print("="*80)
        print("\nChecking common locations...\n")
        
        locations = [
            "data/gtex",
            str(Path.home() / "Downloads"),
            "."
        ]
        
        for loc in locations:
            dir_path = Path(loc)
            if dir_path.exists():
                files = [f for f in dir_path.iterdir() 
                        if f.is_file() 
                        and f.suffix in ['.gz', '.gct', '.txt', '.tsv', '.csv', '.json']
                        and not f.name.endswith('.md')]
                if files:
                    print(f"Found {len(files)} file(s) in: {loc}")
                    check_files_in_directory(loc)
                    break
        else:
            print("❌ No GTEx files found in common locations")
            print("\nUsage: python check_gtex_downloads.py <directory_path>")
            print("Example: python check_gtex_downloads.py C:\\Users\\YourName\\Downloads")






