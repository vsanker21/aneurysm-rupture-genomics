"""
Download GSE186927 using GEOquery-like approach
"""

import requests
from pathlib import Path
import sys

def print_flush(*args, **kwargs):
    print(*args, **kwargs)
    sys.stdout.flush()

gse_dir = Path("data/sc_rnaseq/GSE186927")
gse_dir.mkdir(parents=True, exist_ok=True)

print_flush("Attempting to download GSE186927 files...")

# Try different GEO download methods
urls_to_try = [
    # Method 1: Direct GEO download
    ("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE186927&format=file&file=GSE186927%5Fseries%5Fmatrix%2Etxt%2Egz",
     gse_dir / "GSE186927_series_matrix.txt.gz"),
    
    # Method 2: FTP matrix file
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186927/matrix/GSE186927_series_matrix.txt.gz",
     gse_dir / "GSE186927_series_matrix.txt.gz"),
    
    # Method 3: Try to get file list
    ("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186927&targ=gsm&form=text&view=brief",
     None),  # Just to check
]

for url, output_path in urls_to_try:
    if output_path is None:
        continue
    if output_path.exists():
        print_flush(f"  File already exists: {output_path}")
        continue
        
    print_flush(f"\nTrying: {url}")
    try:
        response = requests.get(url, stream=True, timeout=60, allow_redirects=True)
        if response.status_code == 200:
            total_size = int(response.headers.get('content-length', 0))
            downloaded = 0
            
            with open(output_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        downloaded += len(chunk)
                        if total_size > 0:
                            percent = (downloaded / total_size) * 100
                            print_flush(f"\r  Progress: {percent:.1f}%", end='')
            
            print_flush(f"\n  ✓ Downloaded: {output_path} ({output_path.stat().st_size:,} bytes)")
        else:
            print_flush(f"  ✗ Failed: HTTP {response.status_code}")
    except Exception as e:
        print_flush(f"  ✗ Error: {e}")

print_flush("\nNote: GSE186927 raw counts may require:")
print_flush("  1. Access through EGA (European Genome-Phenome Archive)")
print_flush("  2. Using GEOquery R package: library(GEOquery); getGEO('GSE186927')")
print_flush("  3. Direct contact with authors")


