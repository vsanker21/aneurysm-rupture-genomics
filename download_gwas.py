"""
Download GWAS summary statistics from public sources
"""

import requests
from pathlib import Path
import sys

def print_flush(*args, **kwargs):
    print(*args, **kwargs)
    sys.stdout.flush()

gwas_dir = Path("data/gwas")
gwas_dir.mkdir(parents=True, exist_ok=True)

print_flush("="*70)
print_flush("Downloading GWAS Summary Statistics")
print_flush("="*70)

# Try GWAS Catalog API
print_flush("\n1. Checking GWAS Catalog for Bakker et al. 2020 data...")

# GWAS Catalog REST API
gwas_catalog_url = "https://www.ebi.ac.uk/gwas/api/search/downloads"
params = {
    'efo': 'EFO_0004709',  # Intracranial aneurysm
    'pvalfilter': '5e-8',
}

try:
    response = requests.get(gwas_catalog_url, params=params, timeout=30)
    if response.status_code == 200:
        output_file = gwas_dir / "gwas_catalog_aneurysm.tsv"
        with open(output_file, 'wb') as f:
            f.write(response.content)
        print_flush(f"  ✓ Downloaded GWAS Catalog results: {output_file}")
        print_flush(f"  File size: {output_file.stat().st_size:,} bytes")
    else:
        print_flush(f"  ✗ GWAS Catalog returned status {response.status_code}")
except Exception as e:
    print_flush(f"  ✗ Error accessing GWAS Catalog: {e}")

# Create download instructions for Bakker et al. 2020
bakker_dir = gwas_dir / "bakker_2020"
bakker_dir.mkdir(parents=True, exist_ok=True)

print_flush("\n2. Creating download instructions for Bakker et al. 2020...")
print_flush("   Note: Full summary statistics may require:")
print_flush("   - Journal supplementary materials")
print_flush("   - Direct contact with authors")
print_flush("   - Institutional repositories")

print_flush("\n" + "="*70)
print_flush("GWAS Download Complete")
print_flush("="*70)


