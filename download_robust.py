"""
Robust programmatic download with better error handling
"""

import requests
import json
import time
from pathlib import Path
import sys

# Force unbuffered output
sys.stdout.reconfigure(line_buffering=True)
sys.stderr.reconfigure(line_buffering=True)

def print_now(msg):
    """Print with immediate flush"""
    # Remove Unicode characters that cause issues in Windows console
    msg = msg.replace('✓', '[OK]').replace('✗', '[FAIL]').replace('📥', '[DL]').replace('📄', '[PAGE]').replace('ℹ️', '[INFO]')
    print(msg, flush=True)

BASE_DIR = Path("data")
BASE_DIR.mkdir(exist_ok=True)

def safe_download(url, output_path, description="file", timeout=60):
    """Download with timeout and error handling"""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    if output_path.exists() and output_path.stat().st_size > 0:
        print_now(f"  ✓ Exists: {output_path.name}")
        return True
    
    try:
        print_now(f"  📥 {description}")
        response = requests.get(url, stream=True, timeout=timeout, allow_redirects=True)
        
        if response.status_code == 404:
            print_now(f"  ✗ Not found (404)")
            return False
            
        response.raise_for_status()
        
        total_size = int(response.headers.get('content-length', 0))
        downloaded = 0
        
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total_size > 0 and downloaded % (1024*1024) == 0:  # Print every MB
                        percent = (downloaded / total_size) * 100
                        print_now(f"    {percent:.1f}% ({downloaded//1024//1024} MB)")
        
        size = output_path.stat().st_size
        print_now(f"  ✅ {size:,} bytes")
        return True
        
    except requests.exceptions.Timeout:
        print_now(f"  ✗ Timeout")
        if output_path.exists():
            output_path.unlink()
        return False
    except Exception as e:
        print_now(f"  ✗ Error: {str(e)[:80]}")
        if output_path.exists():
            output_path.unlink()
        return False

print_now("="*80)
print_now("PROGRAMMATIC DATASET DOWNLOAD")
print_now("="*80)

# ============================================================================
# 1. GSE186927 - Try multiple file patterns
# ============================================================================
print_now("\n[1/4] GSE186927 - Single-cell RNA-seq")
print_now("-" * 80)

gse_dir = BASE_DIR / "sc_rnaseq" / "GSE186927"
gse_dir.mkdir(parents=True, exist_ok=True)

# Known working: series matrix
safe_download(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186927/matrix/GSE186927_series_matrix.txt.gz",
    gse_dir / "GSE186927_series_matrix.txt.gz",
    "Series matrix"
)

# Try common file name patterns
patterns = [
    "GSE186927_raw_counts.tsv.gz",
    "GSE186927_counts.tsv.gz", 
    "GSE186927_expression.tsv.gz",
    "GSE186927_processed.h5ad",
    "GSE186927_raw_counts.txt.gz",
]

for pattern in patterns:
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186927/suppl/{pattern}"
    safe_download(url, gse_dir / pattern, f"Trying {pattern}", timeout=30)

# ============================================================================
# 2. GWAS Catalog - Intracranial Aneurysm
# ============================================================================
print_now("\n[2/4] GWAS Summary Statistics")
print_now("-" * 80)

gwas_dir = BASE_DIR / "gwas"
gwas_dir.mkdir(parents=True, exist_ok=True)

# GWAS Catalog downloads
gwas_urls = [
    ("https://www.ebi.ac.uk/gwas/api/search/downloads?efo=EFO_0004709&pvalfilter=5e-8",
     gwas_dir / "gwas_catalog_aneurysm_5e8.tsv"),
    ("https://www.ebi.ac.uk/gwas/api/search/downloads?efo=EFO_0004709",
     gwas_dir / "gwas_catalog_aneurysm_all.tsv"),
]

for url, output in gwas_urls:
    safe_download(url, output, f"GWAS Catalog: {output.name}", timeout=120)

# ============================================================================
# 3. GTEx - Sample attributes and metadata (skip large expression file)
# ============================================================================
print_now("\n[3/4] GTEx - Vessel Reference Data")
print_now("-" * 80)

gtex_dir = BASE_DIR / "gtex"
gtex_dir.mkdir(parents=True, exist_ok=True)

# GTEx metadata files (small, fast)
gtex_files = [
    ("https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
     gtex_dir / "GTEx_v8_sample_attributes.txt"),
    ("https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
     gtex_dir / "GTEx_v8_subject_phenotypes.txt"),
]

for url, output in gtex_files:
    safe_download(url, output, f"GTEx: {output.name}", timeout=120)

print_now("\n  ℹ️  GTEx expression matrix (~2GB) skipped - download manually if needed:")
print_now("     https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")

# ============================================================================
# 4. UK Biobank - Save pages for manual inspection
# ============================================================================
print_now("\n[4/4] UK Biobank - Checking available resources")
print_now("-" * 80)

ukbb_dir = BASE_DIR / "uk_biobank"
ukbb_dir.mkdir(parents=True, exist_ok=True)

# Save Neale Lab page
try:
    print_now("  📄 Fetching Neale Lab page...")
    response = requests.get("http://www.nealelab.is/uk-biobank", timeout=30)
    if response.status_code == 200:
        (ukbb_dir / "neale_lab_page.html").write_bytes(response.content)
        print_now(f"  ✅ Saved: neale_lab_page.html ({len(response.content):,} bytes)")
except Exception as e:
    print_now(f"  ✗ Failed: {e}")

# ============================================================================
# Summary
# ============================================================================
print_now("\n" + "="*80)
print_now("DOWNLOAD COMPLETE")
print_now("="*80)

# Count files
data_files = [f for f in BASE_DIR.rglob("*") 
              if f.is_file() and f.suffix in ['.gz', '.tsv', '.txt', '.h5ad', '.json', '.html']
              and not f.name.endswith('README.md')]

print_now(f"\n✅ Downloaded {len(data_files)} data files")
print_now("\nFiles by dataset:")
for subdir in ['blood_rnaseq', 'sc_rnaseq', 'gwas', 'gtex', 'uk_biobank']:
    dir_path = BASE_DIR / subdir
    if dir_path.exists():
        files = [f for f in dir_path.rglob("*") 
                if f.is_file() and f.suffix in ['.gz', '.tsv', '.txt', '.h5ad', '.json', '.html']]
        if files:
            print_now(f"\n  {subdir}/")
            for f in files:
                size = f.stat().st_size
                print_now(f"    {f.name} ({size:,} bytes)")

print_now("\n" + "="*80)

