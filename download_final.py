"""
Final robust download script - handles all errors gracefully
"""

import requests
from pathlib import Path
import sys
import time

# Simple print function
def p(msg):
    try:
        print(msg)
        sys.stdout.flush()
    except:
        pass

BASE_DIR = Path("data")
BASE_DIR.mkdir(exist_ok=True)

def download(url, output, desc="file", timeout=30):
    """Download with error handling"""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    
    if output.exists() and output.stat().st_size > 100:
        p(f"  [SKIP] {desc} (exists)")
        return True
    
    try:
        p(f"  [DL] {desc}")
        r = requests.get(url, stream=True, timeout=timeout, allow_redirects=True)
        
        if r.status_code == 404:
            p(f"  [404] Not found")
            return False
        
        r.raise_for_status()
        
        with open(output, 'wb') as f:
            for chunk in r.iter_content(8192):
                if chunk:
                    f.write(chunk)
        
        size = output.stat().st_size
        p(f"  [OK] {size:,} bytes")
        return True
        
    except requests.exceptions.Timeout:
        p(f"  [TIMEOUT]")
        if output.exists():
            output.unlink()
        return False
    except Exception as e:
        p(f"  [ERROR] {str(e)[:50]}")
        if output.exists():
            output.unlink()
        return False

p("="*70)
p("DOWNLOADING AVAILABLE DATASETS")
p("="*70)

# 1. GSE186927
p("\n[1] GSE186927 - Single-cell RNA-seq")
p("-"*70)
gse_dir = BASE_DIR / "sc_rnaseq" / "GSE186927"
gse_dir.mkdir(parents=True, exist_ok=True)

download(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186927/matrix/GSE186927_series_matrix.txt.gz",
    gse_dir / "GSE186927_series_matrix.txt.gz",
    "Series matrix"
)

# Try different file patterns
patterns = [
    "GSE186927_raw_counts.tsv.gz",
    "GSE186927_counts.tsv.gz",
    "GSE186927_expression.tsv.gz",
]
for pattern in patterns:
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186927/suppl/{pattern}"
    download(url, gse_dir / pattern, pattern, timeout=20)

# 2. GWAS Catalog
p("\n[2] GWAS Catalog - Intracranial Aneurysm")
p("-"*70)
gwas_dir = BASE_DIR / "gwas"
gwas_dir.mkdir(parents=True, exist_ok=True)

download(
    "https://www.ebi.ac.uk/gwas/api/search/downloads?efo=EFO_0004709&pvalfilter=5e-8",
    gwas_dir / "gwas_catalog_aneurysm.tsv",
    "GWAS Catalog results",
    timeout=120
)

# 3. GTEx metadata
p("\n[3] GTEx - Metadata files")
p("-"*70)
gtex_dir = BASE_DIR / "gtex"
gtex_dir.mkdir(parents=True, exist_ok=True)

download(
    "https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
    gtex_dir / "GTEx_v8_sample_attributes.txt",
    "Sample attributes",
    timeout=120
)

download(
    "https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
    gtex_dir / "GTEx_v8_subject_phenotypes.txt",
    "Subject phenotypes",
    timeout=120
)

p("\n  [INFO] GTEx expression matrix (2GB) - download manually:")
p("        https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")

# 4. UK Biobank pages
p("\n[4] UK Biobank - Resource pages")
p("-"*70)
ukbb_dir = BASE_DIR / "uk_biobank"
ukbb_dir.mkdir(parents=True, exist_ok=True)

try:
    p("  [DL] Neale Lab page")
    r = requests.get("http://www.nealelab.is/uk-biobank", timeout=30)
    if r.status_code == 200:
        (ukbb_dir / "neale_lab_page.html").write_bytes(r.content)
        p(f"  [OK] {len(r.content):,} bytes")
except Exception as e:
    p(f"  [ERROR] {e}")

# Summary
p("\n" + "="*70)
p("SUMMARY")
p("="*70)

files = []
for ext in ['.gz', '.tsv', '.txt', '.html']:
    files.extend(list(BASE_DIR.rglob(f"*{ext}")))

# Filter out READMEs
data_files = [f for f in files if 'README' not in f.name and f.stat().st_size > 0]

p(f"\nDownloaded {len(data_files)} data files:\n")
for f in sorted(data_files):
    size = f.stat().st_size
    rel_path = f.relative_to(BASE_DIR.parent)
    p(f"  {rel_path} ({size:,} bytes)")

p("\n" + "="*70)


