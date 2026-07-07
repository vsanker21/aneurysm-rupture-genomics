"""
Final download script with corrected URLs and alternative methods
"""

import requests
from pathlib import Path
import sys
import json
import time

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
        p(f"  [ERROR] {str(e)[:60]}")
        if output.exists():
            output.unlink()
        return False

p("="*70)
p("PROGRAMMATIC DATASET DOWNLOAD - FINAL ATTEMPT")
p("="*70)

# 1. GSE186927 - Already have series matrix, try GEOquery approach info
p("\n[1] GSE186927 - Single-cell RNA-seq")
p("-"*70)
gse_dir = BASE_DIR / "sc_rnaseq" / "GSE186927"
gse_dir.mkdir(parents=True, exist_ok=True)

p("  [INFO] Series matrix already downloaded")
p("  [INFO] Raw counts not available via direct FTP")
p("  [INFO] Use R GEOquery: library(GEOquery); getGEO('GSE186927')")

# 2. GWAS Catalog - Try different API endpoints
p("\n[2] GWAS Catalog - Intracranial Aneurysm")
p("-"*70)
gwas_dir = BASE_DIR / "gwas"
gwas_dir.mkdir(parents=True, exist_ok=True)

# Try simpler GWAS Catalog download
gwas_urls = [
    ("https://www.ebi.ac.uk/gwas/api/search/downloads?efo=EFO_0004709",
     gwas_dir / "gwas_catalog_aneurysm_all.tsv"),
    # Try without EFO filter
    ("https://www.ebi.ac.uk/gwas/api/search/downloads?q=intracranial+aneurysm",
     gwas_dir / "gwas_catalog_aneurysm_search.tsv"),
]

for url, output in gwas_urls:
    download(url, output, f"GWAS: {output.name}", timeout=120)

# Try to get study list
try:
    p("  [DL] GWAS Catalog study list")
    r = requests.get("https://www.ebi.ac.uk/gwas/api/search/studies?q=intracranial+aneurysm&size=20", timeout=60)
    if r.status_code == 200:
        data = r.json()
        output = gwas_dir / "gwas_catalog_studies.json"
        with open(output, 'w') as f:
            json.dump(data, f, indent=2)
        p(f"  [OK] Found {len(data.get('_embedded', {}).get('studies', []))} studies")
except Exception as e:
    p(f"  [ERROR] {e}")

# 3. GTEx - Try updated URLs
p("\n[3] GTEx - Metadata files")
p("-"*70)
gtex_dir = BASE_DIR / "gtex"
gtex_dir.mkdir(parents=True, exist_ok=True)

# Try v8 URLs
gtex_urls = [
    ("https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
     gtex_dir / "GTEx_v8_sample_attributes.txt"),
    ("https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
     gtex_dir / "GTEx_v8_subject_phenotypes.txt"),
    # Try v9 if v8 doesn't work
    ("https://storage.googleapis.com/gtex_analysis_v9/annotations/GTEx_Analysis_v9_Annotations_SampleAttributesDS.txt",
     gtex_dir / "GTEx_v9_sample_attributes.txt"),
]

for url, output in gtex_urls:
    if not download(url, output, f"GTEx: {output.name}", timeout=120):
        time.sleep(1)  # Brief pause between attempts

# Try GTEx Portal API for vessel tissue data
p("\n  [INFO] Attempting GTEx API for vessel tissue expression...")
try:
    # Artery - Aorta
    api_url = "https://gtexportal.org/api/v2/expression/medianGeneExpression"
    params = {"tissueSiteDetailId": "Artery_Aorta", "datasetId": "gtex_v8"}
    r = requests.get(api_url, params=params, timeout=60)
    if r.status_code == 200:
        data = r.json()
        output = gtex_dir / "GTEx_v8_Artery_Aorta_expression.json"
        with open(output, 'w') as f:
            json.dump(data, f, indent=2)
        p(f"  [OK] Artery_Aorta: {len(data.get('data', []))} genes")
except Exception as e:
    p(f"  [ERROR] GTEx API: {e}")

p("\n  [INFO] GTEx full expression matrix (2GB):")
p("        https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")

# 4. UK Biobank - Already downloaded Neale Lab page
p("\n[4] UK Biobank")
p("-"*70)
p("  [OK] Neale Lab page already downloaded")
p("  [INFO] Check neale_lab_page.html for download links")

# Try Pan-UK Biobank
ukbb_dir = BASE_DIR / "uk_biobank"
try:
    p("  [DL] Pan-UK Biobank page")
    r = requests.get("https://pan.ukbb.broadinstitute.org/downloads", timeout=30)
    if r.status_code == 200:
        (ukbb_dir / "pan_uk_biobank_page.html").write_bytes(r.content)
        p(f"  [OK] {len(r.content):,} bytes")
except Exception as e:
    p(f"  [ERROR] {e}")

# Summary
p("\n" + "="*70)
p("FINAL SUMMARY")
p("="*70)

files = []
for ext in ['.gz', '.tsv', '.txt', '.html', '.json']:
    files.extend(list(BASE_DIR.rglob(f"*{ext}")))

data_files = [f for f in files if 'README' not in f.name and f.stat().st_size > 0]

p(f"\nTotal data files: {len(data_files)}\n")
for subdir in ['blood_rnaseq', 'sc_rnaseq', 'gwas', 'gtex', 'uk_biobank']:
    dir_files = [f for f in data_files if subdir in str(f)]
    if dir_files:
        p(f"\n{subdir}/:")
        for f in sorted(dir_files):
            size = f.stat().st_size
            rel_path = f.relative_to(BASE_DIR.parent)
            p(f"  {rel_path.name} ({size:,} bytes)")

p("\n" + "="*70)
p("MANUAL DOWNLOADS NEEDED:")
p("="*70)
p("1. GSE186927 raw counts - Use R GEOquery or contact authors")
p("2. Bakker 2020 GWAS - Download from journal or GWAS Catalog website")
p("3. GTEx expression matrix (2GB) - Download from Google Cloud Storage")
p("4. UK Biobank summary stats - Download from Neale Lab or Pan-UK Biobank")
p("="*70)

