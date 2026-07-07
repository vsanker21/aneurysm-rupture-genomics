"""
Aggressive programmatic download of all available datasets
Tries multiple methods and sources for each dataset
"""

import requests
import json
import time
from pathlib import Path
import sys
import gzip
import shutil

def print_flush(*args, **kwargs):
    print(*args, **kwargs)
    sys.stdout.flush()

BASE_DIR = Path("data")
BASE_DIR.mkdir(exist_ok=True)

def download_file(url, output_path, description="file", timeout=300):
    """Download a file with progress tracking"""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    if output_path.exists():
        size = output_path.stat().st_size
        if size > 0:
            print_flush(f"  ⏭️  Already exists: {output_path} ({size:,} bytes)")
            return True
    
    try:
        print_flush(f"  📥 Downloading: {description}")
        response = requests.get(url, stream=True, timeout=timeout, allow_redirects=True)
        response.raise_for_status()
        
        total_size = int(response.headers.get('content-length', 0))
        downloaded = 0
        
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total_size > 0:
                        percent = (downloaded / total_size) * 100
                        print_flush(f"\r    Progress: {percent:.1f}% ({downloaded:,} / {total_size:,} bytes)", end='')
        
        print_flush(f"\n  ✅ Success: {output_path} ({output_path.stat().st_size:,} bytes)")
        return True
        
    except Exception as e:
        print_flush(f"\n  ❌ Failed: {str(e)[:100]}")
        if output_path.exists():
            output_path.unlink()  # Remove partial file
        return False

print_flush("="*80)
print_flush("AGGRESSIVE PROGRAMMATIC DATASET DOWNLOAD")
print_flush("="*80)

# ============================================================================
# 1. GSE186927 - Single-cell RNA-seq (Try multiple methods)
# ============================================================================
print_flush("\n" + "="*80)
print_flush("1. GSE186927 - Single-cell RNA-seq")
print_flush("="*80)

gse_dir = BASE_DIR / "sc_rnaseq" / "GSE186927"
gse_dir.mkdir(parents=True, exist_ok=True)

# Try multiple GEO download URLs
geo_urls = [
    # Series matrix (metadata)
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186927/matrix/GSE186927_series_matrix.txt.gz",
     gse_dir / "GSE186927_series_matrix.txt.gz"),
    
    # Try different file name patterns
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186927/suppl/GSE186927_raw_counts.tsv.gz",
     gse_dir / "GSE186927_raw_counts.tsv.gz"),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186927/suppl/GSE186927_raw_counts.txt.gz",
     gse_dir / "GSE186927_raw_counts.txt.gz"),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186927/suppl/GSE186927_counts.tsv.gz",
     gse_dir / "GSE186927_counts.tsv.gz"),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186927/suppl/GSE186927_expression_matrix.tsv.gz",
     gse_dir / "GSE186927_expression_matrix.tsv.gz"),
    
    # Try h5ad format (processed)
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186927/suppl/GSE186927_processed.h5ad",
     gse_dir / "GSE186927_processed.h5ad"),
    
    # Try GEO download API
    ("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE186927&format=file&file=GSE186927%5Fraw%5Fcounts%2Etsv%2Egz",
     gse_dir / "GSE186927_raw_counts_via_api.tsv.gz"),
]

for url, output_path in geo_urls:
    download_file(url, output_path, f"GSE186927 {output_path.name}")

# Try to get file list from GEO
print_flush("\n  🔍 Checking GEO for available files...")
try:
    geo_page = requests.get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186927&targ=self&form=text&view=full", timeout=30)
    if geo_page.status_code == 200:
        # Save the page for manual inspection
        with open(gse_dir / "GSE186927_geo_page.txt", 'w', encoding='utf-8') as f:
            f.write(geo_page.text)
        print_flush(f"  ✅ Saved GEO page: {gse_dir / 'GSE186927_geo_page.txt'}")
        # Look for file URLs in the page
        if 'suppl' in geo_page.text.lower() or 'ftp' in geo_page.text.lower():
            print_flush("  ℹ️  Check GEO page for additional file links")
except Exception as e:
    print_flush(f"  ⚠️  Could not fetch GEO page: {e}")

# ============================================================================
# 2. GWAS Data - Bakker et al. 2020 (Try GWAS Catalog API)
# ============================================================================
print_flush("\n" + "="*80)
print_flush("2. GWAS Summary Statistics - Bakker et al. 2020")
print_flush("="*80)

gwas_dir = BASE_DIR / "gwas"
gwas_dir.mkdir(parents=True, exist_ok=True)

# Try GWAS Catalog API
print_flush("\n  🔍 Querying GWAS Catalog API...")
gwas_urls = [
    # GWAS Catalog search for intracranial aneurysm
    ("https://www.ebi.ac.uk/gwas/api/search/downloads?efo=EFO_0004709&pvalfilter=5e-8&genomiccontexts=true",
     gwas_dir / "gwas_catalog_aneurysm_all.tsv"),
    
    # Try without p-value filter
    ("https://www.ebi.ac.uk/gwas/api/search/downloads?efo=EFO_0004709",
     gwas_dir / "gwas_catalog_aneurysm_no_filter.tsv"),
    
    # Try with different EFO terms
    ("https://www.ebi.ac.uk/gwas/api/search/downloads?efo=EFO_0004709&pvalfilter=1e-5",
     gwas_dir / "gwas_catalog_aneurysm_1e5.tsv"),
]

for url, output_path in gwas_urls:
    download_file(url, output_path, f"GWAS Catalog {output_path.name}")

# Try to find Bakker 2020 specific study
print_flush("\n  🔍 Searching for Bakker 2020 study...")
try:
    # Search for studies
    search_url = "https://www.ebi.ac.uk/gwas/api/search/studies"
    params = {"q": "Bakker aneurysm", "size": 10}
    response = requests.get(search_url, params=params, timeout=30)
    if response.status_code == 200:
        studies = response.json()
        if studies.get('_embedded', {}).get('studies'):
            study_file = gwas_dir / "bakker_2020" / "gwas_catalog_studies.json"
            study_file.parent.mkdir(parents=True, exist_ok=True)
            with open(study_file, 'w') as f:
                json.dump(studies, f, indent=2)
            print_flush(f"  ✅ Found studies, saved to: {study_file}")
            
            # Try to download summary stats for each study
            for study in studies['_embedded']['studies']:
                study_id = study.get('studyId', '')
                if 'aneurysm' in study.get('diseaseTrait', {}).get('trait', '').lower():
                    print_flush(f"  📋 Found study: {study_id}")
except Exception as e:
    print_flush(f"  ⚠️  GWAS Catalog search failed: {e}")

# ============================================================================
# 3. GTEx Data - Try API and direct downloads
# ============================================================================
print_flush("\n" + "="*80)
print_flush("3. GTEx - Normal Vessel Reference Transcriptomes")
print_flush("="*80)

gtex_dir = BASE_DIR / "gtex"
gtex_dir.mkdir(parents=True, exist_ok=True)

# GTEx public data URLs
gtex_urls = [
    # GTEx v8 gene expression (very large ~2GB, but we'll try)
    ("https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz",
     gtex_dir / "GTEx_v8_gene_tpm.gct.gz"),
    
    # GTEx v8 sample attributes
    ("https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
     gtex_dir / "GTEx_v8_sample_attributes.txt"),
    
    # GTEx v8 subject phenotypes
    ("https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
     gtex_dir / "GTEx_v8_subject_phenotypes.txt"),
]

print_flush("\n  ⚠️  Note: GTEx expression file is ~2GB. Download may take time.")
for url, output_path in gtex_urls:
    # Skip the large expression file for now, or download it
    if 'gene_tpm' in str(output_path):
        print_flush(f"  ⏭️  Skipping large file (2GB): {output_path.name}")
        print_flush(f"     Download manually if needed: {url}")
        continue
    download_file(url, output_path, f"GTEx {output_path.name}", timeout=600)

# Try GTEx API for vessel-specific data
print_flush("\n  🔍 Querying GTEx API for vessel tissue data...")
try:
    # GTEx REST API endpoint for median expression
    gtex_api_url = "https://gtexportal.org/api/v2/expression/medianGeneExpression"
    # Vessel tissues: Artery - Aorta, Artery - Coronary, Artery - Tibial
    vessel_tissues = ["Artery_Aorta", "Artery_Coronary", "Artery_Tibial"]
    
    for tissue in vessel_tissues:
        params = {"tissueSiteDetailId": tissue, "datasetId": "gtex_v8"}
        try:
            response = requests.get(gtex_api_url, params=params, timeout=60)
            if response.status_code == 200:
                data = response.json()
                output_file = gtex_dir / f"GTEx_v8_{tissue}_median_expression.json"
                with open(output_file, 'w') as f:
                    json.dump(data, f, indent=2)
                print_flush(f"  ✅ Downloaded: {tissue} ({len(data.get('data', []))} genes)")
            else:
                print_flush(f"  ⚠️  API returned {response.status_code} for {tissue}")
        except Exception as e:
            print_flush(f"  ⚠️  API error for {tissue}: {e}")
        time.sleep(1)  # Be nice to the API
except Exception as e:
    print_flush(f"  ⚠️  GTEx API query failed: {e}")

# ============================================================================
# 4. UK Biobank - Try Neale Lab and Pan-UK Biobank
# ============================================================================
print_flush("\n" + "="*80)
print_flush("4. UK Biobank GWAS Summary Statistics")
print_flush("="*80)

ukbb_dir = BASE_DIR / "uk_biobank"
ukbb_dir.mkdir(parents=True, exist_ok=True)

# Try Neale Lab UK Biobank
print_flush("\n  🔍 Checking Neale Lab UK Biobank...")
try:
    # Neale Lab provides phenotype list
    neale_url = "http://www.nealelab.is/uk-biobank"
    response = requests.get(neale_url, timeout=30)
    if response.status_code == 200:
        # Save the page to check for download links
        with open(ukbb_dir / "neale_lab_page.html", 'w', encoding='utf-8') as f:
            f.write(response.text)
        print_flush(f"  ✅ Saved Neale Lab page: {ukbb_dir / 'neale_lab_page.html'}")
        print_flush("  ℹ️  Check page for aneurysm-related phenotype download links")
except Exception as e:
    print_flush(f"  ⚠️  Neale Lab page access failed: {e}")

# Try Pan-UK Biobank
print_flush("\n  🔍 Checking Pan-UK Biobank...")
try:
    pan_uk_url = "https://pan.ukbb.broadinstitute.org/downloads"
    response = requests.get(pan_uk_url, timeout=30)
    if response.status_code == 200:
        with open(ukbb_dir / "pan_uk_biobank_page.html", 'w', encoding='utf-8') as f:
            f.write(response.text)
        print_flush(f"  ✅ Saved Pan-UK Biobank page: {ukbb_dir / 'pan_uk_biobank_page.html'}")
except Exception as e:
    print_flush(f"  ⚠️  Pan-UK Biobank page access failed: {e}")

# ============================================================================
# Summary
# ============================================================================
print_flush("\n" + "="*80)
print_flush("DOWNLOAD SUMMARY")
print_flush("="*80)

# Count downloaded files
all_files = list(BASE_DIR.rglob("*"))
files = [f for f in all_files if f.is_file() and not f.name.endswith('.md') and not f.name.endswith('.html') and not f.name.endswith('.txt')]

print_flush(f"\n✅ Total data files downloaded: {len(files)}")
print_flush("\nFiles by directory:")
for subdir in ['blood_rnaseq', 'sc_rnaseq', 'gwas', 'gtex', 'uk_biobank']:
    dir_path = BASE_DIR / subdir
    if dir_path.exists():
        files_in_dir = [f for f in dir_path.rglob("*") if f.is_file() and f.suffix not in ['.md', '.html', '.txt']]
        if files_in_dir:
            print_flush(f"  {subdir}: {len(files_in_dir)} files")
            for f in files_in_dir[:5]:  # Show first 5
                size = f.stat().st_size
                print_flush(f"    - {f.name} ({size:,} bytes)")

print_flush("\n" + "="*80)
print_flush("Next Steps:")
print_flush("="*80)
print_flush("1. Check downloaded files in data/ directories")
print_flush("2. For large files (GTEx 2GB), download manually if needed")
print_flush("3. Review HTML/text files saved for additional download links")
print_flush("4. Use R GEOquery for GSE186927 if raw counts still missing")
print_flush("="*80)


