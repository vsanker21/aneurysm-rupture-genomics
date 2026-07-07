"""
Direct Dataset Download Script
"""

import os
import sys
import requests
from pathlib import Path

# Force output flushing
def print_flush(*args, **kwargs):
    print(*args, **kwargs)
    sys.stdout.flush()

BASE_DIR = Path("data")
BASE_DIR.mkdir(exist_ok=True)

def download_file(url, output_path, description="file"):
    """Download a file with progress tracking"""
    print_flush(f"\nDownloading {description}...")
    print_flush(f"  URL: {url}")
    print_flush(f"  Destination: {output_path}")
    
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        response = requests.get(url, stream=True, timeout=300)
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
                        print_flush(f"\r  Progress: {percent:.1f}% ({downloaded:,} / {total_size:,} bytes)", end='')
        
        print_flush(f"\n  ✓ Download complete: {output_path}")
        print_flush(f"  File size: {output_path.stat().st_size:,} bytes")
        return True
        
    except Exception as e:
        print_flush(f"\n  ✗ Download failed: {e}")
        return False

print_flush("="*70)
print_flush("Cerebral Aneurysm Genomics Project - Dataset Download")
print_flush("="*70)

# 1. Download GSE186927 (single-cell RNA-seq)
print_flush("\n" + "="*70)
print_flush("1. Downloading GSE186927: Single-cell RNA-seq")
print_flush("="*70)

gse_dir = BASE_DIR / "sc_rnaseq" / "GSE186927"
gse_dir.mkdir(parents=True, exist_ok=True)

# Try direct FTP download from GEO
geo_urls = [
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186927/matrix/GSE186927_series_matrix.txt.gz",
     gse_dir / "GSE186927_series_matrix.txt.gz"),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186927/suppl/GSE186927_raw_counts.tsv.gz",
     gse_dir / "GSE186927_raw_counts.tsv.gz"),
]

for url, output_path in geo_urls:
    if not output_path.exists():
        download_file(url, output_path, f"GSE186927 {output_path.name}")
    else:
        print_flush(f"  File already exists: {output_path}")

# 2. Create README files for datasets that need manual download
print_flush("\n" + "="*70)
print_flush("2. Creating documentation for manual downloads")
print_flush("="*70)

# GWAS README
gwas_dir = BASE_DIR / "gwas" / "bakker_2020"
gwas_dir.mkdir(parents=True, exist_ok=True)
readme_gwas = gwas_dir / "README.md"
readme_gwas.write_text("""# Bakker et al. 2020 GWAS Summary Statistics

## Publication
Bakker et al. (2020) - Meta-analysis of intracranial aneurysm GWAS
- N = 10,754 cases + 306,882 controls
- 17 risk loci identified
- Published in European Heart Journal

## Download Instructions

### Option 1: GWAS Catalog
1. Visit: https://www.ebi.ac.uk/gwas/
2. Search for "intracranial aneurysm" or "Bakker 2020"
3. Download summary statistics

### Option 2: Journal Supplementary Materials
1. Access European Heart Journal article
2. Download supplementary data files
3. Look for summary statistics files

### Option 3: Authors' Repository
- Contact authors for direct access
- Check GitHub or institutional repositories

## Expected Files
- Summary statistics file (usually .txt or .gz format)
- Should contain: SNP, CHR, BP, A1, A2, OR, SE, P-value, etc.
""")
print_flush(f"  ✓ Created: {readme_gwas}")

# GTEx README
gtex_dir = BASE_DIR / "gtex"
gtex_dir.mkdir(parents=True, exist_ok=True)
readme_gtex = gtex_dir / "README.md"
readme_gtex.write_text("""# GTEx Normal Vessel Reference Transcriptomes

## Data Source
GTEx (Genotype-Tissue Expression) Project
- Portal: https://gtexportal.org/
- Version: GTEx v8 (recommended)

## Download Instructions

### Option 1: GTEx Portal (Public Data)
1. Visit: https://gtexportal.org/home/datasets
2. Navigate to "Gene Expression" section
3. Download vessel tissue expression data:
   - Artery - Aorta
   - Artery - Coronary
   - Artery - Tibial
4. Select format: TPM (Transcripts Per Million) or raw counts

### Option 2: GTEx API
Use GTEx REST API to query vessel tissue expression:
- API docs: https://gtexportal.org/api/v2/redoc
- Endpoint: /api/v2/expression/medianGeneExpression

### Option 3: dbGaP (Full Access)
1. Apply for dbGaP access: https://www.ncbi.nlm.nih.gov/gap/
2. Download full GTEx dataset
3. Extract vessel tissue samples

## Expected Files
- Gene expression matrix (vessel tissues)
- Sample metadata
- Gene annotation files

## Tissue Types Needed
- Artery - Aorta
- Artery - Coronary  
- Artery - Tibial
- (Other vascular tissues as available)
""")
print_flush(f"  ✓ Created: {readme_gtex}")

# UK Biobank README
ukbb_dir = BASE_DIR / "uk_biobank"
ukbb_dir.mkdir(parents=True, exist_ok=True)
readme_ukbb = ukbb_dir / "README.md"
readme_ukbb.write_text("""# UK Biobank GWAS Summary Statistics

## Data Source
UK Biobank - Large-scale genetic and health database
- Website: https://www.ukbiobank.ac.uk/
- Neale Lab: http://www.nealelab.is/uk-biobank/

## Download Instructions

### Option 1: Neale Lab (Easiest - Processed Summary Stats)
1. Visit: http://www.nealelab.is/uk-biobank/
2. Browse available phenotypes
3. Search for "intracranial aneurysm" or related vascular phenotypes
4. Download summary statistics files

Direct link: http://www.nealelab.is/uk-biobank/

### Option 2: UK Biobank Official
1. Visit: https://www.ukbiobank.ac.uk/
2. Apply for access (if needed)
3. Access data through UK Biobank Showcase
4. Download GWAS summary statistics

### Option 3: GWAS Catalog
1. Visit: https://www.ebi.ac.uk/gwas/
2. Search for UK Biobank studies
3. Filter by phenotype: intracranial aneurysm

## Expected Files
- Summary statistics (usually .txt or .tsv format)
- Should contain: variant ID, chromosome, position, alleles, effect sizes, p-values

## Alternative: Use Bakker et al. 2020 GWAS
The Bakker et al. 2020 meta-analysis already includes UK Biobank data,
so you may not need separate UK Biobank files if you have the Bakker summary stats.
""")
print_flush(f"  ✓ Created: {readme_ukbb}")

print_flush("\n" + "="*70)
print_flush("Download Summary")
print_flush("="*70)
print_flush("✓ GSE186927: Attempted download (check files)")
print_flush("✓ Documentation created for manual downloads")
print_flush("="*70)


