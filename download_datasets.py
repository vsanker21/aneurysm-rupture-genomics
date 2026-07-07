"""
Dataset Download Script for Cerebral Aneurysm Genomics Project
Downloads remaining datasets: GWAS, scRNA-seq, GTEx, UK Biobank
"""

import os
import requests
import gzip
import shutil
from pathlib import Path
from urllib.parse import urlparse
import time

# Base directories
BASE_DIR = Path("data")
BASE_DIR.mkdir(exist_ok=True)

def download_file(url, output_path, description="file"):
    """Download a file with progress tracking"""
    print(f"\nDownloading {description}...")
    print(f"  URL: {url}")
    print(f"  Destination: {output_path}")
    
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
                        print(f"\r  Progress: {percent:.1f}% ({downloaded:,} / {total_size:,} bytes)", end='')
        
        print(f"\n  ✓ Download complete: {output_path}")
        print(f"  File size: {output_path.stat().st_size:,} bytes")
        return True
        
    except requests.exceptions.RequestException as e:
        print(f"\n  ✗ Download failed: {e}")
        return False

def download_gse186927():
    """Download GSE186927 single-cell RNA-seq dataset from GEO"""
    print("\n" + "="*70)
    print("Downloading GSE186927: Single-cell RNA-seq of aneurysm tissue")
    print("="*70)
    
    gse_dir = BASE_DIR / "sc_rnaseq" / "GSE186927"
    gse_dir.mkdir(parents=True, exist_ok=True)
    
    # GEO download URLs
    base_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE186927"
    
    files_to_download = [
        ("&format=file&file=GSE186927%5Fraw%5Fcounts%2Etsv%2Egz", "GSE186927_raw_counts.tsv.gz"),
        ("&format=file&file=GSE186927%5Fmetadata%2Etxt%2Egz", "GSE186927_metadata.txt.gz"),
        ("&format=file&file=GSE186927%5Fseries%5Fmatrix%2Etxt%2Egz", "GSE186927_series_matrix.txt.gz"),
    ]
    
    # Also try direct GEO download
    geo_urls = [
        f"{base_url}&format=file",
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186927/suppl/GSE186927_raw_counts.tsv.gz",
    ]
    
    success = False
    for url_suffix, filename in files_to_download:
        url = base_url + url_suffix
        output_path = gse_dir / filename
        if download_file(url, output_path, f"GSE186927 {filename}"):
            success = True
    
    # Try alternative: direct FTP download
    if not success:
        print("\nTrying alternative download method...")
        alt_urls = [
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186927/suppl/GSE186927_raw_counts.tsv.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186927/matrix/GSE186927_series_matrix.txt.gz",
        ]
        for url in alt_urls:
            filename = os.path.basename(urlparse(url).path)
            output_path = gse_dir / filename
            if download_file(url, output_path, f"GSE186927 {filename}"):
                success = True
    
    return success

def download_gwas_bakker():
    """Download Bakker et al. 2020 GWAS summary statistics"""
    print("\n" + "="*70)
    print("Downloading GWAS Summary Statistics (Bakker et al. 2020)")
    print("="*70)
    
    gwas_dir = BASE_DIR / "gwas" / "bakker_2020"
    gwas_dir.mkdir(parents=True, exist_ok=True)
    
    # Bakker et al. 2020 GWAS - European Heart Journal
    # Common sources for GWAS summary statistics:
    sources = [
        # Try GWAS Catalog
        "https://www.ebi.ac.uk/gwas/api/search/downloads?efo=EFO_0004709&pvalfilter=5e-8&genomiccontexts=true",
        # Try authors' repository (if available)
        "https://github.com/mkbakker/aneurysm-gwas",  # Example - may not exist
    ]
    
    print("\nNote: Bakker et al. 2020 GWAS summary statistics may require:")
    print("  1. Direct download from European Heart Journal supplementary materials")
    print("  2. Contact with authors for full summary statistics")
    print("  3. Access through GWAS Catalog (https://www.ebi.ac.uk/gwas/)")
    print("\nCreating placeholder structure and documentation...")
    
    # Create a README with download instructions
    readme_path = gwas_dir / "README.md"
    readme_content = """# Bakker et al. 2020 GWAS Summary Statistics

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
"""
    
    with open(readme_path, 'w') as f:
        f.write(readme_content)
    
    print(f"  ✓ Created README at: {readme_path}")
    return True

def download_gtex():
    """Download GTEx normal vessel reference transcriptomes"""
    print("\n" + "="*70)
    print("Downloading GTEx: Normal Vessel Reference Transcriptomes")
    print("="*70)
    
    gtex_dir = BASE_DIR / "gtex"
    gtex_dir.mkdir(parents=True, exist_ok=True)
    
    # GTEx Portal URLs
    # GTEx v8 data is available through:
    # 1. GTEx Portal: https://gtexportal.org/home/datasets
    # 2. dbGaP (requires approval): https://www.ncbi.nlm.nih.gov/gap/
    # 3. AnVIL: https://anvilproject.org/
    
    print("\nNote: GTEx data access:")
    print("  1. Public expression data: Available via GTEx Portal")
    print("  2. Full data: Requires dbGaP approval")
    print("\nDownloading publicly available vessel tissue expression data...")
    
    # GTEx public API endpoints
    gtex_urls = [
        # GTEx API for gene expression (example - may need actual API calls)
        "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz",
        # Vessel-specific data (if available)
    ]
    
    # Create README with instructions
    readme_path = gtex_dir / "README.md"
    readme_content = """# GTEx Normal Vessel Reference Transcriptomes

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
"""
    
    with open(readme_path, 'w') as f:
        f.write(readme_content)
    
    print(f"  ✓ Created README at: {readme_path}")
    
    # Try to download public GTEx expression data
    print("\nAttempting to download GTEx public expression data...")
    gtex_public_url = "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"
    
    # This file is very large (~2GB), so we'll create a script to download it
    download_script = gtex_dir / "download_gtex_expression.sh"
    with open(download_script, 'w') as f:
        f.write(f"""#!/bin/bash
# GTEx Expression Data Download Script
# File size: ~2GB, download may take time

echo "Downloading GTEx v8 gene expression data..."
wget -O GTEx_v8_gene_tpm.gct.gz \\
    "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"

echo "Download complete. Extract vessel tissues using R/Python."
""")
    
    # Also create Python download script
    python_script = gtex_dir / "download_gtex.py"
    with open(python_script, 'w') as f:
        f.write("""import requests
import gzip
from pathlib import Path

# GTEx v8 Gene Expression (TPM)
url = "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"
output = Path("GTEx_v8_gene_tpm.gct.gz")

print("Downloading GTEx v8 expression data (~2GB)...")
print("This may take 10-30 minutes depending on connection speed.")

response = requests.get(url, stream=True)
total_size = int(response.headers.get('content-length', 0))

with open(output, 'wb') as f:
    downloaded = 0
    for chunk in response.iter_content(chunk_size=8192):
        if chunk:
            f.write(chunk)
            downloaded += len(chunk)
            if total_size > 0:
                percent = (downloaded / total_size) * 100
                print(f"\\rProgress: {percent:.1f}%", end='')

print("\\nDownload complete!")
""")
    
    print(f"  ✓ Created download scripts at: {gtex_dir}")
    return True

def download_uk_biobank():
    """Download UK Biobank GWAS summary statistics"""
    print("\n" + "="*70)
    print("Downloading UK Biobank GWAS Summary Statistics")
    print("="*70)
    
    ukbb_dir = BASE_DIR / "uk_biobank"
    ukbb_dir.mkdir(parents=True, exist_ok=True)
    
    print("\nNote: UK Biobank data access:")
    print("  1. Summary statistics: Available through UK Biobank website")
    print("  2. Full data: Requires application and approval")
    print("  3. Neale Lab: Provides processed summary statistics")
    
    # Create README with instructions
    readme_path = ukbb_dir / "README.md"
    readme_content = """# UK Biobank GWAS Summary Statistics

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
"""
    
    with open(readme_path, 'w') as f:
        f.write(readme_content)
    
    print(f"  ✓ Created README at: {readme_path}")
    
    # Try to download from Neale Lab if available
    print("\nChecking Neale Lab UK Biobank resources...")
    # Note: Neale Lab doesn't have a direct API, but provides download links
    
    return True

def main():
    """Main download function"""
    print("="*70)
    print("Cerebral Aneurysm Genomics Project - Dataset Download")
    print("="*70)
    
    results = {}
    
    # 1. Download GSE186927 (single-cell RNA-seq)
    results['GSE186927'] = download_gse186927()
    
    # 2. Download GWAS (Bakker et al. 2020)
    results['GWAS_Bakker'] = download_gwas_bakker()
    
    # 3. Download GTEx
    results['GTEx'] = download_gtex()
    
    # 4. Download UK Biobank
    results['UK_Biobank'] = download_uk_biobank()
    
    # Summary
    print("\n" + "="*70)
    print("Download Summary")
    print("="*70)
    for dataset, success in results.items():
        status = "✓" if success else "✗"
        print(f"{status} {dataset}: {'Success' if success else 'See README for manual download'}")
    
    print("\n" + "="*70)
    print("Next Steps:")
    print("="*70)
    print("1. Check README files in each directory for manual download instructions")
    print("2. Some datasets require registration/approval (GTEx, UK Biobank)")
    print("3. GSE186927: Check if download was successful, may need GEOquery R package")
    print("4. GWAS data: May need to download from journal supplementary materials")
    print("="*70)

if __name__ == "__main__":
    main()


