"""
Download remaining datasets that require manual steps or alternative methods
"""

import requests
from pathlib import Path
import sys
import json
import time
import gzip
import shutil

def p(msg):
    """Print with flush"""
    try:
        print(msg)
        sys.stdout.flush()
    except:
        pass

BASE_DIR = Path("data")
BASE_DIR.mkdir(exist_ok=True)

def download_file(url, output_path, desc="file", timeout=300, chunk_size=8192):
    """Download file with progress tracking"""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Skip if file exists and has reasonable size
    if output_path.exists():
        size = output_path.stat().st_size
        if size > 1000:  # At least 1KB
            p(f"  [SKIP] {desc} (exists, {size:,} bytes)")
            return True
    
    try:
        p(f"  [DL] {desc}")
        p(f"        URL: {url[:80]}...")
        
        response = requests.get(url, stream=True, timeout=timeout, allow_redirects=True)
        
        if response.status_code == 404:
            p(f"  [404] Not found")
            return False
        
        response.raise_for_status()
        
        total_size = int(response.headers.get('content-length', 0))
        downloaded = 0
        
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total_size > 0:
                        percent = (downloaded / total_size) * 100
                        if downloaded % (chunk_size * 100) == 0:  # Print every 100 chunks
                            p(f"        Progress: {percent:.1f}% ({downloaded:,}/{total_size:,} bytes)")
        
        final_size = output_path.stat().st_size
        p(f"  [OK] Downloaded {final_size:,} bytes")
        return True
        
    except requests.exceptions.Timeout:
        p(f"  [TIMEOUT] Download took too long")
        if output_path.exists():
            output_path.unlink()
        return False
    except requests.exceptions.RequestException as e:
        p(f"  [ERROR] {str(e)[:100]}")
        if output_path.exists():
            output_path.unlink()
        return False
    except Exception as e:
        p(f"  [ERROR] {type(e).__name__}: {str(e)[:100]}")
        if output_path.exists():
            output_path.unlink()
        return False

def download_gwas_bakker():
    """Download Bakker et al. 2020 GWAS summary statistics"""
    p("\n" + "="*80)
    p("[1] GWAS Summary Statistics - Bakker et al. 2020")
    p("="*80)
    
    gwas_dir = BASE_DIR / "gwas" / "bakker_2020"
    gwas_dir.mkdir(parents=True, exist_ok=True)
    
    # Try Knowledge Portal Network (from web search)
    kp4cd_url = "https://www.kp4cd.org/node/451"
    p(f"\n  [INFO] Checking Knowledge Portal Network: {kp4cd_url}")
    
    # Try to get the page and find download links
    try:
        response = requests.get(kp4cd_url, timeout=30)
        if response.status_code == 200:
            html_content = response.text
            # Save the page for reference
            (gwas_dir / "kp4cd_page.html").write_text(html_content, encoding='utf-8')
            p(f"  [OK] Saved page for reference")
            
            # Look for download links in the HTML
            import re
            # Common patterns for download links
            download_patterns = [
                r'href=["\']([^"\']*\.(?:txt|tsv|csv|gz|zip)[^"\']*)["\']',
                r'https://[^"\s]+\.(?:txt|tsv|csv|gz|zip)',
            ]
            
            found_urls = []
            for pattern in download_patterns:
                matches = re.findall(pattern, html_content, re.IGNORECASE)
                found_urls.extend(matches)
            
            if found_urls:
                p(f"  [INFO] Found {len(found_urls)} potential download links")
                for url in found_urls[:5]:  # Try first 5
                    if 'summary' in url.lower() or 'statistics' in url.lower() or 'gwas' in url.lower():
                        output_file = gwas_dir / "bakker_2020_summary_stats.txt.gz"
                        if url.startswith('/'):
                            url = "https://www.kp4cd.org" + url
                        if download_file(url, output_file, "GWAS summary statistics", timeout=300):
                            return True
    except Exception as e:
        p(f"  [WARN] Could not parse KP4CD page: {e}")
    
    # Try GWAS Catalog API
    p("\n  [INFO] Trying GWAS Catalog API...")
    gwas_catalog_urls = [
        "https://www.ebi.ac.uk/gwas/api/search/downloads?q=intracranial+aneurysm+Bakker",
        "https://www.ebi.ac.uk/gwas/api/search/downloads?efo=EFO_0004709",
    ]
    
    for url in gwas_catalog_urls:
        output_file = gwas_dir / "gwas_catalog_aneurysm.tsv"
        if download_file(url, output_file, "GWAS Catalog search results", timeout=120):
            # Check if file contains Bakker study
            try:
                content = output_file.read_text(encoding='utf-8', errors='ignore')
                if 'bakker' in content.lower() or '2020' in content:
                    p(f"  [OK] Found Bakker study in GWAS Catalog results")
                    return True
            except:
                pass
    
    # Try European Heart Journal supplementary materials (common pattern)
    p("\n  [INFO] Note: Summary statistics may be in journal supplementary materials")
    p("        Check: https://academic.oup.com/eurheartj/article-lookup/doi/10.1093/eurheartj/ehaa520")
    
    return False

def download_gtex():
    """Download GTEx vessel tissue expression data"""
    p("\n" + "="*80)
    p("[2] GTEx Vessel Tissue Expression Data")
    p("="*80)
    
    gtex_dir = BASE_DIR / "gtex"
    gtex_dir.mkdir(parents=True, exist_ok=True)
    
    # Download sample attributes first (small file)
    p("\n  [INFO] Downloading GTEx metadata files...")
    
    metadata_urls = [
        ("https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
         gtex_dir / "GTEx_v8_sample_attributes.txt"),
        ("https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
         gtex_dir / "GTEx_v8_subject_phenotypes.txt"),
    ]
    
    for url, output in metadata_urls:
        download_file(url, output, f"GTEx: {output.name}", timeout=120)
    
    # Try GTEx Portal API for vessel tissues
    p("\n  [INFO] Attempting GTEx Portal API for vessel tissue expression...")
    
    vessel_tissues = ["Artery_Aorta", "Artery_Coronary", "Artery_Tibial"]
    
    for tissue in vessel_tissues:
        try:
            api_url = "https://gtexportal.org/api/v2/expression/medianGeneExpression"
            params = {
                "tissueSiteDetailId": tissue,
                "datasetId": "gtex_v8"
            }
            p(f"\n  [DL] {tissue} expression via API...")
            response = requests.get(api_url, params=params, timeout=120)
            
            if response.status_code == 200:
                data = response.json()
                output_file = gtex_dir / f"GTEx_v8_{tissue}_expression.json"
                with open(output_file, 'w', encoding='utf-8') as f:
                    json.dump(data, f, indent=2)
                
                gene_count = len(data.get('data', []))
                p(f"  [OK] {tissue}: {gene_count} genes")
            else:
                p(f"  [WARN] API returned status {response.status_code}")
        except Exception as e:
            p(f"  [ERROR] {tissue}: {e}")
    
    # Full expression matrix (large file - 2GB)
    p("\n  [INFO] Full GTEx expression matrix available (2GB):")
    full_matrix_url = "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"
    p(f"        URL: {full_matrix_url}")
    p("        This is a large file. Download manually if needed.")
    
    # Ask user if they want to download the full matrix
    response = input("\n  Download full GTEx expression matrix (2GB)? [y/N]: ").strip().lower()
    if response == 'y':
        output_file = gtex_dir / "GTEx_v8_full_expression_matrix.gct.gz"
        download_file(full_matrix_url, output_file, "GTEx full expression matrix", timeout=3600, chunk_size=65536)
    
    return True

def download_gse186927():
    """Attempt to download GSE186927 raw counts using alternative methods"""
    p("\n" + "="*80)
    p("[3] GSE186927 - Single-cell RNA-seq Raw Counts")
    p("="*80)
    
    gse_dir = BASE_DIR / "sc_rnaseq" / "GSE186927"
    gse_dir.mkdir(parents=True, exist_ok=True)
    
    # Try alternative GEO download URLs
    p("\n  [INFO] Trying alternative GEO download methods...")
    
    # Method 1: Try GEO soft format
    geo_soft_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE186927&format=soft&file=GSE186927_family.soft.gz"
    output_file = gse_dir / "GSE186927_family.soft.gz"
    if download_file(geo_soft_url, output_file, "GEO SOFT format", timeout=300):
        p("  [OK] Downloaded SOFT format (contains metadata and some expression data)")
    
    # Method 2: Try GEO supplementary files
    p("\n  [INFO] Checking GEO supplementary files...")
    geo_supp_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186927"
    try:
        response = requests.get(geo_supp_url, timeout=30)
        if response.status_code == 200:
            html = response.text
            (gse_dir / "GSE186927_geo_page_full.html").write_text(html, encoding='utf-8')
            
            # Look for supplementary file links
            import re
            supp_pattern = r'href=["\']([^"\']*suppl[^"\']*\.(?:txt|tsv|csv|h5|h5ad|gz|zip)[^"\']*)["\']'
            supp_links = re.findall(supp_pattern, html, re.IGNORECASE)
            
            if supp_links:
                p(f"  [INFO] Found {len(supp_links)} potential supplementary file links")
                for link in supp_links[:3]:  # Try first 3
                    if not link.startswith('http'):
                        link = "https://www.ncbi.nlm.nih.gov" + link
                    filename = link.split('/')[-1].split('?')[0]
                    output_file = gse_dir / filename
                    if download_file(link, output_file, f"Supplementary: {filename}", timeout=300):
                        p(f"  [OK] Downloaded supplementary file")
    except Exception as e:
        p(f"  [WARN] Could not check supplementary files: {e}")
    
    # Method 3: Create R script for GEOquery
    p("\n  [INFO] Creating R script for GEOquery download...")
    r_script = """# R script to download GSE186927 using GEOquery
# Run this in R: source("download_gse186927.R")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("GEOquery", quietly = TRUE))
    BiocManager::install("GEOquery")

library(GEOquery)

# Download GSE186927
cat("Downloading GSE186927...\\n")
gse <- getGEO("GSE186927", GSEMatrix = TRUE, getGPL = FALSE)

# Extract expression data
if (length(gse) > 0) {
    expr_data <- exprs(gse[[1]])
    pdata <- pData(gse[[1]])
    
    # Save expression matrix
    write.csv(expr_data, "data/sc_rnaseq/GSE186927/GSE186927_expression_matrix.csv", quote = FALSE)
    write.csv(pdata, "data/sc_rnaseq/GSE186927/GSE186927_phenotype_data.csv", quote = FALSE)
    
    cat("Expression data saved to GSE186927_expression_matrix.csv\\n")
    cat("Phenotype data saved to GSE186927_phenotype_data.csv\\n")
    cat("Dimensions:", dim(expr_data), "\\n")
} else {
    cat("Error: Could not download GSE186927\\n")
}
"""
    
    r_script_path = gse_dir / "download_gse186927.R"
    r_script_path.write_text(r_script, encoding='utf-8')
    p(f"  [OK] Created R script: {r_script_path}")
    p("        Run this script in R to download raw expression data")
    
    return True

def main():
    """Main download function"""
    p("="*80)
    p("DOWNLOADING REMAINING DATASETS")
    p("="*80)
    p(f"Output directory: {BASE_DIR.absolute()}")
    
    results = {
        'gwas': download_gwas_bakker(),
        'gtex': download_gtex(),
        'gse186927': download_gse186927(),
    }
    
    # Summary
    p("\n" + "="*80)
    p("DOWNLOAD SUMMARY")
    p("="*80)
    
    for dataset, success in results.items():
        status = "✓" if success else "⚠"
        p(f"{status} {dataset.upper()}: {'Downloaded' if success else 'Partial/Manual steps needed'}")
    
    p("\n" + "="*80)
    p("NEXT STEPS:")
    p("="*80)
    p("1. GWAS: Check data/gwas/bakker_2020/ for downloaded files")
    p("2. GTEx: Check data/gtex/ for expression data")
    p("3. GSE186927: Run download_gse186927.R in R, or check supplementary files")
    p("="*80)

if __name__ == "__main__":
    main()







