"""Check GTEx URL availability and find working alternatives"""

import requests
from pathlib import Path

def check_url(url, desc=""):
    """Check if URL is accessible"""
    print(f"\nChecking: {desc}")
    print(f"URL: {url[:80]}...")
    
    try:
        # Try HEAD request first (faster)
        response = requests.head(url, timeout=10, allow_redirects=True)
        status = response.status_code
        
        if status == 200:
            size = response.headers.get('Content-Length', 'Unknown')
            content_type = response.headers.get('Content-Type', 'Unknown')
            print(f"  ✅ Status: {status} OK")
            print(f"  Size: {size} bytes" if size != 'Unknown' else "  Size: Unknown")
            print(f"  Type: {content_type}")
            return True
        else:
            print(f"  ❌ Status: {status}")
            return False
            
    except requests.exceptions.RequestException as e:
        print(f"  ❌ Error: {e}")
        return False

# Test the original URL
original_url = "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"
print("="*80)
print("GTEx URL Verification")
print("="*80)

original_works = check_url(original_url, "Original URL (v8, 2017)")

# Try alternative URLs
alternative_urls = [
    ("GTEx v8 alternative path", "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"),
    ("GTEx v8 without date", "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"),
    ("GTEx v9", "https://storage.googleapis.com/gtex_analysis_v9/rna_seq_data/GTEx_Analysis_v9_RNASeQCv1.1.9_gene_tpm.gct.gz"),
    ("GTEx v8 gene read counts", "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"),
]

print("\n" + "="*80)
print("Testing Alternative URLs")
print("="*80)

working_urls = []
for desc, url in alternative_urls:
    if check_url(url, desc):
        working_urls.append((desc, url))

# Check GTEx Portal API
print("\n" + "="*80)
print("Checking GTEx Portal API")
print("="*80)

api_base = "https://gtexportal.org/api/v2"
api_endpoints = [
    ("API base", f"{api_base}/"),
    ("Expression endpoint", f"{api_base}/expression/medianGeneExpression"),
]

for desc, url in api_endpoints:
    try:
        response = requests.get(url, timeout=10)
        print(f"\n{desc}: {url}")
        print(f"  Status: {response.status_code}")
        if response.status_code == 200:
            print(f"  ✅ Accessible")
        else:
            print(f"  Response: {response.text[:200]}")
    except Exception as e:
        print(f"\n{desc}: {url}")
        print(f"  ❌ Error: {e}")

# Summary
print("\n" + "="*80)
print("SUMMARY")
print("="*80)

if original_works:
    print("✅ Original URL is WORKING")
    print(f"   URL: {original_url}")
else:
    print("❌ Original URL is BROKEN")
    
if working_urls:
    print(f"\n✅ Found {len(working_urls)} working alternative URL(s):")
    for desc, url in working_urls:
        print(f"   {desc}: {url}")
else:
    print("\n❌ No working alternative URLs found")
    print("\nRecommendations:")
    print("1. Visit GTEx Portal: https://gtexportal.org/home/datasets")
    print("2. Check GTEx GitHub: https://github.com/broadinstitute/gtex-pipeline")
    print("3. Use GTEx Portal API with proper authentication")
    print("4. Access via dbGaP (requires approval): https://www.ncbi.nlm.nih.gov/gap/")






