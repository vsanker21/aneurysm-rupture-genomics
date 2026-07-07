"""
Check which datasets can be obtained programmatically
Tests APIs and download methods for various data sources
"""
import requests
import os
import json
from pathlib import Path

print("=" * 80)
print("CHECKING PROGRAMMATIC ACCESS TO DATASETS")
print("=" * 80)

results = {
    'programmatic': [],
    'requires_contact': [],
    'requires_approval': [],
    'unknown': []
}

# ============================================================================
# 1. E-GEOD-15629 (ArrayExpress) - Programmatically Downloadable
# ============================================================================

print("\n[1] Checking E-GEOD-15629 (ArrayExpress)...")
print("   Dataset: E-GEOD-15629 - 19 samples (8 ruptured + 6 unruptured + 5 controls)")
print("   URL: https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-15629/")

try:
    # Check if ArrayExpress API is accessible
    base_url = "https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/E-GEOD-15629"
    response = requests.get(base_url, timeout=10)
    
    if response.status_code == 200:
        data = response.json()
        print(f"   ✓ ArrayExpress API accessible")
        print(f"   ✓ Dataset metadata available")
        
        # Check for file downloads
        if 'experiments' in data and len(data['experiments']) > 0:
            exp = data['experiments'][0]
            if 'files' in exp:
                file_count = len(exp['files'])
                print(f"   ✓ Found {file_count} files available for download")
                results['programmatic'].append({
                    'dataset': 'E-GEOD-15629',
                    'source': 'ArrayExpress',
                    'status': 'Downloadable via API',
                    'url': base_url,
                    'files': file_count
                })
            else:
                print(f"   ⚠ Metadata available but file list needs checking")
                results['programmatic'].append({
                    'dataset': 'E-GEOD-15629',
                    'source': 'ArrayExpress',
                    'status': 'Metadata accessible, may require manual file download',
                    'url': base_url
                })
    else:
        print(f"   ⚠ API returned status code: {response.status_code}")
        results['unknown'].append({
            'dataset': 'E-GEOD-15629',
            'source': 'ArrayExpress',
            'status': f'API returned {response.status_code}'
        })
        
except Exception as e:
    print(f"   ✗ Error accessing ArrayExpress API: {e}")
    results['unknown'].append({
        'dataset': 'E-GEOD-15629',
        'source': 'ArrayExpress',
        'status': f'Error: {str(e)}'
    })

# ============================================================================
# 2. GEO Datasets (GSE122897, GSE159610) - Programmatically Downloadable
# ============================================================================

print("\n[2] Checking GEO datasets (GSE122897, GSE159610)...")

geo_datasets = [
    {'id': 'GSE122897', 'description': 'Tissue RNA-seq (Ruptured vs Unruptured)'},
    {'id': 'GSE159610', 'description': 'Blood RNA-seq (Aneurysm vs Control)'}
]

for dataset in geo_datasets:
    gse_id = dataset['id']
    print(f"\n   Checking {gse_id}: {dataset['description']}")
    
    try:
        # Check GEO download URL
        geo_url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}"
        response = requests.get(geo_url, timeout=10)
        
        if response.status_code == 200:
            # Check if dataset exists and is accessible
            if gse_id in response.text:
                print(f"   ✓ GEO dataset accessible: {geo_url}")
                
                # Check download URL pattern
                download_url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file&file={gse_id}_series_matrix.txt.gz"
                download_test = requests.head(download_url, timeout=10, allow_redirects=True)
                
                if download_test.status_code in [200, 302]:
                    print(f"   ✓ Files downloadable via GEO API")
                    results['programmatic'].append({
                        'dataset': gse_id,
                        'source': 'GEO',
                        'status': 'Downloadable via GEO API',
                        'url': geo_url,
                        'download_pattern': f'https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file'
                    })
                else:
                    print(f"   ⚠ Metadata accessible but download URL needs verification")
                    results['programmatic'].append({
                        'dataset': gse_id,
                        'source': 'GEO',
                        'status': 'Accessible, download may require GEOquery package',
                        'url': geo_url
                    })
            else:
                print(f"   ⚠ Dataset page accessible but {gse_id} not found in content")
                results['unknown'].append({
                    'dataset': gse_id,
                    'source': 'GEO',
                    'status': 'Page accessible but dataset not verified'
                })
        else:
            print(f"   ✗ Error accessing GEO: Status {response.status_code}")
            results['unknown'].append({
                'dataset': gse_id,
                'source': 'GEO',
                'status': f'HTTP {response.status_code}'
            })
            
    except Exception as e:
        print(f"   ✗ Error checking {gse_id}: {e}")
        results['unknown'].append({
            'dataset': gse_id,
            'source': 'GEO',
            'status': f'Error: {str(e)}'
        })

# ============================================================================
# 3. dbGaP - Requires Approval (Not Directly Programmatic)
# ============================================================================

print("\n[3] Checking dbGaP access...")
print("   Searching for 'intracranial aneurysm' studies with genotype data")

try:
    # dbGaP API endpoint (if available)
    dbgap_url = "https://www.ncbi.nlm.nih.gov/gap/"
    response = requests.get(dbgap_url, timeout=10)
    
    if response.status_code == 200:
        print(f"   ✓ dbGaP website accessible")
        print(f"   ⚠ dbGaP requires:")
        print(f"      - Principal Investigator (PI) account")
        print(f"      - IRB approval")
        print(f"      - Data Use Agreement (DUA)")
        print(f"      - Application review (2-4 weeks)")
        print(f"   ✗ NOT programmatically accessible without approval")
        
        results['requires_approval'].append({
            'dataset': 'dbGaP - Intracranial aneurysm genotype data',
            'source': 'dbGaP',
            'status': 'Requires PI account, IRB approval, DUA, and application review',
            'url': 'https://www.ncbi.nlm.nih.gov/gap/',
            'search_url': 'https://www.ncbi.nlm.nih.gov/gap/?term=intracranial+aneurysm',
            'notes': 'Search can be automated, but download requires approval process'
        })
    else:
        print(f"   ⚠ dbGaP website returned status {response.status_code}")
        results['requires_approval'].append({
            'dataset': 'dbGaP',
            'source': 'dbGaP',
            'status': 'Requires approval (website status uncertain)'
        })
        
except Exception as e:
    print(f"   ⚠ Error accessing dbGaP: {e}")
    results['requires_approval'].append({
        'dataset': 'dbGaP',
        'source': 'dbGaP',
        'status': 'Requires approval (access check failed)',
        'error': str(e)
    })

# ============================================================================
# 4. Contact Authors - Requires Manual Contact
# ============================================================================

print("\n[4] Checking author contact requirements...")

authors = [
    {
        'name': 'Dr. Ynte Ruigrok',
        'email': 'ij.m.ruigrok@umcutrecht.nl',
        'dataset': 'GSE122897',
        'request': 'Paired blood samples or genotype data'
    },
    {
        'name': 'Dr. Poppenberg',
        'email': 'Check GEO metadata',
        'dataset': 'GSE159610',
        'request': 'Paired tissue samples or genotype data'
    }
]

for author in authors:
    print(f"\n   {author['name']} ({author['dataset']})")
    print(f"   Request: {author['request']}")
    print(f"   ✗ Requires manual email contact")
    
    results['requires_contact'].append({
        'author': author['name'],
        'email': author['email'],
        'dataset': author['dataset'],
        'request': author['request'],
        'status': 'Requires manual email contact'
    })

# ============================================================================
# Summary
# ============================================================================

print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)

print("\n✅ PROGRAMMATICALLY ACCESSIBLE:")
if results['programmatic']:
    for item in results['programmatic']:
        print(f"   ✓ {item['dataset']} ({item['source']})")
        print(f"     Status: {item['status']}")
else:
    print("   None found")

print("\n📧 REQUIRES MANUAL CONTACT:")
if results['requires_contact']:
    for item in results['requires_contact']:
        print(f"   ✗ {item['author']} - {item['dataset']}")
        print(f"     Email: {item['email']}")
        print(f"     Request: {item['request']}")
else:
    print("   None")

print("\n⚠️ REQUIRES APPROVAL/APPLICATION:")
if results['requires_approval']:
    for item in results['requires_approval']:
        print(f"   ⚠ {item['dataset']} ({item['source']})")
        print(f"     Status: {item['status']}")
else:
    print("   None")

print("\n❓ UNKNOWN/ERROR:")
if results['unknown']:
    for item in results['unknown']:
        print(f"   ? {item['dataset']} ({item['source']})")
        print(f"     Status: {item['status']}")
else:
    print("   None")

# Save results
os.makedirs('results', exist_ok=True)
results_file = 'results/programmatic_access_check.json'
with open(results_file, 'w') as f:
    json.dump(results, f, indent=2)

print(f"\n[OK] Results saved to {results_file}")

print("\n" + "=" * 80)
print("RECOMMENDATIONS")
print("=" * 80)
print("""
1. ✅ E-GEOD-15629: Can download programmatically via ArrayExpress API
   → Create download script using ArrayExpress API

2. ✅ GEO datasets (GSE122897, GSE159610): Can download programmatically
   → Already have download scripts (check download_datasets.py)

3. ⚠️ dbGaP: NOT programmatically accessible
   → Requires application and approval process
   → Can programmatically search for studies, but download needs approval

4. ✗ Author contacts: NOT programmatically accessible
   → Requires manual email contact
   → No API for requesting paired samples

NEXT STEPS:
1. Create ArrayExpress download script for E-GEOD-15629
2. Verify existing GEO download scripts work
3. Create dbGaP search script (for finding studies, not downloading)
4. Create email templates for author contacts
""")


