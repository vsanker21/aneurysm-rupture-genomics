"""
Verify potential intracranial aneurysm single-cell RNA-seq datasets
"""

import requests
from bs4 import BeautifulSoup
import re

def check_geo_dataset(accession):
    """Check a GEO dataset and extract key information"""
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"
    
    try:
        response = requests.get(url, timeout=30)
        if response.status_code != 200:
            return None
        
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Extract title
        title_elem = soup.find('td', string=re.compile('Title', re.I))
        title = title_elem.find_next_sibling('td').text.strip() if title_elem else "Not found"
        
        # Extract organism
        org_elem = soup.find('td', string=re.compile('Organism', re.I))
        organism = org_elem.find_next_sibling('td').text.strip() if org_elem else "Not found"
        
        # Extract experiment type
        exp_elem = soup.find('td', string=re.compile('Experiment type', re.I))
        exp_type = exp_elem.find_next_sibling('td').text.strip() if exp_elem else "Not found"
        
        # Extract summary
        summary_elem = soup.find('td', string=re.compile('Summary', re.I))
        summary = summary_elem.find_next_sibling('td').text.strip()[:200] if summary_elem else "Not found"
        
        return {
            'accession': accession,
            'title': title,
            'organism': organism,
            'exp_type': exp_type,
            'summary': summary,
            'url': url
        }
    except Exception as e:
        return {'accession': accession, 'error': str(e)}

# Candidate datasets to check
candidates = [
    'GSE144678',  # Mentioned as "Single-cell RNA sequencing of human intracranial aneurysm"
    'GSE193533',  # Mouse cerebral aneurysm model
    'GSE158558',  # Intracranial aneurysm sequencing data
    'GSE6551',    # Expression data from intracranial arteries and aneurysms
]

print("="*80)
print("VERIFYING INTRACRANIAL ANEURYSM SINGLE-CELL DATASETS")
print("="*80)

results = []
for acc in candidates:
    print(f"\nChecking {acc}...")
    result = check_geo_dataset(acc)
    if result:
        results.append(result)
        print(f"  Title: {result.get('title', 'N/A')}")
        print(f"  Organism: {result.get('organism', 'N/A')}")
        print(f"  Type: {result.get('exp_type', 'N/A')}")
    else:
        print(f"  ❌ Could not retrieve information")

print("\n" + "="*80)
print("SUMMARY")
print("="*80)

# Filter for human single-cell datasets
human_sc_datasets = [r for r in results if 'homo sapiens' in r.get('organism', '').lower() 
                     and ('single' in r.get('exp_type', '').lower() or 'scRNA' in r.get('exp_type', ''))]

if human_sc_datasets:
    print("\n✅ POTENTIAL CORRECT DATASETS (Human Single-Cell):")
    for r in human_sc_datasets:
        print(f"\n  {r['accession']}: {r['title']}")
        print(f"    Organism: {r['organism']}")
        print(f"    Type: {r['exp_type']}")
        print(f"    URL: {r['url']}")
else:
    print("\n⚠️  No human single-cell datasets found in candidates")
    print("\nAll checked datasets:")
    for r in results:
        print(f"\n  {r['accession']}: {r['title']}")
        print(f"    Organism: {r['organism']}")
        print(f"    Type: {r['exp_type']}")






