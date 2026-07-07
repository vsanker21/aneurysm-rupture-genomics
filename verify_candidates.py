"""
Verify candidate datasets for human intracranial aneurysm single-cell RNA-seq
"""

import requests
from bs4 import BeautifulSoup
import re

def check_dataset(accession):
    """Check a GEO dataset"""
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"
    
    try:
        response = requests.get(url, timeout=30)
        if response.status_code != 200:
            return None
        
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Extract information
        info = {'accession': accession, 'url': url}
        
        # Title
        title_elem = soup.find('td', string=re.compile('Title', re.I))
        info['title'] = title_elem.find_next_sibling('td').text.strip() if title_elem else "Not found"
        
        # Organism
        org_elem = soup.find('td', string=re.compile('Organism', re.I))
        info['organism'] = org_elem.find_next_sibling('td').text.strip() if org_elem else "Not found"
        
        # Experiment type
        exp_elem = soup.find('td', string=re.compile('Experiment type', re.I))
        info['exp_type'] = exp_elem.find_next_sibling('td').text.strip() if exp_elem else "Not found"
        
        # Summary
        summary_elem = soup.find('td', string=re.compile('Summary', re.I))
        if summary_elem:
            summary_text = summary_elem.find_next_sibling('td').text.strip()
            info['summary'] = summary_text[:300] if len(summary_text) > 300 else summary_text
        else:
            info['summary'] = "Not found"
        
        # Check for single-cell keywords
        text_content = (info['title'] + ' ' + info['summary']).lower()
        info['is_single_cell'] = any(kw in text_content for kw in ['single cell', 'single-cell', 'scrna', 'scRNA'])
        info['is_human'] = 'homo sapiens' in info['organism'].lower() or 'human' in info['organism'].lower()
        info['is_aneurysm'] = 'aneurysm' in text_content
        
        return info
    except Exception as e:
        return {'accession': accession, 'error': str(e)}

# Candidates from search
candidates = [
    'GSE115304',
    'GSE107196', 
    'GSE165824',
    'GSE172230',
    'GSE287300',
    'GSE288785',
    'GSE296197',
    'GSE72180'
]

print("="*80)
print("VERIFYING CANDIDATE DATASETS")
print("="*80)

results = []
for acc in candidates:
    print(f"\nChecking {acc}...")
    result = check_dataset(acc)
    if result and 'error' not in result:
        results.append(result)
        print(f"  Title: {result['title'][:80]}...")
        print(f"  Organism: {result['organism']}")
        print(f"  Type: {result['exp_type']}")
        print(f"  Single-cell: {result.get('is_single_cell', False)}")
        print(f"  Human: {result.get('is_human', False)}")
        print(f"  Aneurysm: {result.get('is_aneurysm', False)}")
    else:
        print(f"  ❌ Error or not found")

print("\n" + "="*80)
print("POTENTIAL CORRECT DATASETS")
print("="*80)

# Filter for human, single-cell, aneurysm datasets
correct_candidates = [r for r in results 
                      if r.get('is_human') and r.get('is_single_cell') and r.get('is_aneurysm')]

if correct_candidates:
    print("\n✅ LIKELY CORRECT DATASETS:")
    for r in correct_candidates:
        print(f"\n  {r['accession']}: {r['title']}")
        print(f"    Organism: {r['organism']}")
        print(f"    Type: {r['exp_type']}")
        print(f"    URL: {r['url']}")
        print(f"    Summary: {r['summary'][:150]}...")
else:
    print("\n⚠️  No perfect matches found. Showing closest matches:")
    # Show human single-cell datasets (even if not explicitly aneurysm)
    human_sc = [r for r in results if r.get('is_human') and r.get('is_single_cell')]
    for r in human_sc[:5]:
        print(f"\n  {r['accession']}: {r['title']}")
        print(f"    Aneurysm-related: {r.get('is_aneurysm', False)}")






