"""
Verify GSE66238 dataset and check if it's relevant to aneurysm project
"""

import requests
from bs4 import BeautifulSoup
import re
from pathlib import Path

def check_geo_dataset(accession):
    """Check GEO dataset information"""
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"
    
    try:
        response = requests.get(url, timeout=30)
        if response.status_code != 200:
            return None
        
        soup = BeautifulSoup(response.text, 'html.parser')
        
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
            info['summary'] = summary_text[:500] if len(summary_text) > 500 else summary_text
        else:
            info['summary'] = "Not found"
        
        # Check relevance
        text_content = (info['title'] + ' ' + info['summary']).lower()
        info['is_human'] = 'homo sapiens' in info['organism'].lower() or 'human' in info['organism'].lower()
        info['is_aneurysm'] = 'aneurysm' in text_content
        info['is_rnaseq'] = 'rna' in text_content and 'seq' in text_content
        info['is_vascular'] = any(kw in text_content for kw in ['vascular', 'artery', 'vessel', 'endothelial', 'cerebral', 'intracranial'])
        
        return info
    except Exception as e:
        return {'accession': accession, 'error': str(e)}

# Check GSE66238
print("="*80)
print("VERIFYING GSE66238 DATASET")
print("="*80)

result = check_geo_dataset('GSE66238')

if result and 'error' not in result:
    print(f"\n✅ Dataset Found: {result['accession']}")
    print(f"\nTitle: {result['title']}")
    print(f"Organism: {result['organism']}")
    print(f"Type: {result['exp_type']}")
    print(f"\nSummary:")
    print(result['summary'][:400] + "..." if len(result['summary']) > 400 else result['summary'])
    
    print("\n" + "="*80)
    print("RELEVANCE CHECK")
    print("="*80)
    
    print(f"\nHuman: {result.get('is_human', False)}")
    print(f"Aneurysm: {result.get('is_aneurysm', False)}")
    print(f"RNA-seq: {result.get('is_rnaseq', False)}")
    print(f"Vascular: {result.get('is_vascular', False)}")
    
    relevance_score = sum([
        result.get('is_human', False),
        result.get('is_aneurysm', False),
        result.get('is_rnaseq', False),
        result.get('is_vascular', False)
    ])
    
    if relevance_score >= 2:
        print(f"\n✅ RELEVANT (Score: {relevance_score}/4) - This dataset can be included!")
    elif relevance_score == 1:
        print(f"\n⚠️  PARTIALLY RELEVANT (Score: {relevance_score}/4) - May be useful")
    else:
        print(f"\n❌ NOT RELEVANT (Score: {relevance_score}/4)")
    
    print(f"\nURL: {result['url']}")
else:
    print("\n❌ Could not retrieve dataset information")

# Check the downloaded files
print("\n" + "="*80)
print("CHECKING DOWNLOADED FILES")
print("="*80)

downloads_dir = Path("C:/Users/vsanker/Downloads")
gse66238_files = list(downloads_dir.glob("GSE66238*"))

if gse66238_files:
    print(f"\n✅ Found {len(gse66238_files)} file(s) for GSE66238:")
    total_size = 0
    for f in sorted(gse66238_files):
        size_mb = f.stat().st_size / (1024 * 1024)
        total_size += f.stat().st_size
        print(f"  - {f.name} ({size_mb:.2f} MB)")
    
    total_size_mb = total_size / (1024 * 1024)
    print(f"\nTotal size: {total_size_mb:.2f} MB")
    
    # Check file types
    print("\n" + "="*80)
    print("FILE TYPE ANALYSIS")
    print("="*80)
    
    file_types = {
        'raw_counts': [f for f in gse66238_files if 'raw_counts' in f.name.lower()],
        'normalized': [f for f in gse66238_files if 'norm_counts' in f.name.lower() or 'tpm' in f.name.lower() or 'fpkm' in f.name.lower()],
        'differential': [f for f in gse66238_files if 'diff' in f.name.lower() or 'differential' in f.name.lower()],
        'isoform': [f for f in gse66238_files if 'isoform' in f.name.lower()],
    }
    
    for file_type, matches in file_types.items():
        if matches:
            print(f"\n{file_type.upper().replace('_', ' ')} files ({len(matches)}):")
            for f in matches:
                size_mb = f.stat().st_size / (1024 * 1024)
                print(f"  - {f.name} ({size_mb:.2f} MB)")
else:
    print(f"\n❌ No GSE66238 files found in Downloads")

print("\n" + "="*80)
print("RECOMMENDATION")
print("="*80)

if result and relevance_score >= 2:
    print("\n✅ RECOMMENDED: Include this dataset in the project")
    print("\nNext steps:")
    print("1. Move files to appropriate project directory")
    print("2. Create README with dataset information")
    print("3. Update project documentation")
elif result and relevance_score == 1:
    print("\n⚠️  PARTIALLY RELEVANT - Review dataset details")
    print("   May be useful depending on specific research needs")
else:
    print("\n❌ NOT RELEVANT - Dataset does not match project criteria")






