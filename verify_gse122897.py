"""
Verify GSE122897 dataset (Human Intracranial Aneurysm RNA-seq)
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
            info['summary'] = summary_text[:400] if len(summary_text) > 400 else summary_text
        else:
            info['summary'] = "Not found"
        
        # Check relevance
        text_content = (info['title'] + ' ' + info['summary']).lower()
        info['is_human'] = 'homo sapiens' in info['organism'].lower() or 'human' in info['organism'].lower()
        info['is_aneurysm'] = 'aneurysm' in text_content
        info['is_rnaseq'] = 'rna' in text_content and 'seq' in text_content
        
        return info
    except Exception as e:
        return {'accession': accession, 'error': str(e)}

# Check GSE122897
print("="*80)
print("VERIFYING GSE122897 DATASET")
print("="*80)

result = check_geo_dataset('GSE122897')

if result and 'error' not in result:
    print(f"\n✅ Dataset Found: {result['accession']}")
    print(f"\nTitle: {result['title']}")
    print(f"Organism: {result['organism']}")
    print(f"Type: {result['exp_type']}")
    print(f"\nSummary:")
    print(result['summary'][:300] + "..." if len(result['summary']) > 300 else result['summary'])
    
    print("\n" + "="*80)
    print("RELEVANCE CHECK")
    print("="*80)
    
    print(f"\nHuman: {result.get('is_human', False)}")
    print(f"Aneurysm: {result.get('is_aneurysm', False)}")
    print(f"RNA-seq: {result.get('is_rnaseq', False)}")
    
    if result.get('is_human') and result.get('is_aneurysm') and result.get('is_rnaseq'):
        print("\n✅ HIGHLY RELEVANT - This dataset can be included!")
    elif result.get('is_human') and result.get('is_aneurysm'):
        print("\n✅ RELEVANT - Human aneurysm dataset")
    else:
        print("\n⚠️  May not be directly relevant")
    
    print(f"\nURL: {result['url']}")
else:
    print("\n❌ Could not retrieve dataset information")

# Check files in Downloads
print("\n" + "="*80)
print("FILES FOUND IN DOWNLOADS")
print("="*80)

downloads_dir = Path("C:/Users/vsanker/Downloads")
gse122897_files = [
    "GSE122897_readCounts_raw.txt.gz",
    "GSE122897_raw_counts_GRCh38.p13_NCBI.tsv.gz",
    "GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz",
    "GSE122897_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz",
    "Human.GRCh38.p13.annot.tsv.gz"
]

found_files = []
for filename in gse122897_files:
    source = downloads_dir / filename
    if source.exists():
        size_mb = source.stat().st_size / (1024 * 1024)
        found_files.append((filename, size_mb))
        print(f"✅ {filename} ({size_mb:.2f} MB)")
    else:
        print(f"❌ {filename} - Not found")

if found_files:
    total_size = sum(size for _, size in found_files)
    print(f"\nTotal: {len(found_files)} files, {total_size:.2f} MB")

print("\n" + "="*80)
print("COMPARISON WITH GSE158558")
print("="*80)

print("\nGSE158558 (already in project):")
print("  - Intracranial aneurysm sequencing data")
print("  - Bulk RNA-seq")
print("  - mRNA and lncRNA expression")

print("\nGSE122897 (newly downloaded):")
if result:
    print(f"  - {result['title']}")
    print(f"  - {result['exp_type']}")
    print("  - Similar file structure to GSE158558")

print("\n✅ Both datasets can be used together for:")
print("  - Validation across multiple studies")
print("  - Meta-analysis of aneurysm expression")
print("  - Increased sample size")
print("  - Cross-study validation of findings")






