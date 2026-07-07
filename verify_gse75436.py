"""
Verify GSE75436 dataset and check if it's relevant to aneurysm project
"""

import requests
from bs4 import BeautifulSoup
import re
from pathlib import Path
import tarfile

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
        info['is_vascular'] = any(kw in text_content for kw in ['vascular', 'artery', 'vessel', 'endothelial'])
        
        return info
    except Exception as e:
        return {'accession': accession, 'error': str(e)}

# Check GSE75436
print("="*80)
print("VERIFYING GSE75436 DATASET")
print("="*80)

result = check_geo_dataset('GSE75436')

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

# Check the tar file
print("\n" + "="*80)
print("CHECKING TAR FILE CONTENTS")
print("="*80)

downloads_dir = Path("C:/Users/vsanker/Downloads")
tar_file = downloads_dir / "GSE75436_RAW.tar"

if tar_file.exists():
    size_mb = tar_file.stat().st_size / (1024 * 1024)
    print(f"\nFile: {tar_file.name}")
    print(f"Size: {size_mb:.2f} MB ({tar_file.stat().st_size:,} bytes)")
    
    try:
        print("\nAttempting to read tar file contents...")
        with tarfile.open(tar_file, 'r') as tar:
            members = tar.getmembers()
            print(f"\n✅ Tar file contains {len(members)} file(s)")
            
            # Show first 20 files
            print("\nFirst 20 files in archive:")
            for i, member in enumerate(members[:20], 1):
                size_kb = member.size / 1024
                print(f"  {i}. {member.name} ({size_kb:.1f} KB)")
            
            if len(members) > 20:
                print(f"  ... and {len(members) - 20} more files")
            
            # Check for common data file patterns
            data_patterns = {
                'expression': [m for m in members if any(kw in m.name.lower() for kw in ['expr', 'count', 'tpm', 'fpkm'])],
                'metadata': [m for m in members if any(kw in m.name.lower() for kw in ['meta', 'pheno', 'sample', 'info'])],
                'raw_data': [m for m in members if any(kw in m.name.lower() for kw in ['raw', 'fastq', 'bam', 'sam'])],
            }
            
            print("\n" + "="*80)
            print("FILE TYPE ANALYSIS")
            print("="*80)
            
            for pattern, matches in data_patterns.items():
                if matches:
                    print(f"\n{pattern.upper()} files ({len(matches)}):")
                    for m in matches[:5]:
                        print(f"  - {m.name}")
                    if len(matches) > 5:
                        print(f"  ... and {len(matches) - 5} more")
            
    except Exception as e:
        print(f"\n⚠️  Could not read tar file: {e}")
        print("   File may be corrupted or in unexpected format")
else:
    print(f"\n❌ File not found: {tar_file}")

print("\n" + "="*80)
print("RECOMMENDATION")
print("="*80)

if result and relevance_score >= 2:
    print("\n✅ RECOMMENDED: Include this dataset in the project")
    print("\nNext steps:")
    print("1. Extract tar file to examine contents")
    print("2. Move to appropriate project directory")
    print("3. Update documentation")
elif result and relevance_score == 1:
    print("\n⚠️  PARTIALLY RELEVANT - Review dataset details")
    print("   May be useful depending on specific research needs")
else:
    print("\n❌ NOT RELEVANT - Dataset does not match project criteria")






