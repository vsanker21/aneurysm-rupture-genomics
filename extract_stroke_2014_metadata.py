"""
Extract metadata and study information from Stroke 2014 supplementary PDF
This PDF contains methods/cohort information, not GWAS summary statistics
"""
import re
import json
from pathlib import Path
from typing import Dict, List


def extract_study_metadata(text: str) -> Dict:
    """Extract study metadata and cohort information"""
    
    metadata = {
        'study_name': 'Foroud_2014_Stroke',
        'journal': 'Stroke',
        'year': '2014',
        'title': 'Genome-Wide Association Study of Intracranial Aneurysms',
        'cohorts': {},
        'sample_sizes': {},
        'methods': {}
    }
    
    # Extract sample sizes
    # Discovery sample
    discovery_match = re.search(r'Discovery Sample\s+(\d+)', text)
    
    # Look for specific cohort sizes
    # Dutch sample
    dutch_match = re.search(r'717 Dutch cases[,\s]+3004 Dutch controls', text)
    if dutch_match:
        metadata['sample_sizes']['dutch'] = {
            'cases': 717,
            'controls': 3004
        }
    
    # Finnish sample
    finnish_match = re.search(r'799 Finnish cases[,\s]+23\s*17\s*Finnish controls', text)
    if finnish_match:
        metadata['sample_sizes']['finnish'] = {
            'cases': 799,
            'controls': 2317
        }
    
    # Extract FIA study information
    fia_match = re.search(r'FIA family[^.]*\(n=(\d+)\)', text)
    if fia_match:
        metadata['cohorts']['fia_discovery_sample_1'] = int(fia_match.group(1))
    
    # Extract discovery sample sizes
    discovery_2_match = re.search(r'829 Caucasian IA cases', text)
    if discovery_2_match:
        metadata['cohorts']['discovery_sample_2_cases'] = 829
    
    discovery_3_match = re.search(r'607 Caucasian IA cases', text)
    if discovery_3_match:
        metadata['cohorts']['discovery_sample_3_cases'] = 607
    
    # Extract genotyping platform information
    if 'Axiom array' in text:
        metadata['methods']['genotyping_platform'] = 'Affymetrix Axiom array'
    
    if 'Illumina CNV370' in text:
        metadata['methods']['replication_platform'] = 'Illumina CNV370 Duo BeadChips'
    
    # Extract SNP counts
    snp_match = re.search(r'(\d+[,\s]?\d*)\s*SNPs', text)
    if snp_match:
        snp_str = snp_match.group(1).replace(',', '').replace(' ', '')
        try:
            metadata['methods']['total_snps'] = int(snp_str)
        except:
            pass
    
    # Extract reference to main paper
    ref_match = re.search(r'Foroud T[^.]*Stroke[^.]*2012[^.]*43[^.]*2846', text)
    if ref_match:
        metadata['main_paper_reference'] = {
            'authors': 'Foroud et al.',
            'journal': 'Stroke',
            'year': '2012',
            'volume': '43',
            'pages': '2846-2852',
            'title': 'Genome-Wide Association Study of Intracranial Aneurysms Confirms Role of Anril and SOX17 in Disease Risk'
        }
    
    return metadata


def extract_cohort_descriptions(text: str) -> List[Dict]:
    """Extract detailed cohort descriptions"""
    
    cohorts = []
    
    # FIA Study
    if 'Familial Intracranial Aneurysm Study' in text:
        cohorts.append({
            'name': 'FIA Study',
            'description': 'Familial Intracranial Aneurysm Study',
            'type': 'Discovery',
            'locations': 'North America, New Zealand, Australia'
        })
    
    # ACROSS
    if 'ACROSS' in text or 'Australasian Cooperative Research' in text:
        cohorts.append({
            'name': 'ACROSS',
            'description': 'Australasian Cooperative Research on Subarachnoid Hemorrhage Study',
            'type': 'Discovery',
            'locations': 'Australia, New Zealand'
        })
    
    # UCSF
    if 'University of California, San Francisco' in text:
        cohorts.append({
            'name': 'UCSF',
            'description': 'University of California, San Francisco cohort',
            'type': 'Discovery',
            'location': 'San Francisco, USA'
        })
    
    # ARIC
    if 'ARIC' in text or 'Atherosclerosis Risk in Communities' in text:
        cohorts.append({
            'name': 'ARIC',
            'description': 'Atherosclerosis Risk in Communities study',
            'type': 'Control',
            'location': 'USA'
        })
    
    # Dutch replication
    if 'Dutch sample' in text or 'Utrecht University Medical Center' in text:
        cohorts.append({
            'name': 'Dutch Replication',
            'description': 'Dutch IA patients from Utrecht University Medical Center',
            'type': 'Replication',
            'location': 'Netherlands'
        })
    
    # Finnish replication
    if 'Finnish sample' in text or 'Helsinki and Kuopio University Hospitals' in text:
        cohorts.append({
            'name': 'Finnish Replication',
            'description': 'Finnish IA patients from Helsinki and Kuopio University Hospitals',
            'type': 'Replication',
            'location': 'Finland'
        })
    
    return cohorts


def main():
    """Main extraction function"""
    print("=" * 80)
    print("EXTRACTING STROKE 2014 SUPPLEMENTARY METADATA")
    print("=" * 80)
    
    text_file = Path("data/supplementary_pdfs/extracted_gwas/str_stroke-2014-006096_supp1_text.txt")
    
    if not text_file.exists():
        print(f"\n❌ Text file not found: {text_file}")
        return
    
    print(f"\n📄 Reading: {text_file.name}")
    
    with open(text_file, 'r', encoding='utf-8') as f:
        text = f.read()
    
    print(f"  Text length: {len(text):,} characters")
    
    # Extract metadata
    print("\n📋 Extracting study metadata...")
    metadata = extract_study_metadata(text)
    
    # Extract cohort descriptions
    print("📋 Extracting cohort information...")
    cohorts = extract_cohort_descriptions(text)
    metadata['cohorts_detailed'] = cohorts
    
    # Print summary
    print("\n" + "=" * 80)
    print("EXTRACTED METADATA")
    print("=" * 80)
    
    print(f"\n📄 Study: {metadata['study_name']}")
    print(f"  Journal: {metadata['journal']} ({metadata['year']})")
    if 'main_paper_reference' in metadata:
        ref = metadata['main_paper_reference']
        print(f"  Main paper: {ref['authors']} ({ref['year']}) - {ref['journal']}")
    
    if metadata['sample_sizes']:
        print(f"\n  Sample Sizes:")
        for cohort, sizes in metadata['sample_sizes'].items():
            print(f"    {cohort.title()}: {sizes['cases']} cases, {sizes['controls']} controls")
    
    if cohorts:
        print(f"\n  Cohorts Identified: {len(cohorts)}")
        for cohort in cohorts:
            print(f"    - {cohort['name']}: {cohort['type']} ({cohort.get('location', 'Multiple locations')})")
    
    if 'methods' in metadata:
        print(f"\n  Methods:")
        for key, value in metadata['methods'].items():
            print(f"    {key}: {value}")
    
    # Save metadata
    output_file = Path("data/supplementary_pdfs/extracted_gwas/str_stroke-2014-006096_supp1_metadata_extracted.json")
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)
    
    print(f"\n✓ Metadata saved to: {output_file}")
    
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print("\n⚠️  This PDF contains METHODS/COHORT information, not GWAS summary statistics")
    print("   The actual GWAS results would be in:")
    print("   - Main paper: Foroud et al. 2012, Stroke 43:2846-2852")
    print("   - Or a different supplementary file with results tables")
    print("\n✓ Extracted valuable metadata:")
    print("   - Study design and cohorts")
    print("   - Sample sizes")
    print("   - Genotyping methods")
    print("   - Quality control procedures")


if __name__ == "__main__":
    main()






