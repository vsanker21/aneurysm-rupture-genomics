"""
Integrate all new datasets: ZIP files, TAR, PDFs
Extract metadata and prepare for analysis
"""
import pandas as pd
import numpy as np
from pathlib import Path
import zipfile
import json
import os

def check_gse26969_samples():
    """Check GSE26969 extracted samples"""
    extract_dir = Path('results/gse26969_extracted')
    if extract_dir.exists():
        cel_files = list(extract_dir.glob('*.CEL.gz'))
        sample_ids = [f.stem.replace('.CEL', '') for f in cel_files]
        return sample_ids
    return []

def create_gse26969_metadata():
    """
    Create GSE26969 metadata
    Based on web search: GSE26969 contains IA samples and STA controls
    We'll need to infer or look up which samples are which
    """
    print("=" * 80)
    print("CREATING GSE26969 METADATA")
    print("=" * 80)
    
    sample_ids = check_gse26969_samples()
    
    if not sample_ids:
        print("[ERROR] No GSE26969 samples found")
        return None
    
    print(f"\nFound {len(sample_ids)} samples: {sample_ids}")
    
    # Based on typical GEO patterns and the dataset description:
    # GSE26969 has IA (intracranial aneurysm) and STA (superficial temporal artery) samples
    # We'll create a template that can be updated with actual annotations
    # For now, we'll use a reasonable split (assuming first 3 are one group, last 3 are another)
    # This should be verified from GEO or publication
    
    metadata = pd.DataFrame({
        'Sample': sample_ids,
        'Group': ['IA', 'IA', 'IA', 'STA', 'STA', 'STA'],  # Tentative - needs verification
        'Type': ['Aneurysm', 'Aneurysm', 'Aneurysm', 'Control', 'Control', 'Control'],
        'Tissue': ['Intracranial Aneurysm', 'Intracranial Aneurysm', 'Intracranial Aneurysm',
                   'Superficial Temporal Artery', 'Superficial Temporal Artery', 'Superficial Temporal Artery'],
        'Rupture_Status': ['Unknown', 'Unknown', 'Unknown', 'N/A', 'N/A', 'N/A'],
        'Note': ['Verify from GEO', 'Verify from GEO', 'Verify from GEO', 
                'Control tissue', 'Control tissue', 'Control tissue']
    })
    
    output_file = Path('results/gse26969_metadata.tsv')
    metadata.to_csv(output_file, sep='\t', index=False)
    
    print(f"\n[OK] Created metadata: {output_file}")
    print("\n[NOTE] Group assignments are tentative and should be verified from:")
    print("  - GEO series matrix file")
    print("  - Publication metadata")
    print("  - Sample descriptions")
    
    print("\nCurrent assignments:")
    print(metadata[['Sample', 'Group', 'Type']].to_string())
    
    return metadata

def extract_supplementary_data_info():
    """Extract and organize information from supplementary data ZIP files"""
    print("\n" + "=" * 80)
    print("EXTRACTING SUPPLEMENTARY DATA INFORMATION")
    print("=" * 80)
    
    # Check cvae186 supplementary data (already extracted)
    supp_dir = Path('results/supplementary_data_check')
    if supp_dir.exists():
        print(f"\n[OK] Supplementary data already extracted: {supp_dir}")
        
        # Check Excel file
        excel_file = supp_dir / 'Supplemental_Tables.xlsx'
        if excel_file.exists():
            print(f"\n[OK] Found supplemental tables: {excel_file}")
            try:
                xl = pd.ExcelFile(excel_file)
                print(f"  Sheets: {len(xl.sheet_names)}")
                print(f"  Sheet names: {xl.sheet_names[:10]}")
                xl.close()
            except Exception as e:
                print(f"  [WARNING] Could not read Excel file: {e}")
    
    # Check jpm supplementary
    jpm_dir = Path('results/jpm_supplementary')
    if jpm_dir.exists():
        print(f"\n[OK] JPM supplementary data extracted: {jpm_dir}")
        files = list(jpm_dir.glob('*'))
        print(f"  Files: {[f.name for f in files]}")

def extract_pdf_information():
    """Extract and organize information from PDFs"""
    print("\n" + "=" * 80)
    print("ORGANIZING PDF INFORMATION")
    print("=" * 80)
    
    # Check extracted PDF data
    pdf_data_file = Path('results/extracted_pdf_data.json')
    if pdf_data_file.exists():
        print(f"\n[OK] PDF data already extracted: {pdf_data_file}")
        try:
            with open(pdf_data_file, 'r', encoding='utf-8') as f:
                pdf_data = json.load(f)
            print(f"  Title: {pdf_data.get('title', 'N/A')[:80]}")
            print(f"  Authors: {', '.join(pdf_data.get('authors', []))}")
        except Exception as e:
            print(f"  [WARNING] Could not read PDF data: {e}")
    
    # Check summary document
    summary_file = Path('results/EXTRACTED_DATA_SUMMARY.md')
    if summary_file.exists():
        print(f"\n[OK] Summary document exists: {summary_file}")

def create_dataset_inventory():
    """Create comprehensive inventory of all datasets"""
    print("\n" + "=" * 80)
    print("CREATING DATASET INVENTORY")
    print("=" * 80)
    
    inventory = {
        'datasets': {
            'GSE122897': {
                'type': 'RNA-seq',
                'status': 'Processed',
                'location': 'data/sc_rnaseq/GSE122897/',
                'description': 'Ruptured vs Unruptured aneurysm RNA-seq',
                'samples': 'Multiple (check metadata)'
            },
            'GSE26969': {
                'type': 'Microarray (Affymetrix)',
                'status': 'Extracted, needs normalization',
                'location': 'results/gse26969_extracted/',
                'description': 'IA samples vs STA controls',
                'samples': check_gse26969_samples()
            }
        },
        'supplementary_data': {
            'cvae186': {
                'type': 'Excel/Word supplementary tables',
                'status': 'Extracted',
                'location': 'results/supplementary_data_check/',
                'description': 'Clinical characteristics and sequencing samples'
            },
            'jpm_13_00266': {
                'type': 'PDF supplementary',
                'status': 'Extracted',
                'location': 'results/jpm_supplementary/',
                'description': 'Supplementary PDF from JPM paper'
            }
        },
        'pdfs': {
            'maas_2025': {
                'file': '2025.11.06.25339724.full.pdf',
                'status': 'Extracted',
                'location': 'results/extracted_pdf_data.json',
                'description': 'Single-cell atlas of intracranial arteries'
            },
            'foroud_2014': {
                'file': 'foroud-et-al-2014-genome-wide-association-study...pdf',
                'status': 'Available in Downloads',
                'description': 'GWAS study of intracranial aneurysms'
            }
        }
    }
    
    inventory_file = Path('results/dataset_inventory.json')
    with open(inventory_file, 'w') as f:
        json.dump(inventory, f, indent=2)
    
    print(f"\n[OK] Created dataset inventory: {inventory_file}")
    
    return inventory

def main():
    """Main integration function"""
    print("=" * 80)
    print("INTEGRATING ALL NEW DATASETS")
    print("=" * 80)
    
    # Step 1: Create GSE26969 metadata
    metadata = create_gse26969_metadata()
    
    # Step 2: Extract supplementary data info
    extract_supplementary_data_info()
    
    # Step 3: Organize PDF information
    extract_pdf_information()
    
    # Step 4: Create dataset inventory
    inventory = create_dataset_inventory()
    
    print("\n" + "=" * 80)
    print("INTEGRATION COMPLETE")
    print("=" * 80)
    print("\nNext steps:")
    print("1. Verify GSE26969 metadata (update Group column if needed)")
    print("2. Run preprocessing: python preprocess_gse26969.py")
    print("3. Run integrated analysis: python analyze_integrated_datasets.py")
    print("4. Review results in results/ directory")

if __name__ == "__main__":
    main()

