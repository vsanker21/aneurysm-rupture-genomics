"""
Comprehensive Dataset Download Script
Downloads ALL available public datasets for aneurysm genomics project
"""
import os
import requests
import pandas as pd
import time
from pathlib import Path
import gzip
import shutil
from urllib.parse import urlparse
import json

class ComprehensiveDatasetDownloader:
    """Download all available public datasets for aneurysm genomics"""
    
    def __init__(self):
        self.base_dir = Path("data")
        self.base_dir.mkdir(exist_ok=True)
        self.downloaded = []
        self.failed = []
        
    def search_geo_for_aneurysm_datasets(self):
        """Search GEO for all aneurysm-related datasets"""
        print("=" * 80)
        print("SEARCHING GEO FOR ANEURYSM DATASETS")
        print("=" * 80)
        
        # Known aneurysm GEO datasets
        aneurysm_datasets = {
            # Blood transcriptomics
            'GSE159610': {
                'name': 'Blood transcriptomics (aneurysm vs control)',
                'type': 'blood_rnaseq',
                'downloaded': os.path.exists('data/blood_rnaseq/GSE159610')
            },
            
            # Tissue transcriptomics
            'GSE122897': {
                'name': 'Tissue RNA-seq (ruptured vs unruptured)',
                'type': 'tissue_rnaseq',
                'downloaded': os.path.exists('data/sc_rnaseq/GSE122897')
            },
            'GSE158558': {
                'name': 'Tissue RNA-seq (intracranial aneurysm)',
                'type': 'tissue_rnaseq',
                'downloaded': os.path.exists('data/sc_rnaseq/GSE158558')
            },
            'GSE66238': {
                'name': 'Tissue RNA-seq (intracranial aneurysm)',
                'type': 'tissue_rnaseq',
                'downloaded': os.path.exists('data/sc_rnaseq/GSE66238')
            },
            
            # Microarray datasets
            'GSE13353': {
                'name': 'Aneurysm microarray',
                'type': 'microarray',
                'downloaded': os.path.exists('data/microarray/GSE13353')
            },
            'GSE54083': {
                'name': 'Aneurysm microarray',
                'type': 'microarray',
                'downloaded': os.path.exists('data/microarray/GSE54083')
            },
            'GSE75436': {
                'name': 'Aneurysm microarray',
                'type': 'microarray',
                'downloaded': os.path.exists('data/microarray/GSE75436')
            },
            
            # Additional datasets to check
            'GSE186927': {
                'name': 'GSE186927 (verify if aneurysm-related)',
                'type': 'tissue_rnaseq',
                'downloaded': os.path.exists('data/sc_rnaseq/GSE186927')
            },
        }
        
        # Check for additional datasets
        additional_to_search = [
            'GSE54083',  # Aneurysm microarray
            'GSE66238',  # Already have
            'GSE75436',  # Aneurysm microarray
            'GSE13353',  # Aneurysm microarray
        ]
        
        print("\nKnown Aneurysm Datasets:")
        for gse_id, info in aneurysm_datasets.items():
            status = "✅ Downloaded" if info['downloaded'] else "❌ Missing"
            print(f"  {gse_id}: {info['name']} - {status}")
        
        return aneurysm_datasets
    
    def download_geo_dataset(self, gse_id, output_dir):
        """Download a GEO dataset"""
        print(f"\nDownloading {gse_id}...")
        
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Try multiple download methods
        urls = [
            f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file&file={gse_id}_series_matrix.txt.gz",
            f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_id[:6]}nnn/{gse_id}/matrix/{gse_id}_series_matrix.txt.gz",
            f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file",
        ]
        
        downloaded_files = []
        
        for url in urls:
            try:
                print(f"  Trying: {url}")
                response = requests.get(url, timeout=60, stream=True)
                
                if response.status_code == 200:
                    # Determine filename
                    if 'series_matrix' in url:
                        filename = f"{gse_id}_series_matrix.txt.gz"
                    else:
                        # Try to get filename from headers
                        content_disposition = response.headers.get('Content-Disposition', '')
                        if 'filename=' in content_disposition:
                            filename = content_disposition.split('filename=')[1].strip('"')
                        else:
                            filename = f"{gse_id}_raw.tar.gz"
                    
                    filepath = output_dir / filename
                    
                    # Download file
                    total_size = int(response.headers.get('content-length', 0))
                    print(f"  Downloading {filename} ({total_size / 1024 / 1024:.1f} MB)...")
                    
                    with open(filepath, 'wb') as f:
                        downloaded = 0
                        for chunk in response.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)
                                downloaded += len(chunk)
                                if total_size > 0:
                                    progress = (downloaded / total_size) * 100
                                    print(f"\r  Progress: {progress:.1f}%", end='', flush=True)
                    
                    print(f"\n  ✅ Downloaded: {filename}")
                    downloaded_files.append(filepath)
                    
                    # Try to download supplementary files
                    supp_url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file&file={gse_id}_RAW.tar"
                    try:
                        supp_response = requests.get(supp_url, timeout=60, stream=True)
                        if supp_response.status_code == 200:
                            supp_file = output_dir / f"{gse_id}_RAW.tar"
                            with open(supp_file, 'wb') as f:
                                for chunk in supp_response.iter_content(chunk_size=8192):
                                    if chunk:
                                        f.write(chunk)
                            print(f"  ✅ Downloaded supplementary files: {gse_id}_RAW.tar")
                            downloaded_files.append(supp_file)
                    except:
                        pass
                    
                    break  # Success, move to next
                    
            except Exception as e:
                print(f"    Error: {e}")
                continue
        
        if downloaded_files:
            return downloaded_files
        else:
            print(f"  ❌ Failed to download {gse_id}")
            return None
    
    def download_additional_gwas_datasets(self):
        """Download additional GWAS datasets"""
        print("\n" + "=" * 80)
        print("DOWNLOADING ADDITIONAL GWAS DATASETS")
        print("=" * 80)
        
        gwas_dir = self.base_dir / "gwas"
        gwas_dir.mkdir(exist_ok=True)
        
        # Check GWAS Catalog for aneurysm GWAS
        print("\nChecking GWAS Catalog for aneurysm studies...")
        
        # Try to download from GWAS Catalog
        gwas_catalog_url = "https://www.ebi.ac.uk/gwas/api/search/downloads/alternative"
        
        try:
            # Search for intracranial aneurysm
            search_term = "intracranial aneurysm"
            # Note: GWAS Catalog API requires proper query format
            print("  Note: GWAS Catalog requires manual download or API access")
            print("  Available studies:")
            print("    - Bakker et al. 2020 (already downloaded)")
            print("    - Foroud et al. 2014 (check journal supplement)")
            print("    - Bilguvar et al. 2008 (check journal supplement)")
        except:
            pass
    
    def download_1000genomes_data(self):
        """Download 1000 Genomes reference data for PRS calculation"""
        print("\n" + "=" * 80)
        print("DOWNLOADING 1000 GENOMES REFERENCE DATA")
        print("=" * 80)
        
        kg_dir = self.base_dir / "1000genomes"
        kg_dir.mkdir(exist_ok=True)
        
        print("\n1000 Genomes Phase 3 data available at:")
        print("  https://www.internationalgenome.org/data")
        print("\nNote: Full dataset is very large (several TB)")
        print("For PRS calculation, typically only need:")
        print("  - VCF files for specific chromosomes")
        print("  - Population allele frequencies")
        
        # Download population frequency file (smaller)
        freq_url = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
        # Note: This is just an example - would need to download chromosome by chromosome
    
    def download_gtex_additional_tissues(self):
        """Download additional GTEx tissues if needed"""
        print("\n" + "=" * 80)
        print("CHECKING ADDITIONAL GTEx TISSUES")
        print("=" * 80)
        
        gtex_dir = self.base_dir / "gtex"
        gtex_dir.mkdir(exist_ok=True)
        
        # Already have vessel tissues
        vessel_tissues = ['artery_aorta', 'artery_coronary', 'artery_tibial']
        
        print("\nCurrent GTEx vessel tissues:")
        for tissue in vessel_tissues:
            tpm_file = gtex_dir / f"gene_tpm_v10_artery_{tissue.split('_')[1]}.gct.gz"
            if tpm_file.exists():
                print(f"  ✅ {tissue}")
            else:
                print(f"  ❌ {tissue} - Missing")
    
    def download_uk_biobank_summary_stats(self):
        """Attempt to download UK Biobank summary statistics"""
        print("\n" + "=" * 80)
        print("UK BIOBANK SUMMARY STATISTICS")
        print("=" * 80)
        
        ukbb_dir = self.base_dir / "uk_biobank"
        ukbb_dir.mkdir(exist_ok=True)
        
        print("\nUK Biobank summary statistics available from:")
        print("  1. Neale Lab: http://www.nealelab.is/uk-biobank/")
        print("  2. Pan-UK Biobank: https://pan.ukbb.broadinstitute.org/")
        print("\nNote: Requires checking for 'intracranial aneurysm' phenotype")
        print("Full access requires application and approval")
        
        # Try to download from Neale Lab if available
        # Note: This requires knowing the exact phenotype code
        print("\nTo download UK Biobank data:")
        print("  1. Visit Neale Lab website")
        print("  2. Search for aneurysm-related phenotypes")
        print("  3. Download summary statistics files")
    
    def download_public_genotype_reference(self):
        """Download public genotype reference data"""
        print("\n" + "=" * 80)
        print("PUBLIC GENOTYPE REFERENCE DATA")
        print("=" * 80)
        
        print("\nAvailable public genotype reference datasets:")
        print("  1. 1000 Genomes Phase 3 - Population reference")
        print("  2. HapMap - Reference panel")
        print("  3. TOPMed - Diverse population reference")
        
        print("\nNote: Patient-level genotype data typically requires:")
        print("  - dbGaP approval (for restricted access)")
        print("  - Institutional agreements")
        print("  - Cannot download without proper authorization")
        
        print("\nFor PRS calculation without patient data, can use:")
        print("  - Summary statistics (already have)")
        print("  - Population allele frequencies")
        print("  - Reference panel for imputation")
    
    def check_dbgap_datasets(self):
        """Check dbGaP for available aneurysm datasets"""
        print("\n" + "=" * 80)
        print("dbGaP DATASET AVAILABILITY")
        print("=" * 80)
        
        print("\ndbGaP datasets (require approval):")
        print("  1. Search dbGaP: https://www.ncbi.nlm.nih.gov/gap/")
        print("  2. Look for aneurysm/genetics studies")
        print("  3. Apply for data access")
        
        print("\nNote: dbGaP access typically requires:")
        print("  - Research proposal")
        print("  - Institutional approval")
        print("  - Data Use Agreement")
        print("  - Review and approval process (2-4 weeks)")
    
    def create_dataset_inventory(self):
        """Create comprehensive inventory of all datasets"""
        print("\n" + "=" * 80)
        print("CREATING DATASET INVENTORY")
        print("=" * 80)
        
        inventory = {
            'blood_rnaseq': [],
            'tissue_rnaseq': [],
            'microarray': [],
            'gwas': [],
            'gtex': [],
            'reference_data': [],
            'missing_but_available': [],
            'restricted_access': []
        }
        
        # Check what's downloaded
        data_dir = Path("data")
        
        # Blood RNA-seq
        blood_dir = data_dir / "blood_rnaseq"
        if blood_dir.exists():
            for subdir in blood_dir.iterdir():
                if subdir.is_dir():
                    inventory['blood_rnaseq'].append(subdir.name)
        
        # Tissue RNA-seq
        tissue_dir = data_dir / "sc_rnaseq"
        if tissue_dir.exists():
            for subdir in tissue_dir.iterdir():
                if subdir.is_dir():
                    inventory['tissue_rnaseq'].append(subdir.name)
        
        # Microarray
        microarray_dir = data_dir / "microarray"
        if microarray_dir.exists():
            for subdir in microarray_dir.iterdir():
                if subdir.is_dir():
                    inventory['microarray'].append(subdir.name)
        
        # GWAS
        gwas_dir = data_dir / "gwas"
        if gwas_dir.exists():
            for subdir in gwas_dir.iterdir():
                if subdir.is_dir():
                    inventory['gwas'].append(subdir.name)
        
        # GTEx
        gtex_dir = data_dir / "gtex"
        if gtex_dir.exists():
            files = list(gtex_dir.glob("*.gct.gz"))
            if files:
                inventory['gtex'] = [f.name for f in files]
        
        # Save inventory
        inventory_file = self.base_dir / "dataset_inventory.json"
        with open(inventory_file, 'w') as f:
            json.dump(inventory, f, indent=2)
        
        print("\nDataset Inventory:")
        for category, datasets in inventory.items():
            if datasets:
                print(f"\n{category.replace('_', ' ').title()}:")
                for ds in datasets:
                    print(f"  ✅ {ds}")
        
        return inventory
    
    def download_missing_datasets(self):
        """Download any missing datasets"""
        print("\n" + "=" * 80)
        print("DOWNLOADING MISSING DATASETS")
        print("=" * 80)
        
        aneurysm_datasets = self.search_geo_for_aneurysm_datasets()
        
        for gse_id, info in aneurysm_datasets.items():
            if not info['downloaded']:
                print(f"\nAttempting to download missing dataset: {gse_id}")
                
                if info['type'] == 'blood_rnaseq':
                    output_dir = self.base_dir / "blood_rnaseq" / gse_id
                elif info['type'] == 'tissue_rnaseq':
                    output_dir = self.base_dir / "sc_rnaseq" / gse_id
                elif info['type'] == 'microarray':
                    output_dir = self.base_dir / "microarray" / gse_id
                else:
                    output_dir = self.base_dir / "other" / gse_id
                
                result = self.download_geo_dataset(gse_id, output_dir)
                if result:
                    self.downloaded.append((gse_id, result))
                else:
                    self.failed.append(gse_id)
                
                time.sleep(2)  # Rate limiting
    
    def generate_download_report(self):
        """Generate comprehensive download report"""
        print("\n" + "=" * 80)
        print("DOWNLOAD REPORT")
        print("=" * 80)
        
        print(f"\nSuccessfully downloaded: {len(self.downloaded)} datasets")
        for gse_id, files in self.downloaded:
            print(f"  ✅ {gse_id}: {len(files)} file(s)")
        
        if self.failed:
            print(f"\nFailed to download: {len(self.failed)} datasets")
            for gse_id in self.failed:
                print(f"  ❌ {gse_id}")
        
        # Create inventory
        inventory = self.create_dataset_inventory()
        
        # Generate recommendations
        print("\n" + "=" * 80)
        print("RECOMMENDATIONS FOR ADDITIONAL DATA")
        print("=" * 80)
        
        print("\nFor Patient-Level PRS Calculation:")
        print("  1. Apply for dbGaP access to aneurysm studies")
        print("  2. Check if authors provide individual-level data")
        print("  3. Use population-level frequencies as alternative")
        
        print("\nFor Predictive Model Training:")
        print("  1. Contact dataset authors for clinical outcome data")
        print("  2. Look for published supplementary tables with patient data")
        print("  3. Check if metadata includes rupture outcomes")
        
        print("\nFor Exosomal Biomarker Discovery:")
        print("  1. Search GEO for exosomal RNA datasets")
        print("  2. Check for aneurysm-specific exosomal studies")
        print("  3. May need to analyze existing datasets for exosomal markers")

def main():
    """Main download pipeline"""
    print("=" * 80)
    print("COMPREHENSIVE DATASET DOWNLOAD")
    print("Downloading ALL available public datasets")
    print("=" * 80)
    
    downloader = ComprehensiveDatasetDownloader()
    
    # 1. Search and download GEO datasets
    downloader.search_geo_for_aneurysm_datasets()
    downloader.download_missing_datasets()
    
    # 2. Check additional GWAS datasets
    downloader.download_additional_gwas_datasets()
    
    # 3. Check GTEx additional tissues
    downloader.download_gtex_additional_tissues()
    
    # 4. UK Biobank info
    downloader.download_uk_biobank_summary_stats()
    
    # 5. Reference data info
    downloader.download_public_genotype_reference()
    
    # 6. dbGaP info
    downloader.check_dbgap_datasets()
    
    # 7. Generate report
    downloader.generate_download_report()
    
    print("\n" + "=" * 80)
    print("DOWNLOAD PIPELINE COMPLETE!")
    print("=" * 80)

if __name__ == "__main__":
    main()






