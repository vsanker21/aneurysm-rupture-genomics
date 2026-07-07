"""
Comprehensive Dataset Download - Downloads ALL available public datasets
Searches and downloads all missing aneurysm-related datasets
"""
import os
import requests
import pandas as pd
import time
from pathlib import Path
import json
from typing import List, Dict

class AneurysmDatasetDownloader:
    """Comprehensive downloader for all aneurysm genomics datasets"""
    
    def __init__(self):
        self.base_dir = Path("data")
        self.base_dir.mkdir(exist_ok=True)
        self.session = requests.Session()
        self.session.headers.update({'User-Agent': 'Mozilla/5.0'})
        
    def download_geo_series_matrix(self, gse_id: str, output_dir: Path):
        """Download GEO series matrix file"""
        output_dir.mkdir(parents=True, exist_ok=True)
        
        url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file&file={gse_id}_series_matrix.txt.gz"
        output_file = output_dir / f"{gse_id}_series_matrix.txt.gz"
        
        if output_file.exists():
            print(f"  ✅ Already exists: {output_file.name}")
            return output_file
        
        try:
            print(f"  Downloading series matrix for {gse_id}...")
            response = self.session.get(url, timeout=120, stream=True)
            
            if response.status_code == 200:
                with open(output_file, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                
                size_mb = output_file.stat().st_size / (1024 * 1024)
                print(f"  ✅ Downloaded: {output_file.name} ({size_mb:.1f} MB)")
                return output_file
            else:
                print(f"  ❌ Failed: HTTP {response.status_code}")
                return None
                
        except Exception as e:
            print(f"  ❌ Error: {e}")
            return None
    
    def download_geo_raw_files(self, gse_id: str, output_dir: Path):
        """Download GEO raw data files"""
        url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file&file={gse_id}_RAW.tar"
        output_file = output_dir / f"{gse_id}_RAW.tar"
        
        if output_file.exists():
            print(f"  ✅ Already exists: {output_file.name}")
            return output_file
        
        try:
            print(f"  Downloading raw files for {gse_id}...")
            response = self.session.get(url, timeout=300, stream=True)
            
            if response.status_code == 200:
                total_size = int(response.headers.get('content-length', 0))
                print(f"    Size: {total_size / 1024 / 1024:.1f} MB")
                
                with open(output_file, 'wb') as f:
                    downloaded = 0
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                            downloaded += len(chunk)
                            if total_size > 0 and downloaded % (10 * 1024 * 1024) == 0:
                                progress = (downloaded / total_size) * 100
                                print(f"    Progress: {progress:.1f}%", end='\r')
                
                print(f"\n  ✅ Downloaded: {output_file.name}")
                return output_file
            else:
                print(f"  ⚠️ Raw files not available (HTTP {response.status_code})")
                return None
                
        except Exception as e:
            print(f"  ⚠️ Error downloading raw files: {e}")
            return None
    
    def download_geo_dataset(self, gse_id: str, dataset_type: str = "unknown"):
        """Download a complete GEO dataset"""
        print(f"\n{'='*80}")
        print(f"Processing {gse_id} ({dataset_type})")
        print(f"{'='*80}")
        
        # Determine output directory
        if dataset_type == "blood_rnaseq":
            output_dir = self.base_dir / "blood_rnaseq" / gse_id
        elif dataset_type == "tissue_rnaseq":
            output_dir = self.base_dir / "sc_rnaseq" / gse_id
        elif dataset_type == "microarray":
            output_dir = self.base_dir / "microarray" / gse_id
        else:
            output_dir = self.base_dir / "other" / gse_id
        
        downloaded_files = []
        
        # Download series matrix
        matrix_file = self.download_geo_series_matrix(gse_id, output_dir)
        if matrix_file:
            downloaded_files.append(matrix_file)
        
        # Download raw files
        raw_file = self.download_geo_raw_files(gse_id, output_dir)
        if raw_file:
            downloaded_files.append(raw_file)
        
        time.sleep(1)  # Rate limiting
        
        return downloaded_files
    
    def download_all_aneurysm_datasets(self):
        """Download all known aneurysm datasets"""
        print("=" * 80)
        print("DOWNLOADING ALL ANEURYSM DATASETS")
        print("=" * 80)
        
        # Known aneurysm datasets
        datasets_to_download = {
            # Blood RNA-seq
            'GSE159610': {
                'type': 'blood_rnaseq',
                'name': 'Blood transcriptomics (aneurysm vs control)',
                'priority': 'high'
            },
            
            # Tissue RNA-seq
            'GSE122897': {
                'type': 'tissue_rnaseq',
                'name': 'Tissue RNA-seq (ruptured vs unruptured)',
                'priority': 'high'
            },
            'GSE158558': {
                'type': 'tissue_rnaseq',
                'name': 'Tissue RNA-seq (intracranial aneurysm)',
                'priority': 'high'
            },
            'GSE66238': {
                'type': 'tissue_rnaseq',
                'name': 'Tissue RNA-seq (intracranial aneurysm)',
                'priority': 'medium'
            },
            
            # Microarray
            'GSE13353': {
                'type': 'microarray',
                'name': 'Aneurysm microarray',
                'priority': 'medium'
            },
            'GSE54083': {
                'type': 'microarray',
                'name': 'Aneurysm microarray',
                'priority': 'medium'
            },
            'GSE75436': {
                'type': 'microarray',
                'name': 'Aneurysm microarray',
                'priority': 'medium'
            },
            
            # Verify these
            'GSE186927': {
                'type': 'tissue_rnaseq',
                'name': 'GSE186927 (verify if aneurysm-related)',
                'priority': 'low'
            },
        }
        
        downloaded = {}
        failed = {}
        
        for gse_id, info in datasets_to_download.items():
            files = self.download_geo_dataset(gse_id, info['type'])
            
            if files:
                downloaded[gse_id] = {
                    'files': [str(f) for f in files],
                    'count': len(files)
                }
            else:
                failed[gse_id] = info['name']
        
        return downloaded, failed
    
    def search_additional_gwas_studies(self):
        """Search for additional GWAS studies"""
        print("\n" + "=" * 80)
        print("ADDITIONAL GWAS STUDIES")
        print("=" * 80)
        
        additional_gwas = {
            'Foroud_2014': {
                'study': 'Foroud et al. 2014 (Japanese cohort)',
                'status': 'Check journal supplement',
                'url': 'https://pubmed.ncbi.nlm.nih.gov/25355837/'
            },
            'Bilguvar_2008': {
                'study': 'Bilguvar et al. 2008 (Nature Genetics)',
                'status': 'Check journal supplement',
                'url': 'https://pubmed.ncbi.nlm.nih.gov/18758463/'
            }
        }
        
        print("\nAdditional GWAS studies to check:")
        for study_id, info in additional_gwas.items():
            print(f"\n  {study_id}:")
            print(f"    Study: {info['study']}")
            print(f"    Status: {info['status']}")
            print(f"    URL: {info['url']}")
            print(f"    Action: Download from journal supplement or contact authors")
    
    def check_gtex_additional_data(self):
        """Check for additional GTEx data"""
        print("\n" + "=" * 80)
        print("GTEx ADDITIONAL DATA")
        print("=" * 80)
        
        gtex_dir = self.base_dir / "gtex"
        existing_files = list(gtex_dir.glob("*.gct.gz")) if gtex_dir.exists() else []
        
        print(f"\nCurrent GTEx files: {len(existing_files)}")
        
        if len(existing_files) < 6:
            print("\n  ⚠️ Some GTEx vessel tissue files may be missing")
            print("  Required: TPM and read count files for:")
            print("    - Artery Aorta")
            print("    - Artery Coronary")
            print("    - Artery Tibial")
    
    def generate_download_report(self, downloaded: Dict, failed: Dict):
        """Generate comprehensive download report"""
        print("\n" + "=" * 80)
        print("DOWNLOAD SUMMARY REPORT")
        print("=" * 80)
        
        print(f"\n✅ Successfully downloaded: {len(downloaded)} datasets")
        total_files = 0
        for gse_id, info in downloaded.items():
            print(f"\n  {gse_id}:")
            print(f"    Files: {info['count']}")
            for filepath in info['files']:
                file_path = Path(filepath)
                if file_path.exists():
                    size_mb = file_path.stat().st_size / (1024 * 1024)
                    print(f"      ✅ {file_path.name} ({size_mb:.1f} MB)")
                    total_files += 1
        
        if failed:
            print(f"\n❌ Failed to download: {len(failed)} datasets")
            for gse_id, name in failed.items():
                print(f"  {gse_id}: {name}")
        
        print(f"\nTotal files downloaded: {total_files}")
        
        # Save report
        report = {
            'downloaded': downloaded,
            'failed': failed,
            'total_datasets': len(downloaded) + len(failed),
            'total_files': total_files
        }
        
        report_file = self.base_dir / "download_report.json"
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        print(f"\n📄 Report saved to: {report_file}")
    
    def create_data_requirements_document(self):
        """Create document listing data requirements vs availability"""
        print("\n" + "=" * 80)
        print("DATA REQUIREMENTS ASSESSMENT")
        print("=" * 80)
        
        requirements = {
            'available_publicly': {
                'blood_rnaseq': ['GSE159610'],
                'tissue_rnaseq': ['GSE122897', 'GSE158558', 'GSE66238'],
                'microarray': ['GSE13353', 'GSE54083', 'GSE75436'],
                'gwas_summary': ['Bakker_2020'],
                'gtex_reference': ['Artery Aorta', 'Artery Coronary', 'Artery Tibial']
            },
            'requires_approval': {
                'patient_genotypes': ['dbGaP', 'UK Biobank', 'EGA'],
                'patient_outcomes': ['Contact authors', 'Journal supplements'],
                'detailed_clinical': ['PHASES scores', 'Rupture dates']
            },
            'missing_but_needed': {
                'individual_prs_data': 'Patient-level genotype data for PRS calculation',
                'outcome_data': 'Patient-level rupture outcomes for model training',
                'exosomal_rna': 'Exosomal RNA sequencing data',
                'scrnaseq_aneurysm': 'Single-cell RNA-seq of aneurysm tissue'
            }
        }
        
        print("\n📊 Data Availability Summary:")
        print("\n✅ Available Publicly:")
        for category, datasets in requirements['available_publicly'].items():
            print(f"\n  {category.replace('_', ' ').title()}:")
            for ds in datasets:
                print(f"    ✅ {ds}")
        
        print("\n⚠️ Requires Approval/Access:")
        for category, sources in requirements['requires_approval'].items():
            print(f"\n  {category.replace('_', ' ').title()}:")
            for source in sources:
                print(f"    ⚠️ {source}")
        
        print("\n❌ Missing But Needed:")
        for item, description in requirements['missing_but_needed'].items():
            print(f"  ❌ {item.replace('_', ' ').title()}: {description}")
        
        # Save requirements
        req_file = self.base_dir / "data_requirements.json"
        with open(req_file, 'w') as f:
            json.dump(requirements, f, indent=2)
        
        print(f"\n📄 Requirements document saved to: {req_file}")

def main():
    """Main download pipeline"""
    print("=" * 80)
    print("COMPREHENSIVE ANEURYSM DATASET DOWNLOAD")
    print("Downloading ALL available public datasets")
    print("=" * 80)
    
    downloader = AneurysmDatasetDownloader()
    
    # 1. Download all GEO datasets
    downloaded, failed = downloader.download_all_aneurysm_datasets()
    
    # 2. Check additional GWAS studies
    downloader.search_additional_gwas_studies()
    
    # 3. Check GTEx data
    downloader.check_gtex_additional_data()
    
    # 4. Generate reports
    downloader.generate_download_report(downloaded, failed)
    
    # 5. Create requirements assessment
    downloader.create_data_requirements_document()
    
    print("\n" + "=" * 80)
    print("DOWNLOAD PIPELINE COMPLETE!")
    print("=" * 80)
    print("\nNext Steps:")
    print("  1. Review download report: data/download_report.json")
    print("  2. Review data requirements: data/data_requirements.json")
    print("  3. For missing data requiring approval, follow instructions in requirements document")

if __name__ == "__main__":
    main()






