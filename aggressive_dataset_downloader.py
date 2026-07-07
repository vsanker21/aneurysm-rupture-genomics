"""
AGGRESSIVE COMPREHENSIVE DATASET DOWNLOADER
Searches and downloads ALL publicly available aneurysm genomics datasets
"""
import os
import requests
import time
from pathlib import Path
import json
import re
from typing import List, Dict, Optional
from urllib.parse import quote
import xml.etree.ElementTree as ET

class AggressiveDatasetDownloader:
    """Aggressively search and download ALL publicly available datasets"""
    
    def __init__(self):
        self.base_dir = Path("data")
        self.base_dir.mkdir(exist_ok=True)
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
        })
        self.downloaded = []
        self.failed = []
        self.found_datasets = []
        
    def search_geo_database(self, search_terms: List[str]) -> List[Dict]:
        """Search GEO database for aneurysm datasets"""
        print("\n" + "=" * 80)
        print("SEARCHING GEO DATABASE")
        print("=" * 80)
        
        found = []
        
        for term in search_terms:
            print(f"\n🔍 Searching GEO for: '{term}'...")
            try:
                # GEO search API
                url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
                params = {
                    'db': 'gds',
                    'term': f'{term}[All Fields] AND "Homo sapiens"[Organism]',
                    'retmax': 100,
                    'retmode': 'json'
                }
                
                response = self.session.get(url, params=params, timeout=30)
                if response.status_code == 200:
                    data = response.json()
                    if 'esearchresult' in data and 'idlist' in data['esearchresult']:
                        ids = data['esearchresult']['idlist']
                        print(f"  Found {len(ids)} potential datasets")
                        
                        # Get details for each (limit to avoid timeout)
                        for gds_id in ids[:10]:  # Limit to first 10 per search
                            try:
                                summary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                                summary_params = {
                                    'db': 'gds',
                                    'id': gds_id,
                                    'retmode': 'json'
                                }
                                summary_resp = self.session.get(summary_url, params=summary_params, timeout=15)
                                if summary_resp.status_code == 200:
                                    summary_data = summary_resp.json()
                                    if 'result' in summary_data and gds_id in summary_data['result']:
                                        dataset_info = summary_data['result'][gds_id]
                                        title = dataset_info.get('title', 'Unknown')
                                        gse_id = dataset_info.get('gse', '')
                                        
                                        if gse_id and gse_id.startswith('GSE'):
                                            found.append({
                                                'gse_id': gse_id,
                                                'title': title,
                                                'gds_id': gds_id,
                                                'search_term': term
                                            })
                                            print(f"    ✅ {gse_id}: {title[:60]}")
                            except Exception as e:
                                continue  # Skip errors silently to speed up
                
                time.sleep(1)  # Be nice to NCBI
                
            except Exception as e:
                print(f"  ⚠️ Error searching for '{term}': {e}")
                continue
        
        # Remove duplicates
        seen = set()
        unique_found = []
        for item in found:
            if item['gse_id'] not in seen:
                seen.add(item['gse_id'])
                unique_found.append(item)
        
        self.found_datasets.extend(unique_found)
        return unique_found
    
    def download_geo_dataset(self, gse_id: str, output_dir: Path) -> bool:
        """Download a GEO dataset"""
        output_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"\n📥 Downloading {gse_id}...")
        
        # Try multiple download methods
        success = False
        
        # Method 1: Series matrix file
        matrix_urls = [
            f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_id[:6]}nnn/{gse_id}/matrix/{gse_id}_series_matrix.txt.gz",
            f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file&file={gse_id}_series_matrix.txt.gz",
        ]
        
        for url in matrix_urls:
            try:
                filename = f"{gse_id}_series_matrix.txt.gz"
                output_file = output_dir / filename
                
                if output_file.exists():
                    print(f"  ✅ Already exists: {filename}")
                    success = True
                    break
                
                print(f"  Trying: {url}")
                response = self.session.get(url, timeout=120, stream=True)
                
                if response.status_code == 200:
                    with open(output_file, 'wb') as f:
                        for chunk in response.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)
                    
                    size_mb = output_file.stat().st_size / (1024 * 1024)
                    print(f"  ✅ Downloaded: {filename} ({size_mb:.1f} MB)")
                    self.downloaded.append(str(output_file))
                    success = True
                    break
            except Exception as e:
                print(f"  ⚠️ Error: {e}")
                continue
        
        # Method 2: RAW files
        if not success:
            raw_urls = [
                f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_id[:6]}nnn/{gse_id}/suppl/{gse_id}_RAW.tar",
                f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file&file={gse_id}_RAW.tar",
            ]
            
            for url in raw_urls:
                try:
                    filename = f"{gse_id}_RAW.tar"
                    output_file = output_dir / filename
                    
                    if output_file.exists():
                        print(f"  ✅ Already exists: {filename}")
                        success = True
                        break
                    
                    print(f"  Trying RAW files: {url}")
                    response = self.session.get(url, timeout=300, stream=True)
                    
                    if response.status_code == 200:
                        with open(output_file, 'wb') as f:
                            for chunk in response.iter_content(chunk_size=8192):
                                if chunk:
                                    f.write(chunk)
                        
                        size_mb = output_file.stat().st_size / (1024 * 1024)
                        print(f"  ✅ Downloaded: {filename} ({size_mb:.1f} MB)")
                        self.downloaded.append(str(output_file))
                        success = True
                        break
                except Exception as e:
                    print(f"  ⚠️ Error: {e}")
                    continue
        
        # Method 3: Soft file
        if not success:
            soft_urls = [
                f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_id[:6]}nnn/{gse_id}/soft/{gse_id}_family.soft.gz",
                f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file&file={gse_id}_family.soft.gz",
            ]
            
            for url in soft_urls:
                try:
                    filename = f"{gse_id}_family.soft.gz"
                    output_file = output_dir / filename
                    
                    if output_file.exists():
                        print(f"  ✅ Already exists: {filename}")
                        success = True
                        break
                    
                    print(f"  Trying SOFT file: {url}")
                    response = self.session.get(url, timeout=120, stream=True)
                    
                    if response.status_code == 200:
                        with open(output_file, 'wb') as f:
                            for chunk in response.iter_content(chunk_size=8192):
                                if chunk:
                                    f.write(chunk)
                        
                        size_mb = output_file.stat().st_size / (1024 * 1024)
                        print(f"  ✅ Downloaded: {filename} ({size_mb:.1f} MB)")
                        self.downloaded.append(str(output_file))
                        success = True
                        break
                except Exception as e:
                    print(f"  ⚠️ Error: {e}")
                    continue
        
        if not success:
            print(f"  ❌ Failed to download {gse_id}")
            self.failed.append(gse_id)
        
        return success
    
    def search_arrayexpress(self, search_terms: List[str]) -> List[Dict]:
        """Search ArrayExpress for aneurysm datasets"""
        print("\n" + "=" * 80)
        print("SEARCHING ARRAYEXPRESS DATABASE")
        print("=" * 80)
        
        found = []
        
        for term in search_terms:
            print(f"\n🔍 Searching ArrayExpress for: '{term}'...")
            try:
                url = "https://www.ebi.ac.uk/arrayexpress/json/v3/experiments"
                params = {
                    'keywords': term,
                    'species': 'homo sapiens',
                    'pageSize': 50
                }
                
                response = self.session.get(url, params=params, timeout=30)
                if response.status_code == 200:
                    data = response.json()
                    if 'experiments' in data and 'experiment' in data['experiments']:
                        experiments = data['experiments']['experiment']
                        print(f"  Found {len(experiments)} potential datasets")
                        
                        for exp in experiments:
                            accession = exp.get('accession', '')
                            title = exp.get('name', 'Unknown')
                            
                            if accession:
                                found.append({
                                    'accession': accession,
                                    'title': title,
                                    'search_term': term,
                                    'source': 'ArrayExpress'
                                })
                                print(f"    ✅ {accession}: {title[:80]}")
                
                time.sleep(1)
                
            except Exception as e:
                print(f"  ⚠️ Error searching for '{term}': {e}")
                continue
        
        return found
    
    def download_gwas_summary_stats(self):
        """Download GWAS summary statistics from various sources"""
        print("\n" + "=" * 80)
        print("DOWNLOADING GWAS SUMMARY STATISTICS")
        print("=" * 80)
        
        gwas_dir = self.base_dir / "gwas"
        gwas_dir.mkdir(exist_ok=True)
        
        # Known GWAS sources
        gwas_sources = [
            {
                'name': 'Bakker_2020',
                'urls': [
                    'https://www.ebi.ac.uk/gwas/api/search/downloads/field?query=ieu-a-1125',
                    'https://gwas.mrcieu.ac.uk/files/ieu-a-1125/ieu-a-1125.vcf.gz',
                ],
                'dir': gwas_dir / "bakker_2020"
            },
            {
                'name': 'Foroud_2014',
                'urls': [
                    'https://www.nature.com/articles/ng.3094',
                ],
                'dir': gwas_dir / "foroud_2014",
                'note': 'Check journal supplement'
            },
            {
                'name': 'Bilguvar_2008',
                'urls': [
                    'https://www.nature.com/articles/ng.240',
                ],
                'dir': gwas_dir / "bilguvar_2008",
                'note': 'Check journal supplement'
            }
        ]
        
        for source in gwas_sources:
            source['dir'].mkdir(exist_ok=True)
            print(f"\n📋 {source['name']}:")
            
            if 'note' in source:
                print(f"  ⚠️ {source['note']}")
                continue
            
            for url in source['urls']:
                try:
                    filename = url.split('/')[-1]
                    output_file = source['dir'] / filename
                    
                    if output_file.exists():
                        print(f"  ✅ Already exists: {filename}")
                        continue
                    
                    print(f"  Trying: {url}")
                    response = self.session.get(url, timeout=300, stream=True, allow_redirects=True)
                    
                    if response.status_code == 200:
                        with open(output_file, 'wb') as f:
                            for chunk in response.iter_content(chunk_size=8192):
                                if chunk:
                                    f.write(chunk)
                        
                        size_mb = output_file.stat().st_size / (1024 * 1024)
                        print(f"  ✅ Downloaded: {filename} ({size_mb:.1f} MB)")
                        self.downloaded.append(str(output_file))
                    else:
                        print(f"  ⚠️ HTTP {response.status_code}")
                        
                except Exception as e:
                    print(f"  ⚠️ Error: {e}")
                    continue
    
    def download_all_found_datasets(self):
        """Download all datasets found during search"""
        print("\n" + "=" * 80)
        print("DOWNLOADING ALL FOUND DATASETS")
        print("=" * 80)
        
        # Organize by type
        geo_datasets = [d for d in self.found_datasets if 'gse_id' in d]
        
        for dataset in geo_datasets:
            gse_id = dataset['gse_id']
            
            # Determine output directory based on dataset type
            if 'blood' in dataset['title'].lower() or 'serum' in dataset['title'].lower():
                output_dir = self.base_dir / "blood_rnaseq" / gse_id
            elif 'tissue' in dataset['title'].lower() or 'aneurysm' in dataset['title'].lower():
                output_dir = self.base_dir / "sc_rnaseq" / gse_id
            elif 'microarray' in dataset['title'].lower() or 'array' in dataset['title'].lower():
                output_dir = self.base_dir / "microarray" / gse_id
            else:
                output_dir = self.base_dir / "additional" / gse_id
            
            self.download_geo_dataset(gse_id, output_dir)
            time.sleep(2)  # Be nice to servers
    
    def download_known_additional_datasets(self):
        """Download known additional datasets that might have been missed"""
        print("\n" + "=" * 80)
        print("DOWNLOADING ALL KNOWN DATASETS (ENSURING COMPLETE COVERAGE)")
        print("=" * 80)
        
        # ALL known GEO datasets for aneurysm research
        known_datasets = {
            # Blood RNA-seq
            'GSE159610': {'type': 'blood_rnaseq', 'desc': 'Blood RNA-seq (aneurysm vs control)'},
            
            # Tissue RNA-seq
            'GSE122897': {'type': 'sc_rnaseq', 'desc': 'Tissue RNA-seq (ruptured vs unruptured)'},
            'GSE158558': {'type': 'sc_rnaseq', 'desc': 'Tissue RNA-seq (intracranial aneurysm)'},
            'GSE66238': {'type': 'sc_rnaseq', 'desc': 'Tissue RNA-seq (intracranial aneurysm)'},
            
            # Microarray
            'GSE13353': {'type': 'microarray', 'desc': 'Aneurysm microarray'},
            'GSE54083': {'type': 'microarray', 'desc': 'Aneurysm microarray'},
            'GSE75436': {'type': 'microarray', 'desc': 'Aneurysm microarray'},
            
            # Additional potential datasets (to try)
            'GSE186927': {'type': 'sc_rnaseq', 'desc': 'Potential scRNA-seq (verify)'},
            'GSE122897': {'type': 'sc_rnaseq', 'desc': 'Already listed above'},
        }
        
        print(f"\n📋 Attempting to download {len(known_datasets)} known datasets...")
        
        for gse_id, info in known_datasets.items():
            if gse_id in [d.get('gse_id', '') for d in self.found_datasets]:
                continue  # Already found in search
            
            output_dir = self.base_dir / info['type'] / gse_id
            print(f"\n📥 {gse_id}: {info['desc']}")
            success = self.download_geo_dataset(gse_id, output_dir)
            if success:
                self.found_datasets.append({
                    'gse_id': gse_id,
                    'title': info['desc'],
                    'source': 'known_list'
                })
            time.sleep(2)
    
    def generate_download_report(self):
        """Generate comprehensive download report"""
        print("\n" + "=" * 80)
        print("DOWNLOAD REPORT")
        print("=" * 80)
        
        report = {
            'downloaded_files': self.downloaded,
            'failed_datasets': self.found_datasets,
            'found_datasets': [d['gse_id'] if 'gse_id' in d else d.get('accession', 'unknown') 
                              for d in self.found_datasets],
            'summary': {
                'total_downloaded': len(self.downloaded),
                'total_found': len(self.found_datasets),
                'total_failed': len(self.failed)
            }
        }
        
        report_file = self.base_dir / "aggressive_download_report.json"
        with open(report_file, 'w', encoding='utf-8') as f:
            json.dump(report, f, indent=2)
        
        print(f"\n✅ Report saved to: {report_file}")
        print(f"\n📊 Summary:")
        print(f"  Found datasets: {len(self.found_datasets)}")
        print(f"  Downloaded files: {len(self.downloaded)}")
        print(f"  Failed: {len(self.failed)}")

def main():
    """Main aggressive download pipeline"""
    print("=" * 80)
    print("AGGRESSIVE COMPREHENSIVE DATASET DOWNLOADER")
    print("Searching and downloading ALL publicly available datasets")
    print("=" * 80)
    
    downloader = AggressiveDatasetDownloader()
    
    # Comprehensive search terms
    search_terms = [
        "intracranial aneurysm",
        "cerebral aneurysm",
        "brain aneurysm",
        "aneurysm rupture",
        "aneurysm transcriptome",
        "aneurysm RNA-seq",
        "aneurysm microarray",
        "aneurysm genomics",
        "vascular aneurysm",
        "subarachnoid hemorrhage",
        "SAH aneurysm",
        "unruptured aneurysm",
        "ruptured aneurysm"
    ]
    
    # 1. Download ALL known datasets first (most important)
    print("\n" + "=" * 80)
    print("STEP 1: DOWNLOADING ALL KNOWN DATASETS")
    print("=" * 80)
    downloader.download_known_additional_datasets()
    
    # 2. Download GWAS summary statistics
    print("\n" + "=" * 80)
    print("STEP 2: DOWNLOADING GWAS SUMMARY STATISTICS")
    print("=" * 80)
    downloader.download_gwas_summary_stats()
    
    # 3. Search GEO for additional datasets (limited search)
    print("\n" + "=" * 80)
    print("STEP 3: SEARCHING FOR ADDITIONAL DATASETS")
    print("=" * 80)
    # Use more focused search terms
    focused_terms = [
        "intracranial aneurysm",
        "cerebral aneurysm",
        "aneurysm rupture",
        "unruptured aneurysm"
    ]
    geo_results = downloader.search_geo_database(focused_terms)
    
    # 4. Search ArrayExpress
    arrayexpress_results = downloader.search_arrayexpress(focused_terms)
    downloader.found_datasets.extend(arrayexpress_results)
    
    # 5. Download all newly found datasets
    print("\n" + "=" * 80)
    print("STEP 4: DOWNLOADING NEWLY FOUND DATASETS")
    print("=" * 80)
    downloader.download_all_found_datasets()
    
    # 6. Generate report
    downloader.generate_download_report()
    
    print("\n" + "=" * 80)
    print("AGGRESSIVE DOWNLOAD COMPLETE!")
    print("=" * 80)
    print(f"\n✅ Found {len(downloader.found_datasets)} datasets")
    print(f"✅ Downloaded {len(downloader.downloaded)} files")
    if downloader.failed:
        print(f"⚠️ Failed to download {len(downloader.failed)} datasets")

if __name__ == "__main__":
    main()

