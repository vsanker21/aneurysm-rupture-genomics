"""
FINAL COMPREHENSIVE DATASET DOWNLOADER
Ensures ALL publicly available datasets are downloaded using multiple methods
"""
import os
import requests
import time
from pathlib import Path
import json
import subprocess
import sys

class FinalComprehensiveDownloader:
    """Final comprehensive downloader using all available methods"""
    
    def __init__(self):
        self.base_dir = Path("data")
        self.base_dir.mkdir(exist_ok=True)
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
        })
        self.downloaded = []
        self.failed = []
        
    def download_geo_via_multiple_methods(self, gse_id: str, output_dir: Path):
        """Try downloading using multiple methods"""
        output_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"\n📥 {gse_id}: Trying all download methods...")
        
        # Method 1: Direct FTP series matrix
        methods = [
            {
                'name': 'FTP Series Matrix',
                'url': f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_id[:6]}nnn/{gse_id}/matrix/{gse_id}_series_matrix.txt.gz",
                'filename': f"{gse_id}_series_matrix.txt.gz"
            },
            {
                'name': 'GEO Download API',
                'url': f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file&file={gse_id}_series_matrix.txt.gz",
                'filename': f"{gse_id}_series_matrix.txt.gz"
            },
            {
                'name': 'FTP RAW files',
                'url': f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_id[:6]}nnn/{gse_id}/suppl/{gse_id}_RAW.tar",
                'filename': f"{gse_id}_RAW.tar"
            },
            {
                'name': 'GEO RAW API',
                'url': f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file&file={gse_id}_RAW.tar",
                'filename': f"{gse_id}_RAW.tar"
            },
            {
                'name': 'FTP SOFT file',
                'url': f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_id[:6]}nnn/{gse_id}/soft/{gse_id}_family.soft.gz",
                'filename': f"{gse_id}_family.soft.gz"
            },
            {
                'name': 'GEO SOFT API',
                'url': f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file&file={gse_id}_family.soft.gz",
                'filename': f"{gse_id}_family.soft.gz"
            }
        ]
        
        success = False
        for method in methods:
            output_file = output_dir / method['filename']
            
            # Skip if file exists and has content
            if output_file.exists() and output_file.stat().st_size > 1000:  # > 1KB
                print(f"  ✅ Already exists: {method['filename']} ({output_file.stat().st_size / 1024:.1f} KB)")
                success = True
                self.downloaded.append(str(output_file))
                continue
            
            try:
                print(f"  Trying {method['name']}...")
                response = self.session.get(method['url'], timeout=180, stream=True, allow_redirects=True)
                
                if response.status_code == 200:
                    # Check if we got actual content (not HTML error page)
                    content_length = response.headers.get('Content-Length')
                    if content_length and int(content_length) < 1000:
                        print(f"    ⚠️ File too small, might be error page")
                        continue
                    
                    with open(output_file, 'wb') as f:
                        for chunk in response.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)
                    
                    # Verify file has content
                    if output_file.stat().st_size > 1000:
                        size_mb = output_file.stat().st_size / (1024 * 1024)
                        print(f"    ✅ Downloaded: {method['filename']} ({size_mb:.1f} MB)")
                        self.downloaded.append(str(output_file))
                        success = True
                        break
                    else:
                        print(f"    ⚠️ File too small, deleting...")
                        output_file.unlink()
                        
            except Exception as e:
                print(f"    ⚠️ Error: {str(e)[:50]}")
                continue
        
        if not success:
            print(f"  ❌ Failed to download {gse_id}")
            self.failed.append(gse_id)
        
        return success
    
    def download_all_known_datasets(self):
        """Download ALL known datasets"""
        print("=" * 80)
        print("DOWNLOADING ALL KNOWN ANEURYSM DATASETS")
        print("=" * 80)
        
        # Complete list of ALL known datasets
        all_datasets = {
            # Blood RNA-seq
            'GSE159610': {
                'type': 'blood_rnaseq',
                'desc': 'Blood RNA-seq (aneurysm vs control)',
                'priority': 1
            },
            
            # Tissue RNA-seq
            'GSE122897': {
                'type': 'sc_rnaseq',
                'desc': 'Tissue RNA-seq (ruptured vs unruptured)',
                'priority': 1
            },
            'GSE158558': {
                'type': 'sc_rnaseq',
                'desc': 'Tissue RNA-seq (intracranial aneurysm)',
                'priority': 1
            },
            'GSE66238': {
                'type': 'sc_rnaseq',
                'desc': 'Tissue RNA-seq (intracranial aneurysm)',
                'priority': 1
            },
            
            # Microarray
            'GSE13353': {
                'type': 'microarray',
                'desc': 'Aneurysm microarray',
                'priority': 1
            },
            'GSE54083': {
                'type': 'microarray',
                'desc': 'Aneurysm microarray',
                'priority': 1
            },
            'GSE75436': {
                'type': 'microarray',
                'desc': 'Aneurysm microarray',
                'priority': 1
            },
            
            # Additional potential
            'GSE186927': {
                'type': 'sc_rnaseq',
                'desc': 'Potential scRNA-seq (verify)',
                'priority': 2
            },
        }
        
        # Sort by priority
        sorted_datasets = sorted(all_datasets.items(), key=lambda x: x[1]['priority'])
        
        for gse_id, info in sorted_datasets:
            output_dir = self.base_dir / info['type'] / gse_id
            self.download_geo_via_multiple_methods(gse_id, output_dir)
            time.sleep(1)  # Be nice to servers
    
    def verify_all_downloads(self):
        """Verify all downloads are complete"""
        print("\n" + "=" * 80)
        print("VERIFYING ALL DOWNLOADS")
        print("=" * 80)
        
        verification = {
            'complete': [],
            'incomplete': [],
            'missing': []
        }
        
        expected_datasets = {
            'GSE159610': 'blood_rnaseq',
            'GSE122897': 'sc_rnaseq',
            'GSE158558': 'sc_rnaseq',
            'GSE66238': 'sc_rnaseq',
            'GSE13353': 'microarray',
            'GSE54083': 'microarray',
            'GSE75436': 'microarray',
        }
        
        for gse_id, data_type in expected_datasets.items():
            data_dir = self.base_dir / data_type / gse_id
            
            if not data_dir.exists():
                verification['missing'].append(gse_id)
                print(f"  ❌ {gse_id}: MISSING")
                continue
            
            files = list(data_dir.glob("*"))
            data_files = [f for f in files if f.is_file() and f.stat().st_size > 1000]
            
            if data_files:
                total_size = sum(f.stat().st_size for f in data_files) / (1024 * 1024)
                verification['complete'].append({
                    'gse_id': gse_id,
                    'files': len(data_files),
                    'size_mb': round(total_size, 2)
                })
                print(f"  ✅ {gse_id}: {len(data_files)} files ({total_size:.1f} MB)")
            else:
                verification['incomplete'].append(gse_id)
                print(f"  ⚠️ {gse_id}: INCOMPLETE (no data files)")
        
        return verification
    
    def generate_final_report(self, verification):
        """Generate final comprehensive report"""
        print("\n" + "=" * 80)
        print("FINAL DOWNLOAD REPORT")
        print("=" * 80)
        
        report = {
            'download_summary': {
                'total_downloaded': len(self.downloaded),
                'total_failed': len(self.failed),
                'failed_datasets': self.failed
            },
            'verification': verification,
            'downloaded_files': self.downloaded[:20]  # First 20
        }
        
        report_file = self.base_dir / "final_download_report.json"
        with open(report_file, 'w', encoding='utf-8') as f:
            json.dump(report, f, indent=2)
        
        print(f"\n✅ Report saved: {report_file}")
        print(f"\n📊 Summary:")
        print(f"  Complete datasets: {len(verification['complete'])}")
        print(f"  Incomplete datasets: {len(verification['incomplete'])}")
        print(f"  Missing datasets: {len(verification['missing'])}")
        print(f"  Downloaded files: {len(self.downloaded)}")
        
        if verification['complete']:
            print(f"\n✅ Successfully downloaded datasets:")
            for item in verification['complete']:
                print(f"    - {item['gse_id']}: {item['files']} files ({item['size_mb']} MB)")
        
        if verification['incomplete'] or verification['missing']:
            print(f"\n⚠️ Datasets needing attention:")
            for gse_id in verification['incomplete'] + verification['missing']:
                print(f"    - {gse_id}")

def main():
    """Main download pipeline"""
    print("=" * 80)
    print("FINAL COMPREHENSIVE DATASET DOWNLOADER")
    print("Ensuring ALL publicly available datasets are downloaded")
    print("=" * 80)
    
    downloader = FinalComprehensiveDownloader()
    
    # Download all known datasets
    downloader.download_all_known_datasets()
    
    # Verify downloads
    verification = downloader.verify_all_downloads()
    
    # Generate report
    downloader.generate_final_report(verification)
    
    print("\n" + "=" * 80)
    print("FINAL DOWNLOAD COMPLETE!")
    print("=" * 80)

if __name__ == "__main__":
    main()






