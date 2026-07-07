"""
Analyze New Downloads - PDF and ZIP/TAR files
Analyze 2025.11.06.25339724.full.pdf and GSE26969_RAW.tar
"""
import sys
import os
from pathlib import Path
import json
import re

# Add extraction module to path
sys.path.insert(0, str(Path(__file__).parent))

# Import the PDF extraction class
from extract_pdf_data_comprehensive import PDFDataExtractor, PDF2_AVAILABLE, PDFPLUMBER_AVAILABLE, PYMUPDF_AVAILABLE

def analyze_new_files():
    """Analyze the new PDF and TAR files"""
    
    print("=" * 80)
    print("ANALYZING NEW DOWNLOADS - PDF AND TAR FILES")
    print("=" * 80)
    
    downloads = Path.home() / "Downloads"
    
    # Files to analyze
    pdf_file = downloads / "2025.11.06.25339724.full.pdf"
    tar_file = downloads / "GSE26969_RAW.tar"
    
    files_to_analyze = []
    
    if pdf_file.exists():
        files_to_analyze.append(('pdf', pdf_file))
        print(f"\n[OK] Found PDF: {pdf_file.name} ({pdf_file.stat().st_size / 1024 / 1024:.2f} MB)")
    
    if tar_file.exists():
        files_to_analyze.append(('tar', tar_file))
        print(f"[OK] Found TAR: {tar_file.name} ({tar_file.stat().st_size / 1024 / 1024:.2f} MB)")
    
    if not files_to_analyze:
        print("\n[ERROR] Files not found!")
        return
    
    # Create output directory
    output_dir = Path("results/new_downloads_analysis")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    results = {}
    
    # Analyze PDF
    if ('pdf', pdf_file) in files_to_analyze:
        print("\n" + "=" * 80)
        print(f"Analyzing PDF: {pdf_file.name}")
        print("=" * 80)
        
        extractor = PDFDataExtractor(pdf_file)
        extractor.output_dir = output_dir / "pdf_analysis"
        extractor.output_dir.mkdir(parents=True, exist_ok=True)
        
        extracted = extractor.extract_all()
        results['pdf'] = extracted
        
        text = extracted.get('text', '')
        if text:
            # Check for aneurysm-related keywords
            aneurysm_keywords = [
                'aneurysm', 'rupture', 'cerebral', 'intracranial',
                'gwas', 'genome-wide', 'snp', 'genetic',
                'risk', 'biomarker', 'gene', 'expression',
                'transcriptome', 'rna-seq', 'blood', 'vascular'
            ]
            
            keyword_counts = {}
            for keyword in aneurysm_keywords:
                count = len(re.findall(rf'\b{keyword}\b', text, re.IGNORECASE))
                if count > 0:
                    keyword_counts[keyword] = count
            
            if keyword_counts:
                print(f"\n  [INFO] Relevant keywords found:")
                for keyword, count in sorted(keyword_counts.items(), key=lambda x: x[1], reverse=True)[:15]:
                    print(f"    - {keyword}: {count} occurrences")
            
            # Extract paper information
            paper_info = extracted.get('paper_info', {})
            if paper_info.get('title'):
                print(f"\n  [INFO] Title: {paper_info['title'][:150]}...")
            if paper_info.get('year'):
                print(f"  [INFO] Year: {paper_info['year']}")
            if paper_info.get('doi'):
                print(f"  [INFO] DOI: {paper_info['doi']}")
            if paper_info.get('pubmed_id'):
                print(f"  [INFO] PubMed ID: {paper_info['pubmed_id']}")
            
            # Check for GEO datasets
            geo_pattern = r'GSE\d+'
            geo_matches = re.findall(geo_pattern, text)
            if geo_matches:
                unique_geo = list(set(geo_matches))
                print(f"\n  [INFO] GEO Datasets mentioned: {', '.join(unique_geo)}")
            
            # Check for GWAS data
            gwas_analysis = extracted.get('gwas_analysis', {})
            if gwas_analysis.get('is_gwas'):
                print(f"\n  [OK] GWAS data detected!")
                indicators = gwas_analysis.get('indicators', {})
                print(f"    - RS numbers: {indicators.get('rs_numbers', 0)}")
                print(f"    - P-values: {indicators.get('p_value', 0)}")
    
    # Analyze TAR file
    if ('tar', tar_file) in files_to_analyze:
        print("\n" + "=" * 80)
        print(f"Analyzing TAR file: {tar_file.name}")
        print("=" * 80)
        
        # Extract TAR file
        extract_dir = output_dir / "gse26969_extracted"
        extract_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"\n  Extracting TAR file to: {extract_dir}")
        try:
            import tarfile
            with tarfile.open(tar_file, 'r') as tar:
                tar.extractall(extract_dir)
                members = tar.getnames()
                
            print(f"  [OK] Extracted {len(members)} files/folders")
            
            # List contents
            print(f"\n  [INFO] Contents (first 20):")
            for item in members[:20]:
                size = "DIR" if item.endswith('/') else "FILE"
                print(f"    {size:4s}: {item}")
            if len(members) > 20:
                print(f"    ... and {len(members) - 20} more items")
            
            # Check for GEO dataset info
            print(f"\n  [INFO] GSE26969 is a GEO dataset")
            print(f"  [INFO] Checking dataset type...")
            
            # Look for common GEO file patterns
            ge_files = list(extract_dir.rglob("*.txt"))
            ge_files += list(extract_dir.rglob("*.CEL"))
            ge_files += list(extract_dir.rglob("*.gz"))
            
            if ge_files:
                print(f"  [INFO] Found {len(ge_files)} potential data files")
                print(f"    Types: {', '.join(set([f.suffix for f in ge_files[:10]]))}")
            
            results['tar'] = {
                'extracted': True,
                'n_files': len(members),
                'files': members[:20]
            }
            
        except Exception as e:
            print(f"  [ERROR] Could not extract TAR: {e}")
            results['tar'] = {'error': str(e)}
    
    # Generate summary
    print("\n" + "=" * 80)
    print("ANALYSIS SUMMARY")
    print("=" * 80)
    
    summary = {
        'pdf_file': str(pdf_file.name) if ('pdf', pdf_file) in files_to_analyze else None,
        'tar_file': str(tar_file.name) if ('tar', tar_file) in files_to_analyze else None,
        'results': results
    }
    
    summary_file = output_dir / "analysis_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\n[OK] Summary saved to: {summary_file}")
    
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE!")
    print("=" * 80)
    
    return results

if __name__ == "__main__":
    analyze_new_files()






