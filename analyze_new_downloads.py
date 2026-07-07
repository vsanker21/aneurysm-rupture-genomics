"""
Analyze New Downloads from Downloads Folder
Extract useful information from jpm-13-00266.pdf and jpm-2165218-supplementary.pdf
"""
import sys
import os
from pathlib import Path
import json
import pandas as pd
import re
from typing import Dict, List, Optional

# Add extraction module to path
sys.path.insert(0, str(Path(__file__).parent))

# Import the PDF extraction class
from extract_pdf_data_comprehensive import PDFDataExtractor, PDF2_AVAILABLE, PDFPLUMBER_AVAILABLE, PYMUPDF_AVAILABLE

def get_downloads_folder():
    """Get the Downloads folder path"""
    return Path.home() / "Downloads"

def find_jpm_files():
    """Find the JPM PDF files in Downloads"""
    downloads = get_downloads_folder()
    
    files_to_analyze = []
    
    # Main paper
    main_pdf = downloads / "jpm-13-00266.pdf"
    if main_pdf.exists():
        files_to_analyze.append(main_pdf)
        print(f"[OK] Found main paper: {main_pdf.name}")
    
    # Supplementary PDF from ZIP
    supp_pdf = downloads / "jpm-13-00266-extracted" / "jpm-2165218-supplementary.pdf"
    if supp_pdf.exists():
        files_to_analyze.append(supp_pdf)
        print(f"[OK] Found supplementary PDF: {supp_pdf.name}")
    
    return files_to_analyze

def analyze_jpm_pdfs():
    """Analyze the JPM PDFs for useful data"""
    
    print("=" * 80)
    print("ANALYZING NEW JPM DOWNLOADS")
    print("=" * 80)
    
    # Check available libraries
    print("\n[INFO] Available PDF libraries:")
    print(f"  PyPDF2: {'YES' if PDF2_AVAILABLE else 'NO'}")
    print(f"  pdfplumber: {'YES' if PDFPLUMBER_AVAILABLE else 'NO'}")
    print(f"  pymupdf: {'YES' if PYMUPDF_AVAILABLE else 'NO'}")
    
    if not any([PDF2_AVAILABLE, PDFPLUMBER_AVAILABLE, PYMUPDF_AVAILABLE]):
        print("\n[ERROR] No PDF extraction libraries available!")
        print("   Please install: pip install PyPDF2 pdfplumber pymupdf")
        return
    
    # Find PDFs
    pdf_files = find_jpm_files()
    
    if not pdf_files:
        print("\n[ERROR] No JPM PDF files found in Downloads folder")
        print("   Expected files:")
        print("   - jpm-13-00266.pdf")
        print("   - jpm-13-00266-extracted/jpm-2165218-supplementary.pdf")
        return
    
    print(f"\n[INFO] Found {len(pdf_files)} PDFs to analyze")
    
    # Create output directory
    output_dir = Path("results/jpm_analysis")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Analyze each PDF
    all_results = {}
    
    for pdf_path in pdf_files:
        print(f"\n{'='*80}")
        print(f"Processing: {pdf_path.name}")
        print("=" * 80)
        
        # Create extractor (use output dir in results folder)
        extractor = PDFDataExtractor(pdf_path)
        extractor.output_dir = output_dir / pdf_path.stem
        extractor.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Extract all data
        extracted = extractor.extract_all()
        all_results[pdf_path.name] = extracted
        
        # Check for specific useful data
        text = extracted.get('text', '')
        if text:
            # Check for aneurysm-related keywords
            aneurysm_keywords = [
                'aneurysm', 'rupture', 'cerebral', 'intracranial',
                'gwas', 'genome-wide', 'snp', 'genetic',
                'risk', 'biomarker', 'gene', 'expression',
                'transcriptome', 'rna-seq', 'blood'
            ]
            
            keyword_counts = {}
            for keyword in aneurysm_keywords:
                count = len(re.findall(rf'\b{keyword}\b', text, re.IGNORECASE))
                if count > 0:
                    keyword_counts[keyword] = count
            
            if keyword_counts:
                print(f"\n  [INFO] Relevant keywords found:")
                for keyword, count in sorted(keyword_counts.items(), key=lambda x: x[1], reverse=True)[:10]:
                    print(f"    - {keyword}: {count} occurrences")
            
            # Check for dataset/GEO accession numbers
            geo_pattern = r'GSE\d+'
            geo_matches = re.findall(geo_pattern, text)
            if geo_matches:
                unique_geo = list(set(geo_matches))
                print(f"\n  [INFO] GEO Datasets mentioned: {', '.join(unique_geo)}")
            
            # Check for PubMed IDs
            pubmed_pattern = r'PMID[:\s]+(\d+)|pubmed[:\s]+(\d+)|pmid[:\s]+(\d+)'
            pubmed_matches = re.findall(pubmed_pattern, text, re.IGNORECASE)
            if pubmed_matches:
                pmids = [m[0] or m[1] or m[2] for m in pubmed_matches]
                unique_pmids = list(set(pmids))
                print(f"\n  [INFO] PubMed IDs: {', '.join(unique_pmids[:5])}")
        
        # Check tables for useful data
        tables = extracted.get('tables', [])
        if tables:
            print(f"\n  [INFO] Analyzing {len(tables)} tables for useful data...")
            for i, table in enumerate(tables):
                if table is not None and not table.empty:
                    # Check if table contains GWAS data
                    table_text = table.to_string()
                    if re.search(r'rs\d+|chr\s*\d+|p\s*[<>=]?\s*[\d.e-]+|odds\s*ratio', table_text, re.IGNORECASE):
                        print(f"    [OK] Table {i+1}: Contains GWAS-like data")
                    
                    # Check if table contains gene names
                    if re.search(r'\b[A-Z]{2,8}\b|\bgene\b', table_text, re.IGNORECASE):
                        print(f"    [OK] Table {i+1}: Contains gene names")
        
        # GWAS analysis
        gwas_analysis = extracted.get('gwas_analysis', {})
        if gwas_analysis.get('is_gwas'):
            print(f"\n  [OK] GWAS data detected!")
            indicators = gwas_analysis.get('indicators', {})
            print(f"    - RS numbers: {indicators.get('rs_numbers', 0)}")
            print(f"    - P-values: {indicators.get('p_value', 0)}")
            print(f"    - Odds ratios: {indicators.get('odds_ratio', 0)}")
    
    # Generate summary report
    print("\n" + "=" * 80)
    print("ANALYSIS SUMMARY")
    print("=" * 80)
    
    summary_data = []
    
    for pdf_name, data in all_results.items():
        paper_info = data.get('paper_info', {})
        
        summary = {
            'File': pdf_name,
            'Title': paper_info.get('title', 'N/A')[:100] if paper_info.get('title') else 'N/A',
            'Year': paper_info.get('year', 'N/A'),
            'Journal': paper_info.get('journal', 'N/A'),
            'PubMed_ID': paper_info.get('pubmed_id', 'N/A'),
            'DOI': paper_info.get('doi', 'N/A'),
            'Text_Length': len(data.get('text', '')),
            'Tables': len(data.get('tables', [])),
            'Has_GWAS': 'Yes' if data.get('gwas_analysis', {}).get('is_gwas') else 'No',
        }
        
        summary_data.append(summary)
        
        print(f"\n[FILE] {pdf_name}")
        print(f"  Title: {summary['Title']}")
        print(f"  Year: {summary['Year']}")
        print(f"  Journal: {summary['Journal']}")
        if summary['PubMed_ID'] != 'N/A':
            print(f"  PubMed ID: {summary['PubMed_ID']}")
        print(f"  Text length: {summary['Text_Length']:,} characters")
        print(f"  Tables: {summary['Tables']}")
        print(f"  GWAS data: {summary['Has_GWAS']}")
    
    # Save summary
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_file = output_dir / "jpm_analysis_summary.tsv"
        summary_df.to_csv(summary_file, sep='\t', index=False)
        print(f"\n[OK] Summary saved to: {summary_file}")
    
    # Check for actionable data
    print("\n" + "=" * 80)
    print("ACTIONABLE DATA CHECK")
    print("=" * 80)
    
    actionable_findings = []
    
    for pdf_name, data in all_results.items():
        text = data.get('text', '')
        tables = data.get('tables', [])
        gwas_analysis = data.get('gwas_analysis', {})
        
        findings = []
        
        # Check for GWAS summary statistics
        if gwas_analysis.get('is_gwas'):
            findings.append("[OK] Contains GWAS summary statistics (could extract SNP data)")
        
        # Check for tables with data
        if tables:
            findings.append(f"[OK] Contains {len(tables)} tables (could extract structured data)")
        
        # Check for GEO datasets
        geo_matches = re.findall(r'GSE\d+', text)
        if geo_matches:
            unique_geo = list(set(geo_matches))
            findings.append(f"[OK] Mentions GEO datasets: {', '.join(unique_geo)}")
        
        # Check for gene lists
        gene_list_patterns = [
            r'Table\s+\d+.*genes?',
            r'Supplementary\s+Table\s+\d+.*genes?',
            r'gene\s+list',
            r'differentially\s+expressed\s+genes?'
        ]
        if any(re.search(pattern, text, re.IGNORECASE) for pattern in gene_list_patterns):
            findings.append("[OK] Contains gene lists or differential expression results")
        
        if findings:
            actionable_findings.append({
                'file': pdf_name,
                'findings': findings
            })
    
    if actionable_findings:
        print("\n[OK] Actionable data found:")
        for item in actionable_findings:
            print(f"\n  [FILE] {item['file']}:")
            for finding in item['findings']:
                print(f"    {finding}")
    else:
        print("\n[WARNING] No clearly actionable data found (may need manual review)")
    
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE!")
    print("=" * 80)
    print(f"\n[OK] All extracted data saved to: {output_dir}/")
    print(f"[OK] Check the extracted text files and tables for detailed information")

if __name__ == "__main__":
    analyze_jpm_pdfs()

