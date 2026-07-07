"""
Extract data from JPM paper (jpm-13-00266) and supplementary files
"""
import sys
import os
from pathlib import Path
import json
import re
import zipfile
import pandas as pd
import numpy as np
from datetime import datetime

# Add user site-packages to path for PyPDF2
import site
user_site = site.getusersitepackages()
if user_site and user_site not in sys.path:
    sys.path.insert(0, user_site)

try:
    import PyPDF2
    PDF2_AVAILABLE = True
except ImportError:
    PDF2_AVAILABLE = False
    print("[WARNING] PyPDF2 not available")


def extract_text_from_pdf(pdf_path: Path) -> str:
    """Extract text from PDF using PyPDF2"""
    if not PDF2_AVAILABLE:
        return ""
    
    try:
        text = ""
        with open(pdf_path, 'rb') as pdf_file:
            reader = PyPDF2.PdfReader(pdf_file)
            print(f"  PDF has {len(reader.pages)} pages")
            for i, page in enumerate(reader.pages):
                page_text = page.extract_text()
                if page_text:
                    text += page_text + "\n"
                if (i + 1) % 10 == 0:
                    print(f"  Processed {i + 1} pages...")
        return text
    except Exception as e:
        print(f"  [ERROR] Error extracting text: {e}")
        return ""


def extract_zip_contents(zip_path: Path, output_dir: Path):
    """Extract and analyze ZIP file contents"""
    print(f"\n[INFO] Extracting ZIP file: {zip_path.name}")
    
    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            file_list = zip_ref.namelist()
            print(f"  Found {len(file_list)} files in ZIP")
            
            # Extract all files
            zip_ref.extractall(output_dir)
            print(f"  [OK] Extracted to: {output_dir}")
            
            # List extracted files
            print(f"\n  Extracted files:")
            for f in file_list[:20]:  # Show first 20
                file_path = output_dir / f
                if file_path.exists():
                    size = file_path.stat().st_size
                    print(f"    - {f} ({size:,} bytes)")
            
            return file_list
    except Exception as e:
        print(f"  [ERROR] Error extracting ZIP: {e}")
        return []


def identify_paper_info(text: str) -> dict:
    """Identify paper information from text"""
    info = {
        'title': '',
        'authors': [],
        'journal': '',
        'year': '',
        'doi': '',
        'pubmed_id': ''
    }
    
    # Extract title (usually in first 50 lines)
    lines = text.split('\n')[:50]
    for line in lines:
        line = line.strip()
        if len(line) > 30 and len(line) < 300:
            if any(keyword in line.lower() for keyword in ['aneurysm', 'intracranial', 'cerebral', 'genome', 'association']):
                if not info['title'] or len(line) > len(info['title']):
                    info['title'] = line
    
    # Extract year
    year_match = re.search(r'\b(20\d{2})\b', text[:2000])
    if year_match:
        info['year'] = year_match.group(1)
    
    # Extract DOI
    doi_match = re.search(r'doi[:\s]+([\d.]+/[\w.-]+)', text, re.IGNORECASE)
    if doi_match:
        info['doi'] = doi_match.group(1)
    
    # Extract PubMed ID
    pubmed_match = re.search(r'pubmed[:\s]+(\d+)|pmid[:\s]+(\d+)', text, re.IGNORECASE)
    if pubmed_match:
        info['pubmed_id'] = pubmed_match.group(1) or pubmed_match.group(2)
    
    # Extract journal
    if 'JPM' in text[:500] or 'Journal of Personalized Medicine' in text[:500]:
        info['journal'] = 'Journal of Personalized Medicine'
    
    return info


def identify_data_types(text: str, file_list: list) -> dict:
    """Identify what types of data are present"""
    data_types = {
        'has_gwas': False,
        'has_expression': False,
        'has_clinical': False,
        'has_supplementary_tables': False,
        'indicators': {}
    }
    
    # Check text for indicators
    data_types['indicators'] = {
        'rs_numbers': len(re.findall(r'rs\d+', text, re.IGNORECASE)),
        'p_values': len(re.findall(r'p\s*[<>=]?\s*[\d.e-]+|pvalue|p-value', text, re.IGNORECASE)),
        'gene_names': len(re.findall(r'\b[A-Z]{2,}[0-9]*\b', text[:5000])),
        'tables': len(re.findall(r'Table\s+\d+|Supplementary\s+Table', text, re.IGNORECASE)),
        'samples': len(re.findall(r'\d+\s+(cases|controls|samples|patients)', text, re.IGNORECASE))
    }
    
    data_types['has_gwas'] = (
        data_types['indicators']['rs_numbers'] > 5 or
        (data_types['indicators']['rs_numbers'] > 0 and data_types['indicators']['p_values'] > 10)
    )
    
    data_types['has_expression'] = data_types['indicators']['gene_names'] > 20
    
    data_types['has_clinical'] = data_types['indicators']['samples'] > 0
    
    # Check file list for data files
    for f in file_list:
        f_lower = f.lower()
        if any(ext in f_lower for ext in ['.tsv', '.csv', '.txt', '.xlsx', '.xls']):
            data_types['has_supplementary_tables'] = True
        if 'gwas' in f_lower or 'summary' in f_lower:
            data_types['has_gwas'] = True
        if 'expression' in f_lower or 'counts' in f_lower:
            data_types['has_expression'] = True
    
    return data_types


def analyze_extracted_files(output_dir: Path, file_list: list):
    """Analyze extracted files for useful data"""
    print(f"\n[INFO] Analyzing extracted files...")
    
    useful_files = []
    
    for f in file_list:
        file_path = output_dir / f
        
        if not file_path.exists() or file_path.is_dir():
            continue
        
        ext = file_path.suffix.lower()
        
        # Check if it's a data file
        if ext in ['.tsv', '.csv', '.txt', '.xlsx', '.xls']:
            try:
                size = file_path.stat().st_size
                print(f"\n  [DATA FILE] {f} ({size:,} bytes)")
                
                # Try to read and preview
                if ext == '.csv':
                    try:
                        df = pd.read_csv(file_path, nrows=5)
                        print(f"    Columns: {list(df.columns)}")
                        print(f"    Rows (preview): {len(df)} shown, file may have more")
                        useful_files.append({'file': f, 'type': 'csv', 'columns': list(df.columns)})
                    except:
                        pass
                elif ext == '.tsv':
                    try:
                        df = pd.read_csv(file_path, sep='\t', nrows=5)
                        print(f"    Columns: {list(df.columns)}")
                        print(f"    Rows (preview): {len(df)} shown, file may have more")
                        useful_files.append({'file': f, 'type': 'tsv', 'columns': list(df.columns)})
                    except:
                        pass
                elif ext in ['.xlsx', '.xls']:
                    try:
                        df = pd.read_excel(file_path, nrows=5)
                        print(f"    Columns: {list(df.columns)}")
                        print(f"    Rows (preview): {len(df)} shown, file may have more")
                        useful_files.append({'file': f, 'type': 'excel', 'columns': list(df.columns)})
                    except:
                        pass
            except Exception as e:
                print(f"    [WARNING] Could not analyze {f}: {e}")
    
    return useful_files


def main():
    """Main extraction function"""
    print("=" * 80)
    print("EXTRACTING DATA FROM JPM PAPER (jpm-13-00266)")
    print("=" * 80)
    
    # Find files in Downloads folder
    downloads_path = Path.home() / "Downloads"
    pdf_file = list(downloads_path.glob("jpm-13-00266.pdf"))
    zip_file = list(downloads_path.glob("jpm-13-00266-s001.zip"))
    
    if not pdf_file and not zip_file:
        print(f"\n[ERROR] JPM files not found in Downloads folder")
        print(f"   Searched: {downloads_path}")
        return
    
    output_dir = Path("data/supplementary_pdfs/extracted_gwas/jpm_13_00266")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    all_data = {}
    
    # Process PDF
    if pdf_file:
        pdf_path = pdf_file[0]
        print(f"\n[OK] Found PDF: {pdf_path.name}")
        print(f"   Size: {pdf_path.stat().st_size / (1024*1024):.2f} MB")
        
        print(f"\n[INFO] Extracting text from PDF...")
        text = extract_text_from_pdf(pdf_path)
        
        if text:
            print(f"[OK] Extracted {len(text):,} characters")
            
            # Identify paper info
            print(f"\n[INFO] Identifying paper information...")
            paper_info = identify_paper_info(text)
            print(f"  Title: {paper_info.get('title', 'Not found')[:100]}...")
            if paper_info.get('year'):
                print(f"  Year: {paper_info['year']}")
            if paper_info.get('journal'):
                print(f"  Journal: {paper_info['journal']}")
            if paper_info.get('doi'):
                print(f"  DOI: {paper_info['doi']}")
            if paper_info.get('pubmed_id'):
                print(f"  PubMed ID: {paper_info['pubmed_id']}")
            
            # Save text
            text_file = output_dir / "jpm_13_00266_text.txt"
            with open(text_file, 'w', encoding='utf-8') as f:
                f.write(text)
            print(f"\n[OK] Saved text to: {text_file}")
            
            all_data['pdf_text'] = text
            all_data['paper_info'] = paper_info
        else:
            print("[WARNING] Could not extract text from PDF")
    
    # Process ZIP file
    if zip_file:
        zip_path = zip_file[0]
        print(f"\n[OK] Found ZIP: {zip_path.name}")
        print(f"   Size: {zip_path.stat().st_size / (1024*1024):.2f} MB")
        
        zip_extract_dir = output_dir / "extracted"
        zip_extract_dir.mkdir(exist_ok=True)
        
        file_list = extract_zip_contents(zip_path, zip_extract_dir)
        
        if file_list:
            # Analyze extracted files
            useful_files = analyze_extracted_files(zip_extract_dir, file_list)
            all_data['zip_files'] = file_list
            all_data['useful_files'] = useful_files
    
    # Identify data types
    if 'pdf_text' in all_data:
        print(f"\n[INFO] Identifying data types...")
        file_list = all_data.get('zip_files', [])
        data_types = identify_data_types(all_data['pdf_text'], file_list)
        
        print(f"  GWAS data: {'Yes' if data_types['has_gwas'] else 'No'}")
        print(f"  Expression data: {'Yes' if data_types['has_expression'] else 'No'}")
        print(f"  Clinical data: {'Yes' if data_types['has_clinical'] else 'No'}")
        print(f"  Supplementary tables: {'Yes' if data_types['has_supplementary_tables'] else 'No'}")
        
        if data_types['indicators']:
            print(f"\n  Indicators:")
            print(f"    - RS numbers: {data_types['indicators']['rs_numbers']}")
            print(f"    - P-values: {data_types['indicators']['p_values']}")
            print(f"    - Gene names: {data_types['indicators']['gene_names']}")
            print(f"    - Tables: {data_types['indicators']['tables']}")
            print(f"    - Sample mentions: {data_types['indicators']['samples']}")
        
        all_data['data_types'] = data_types
    
    # Save metadata
    metadata = {
        'extraction_date': datetime.now().isoformat(),
        'paper_info': all_data.get('paper_info', {}),
        'data_types': all_data.get('data_types', {}),
        'useful_files': all_data.get('useful_files', []),
        'zip_files': all_data.get('zip_files', [])
    }
    
    metadata_file = output_dir / "jpm_13_00266_metadata.json"
    with open(metadata_file, 'w', encoding='utf-8') as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)
    print(f"\n[OK] Saved metadata to: {metadata_file}")
    
    print("\n" + "=" * 80)
    print("EXTRACTION COMPLETE")
    print("=" * 80)
    print(f"\n[OK] All extracted data saved to: {output_dir}")


if __name__ == "__main__":
    main()






