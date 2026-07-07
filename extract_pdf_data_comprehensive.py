"""
Comprehensive PDF Data Extraction
Extracts GWAS summary statistics, tables, and other data from supplementary PDFs
"""
import sys
from pathlib import Path
import json
import pandas as pd
import re
from typing import Dict, List, Optional

# Try importing PDF libraries
PDF2_AVAILABLE = False
PDFPLUMBER_AVAILABLE = False
TABULA_AVAILABLE = False
CAMELOT_AVAILABLE = False
PYMUPDF_AVAILABLE = False

try:
    import PyPDF2
    PDF2_AVAILABLE = True
    # Add user site-packages to path if needed
    import site
    import os
    user_site = site.getusersitepackages()
    if user_site and user_site not in sys.path:
        sys.path.insert(0, user_site)
except ImportError:
    PDF2_AVAILABLE = False

try:
    import pdfplumber
    PDFPLUMBER_AVAILABLE = True
except ImportError:
    pass

try:
    import tabula
    TABULA_AVAILABLE = True
except ImportError:
    pass

try:
    import camelot
    CAMELOT_AVAILABLE = True
except ImportError:
    pass

try:
    import fitz  # pymupdf
    PYMUPDF_AVAILABLE = True
except ImportError:
    pass

class PDFDataExtractor:
    """Extract data from PDFs"""
    
    def __init__(self, pdf_path: Path):
        self.pdf_path = pdf_path
        self.output_dir = pdf_path.parent / "extracted"
        self.output_dir.mkdir(exist_ok=True)
        self.extracted_data = {}
        self.text_content = ""
        self.tables = []
        
    def extract_text_pypdf2(self) -> str:
        """Extract text using PyPDF2"""
        if not PDF2_AVAILABLE:
            return ""
        
        try:
            text = ""
            with open(self.pdf_path, 'rb') as pdf_file:
                reader = PyPDF2.PdfReader(pdf_file)
                for page in reader.pages:
                    text += page.extract_text() + "\n"
            return text
        except Exception as e:
            print(f"  ⚠️ PyPDF2 error: {e}")
            return ""
    
    def extract_text_pdfplumber(self) -> str:
        """Extract text using pdfplumber"""
        if not PDFPLUMBER_AVAILABLE:
            return ""
        
        try:
            text = ""
            with pdfplumber.open(self.pdf_path) as pdf:
                for page in pdf.pages:
                    page_text = page.extract_text()
                    if page_text:
                        text += page_text + "\n"
            return text
        except Exception as e:
            print(f"  ⚠️ pdfplumber error: {e}")
            return ""
    
    def extract_text_pymupdf(self) -> str:
        """Extract text using pymupdf (fitz)"""
        if not PYMUPDF_AVAILABLE:
            return ""
        
        try:
            text = ""
            doc = fitz.open(self.pdf_path)
            for page in doc:
                text += page.get_text() + "\n"
            doc.close()
            return text
        except Exception as e:
            print(f"  ⚠️ pymupdf error: {e}")
            return ""
    
    def extract_all_text(self) -> str:
        """Extract text using all available methods"""
        print(f"  Extracting text...")
        
        # Try pdfplumber first (best quality)
        if PDFPLUMBER_AVAILABLE:
            text = self.extract_text_pdfplumber()
            if text and len(text) > 100:
                return text
        
        # Try pymupdf
        if PYMUPDF_AVAILABLE:
            text = self.extract_text_pymupdf()
            if text and len(text) > 100:
                return text
        
        # Try PyPDF2
        if PDF2_AVAILABLE:
            text = self.extract_text_pypdf2()
            if text and len(text) > 100:
                return text
        
        return ""
    
    def extract_tables_tabula(self) -> List[pd.DataFrame]:
        """Extract tables using tabula"""
        if not TABULA_AVAILABLE:
            return []
        
        try:
            tables = tabula.read_pdf(str(self.pdf_path), pages='all', multiple_tables=True)
            return [df for df in tables if df is not None and not df.empty]
        except Exception as e:
            print(f"  ⚠️ tabula error: {e}")
            return []
    
    def extract_tables_camelot(self) -> List[pd.DataFrame]:
        """Extract tables using camelot"""
        if not CAMELOT_AVAILABLE:
            return []
        
        try:
            tables = camelot.read_pdf(str(self.pdf_path), pages='all')
            dfs = []
            for table in tables:
                if table.df is not None and not table.df.empty:
                    dfs.append(table.df)
            return dfs
        except Exception as e:
            print(f"  ⚠️ camelot error: {e}")
            return []
    
    def extract_tables_pdfplumber(self) -> List[pd.DataFrame]:
        """Extract tables using pdfplumber"""
        if not PDFPLUMBER_AVAILABLE:
            return []
        
        try:
            tables = []
            with pdfplumber.open(self.pdf_path) as pdf:
                for page_num, page in enumerate(pdf.pages):
                    page_tables = page.extract_tables()
                    for table in page_tables:
                        if table and len(table) > 1:
                            df = pd.DataFrame(table[1:], columns=table[0])
                            tables.append(df)
            return tables
        except Exception as e:
            print(f"  ⚠️ pdfplumber table extraction error: {e}")
            return []
    
    def extract_all_tables(self) -> List[pd.DataFrame]:
        """Extract tables using all available methods"""
        print(f"  Extracting tables...")
        
        all_tables = []
        
        # Try pdfplumber first (good for tables)
        if PDFPLUMBER_AVAILABLE:
            tables = self.extract_tables_pdfplumber()
            all_tables.extend(tables)
            if tables:
                print(f"    Found {len(tables)} tables with pdfplumber")
        
        # Try tabula
        if TABULA_AVAILABLE:
            tables = self.extract_tables_tabula()
            all_tables.extend(tables)
            if tables:
                print(f"    Found {len(tables)} tables with tabula")
        
        # Try camelot
        if CAMELOT_AVAILABLE:
            tables = self.extract_tables_camelot()
            all_tables.extend(tables)
            if tables:
                print(f"    Found {len(tables)} tables with camelot")
        
        return all_tables
    
    def identify_gwas_data(self, text: str) -> Dict:
        """Identify GWAS summary statistics in text"""
        gwas_indicators = {
            'rs_numbers': len(re.findall(r'rs\d+', text, re.IGNORECASE)),
            'chromosome': len(re.findall(r'chr\s*\d+|chromosome\s*\d+', text, re.IGNORECASE)),
            'p_value': len(re.findall(r'p\s*[<>=]?\s*[\d.e-]+', text, re.IGNORECASE)),
            'odds_ratio': len(re.findall(r'or\s*[=:]\s*[\d.]+|odds\s*ratio', text, re.IGNORECASE)),
            'beta': len(re.findall(r'beta\s*[=:]\s*[\d.-]+', text, re.IGNORECASE)),
            'allele': len(re.findall(r'allele|a[ctg]|g[act]|t[acg]|c[atg]', text, re.IGNORECASE)),
            'snp': len(re.findall(r'snp|single\s*nucleotide', text, re.IGNORECASE)),
        }
        
        # Check if this looks like GWAS data
        gwas_score = sum(gwas_indicators.values())
        is_gwas = (
            gwas_indicators['rs_numbers'] > 5 or
            gwas_indicators['p_value'] > 10 or
            (gwas_indicators['rs_numbers'] > 0 and gwas_indicators['p_value'] > 0)
        )
        
        return {
            'is_gwas': is_gwas,
            'gwas_score': gwas_score,
            'indicators': gwas_indicators
        }
    
    def identify_paper_info(self, text: str) -> Dict:
        """Identify paper information from text"""
        paper_info = {
            'title': '',
            'authors': [],
            'journal': '',
            'year': '',
            'doi': '',
            'pubmed_id': ''
        }
        
        # Extract title (usually in first few lines)
        lines = text.split('\n')[:20]
        for line in lines:
            if len(line) > 20 and len(line) < 200:
                if not paper_info['title']:
                    paper_info['title'] = line.strip()
        
        # Extract year
        year_match = re.search(r'\b(19|20)\d{2}\b', text)
        if year_match:
            paper_info['year'] = year_match.group()
        
        # Extract DOI
        doi_match = re.search(r'doi[:\s]+([\d.]+/[\w.-]+)', text, re.IGNORECASE)
        if doi_match:
            paper_info['doi'] = doi_match.group(1)
        
        # Extract PubMed ID
        pubmed_match = re.search(r'pubmed[:\s]+(\d+)', text, re.IGNORECASE)
        if pubmed_match:
            paper_info['pubmed_id'] = pubmed_match.group(1)
        
        # Extract journal
        journal_keywords = ['Nature', 'Science', 'Stroke', 'Journal', 'Genetics']
        for keyword in journal_keywords:
            if keyword.lower() in text.lower():
                paper_info['journal'] = keyword
                break
        
        return paper_info
    
    def save_extracted_data(self):
        """Save all extracted data"""
        pdf_name = self.pdf_path.stem
        
        # Save text
        if self.text_content:
            text_file = self.output_dir / f"{pdf_name}_extracted_text.txt"
            with open(text_file, 'w', encoding='utf-8') as f:
                f.write(self.text_content)
            print(f"  ✅ Saved text: {text_file}")
        
        # Save tables
        if self.tables:
            for i, table in enumerate(self.tables):
                if table is not None and not table.empty:
                    # Save as TSV
                    tsv_file = self.output_dir / f"{pdf_name}_table_{i+1}.tsv"
                    table.to_csv(tsv_file, sep='\t', index=False, encoding='utf-8')
                    print(f"  ✅ Saved table {i+1}: {tsv_file} ({len(table)} rows, {len(table.columns)} columns)")
                    
                    # Also save as CSV
                    csv_file = self.output_dir / f"{pdf_name}_table_{i+1}.csv"
                    table.to_csv(csv_file, index=False, encoding='utf-8')
        
        # Save metadata
        metadata = {
            'pdf_name': pdf_name,
            'text_length': len(self.text_content),
            'num_tables': len(self.tables),
            'gwas_analysis': self.identify_gwas_data(self.text_content) if self.text_content else {},
            'paper_info': self.identify_paper_info(self.text_content) if self.text_content else {}
        }
        
        metadata_file = self.output_dir / f"{pdf_name}_metadata.json"
        with open(metadata_file, 'w', encoding='utf-8') as f:
            json.dump(metadata, f, indent=2, ensure_ascii=False)
        print(f"  ✅ Saved metadata: {metadata_file}")
    
    def extract_all(self):
        """Extract all data from PDF"""
        print(f"\n📄 Processing: {self.pdf_path.name}")
        print("-" * 80)
        
        # Extract text
        self.text_content = self.extract_all_text()
        if self.text_content:
            print(f"  ✅ Extracted {len(self.text_content)} characters of text")
        else:
            print(f"  ⚠️ Could not extract text")
        
        # Extract tables
        self.tables = self.extract_all_tables()
        if self.tables:
            print(f"  ✅ Found {len(self.tables)} tables")
        else:
            print(f"  ⚠️ No tables found")
        
        # Analyze content
        if self.text_content:
            gwas_analysis = self.identify_gwas_data(self.text_content)
            if gwas_analysis['is_gwas']:
                print(f"  ✅ GWAS data detected (score: {gwas_analysis['gwas_score']})")
                print(f"    - RS numbers: {gwas_analysis['indicators']['rs_numbers']}")
                print(f"    - P-values: {gwas_analysis['indicators']['p_value']}")
                print(f"    - Odds ratios: {gwas_analysis['indicators']['odds_ratio']}")
            else:
                print(f"  ⚠️ GWAS data not clearly detected")
            
            paper_info = self.identify_paper_info(self.text_content)
            if paper_info['title']:
                print(f"  📄 Title: {paper_info['title'][:100]}...")
            if paper_info['year']:
                print(f"  📅 Year: {paper_info['year']}")
            if paper_info['journal']:
                print(f"  📚 Journal: {paper_info['journal']}")
            if paper_info['pubmed_id']:
                print(f"  🔗 PubMed ID: {paper_info['pubmed_id']}")
        
        # Save extracted data
        self.save_extracted_data()
        
        return {
            'text': self.text_content,
            'tables': self.tables,
            'gwas_analysis': self.identify_gwas_data(self.text_content) if self.text_content else {},
            'paper_info': self.identify_paper_info(self.text_content) if self.text_content else {}
        }

def main():
    """Main extraction function"""
    print("=" * 80)
    print("COMPREHENSIVE PDF DATA EXTRACTION")
    print("=" * 80)
    
    # Check available libraries
    print("\n📚 Available PDF libraries:")
    print(f"  PyPDF2: {'✅' if PDF2_AVAILABLE else '❌'}")
    print(f"  pdfplumber: {'✅' if PDFPLUMBER_AVAILABLE else '❌'}")
    print(f"  tabula: {'✅' if TABULA_AVAILABLE else '❌'}")
    print(f"  camelot: {'✅' if CAMELOT_AVAILABLE else '❌'}")
    print(f"  pymupdf: {'✅' if PYMUPDF_AVAILABLE else '❌'}")
    
    if not any([PDF2_AVAILABLE, PDFPLUMBER_AVAILABLE, PYMUPDF_AVAILABLE]):
        print("\n❌ No PDF extraction libraries available!")
        print("   Please install: pip install PyPDF2 pdfplumber pymupdf")
        return
    
    # Find PDFs
    pdf_dir = Path("data/supplementary_pdfs")
    pdfs = list(pdf_dir.glob("*.pdf"))
    
    if not pdfs:
        print(f"\n❌ No PDFs found in {pdf_dir}")
        return
    
    print(f"\n📁 Found {len(pdfs)} PDFs to process")
    
    # Extract data from each PDF
    all_extracted = {}
    for pdf_path in pdfs:
        extractor = PDFDataExtractor(pdf_path)
        extracted = extractor.extract_all()
        all_extracted[pdf_path.name] = extracted
    
    # Generate summary report
    print("\n" + "=" * 80)
    print("EXTRACTION SUMMARY")
    print("=" * 80)
    
    for pdf_name, data in all_extracted.items():
        print(f"\n📄 {pdf_name}")
        print(f"  Text length: {len(data.get('text', ''))} characters")
        print(f"  Tables found: {len(data.get('tables', []))}")
        
        if data.get('gwas_analysis', {}).get('is_gwas'):
            print(f"  ✅ GWAS data detected")
        
        paper_info = data.get('paper_info', {})
        if paper_info.get('title'):
            print(f"  Title: {paper_info['title'][:80]}...")
        if paper_info.get('year'):
            print(f"  Year: {paper_info['year']}")
        if paper_info.get('pubmed_id'):
            print(f"  PubMed ID: {paper_info['pubmed_id']}")
    
    print("\n" + "=" * 80)
    print("EXTRACTION COMPLETE!")
    print("=" * 80)
    print(f"\n✅ All extracted data saved to: data/supplementary_pdfs/extracted/")
    print(f"✅ Tables saved as TSV and CSV files")
    print(f"✅ Metadata saved as JSON files")

if __name__ == "__main__":
    main()

