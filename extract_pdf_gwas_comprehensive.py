"""
Comprehensive PDF GWAS Data Extraction and Verification
Extracts GWAS summary statistics and supplementary tables from PDFs,
converts to standardized format, and verifies content.
"""
import sys
import os
from pathlib import Path
import json
import pandas as pd
import numpy as np
import re
from typing import Dict, List, Optional, Tuple
from datetime import datetime

# Add user site-packages to path for PyPDF2
import site
user_site = site.getusersitepackages()
if user_site and user_site not in sys.path:
    sys.path.insert(0, user_site)

# Try importing PDF libraries
PDF2_AVAILABLE = False
PDFPLUMBER_AVAILABLE = False
TABULA_AVAILABLE = False
CAMELOT_AVAILABLE = False
PYMUPDF_AVAILABLE = False

try:
    import PyPDF2
    PDF2_AVAILABLE = True
except ImportError:
    pass

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


class PDFGWASExtractor:
    """Extract GWAS summary statistics from PDFs"""
    
    def __init__(self, pdf_path: Path):
        self.pdf_path = pdf_path
        self.output_dir = pdf_path.parent / "extracted_gwas"
        self.output_dir.mkdir(exist_ok=True)
        self.text_content = ""
        self.tables = []
        self.gwas_tables = []
        self.paper_info = {}
        
    def extract_text(self) -> str:
        """Extract text using available methods"""
        print(f"  Extracting text from PDF...")
        
        # Try pdfplumber first (best quality)
        if PDFPLUMBER_AVAILABLE:
            try:
                text = ""
                with pdfplumber.open(self.pdf_path) as pdf:
                    for page in pdf.pages:
                        page_text = page.extract_text()
                        if page_text:
                            text += page_text + "\n"
                if text and len(text) > 100:
                    print(f"    ✓ Extracted {len(text):,} characters using pdfplumber")
                    return text
            except Exception as e:
                print(f"    ⚠ pdfplumber error: {e}")
        
        # Try pymupdf
        if PYMUPDF_AVAILABLE:
            try:
                text = ""
                doc = fitz.open(self.pdf_path)
                for page in doc:
                    text += page.get_text() + "\n"
                doc.close()
                if text and len(text) > 100:
                    print(f"    ✓ Extracted {len(text):,} characters using pymupdf")
                    return text
            except Exception as e:
                print(f"    ⚠ pymupdf error: {e}")
        
        # Try PyPDF2
        if PDF2_AVAILABLE:
            try:
                text = ""
                with open(self.pdf_path, 'rb') as pdf_file:
                    reader = PyPDF2.PdfReader(pdf_file)
                    for page in reader.pages:
                        text += page.extract_text() + "\n"
                if text and len(text) > 100:
                    print(f"    ✓ Extracted {len(text):,} characters using PyPDF2")
                    return text
            except Exception as e:
                print(f"    ⚠ PyPDF2 error: {e}")
        
        return ""
    
    def extract_tables(self) -> List[pd.DataFrame]:
        """Extract tables using all available methods"""
        print(f"  Extracting tables from PDF...")
        all_tables = []
        
        # Try pdfplumber first (good for tables)
        if PDFPLUMBER_AVAILABLE:
            try:
                tables = []
                with pdfplumber.open(self.pdf_path) as pdf:
                    for page_num, page in enumerate(pdf.pages):
                        page_tables = page.extract_tables()
                        for table in page_tables:
                            if table and len(table) > 1:
                                # First row as header
                                df = pd.DataFrame(table[1:], columns=table[0])
                                # Clean column names
                                df.columns = [str(col).strip() if col else f"col_{i}" 
                                            for i, col in enumerate(df.columns)]
                                tables.append(df)
                if tables:
                    print(f"    ✓ Found {len(tables)} tables with pdfplumber")
                    all_tables.extend(tables)
            except Exception as e:
                print(f"    ⚠ pdfplumber table extraction error: {e}")
        
        # Try tabula
        if TABULA_AVAILABLE:
            try:
                tables = tabula.read_pdf(str(self.pdf_path), pages='all', 
                                        multiple_tables=True, silent=True)
                valid_tables = [df for df in tables if df is not None and not df.empty]
                if valid_tables:
                    print(f"    ✓ Found {len(valid_tables)} tables with tabula")
                    all_tables.extend(valid_tables)
            except Exception as e:
                print(f"    ⚠ tabula error: {e}")
        
        # Try camelot
        if CAMELOT_AVAILABLE:
            try:
                tables = camelot.read_pdf(str(self.pdf_path), pages='all', flavor='lattice')
                valid_tables = [table.df for table in tables if table.df is not None and not table.df.empty]
                if valid_tables:
                    print(f"    ✓ Found {len(valid_tables)} tables with camelot")
                    all_tables.extend(valid_tables)
            except Exception as e:
                print(f"    ⚠ camelot error: {e}")
        
        return all_tables
    
    def identify_gwas_table(self, table: pd.DataFrame) -> bool:
        """Check if a table contains GWAS summary statistics"""
        if table is None or table.empty:
            return False
        
        # Convert all values to strings for pattern matching
        table_str = table.astype(str).to_string().lower()
        
        # Check for GWAS indicators
        has_rs_numbers = bool(re.search(r'rs\d+', table_str, re.IGNORECASE))
        has_p_values = bool(re.search(r'p\s*[<>=]?\s*[\d.e-]+|pvalue|p-value', table_str, re.IGNORECASE))
        has_odds_ratio = bool(re.search(r'or\s*[=:]\s*[\d.]+|odds\s*ratio', table_str, re.IGNORECASE))
        has_beta = bool(re.search(r'beta\s*[=:]\s*[\d.-]+', table_str, re.IGNORECASE))
        has_chromosome = bool(re.search(r'chr\s*\d+|chromosome\s*\d+', table_str, re.IGNORECASE))
        has_allele = bool(re.search(r'allele|effect\s*allele|a[ctg]|g[act]|t[acg]|c[atg]', table_str, re.IGNORECASE))
        
        # Check column names
        col_names = ' '.join([str(col).lower() for col in table.columns])
        has_snp_col = bool(re.search(r'snp|rsid|variant|marker', col_names))
        has_p_col = bool(re.search(r'p|pvalue|p-value|pval', col_names))
        has_or_col = bool(re.search(r'or|odds|ratio', col_names))
        has_beta_col = bool(re.search(r'beta|effect', col_names))
        has_chr_col = bool(re.search(r'chr|chromosome|chrom', col_names))
        
        # Score the table
        score = sum([
            has_rs_numbers, has_p_values, has_odds_ratio or has_beta,
            has_chromosome, has_allele, has_snp_col, has_p_col,
            has_or_col or has_beta_col, has_chr_col
        ])
        
        # Consider it GWAS if it has multiple indicators
        is_gwas = score >= 3 or (has_rs_numbers and has_p_values)
        
        return is_gwas
    
    def standardize_gwas_table(self, table: pd.DataFrame, table_num: int) -> Optional[pd.DataFrame]:
        """Convert table to standardized GWAS format"""
        if table is None or table.empty:
            return None
        
        # Standard column mappings
        standard_cols = {
            'snp': ['snp', 'rsid', 'rs', 'variant', 'marker', 'snp_id', 'rs_number'],
            'chr': ['chr', 'chromosome', 'chrom', 'ch'],
            'pos': ['pos', 'position', 'bp', 'base_pair', 'coordinate'],
            'ea': ['ea', 'effect_allele', 'a1', 'allele1', 'alt', 'alt_allele', 'effect'],
            'oa': ['oa', 'other_allele', 'a2', 'allele2', 'ref', 'ref_allele', 'non_effect'],
            'or': ['or', 'odds_ratio', 'odds', 'or_95ci'],
            'beta': ['beta', 'effect', 'effect_size', 'b'],
            'se': ['se', 'standard_error', 'std_err'],
            'p': ['p', 'pvalue', 'p-value', 'p_val', 'pval', 'p_value'],
            'maf': ['maf', 'minor_allele_freq', 'freq', 'frequency'],
            'n': ['n', 'sample_size', 'n_total', 'cases', 'controls']
        }
        
        # Find matching columns
        col_mapping = {}
        table_cols_lower = {str(col).lower().strip(): col for col in table.columns}
        
        for std_col, variants in standard_cols.items():
            for variant in variants:
                if variant in table_cols_lower:
                    col_mapping[std_col] = table_cols_lower[variant]
                    break
        
        # Create standardized dataframe
        standardized = pd.DataFrame()
        
        # Map columns
        for std_col, orig_col in col_mapping.items():
            standardized[std_col] = table[orig_col]
        
        # Try to extract SNP IDs from any column
        if 'snp' not in standardized.columns:
            for col in table.columns:
                col_str = table[col].astype(str)
                if col_str.str.contains(r'rs\d+', regex=True, na=False).any():
                    standardized['snp'] = col_str
                    # Extract rs numbers
                    standardized['snp'] = standardized['snp'].str.extract(r'(rs\d+)', expand=False)
                    break
        
        # Try to extract chromosome from any column
        if 'chr' not in standardized.columns:
            for col in table.columns:
                col_str = table[col].astype(str)
                if col_str.str.contains(r'chr|chromosome', regex=True, case=False, na=False).any():
                    standardized['chr'] = col_str.str.extract(r'(\d+)', expand=False)
                    break
        
        # Try to extract p-values
        if 'p' not in standardized.columns:
            for col in table.columns:
                col_str = table[col].astype(str)
                # Look for scientific notation or decimal numbers
                if col_str.str.contains(r'[\d.e-]+', regex=True, na=False).any():
                    # Check if values look like p-values (small numbers)
                    try:
                        numeric_vals = pd.to_numeric(col_str, errors='coerce')
                        if numeric_vals.notna().any():
                            if (numeric_vals < 1).any() and (numeric_vals > 0).any():
                                standardized['p'] = numeric_vals
                                break
                    except:
                        pass
        
        # Clean and convert data types
        for col in ['chr', 'pos', 'n']:
            if col in standardized.columns:
                standardized[col] = pd.to_numeric(standardized[col], errors='coerce')
        
        for col in ['or', 'beta', 'se', 'p', 'maf']:
            if col in standardized.columns:
                standardized[col] = pd.to_numeric(standardized[col], errors='coerce')
        
        # Add metadata
        standardized['source_pdf'] = self.pdf_path.name
        standardized['table_number'] = table_num
        standardized['extraction_date'] = datetime.now().strftime('%Y-%m-%d')
        
        return standardized if not standardized.empty else None
    
    def identify_paper_info(self, text: str) -> Dict:
        """Identify paper information from text"""
        paper_info = {
            'title': '',
            'authors': [],
            'journal': '',
            'year': '',
            'doi': '',
            'pubmed_id': '',
            'study_name': ''
        }
        
        # Extract title (usually in first few lines)
        lines = text.split('\n')[:30]
        for i, line in enumerate(lines):
            line = line.strip()
            if len(line) > 20 and len(line) < 200 and not line.startswith('SUPPLEMENT'):
                if not paper_info['title'] or len(line) > len(paper_info['title']):
                    paper_info['title'] = line
        
        # Extract year
        year_match = re.search(r'\b(19|20)\d{2}\b', text)
        if year_match:
            paper_info['year'] = year_match.group()
        
        # Extract DOI
        doi_match = re.search(r'doi[:\s]+([\d.]+/[\w.-]+)', text, re.IGNORECASE)
        if doi_match:
            paper_info['doi'] = doi_match.group(1)
        
        # Extract PubMed ID
        pubmed_match = re.search(r'pubmed[:\s]+(\d+)|pmid[:\s]+(\d+)', text, re.IGNORECASE)
        if pubmed_match:
            paper_info['pubmed_id'] = pubmed_match.group(1) or pubmed_match.group(2)
        
        # Extract journal
        journal_keywords = ['Nature', 'Science', 'Stroke', 'Journal', 'Genetics', 
                          'Neurology', 'American Journal']
        for keyword in journal_keywords:
            if keyword.lower() in text.lower():
                paper_info['journal'] = keyword
                break
        
        # Identify study name from filename or content
        pdf_name = self.pdf_path.stem
        if 'nihms' in pdf_name.lower():
            paper_info['study_name'] = 'Bilguvar_2008_NatureGenetics'
            if not paper_info['pubmed_id']:
                paper_info['pubmed_id'] = '18758463'
        elif 'stroke' in pdf_name.lower():
            paper_info['study_name'] = 'Stroke_2014'
        
        return paper_info
    
    def verify_content(self) -> Dict:
        """Verify extracted content and check for overlaps"""
        verification = {
            'pdf_name': self.pdf_path.name,
            'has_gwas_data': len(self.gwas_tables) > 0,
            'num_gwas_tables': len(self.gwas_tables),
            'num_total_tables': len(self.tables),
            'paper_info': self.paper_info,
            'gwas_indicators': {},
            'overlap_check': {}
        }
        
        if self.text_content:
            # Count GWAS indicators in text
            verification['gwas_indicators'] = {
                'rs_numbers': len(re.findall(r'rs\d+', self.text_content, re.IGNORECASE)),
                'p_values': len(re.findall(r'p\s*[<>=]?\s*[\d.e-]+', self.text_content, re.IGNORECASE)),
                'odds_ratios': len(re.findall(r'or\s*[=:]\s*[\d.]+|odds\s*ratio', self.text_content, re.IGNORECASE)),
                'chromosomes': len(re.findall(r'chr\s*\d+|chromosome\s*\d+', self.text_content, re.IGNORECASE))
            }
        
        # Check for overlaps with existing datasets
        # This would check against known GWAS datasets
        verification['overlap_check'] = {
            'status': 'pending_manual_review',
            'note': 'Compare extracted SNPs with existing GWAS datasets'
        }
        
        return verification
    
    def save_extracted_data(self):
        """Save all extracted data in standardized format"""
        pdf_name = self.pdf_path.stem
        
        # Save text
        if self.text_content:
            text_file = self.output_dir / f"{pdf_name}_text.txt"
            with open(text_file, 'w', encoding='utf-8') as f:
                f.write(self.text_content)
            print(f"  ✓ Saved text: {text_file}")
        
        # Save all tables
        if self.tables:
            for i, table in enumerate(self.tables):
                if table is not None and not table.empty:
                    tsv_file = self.output_dir / f"{pdf_name}_table_{i+1:02d}.tsv"
                    csv_file = self.output_dir / f"{pdf_name}_table_{i+1:02d}.csv"
                    
                    table.to_csv(tsv_file, sep='\t', index=False, encoding='utf-8')
                    table.to_csv(csv_file, index=False, encoding='utf-8')
                    print(f"  ✓ Saved table {i+1}: {tsv_file} ({len(table)} rows)")
        
        # Save GWAS tables in standardized format
        if self.gwas_tables:
            for i, gwas_table in enumerate(self.gwas_tables):
                if gwas_table is not None and not gwas_table.empty:
                    tsv_file = self.output_dir / f"{pdf_name}_gwas_{i+1:02d}.tsv"
                    csv_file = self.output_dir / f"{pdf_name}_gwas_{i+1:02d}.csv"
                    
                    gwas_table.to_csv(tsv_file, sep='\t', index=False, encoding='utf-8')
                    gwas_table.to_csv(csv_file, index=False, encoding='utf-8')
                    print(f"  ✓ Saved GWAS table {i+1}: {tsv_file} ({len(gwas_table)} SNPs)")
        
        # Save metadata
        metadata = {
            'pdf_name': pdf_name,
            'extraction_date': datetime.now().isoformat(),
            'text_length': len(self.text_content),
            'num_tables': len(self.tables),
            'num_gwas_tables': len(self.gwas_tables),
            'paper_info': self.paper_info,
            'verification': self.verify_content()
        }
        
        metadata_file = self.output_dir / f"{pdf_name}_metadata.json"
        with open(metadata_file, 'w', encoding='utf-8') as f:
            json.dump(metadata, f, indent=2, ensure_ascii=False)
        print(f"  ✓ Saved metadata: {metadata_file}")
    
    def extract_all(self) -> Dict:
        """Extract all data from PDF"""
        print(f"\n{'='*80}")
        print(f"Processing: {self.pdf_path.name}")
        print(f"{'='*80}")
        
        # Extract text
        self.text_content = self.extract_text()
        if not self.text_content:
            print("  ⚠ Could not extract text")
            return {}
        
        # Identify paper info
        self.paper_info = self.identify_paper_info(self.text_content)
        print(f"\n  Paper Information:")
        if self.paper_info.get('title'):
            print(f"    Title: {self.paper_info['title'][:100]}")
        if self.paper_info.get('year'):
            print(f"    Year: {self.paper_info['year']}")
        if self.paper_info.get('journal'):
            print(f"    Journal: {self.paper_info['journal']}")
        if self.paper_info.get('pubmed_id'):
            print(f"    PubMed ID: {self.paper_info['pubmed_id']}")
        if self.paper_info.get('study_name'):
            print(f"    Study: {self.paper_info['study_name']}")
        
        # Extract tables
        self.tables = self.extract_tables()
        print(f"\n  Found {len(self.tables)} total tables")
        
        # Identify GWAS tables
        print(f"\n  Identifying GWAS summary statistics tables...")
        for i, table in enumerate(self.tables):
            if self.identify_gwas_table(table):
                print(f"    ✓ Table {i+1} appears to contain GWAS data")
                standardized = self.standardize_gwas_table(table, i+1)
                if standardized is not None:
                    self.gwas_tables.append(standardized)
                    print(f"      Standardized: {len(standardized)} rows, {len(standardized.columns)} columns")
        
        print(f"\n  ✓ Identified {len(self.gwas_tables)} GWAS tables")
        
        # Save extracted data
        self.save_extracted_data()
        
        # Verification
        verification = self.verify_content()
        
        return {
            'text': self.text_content,
            'tables': self.tables,
            'gwas_tables': self.gwas_tables,
            'paper_info': self.paper_info,
            'verification': verification
        }


def main():
    """Main extraction function"""
    print("=" * 80)
    print("COMPREHENSIVE PDF GWAS DATA EXTRACTION")
    print("=" * 80)
    
    # Check available libraries
    print("\n📚 Available PDF libraries:")
    print(f"  PyPDF2: {'✓' if PDF2_AVAILABLE else '✗'}")
    print(f"  pdfplumber: {'✓' if PDFPLUMBER_AVAILABLE else '✗'}")
    print(f"  tabula: {'✓' if TABULA_AVAILABLE else '✗'}")
    print(f"  camelot: {'✓' if CAMELOT_AVAILABLE else '✗'}")
    print(f"  pymupdf: {'✓' if PYMUPDF_AVAILABLE else '✗'}")
    
    if not any([PDF2_AVAILABLE, PDFPLUMBER_AVAILABLE, PYMUPDF_AVAILABLE]):
        print("\n❌ No PDF extraction libraries available!")
        print("   Please install: pip install --user PyPDF2 pdfplumber pymupdf")
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
        extractor = PDFGWASExtractor(pdf_path)
        extracted = extractor.extract_all()
        all_extracted[pdf_path.name] = extracted
    
    # Generate summary report
    print("\n" + "=" * 80)
    print("EXTRACTION SUMMARY")
    print("=" * 80)
    
    total_gwas_tables = 0
    for pdf_name, data in all_extracted.items():
        print(f"\n📄 {pdf_name}")
        print(f"  Text length: {len(data.get('text', '')):,} characters")
        print(f"  Total tables: {len(data.get('tables', []))}")
        print(f"  GWAS tables: {len(data.get('gwas_tables', []))}")
        total_gwas_tables += len(data.get('gwas_tables', []))
        
        paper_info = data.get('paper_info', {})
        if paper_info.get('title'):
            print(f"  Title: {paper_info['title'][:80]}...")
        if paper_info.get('study_name'):
            print(f"  Study: {paper_info['study_name']}")
        if paper_info.get('pubmed_id'):
            print(f"  PubMed ID: {paper_info['pubmed_id']}")
        
        verification = data.get('verification', {})
        if verification.get('gwas_indicators'):
            indicators = verification['gwas_indicators']
            print(f"  GWAS indicators:")
            print(f"    - RS numbers: {indicators.get('rs_numbers', 0)}")
            print(f"    - P-values: {indicators.get('p_values', 0)}")
            print(f"    - Chromosomes: {indicators.get('chromosomes', 0)}")
    
    print("\n" + "=" * 80)
    print("EXTRACTION COMPLETE!")
    print("=" * 80)
    print(f"\n✓ Processed {len(pdfs)} PDFs")
    print(f"✓ Extracted {total_gwas_tables} GWAS summary statistics tables")
    print(f"\n📁 All extracted data saved to: data/supplementary_pdfs/extracted_gwas/")
    print(f"   - Tables saved as TSV and CSV files")
    print(f"   - GWAS data in standardized format")
    print(f"   - Metadata saved as JSON files")
    print(f"\n📋 Next steps:")
    print(f"   1. Review extracted GWAS tables")
    print(f"   2. Compare with existing GWAS datasets")
    print(f"   3. Integrate into analysis pipeline")


if __name__ == "__main__":
    main()






