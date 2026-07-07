"""
Extract GWAS data from Foroud et al. 2014 main paper PDF
This is the main paper (not the supplement) which may contain GWAS summary statistics
"""
import sys
import os
from pathlib import Path
import json
import re
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
    print("⚠️ PyPDF2 not available")


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
        print(f"  ⚠️ Error extracting text: {e}")
        return ""


def identify_gwas_tables(text: str) -> dict:
    """Identify GWAS summary statistics in text"""
    gwas_indicators = {
        'rs_numbers': len(re.findall(r'rs\d+', text, re.IGNORECASE)),
        'p_values': len(re.findall(r'p\s*[<>=]?\s*[\d.e-]+|pvalue|p-value', text, re.IGNORECASE)),
        'odds_ratios': len(re.findall(r'or\s*[=:]\s*[\d.]+|odds\s*ratio', text, re.IGNORECASE)),
        'chromosomes': len(re.findall(r'chr\s*\d+|chromosome\s*\d+', text, re.IGNORECASE)),
        'tables': len(re.findall(r'Table\s+\d+|Supplementary\s+Table', text, re.IGNORECASE))
    }
    
    is_gwas = (
        gwas_indicators['rs_numbers'] > 5 or
        (gwas_indicators['rs_numbers'] > 0 and gwas_indicators['p_values'] > 10)
    )
    
    return {
        'is_gwas': is_gwas,
        'indicators': gwas_indicators,
        'score': sum(gwas_indicators.values())
    }


def extract_paper_info(text: str) -> dict:
    """Extract paper information"""
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
            if 'genome' in line.lower() and 'association' in line.lower():
                info['title'] = line
                break
    
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
    if 'Stroke' in text[:500]:
        info['journal'] = 'Stroke'
    
    return info


def parse_gwas_data_from_text(text: str) -> pd.DataFrame:
    """Parse GWAS summary statistics from text"""
    lines = text.split('\n')
    gwas_rows = []
    
    # Look for the main results table
    # Pattern: rsID chr pos allele p-value OR (CI) [replication data]
    # Example: rs12669789 18 734 065 C 9.52×10–6 1.42 (1.22–1.67) 0.714 1.05 (0.83–1.32) ...
    
    for i, line in enumerate(lines):
        # Check if line contains rs number and looks like a table row
        if re.search(r'^rs\d+', line.strip()):
            parts = line.split()
            if len(parts) >= 5:
                try:
                    row = {}
                    
                    # First part should be SNP ID
                    snp_match = re.match(r'(rs\d+)\*?', parts[0])
                    if snp_match:
                        row['snp'] = snp_match.group(1)
                    else:
                        continue
                    
                    # Next should be chromosome
                    if len(parts) > 1 and parts[1].isdigit():
                        chr_val = int(parts[1])
                        if 1 <= chr_val <= 22:
                            row['chr'] = chr_val
                    
                    # Next should be position (may be split across multiple parts)
                    if len(parts) > 2:
                        # Position might be split like "734 065" or combined "734065"
                        pos_parts = []
                        for j in range(2, min(5, len(parts))):
                            if parts[j].isdigit():
                                pos_parts.append(parts[j])
                            else:
                                break
                        if pos_parts:
                            pos_str = ''.join(pos_parts)
                            if pos_str.isdigit() and len(pos_str) >= 6:
                                row['pos'] = int(pos_str)
                    
                    # Next should be allele
                    if len(parts) > 3:
                        allele = parts[3]
                        if len(allele) == 1 and allele in 'ACGT':
                            row['effect_allele'] = allele
                    
                    # Find p-values (scientific notation with × or E)
                    for part in parts:
                        # Handle ×10 notation
                        if '×10' in part or 'x10' in part.lower():
                            try:
                                # Extract number before ×10 and exponent
                                p_match = re.search(r'([\d.]+)[×x]10[–−-]?(\d+)', part)
                                if p_match:
                                    base = float(p_match.group(1))
                                    exp = int(p_match.group(2))
                                    p_val = base * (10 ** -exp)
                                    if 0 < p_val < 1:
                                        row['p'] = p_val
                            except:
                                pass
                        # Handle E notation
                        elif ('E' in part.upper() or 'e' in part) and ('-' in part or '+' in part):
                            try:
                                part_clean = part.replace('‐', '-').replace('–', '-').replace('−', '-')
                                part_clean = re.sub(r'[^\d.Ee+-]', '', part_clean)
                                p_val = float(part_clean)
                                if 0 < p_val < 1:
                                    row['p'] = p_val
                            except:
                                pass
                    
                    # Find OR values (typically between 0.5 and 2.5)
                    for part in parts:
                        try:
                            val = float(part.replace(',', ''))
                            if 0.5 <= val <= 2.5:
                                if 'or' not in row:
                                    row['or'] = val
                                elif 'or_meta' not in row:
                                    row['or_meta'] = val
                        except:
                            pass
                    
                    # Extract confidence intervals (pattern: (1.22–1.67))
                    ci_match = re.search(r'\(([\d.]+)[–−-]([\d.]+)\)', line)
                    if ci_match:
                        row['or_ci_lower'] = float(ci_match.group(1))
                        row['or_ci_upper'] = float(ci_match.group(2))
                    
                    if 'snp' in row:
                        gwas_rows.append(row)
                except Exception as e:
                    continue
    
    if not gwas_rows:
        return pd.DataFrame()
    
    df = pd.DataFrame(gwas_rows)
    return df


def main():
    """Main extraction function"""
    print("=" * 80)
    print("EXTRACTING DATA FROM FOROUD 2014 MAIN PAPER")
    print("=" * 80)
    
    # Find PDF in Downloads folder
    downloads_path = Path.home() / "Downloads"
    pdf_files = list(downloads_path.glob("*foroud*.pdf"))
    
    if not pdf_files:
        print(f"\n[ERROR] No Foroud PDF found in Downloads folder")
        print(f"   Searched: {downloads_path}")
        return
    
    pdf_path = pdf_files[0]
    print(f"\n[OK] Found PDF: {pdf_path.name}")
    print(f"   Size: {pdf_path.stat().st_size / (1024*1024):.2f} MB")
    
    # Extract text
    print(f"\n[INFO] Extracting text from PDF...")
    text = extract_text_from_pdf(pdf_path)
    
    if not text:
        print("❌ Could not extract text from PDF")
        return
    
    print(f"✓ Extracted {len(text):,} characters")
    
    # Identify paper info
    print(f"\n[INFO] Identifying paper information...")
    paper_info = extract_paper_info(text)
    print(f"  Title: {paper_info.get('title', 'Not found')[:100]}...")
    if paper_info.get('year'):
        print(f"  Year: {paper_info['year']}")
    if paper_info.get('journal'):
        print(f"  Journal: {paper_info['journal']}")
    if paper_info.get('pubmed_id'):
        print(f"  PubMed ID: {paper_info['pubmed_id']}")
    
    # Check for GWAS data
    print(f"\n[INFO] Checking for GWAS summary statistics...")
    gwas_analysis = identify_gwas_tables(text)
    print(f"  RS numbers found: {gwas_analysis['indicators']['rs_numbers']}")
    print(f"  P-values found: {gwas_analysis['indicators']['p_values']}")
    print(f"  Tables found: {gwas_analysis['indicators']['tables']}")
    
    if gwas_analysis['is_gwas']:
        print(f"  ✓ GWAS data detected!")
        
        # Try to parse GWAS data
        print(f"\n[INFO] Parsing GWAS summary statistics...")
        gwas_df = parse_gwas_data_from_text(text)
        
        if not gwas_df.empty:
            print(f"  [OK] Extracted {len(gwas_df)} SNPs")
            print(f"\n  First 10 rows:")
            print(gwas_df.head(10).to_string())
            
            # Save extracted data
            output_dir = Path("data/supplementary_pdfs/extracted_gwas")
            output_dir.mkdir(exist_ok=True)
            
            # Add metadata
            gwas_df['source_pdf'] = pdf_path.name
            gwas_df['study'] = 'Foroud_2014_Stroke'
            if paper_info.get('pubmed_id'):
                gwas_df['pubmed_id'] = paper_info['pubmed_id']
            gwas_df['extraction_date'] = datetime.now().strftime('%Y-%m-%d')
            
            # Save
            tsv_file = output_dir / "foroud_2014_main_paper_gwas.tsv"
            csv_file = output_dir / "foroud_2014_main_paper_gwas.csv"
            
            gwas_df.to_csv(tsv_file, sep='\t', index=False, encoding='utf-8')
            gwas_df.to_csv(csv_file, index=False, encoding='utf-8')
            
            print(f"\n[OK] Saved to:")
            print(f"  - {tsv_file}")
            print(f"  - {csv_file}")
        else:
            print(f"  [WARNING] Could not parse structured GWAS table from text")
            print(f"     Data may be in figures or require manual extraction")
    else:
        print(f"  [WARNING] Limited GWAS data detected in text")
        print(f"     May be in figures/tables that require image extraction")
    
    # Save full text and metadata
    output_dir = Path("data/supplementary_pdfs/extracted_gwas")
    output_dir.mkdir(exist_ok=True)
    
    text_file = output_dir / "foroud_2014_main_paper_text.txt"
    with open(text_file, 'w', encoding='utf-8') as f:
        f.write(text)
    print(f"\n[OK] Saved full text to: {text_file}")
    
    metadata = {
        'pdf_name': pdf_path.name,
        'paper_info': paper_info,
        'gwas_analysis': gwas_analysis,
        'extraction_date': datetime.now().isoformat()
    }
    
    metadata_file = output_dir / "foroud_2014_main_paper_metadata.json"
    with open(metadata_file, 'w', encoding='utf-8') as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)
    print(f"[OK] Saved metadata to: {metadata_file}")
    
    print("\n" + "=" * 80)
    print("EXTRACTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()

