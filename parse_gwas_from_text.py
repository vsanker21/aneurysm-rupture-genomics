"""
Parse GWAS summary statistics from extracted PDF text
Extracts tabular GWAS data from text format
"""
import re
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional


def parse_gwas_table_from_text(text: str, pdf_name: str) -> Optional[pd.DataFrame]:
    """Parse GWAS summary statistics table from text"""
    
    # Look for table patterns - lines with rs numbers followed by numbers
    lines = text.split('\n')
    
    gwas_rows = []
    header_found = False
    header_line = None
    
    # Pattern to match SNP rows: rs number, chromosome, position, p-value, etc.
    snp_pattern = re.compile(r'rs\d+')
    
    for i, line in enumerate(lines):
        # Check if line contains rs number
        if snp_pattern.search(line):
            # Try to parse the line
            parts = line.split()
            
            # Find rs number
            rs_match = snp_pattern.search(line)
            if rs_match:
                rs_id = rs_match.group()
                
                # Try to extract other fields
                # Common format: rsID chr pos p-value ...
                row_data = {'snp': rs_id}
                
                # Extract chromosome (usually after rs number)
                chr_match = re.search(r'\b([1-9]|1[0-9]|2[0-2]|X|Y|MT)\b', line)
                if chr_match:
                    chr_val = chr_match.group(1)
                    # Convert X, Y, MT to numbers for consistency
                    chr_map = {'X': '23', 'Y': '24', 'MT': '25'}
                    row_data['chr'] = chr_map.get(chr_val, chr_val)
                
                # Extract position (large number, usually 6-9 digits)
                pos_match = re.search(r'\b([1-9]\d{5,9})\b', line)
                if pos_match:
                    row_data['pos'] = pos_match.group(1)
                
                # Extract p-value (scientific notation or decimal)
                p_match = re.search(r'([\d.]+[Ee][+-]?\d+|0\.\d+|1\.\d+e-\d+)', line)
                if p_match:
                    p_val = p_match.group(1)
                    try:
                        row_data['p'] = float(p_val)
                    except:
                        pass
                
                # Extract odds ratio or beta (number between 0.1 and 10, or negative)
                or_match = re.search(r'\b([0-9]\.[0-9]{1,3}|[0-9]\.[0-9]{1,3})\b', line)
                if or_match:
                    or_val = or_match.group(1)
                    try:
                        val = float(or_val)
                        if 0.1 <= val <= 10:
                            row_data['or'] = val
                    except:
                        pass
                
                # Try to extract all numeric values and assign by position
                numbers = re.findall(r'[\d.]+[Ee][+-]?\d+|[\d.]+', line)
                numeric_values = []
                for num_str in numbers:
                    try:
                        numeric_values.append(float(num_str))
                    except:
                        pass
                
                # If we have enough data, add the row
                if len(row_data) >= 3:  # At least SNP, chr, and one other field
                    gwas_rows.append(row_data)
    
    if not gwas_rows:
        return None
    
    # Create DataFrame
    df = pd.DataFrame(gwas_rows)
    
    # Standardize column names
    standard_cols = {
        'snp': 'snp',
        'chr': 'chr',
        'chromosome': 'chr',
        'pos': 'pos',
        'position': 'pos',
        'bp': 'pos',
        'p': 'p',
        'pvalue': 'p',
        'p_value': 'p',
        'pval': 'p',
        'or': 'or',
        'odds_ratio': 'or',
        'beta': 'beta',
        'effect': 'beta'
    }
    
    # Rename columns
    df = df.rename(columns={col: standard_cols.get(col.lower(), col) for col in df.columns})
    
    # Add metadata
    df['source_pdf'] = pdf_name
    df['extraction_method'] = 'text_parsing'
    
    return df


def parse_bilguvar_table(text: str) -> Optional[pd.DataFrame]:
    """Parse the specific Bilguvar 2008 table format"""
    
    lines = text.split('\n')
    gwas_rows = []
    
    # Look for the table section (around line 314 based on inspection)
    in_table = False
    table_end_markers = ['Supplementary Table', '14', '15', '16']
    
    for i, line in enumerate(lines):
        # Check if we're in the table section
        if 'rs911271' in line or ('rs1429412' in line and 'Chr' in lines[max(0, i-3):i]):
            in_table = True
        
        if in_table:
            # Pattern: rsID chr pos p-value [other values]
            # Example: rs911271   1  234,376,647  4.61E‐07  0.002 0.002  0.000  0.84  0.28  0.51
            parts = line.split()
            if len(parts) >= 4:
                # First part should be rs number
                if parts[0].startswith('rs'):
                    try:
                        row = {
                            'snp': parts[0],
                        }
                        
                        # Chromosome (usually second field, but might have spaces)
                        chr_idx = 1
                        if parts[chr_idx].isdigit():
                            row['chr'] = int(parts[chr_idx])
                        else:
                            # Try to find chromosome number
                            for j, part in enumerate(parts[1:4]):
                                if part.isdigit() and 1 <= int(part) <= 22:
                                    row['chr'] = int(part)
                                    chr_idx = j + 1
                                    break
                        
                        # Position (usually after chromosome, may have commas)
                        if chr_idx + 1 < len(parts):
                            pos_str = parts[chr_idx + 1].replace(',', '')
                            if pos_str.replace('.', '').isdigit():
                                row['pos'] = int(float(pos_str))
                        
                        # Find p-value (scientific notation with E or e)
                        # Look for pattern like 4.61E-07 or 1.59E-09
                        for part in parts:
                            if ('E' in part.upper() or 'e' in part):
                                try:
                                    # Handle special dash character (‐) vs regular dash (-)
                                    part_clean = part.replace('‐', '-').replace('–', '-').replace('−', '-')
                                    # Remove any trailing characters
                                    part_clean = re.sub(r'[^\d.Ee+-]', '', part_clean)
                                    p_val = float(part_clean)
                                    if p_val > 0 and p_val < 1:  # Valid p-value range
                                        row['p'] = p_val
                                        break
                                except:
                                    pass
                        
                        # Try to extract OR or beta from numeric values
                        # Look for values that could be OR (typically 0.5-2.0) or beta
                        for part in parts[3:10]:  # Check first few numeric fields after p-value
                            try:
                                val = float(part.replace(',', ''))
                                # OR typically between 0.5 and 2.0
                                if 0.5 <= val <= 2.0 and 'or' not in row:
                                    row['or'] = val
                                # Beta can be negative or small positive
                                elif -2.0 <= val <= 2.0 and 'beta' not in row and val != 0:
                                    row['beta'] = val
                            except:
                                pass
                        
                        if 'snp' in row and row['snp'] and ('chr' in row or 'p' in row):
                            gwas_rows.append(row)
                    except Exception as e:
                        continue
        
        # Stop if we've moved past the table
        if in_table:
            # Check if we've hit a new section
            line_lower = line.strip().lower()
            if any(marker.lower() in line_lower for marker in table_end_markers):
                # Only stop if we've collected a reasonable number of rows
                if len(gwas_rows) >= 10:
                    break
            # Also stop if we see a long gap without rs numbers (but allow more rows)
            if i > 0 and len(gwas_rows) > 0:
                # Check if last 10 lines had no rs numbers
                recent_lines = lines[max(0, i-10):i]
                if not any('rs' in l for l in recent_lines):
                    # Only stop if we have a good number of rows
                    if len(gwas_rows) >= 10:
                        break
    
    if not gwas_rows:
        return None
    
    df = pd.DataFrame(gwas_rows)
    
    # Clean and convert data types
    if 'chr' in df.columns:
        df['chr'] = pd.to_numeric(df['chr'], errors='coerce')
    if 'pos' in df.columns:
        df['pos'] = df['pos'].astype(str).str.replace(',', '').astype(float)
    if 'p' in df.columns:
        df['p'] = pd.to_numeric(df['p'], errors='coerce')
    if 'or' in df.columns:
        df['or'] = pd.to_numeric(df['or'], errors='coerce')
    if 'beta' in df.columns:
        df['beta'] = pd.to_numeric(df['beta'], errors='coerce')
    
    # Add metadata
    df['source_pdf'] = 'NIHMS4348-supplement-s1.pdf'
    df['study'] = 'Bilguvar_2008_NatureGenetics'
    df['pubmed_id'] = '18758463'
    df['extraction_method'] = 'text_parsing'
    df['extraction_date'] = pd.Timestamp.now().strftime('%Y-%m-%d')
    
    return df


def main():
    """Main parsing function"""
    print("=" * 80)
    print("PARSING GWAS DATA FROM EXTRACTED TEXT")
    print("=" * 80)
    
    pdf_dir = Path("data/supplementary_pdfs/extracted_gwas")
    print(f"\nLooking for text files in: {pdf_dir.absolute()}")
    print(f"Directory exists: {pdf_dir.exists()}")
    
    text_files = list(pdf_dir.glob("*_text.txt"))
    print(f"Found {len(text_files)} text files")
    
    if not text_files:
        print(f"\n❌ No text files found in {pdf_dir}")
        # Try alternative location
        alt_dir = Path("data/supplementary_pdfs/extracted")
        alt_files = list(alt_dir.glob("*_text.txt"))
        if alt_files:
            print(f"Found files in alternative location: {alt_dir}")
            text_files = alt_files
        else:
            return
    
    print(f"\n📁 Found {len(text_files)} text files to parse")
    
    all_gwas_data = {}
    
    for text_file in text_files:
        print(f"\n{'='*80}")
        print(f"Parsing: {text_file.name}")
        print(f"{'='*80}")
        
        # Read text
        with open(text_file, 'r', encoding='utf-8') as f:
            text = f.read()
        
        pdf_name = text_file.name.replace('_text.txt', '')
        
        # Try specific parser for Bilguvar
        if 'NIHMS4348' in pdf_name:
            print("  Using Bilguvar-specific parser...")
            gwas_df = parse_bilguvar_table(text)
        else:
            print("  Using general parser...")
            gwas_df = parse_gwas_table_from_text(text, pdf_name)
        
        if gwas_df is not None and not gwas_df.empty:
            print(f"  ✓ Extracted {len(gwas_df)} SNPs")
            print(f"  Columns: {list(gwas_df.columns)}")
            print(f"\n  First 5 rows:")
            print(gwas_df.head().to_string())
            
            # Save to TSV and CSV
            output_base = pdf_dir / pdf_name
            tsv_file = output_base.parent / f"{pdf_name}_gwas_parsed.tsv"
            csv_file = output_base.parent / f"{pdf_name}_gwas_parsed.csv"
            
            gwas_df.to_csv(tsv_file, sep='\t', index=False, encoding='utf-8')
            gwas_df.to_csv(csv_file, index=False, encoding='utf-8')
            
            print(f"\n  ✓ Saved to:")
            print(f"    - {tsv_file}")
            print(f"    - {csv_file}")
            
            all_gwas_data[pdf_name] = gwas_df
        else:
            print("  ⚠ Could not extract GWAS data from text")
    
    # Summary
    print("\n" + "=" * 80)
    print("PARSING SUMMARY")
    print("=" * 80)
    
    total_snps = 0
    for pdf_name, df in all_gwas_data.items():
        print(f"\n📄 {pdf_name}")
        print(f"  SNPs extracted: {len(df)}")
        if 'p' in df.columns:
            sig_count = (df['p'] < 0.05).sum() if df['p'].notna().any() else 0
            print(f"  Significant (p < 0.05): {sig_count}")
        if 'chr' in df.columns:
            unique_chr = df['chr'].nunique()
            print(f"  Unique chromosomes: {unique_chr}")
        total_snps += len(df)
    
    print(f"\n✓ Total SNPs extracted: {total_snps}")
    print(f"✓ All data saved to: data/supplementary_pdfs/extracted_gwas/")


if __name__ == "__main__":
    main()

