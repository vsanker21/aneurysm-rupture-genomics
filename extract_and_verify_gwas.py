"""
Extract and verify Bakker 2020 GWAS summary statistics
"""
import zipfile
import os
import pandas as pd
import numpy as np

def extract_and_verify_gwas():
    """Extract and verify GWAS data format"""
    
    gwas_file = "data/gwas/bakker_2020/bakker_2020_summary_stats.txt.gz"
    output_dir = "data/gwas/bakker_2020/extracted"
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 80)
    print("EXTRACTING AND VERIFYING GWAS DATA")
    print("=" * 80)
    print(f"\nSource file: {gwas_file}")
    print(f"File size: {os.path.getsize(gwas_file) / (1024*1024):.2f} MB")
    
    # Check if it's a ZIP file (magic bytes 'PK')
    with open(gwas_file, 'rb') as f:
        magic = f.read(2)
    
    if magic == b'PK':
        print("\n✓ Detected ZIP archive format")
        
        # Extract ZIP file
        print("\nExtracting ZIP archive...")
        with zipfile.ZipFile(gwas_file, 'r') as zip_ref:
            file_list = zip_ref.namelist()
            print(f"Files in archive: {len(file_list)}")
            for f in file_list:
                print(f"  - {f}")
            
            # Extract all files
            zip_ref.extractall(output_dir)
            print(f"\n✓ Extracted to: {output_dir}")
            
            # Find the main summary statistics file
            summary_files = [f for f in file_list if any(x in f.lower() for x in ['summary', 'stats', '.txt', '.tsv', '.csv'])]
            
            if summary_files:
                main_file = os.path.join(output_dir, summary_files[0])
                print(f"\nMain summary statistics file: {main_file}")
                
                # Try to read first few lines
                print("\nReading file header...")
                try:
                    # Try different encodings and delimiters
                    for encoding in ['utf-8', 'latin-1', 'cp1252']:
                        for sep in ['\t', ',', ' ']:
                            try:
                                df_sample = pd.read_csv(main_file, nrows=5, sep=sep, encoding=encoding, low_memory=False)
                                if len(df_sample.columns) > 3:  # Likely correct format
                                    print(f"\n✓ Successfully read with encoding={encoding}, separator='{sep}'")
                                    print(f"\nColumns ({len(df_sample.columns)}):")
                                    for i, col in enumerate(df_sample.columns):
                                        print(f"  {i+1}. {col}")
                                    
                                    print(f"\nFirst 3 rows:")
                                    print(df_sample.head(3).to_string())
                                    
                                    # Get file info
                                    file_size = os.path.getsize(main_file) / (1024*1024)
                                    print(f"\nFile size: {file_size:.2f} MB")
                                    
                                    # Count total rows (estimate)
                                    print("\nCounting total rows (this may take a moment)...")
                                    total_rows = sum(1 for _ in open(main_file, 'r', encoding=encoding)) - 1
                                    print(f"Total rows (excluding header): {total_rows:,}")
                                    
                                    return main_file, df_sample.columns.tolist(), encoding, sep
                            except:
                                continue
                except Exception as e:
                    print(f"Error reading file: {e}")
                    return None, None, None, None
    else:
        print("\nFile format not recognized. Trying gzip...")
        import gzip
        try:
            with gzip.open(gwas_file, 'rt', encoding='utf-8', errors='ignore') as f:
                lines = [f.readline() for _ in range(5)]
                print("\nFirst 5 lines:")
                for i, line in enumerate(lines):
                    print(f"{i+1}: {line[:200]}")
        except Exception as e:
            print(f"Error: {e}")
            return None, None, None, None
    
    return None, None, None, None

if __name__ == "__main__":
    extract_and_verify_gwas()






