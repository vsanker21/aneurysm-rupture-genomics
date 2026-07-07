"""
Extract Clinical Metadata from GSE122897
Downloads series matrix and extracts clinical variables
"""
import pandas as pd
import numpy as np
import os
import requests
import re
from urllib.parse import urljoin

def download_series_matrix(gse_id="GSE122897"):
    """Download series matrix file from GEO"""
    print(f"\nDownloading series matrix for {gse_id}...")
    
    # GEO series matrix URL
    url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file&file={gse_id}_series_matrix.txt.gz"
    
    output_file = f"data/sc_rnaseq/{gse_id}/{gse_id}_series_matrix.txt.gz"
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    try:
        print(f"  Downloading from: {url}")
        response = requests.get(url, timeout=60, stream=True)
        
        if response.status_code == 200:
            with open(output_file, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            print(f"  [OK] Downloaded to {output_file}")
            return output_file
        else:
            print(f"  [ERROR] HTTP {response.status_code}")
            # Try alternative URL
            alt_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_id[:5]}nnn/{gse_id}/matrix/{gse_id}_series_matrix.txt.gz"
            print(f"  Trying alternative URL: {alt_url}")
            response = requests.get(alt_url, timeout=60, stream=True)
            if response.status_code == 200:
                with open(output_file, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
                print(f"  [OK] Downloaded from alternative URL")
                return output_file
    except Exception as e:
        print(f"  [ERROR] Download failed: {e}")
        return None

def parse_series_matrix(matrix_file):
    """Parse GEO series matrix file"""
    print(f"\nParsing series matrix file...")
    
    import gzip
    
    try:
        # Try to read as gzip
        if matrix_file.endswith('.gz'):
            f = gzip.open(matrix_file, 'rt', encoding='utf-8', errors='ignore')
        else:
            f = open(matrix_file, 'r', encoding='utf-8', errors='ignore')
        
        metadata = {}
        sample_data = {}
        current_sample = None
        in_sample_section = False
        
        for line in f:
            line = line.strip()
            
            if not line:
                continue
            
            # Parse metadata lines (start with !)
            if line.startswith('!'):
                # Remove leading !
                content = line[1:].strip()
                
                if '=' in content:
                    key, value = content.split('=', 1)
                    key = key.strip()
                    value = value.strip().strip('"')
                    
                    # Sample-specific metadata
                    if key.startswith('Sample_'):
                        sample_id = key.replace('Sample_', '').split('_')[0]
                        field = '_'.join(key.split('_')[1:])
                        
                        if sample_id not in sample_data:
                            sample_data[sample_id] = {}
                        
                        sample_data[sample_id][field] = value
                    else:
                        # General metadata
                        metadata[key] = value
        
        f.close()
        
        print(f"  [OK] Parsed metadata for {len(sample_data)} samples")
        
        # Convert to DataFrame
        if sample_data:
            df = pd.DataFrame(sample_data).T
            df.index.name = 'Sample'
            df = df.reset_index()
            return df, metadata
        else:
            return None, metadata
            
    except Exception as e:
        print(f"  [ERROR] Parsing failed: {e}")
        import traceback
        traceback.print_exc()
        return None, None

def extract_clinical_variables(df):
    """Extract and standardize clinical variables"""
    print(f"\nExtracting clinical variables...")
    
    if df is None or len(df) == 0:
        print("  [ERROR] No data to extract from")
        return None
    
    clinical_df = pd.DataFrame()
    clinical_df['Sample'] = df['Sample'] if 'Sample' in df.columns else df.index
    
    # Common clinical variable patterns
    clinical_patterns = {
        'Age': ['age', 'Age', 'AGE', 'years', 'Years'],
        'Sex': ['sex', 'Sex', 'SEX', 'gender', 'Gender', 'GENDER'],
        'Hypertension': ['hypertension', 'Hypertension', 'HTN', 'htn', 'blood_pressure'],
        'Smoking': ['smoking', 'Smoking', 'smoke', 'Smoke', 'tobacco'],
        'Size': ['size', 'Size', 'SIZE', 'diameter', 'Diameter', 'mm'],
        'Location': ['location', 'Location', 'LOCATION', 'site', 'Site'],
        'Status': ['status', 'Status', 'STATUS', 'group', 'Group', 'condition']
    }
    
    found_vars = {}
    
    # Search for clinical variables in column names
    for var_name, patterns in clinical_patterns.items():
        for col in df.columns:
            col_lower = str(col).lower()
            for pattern in patterns:
                if pattern.lower() in col_lower:
                    if var_name not in found_vars:
                        found_vars[var_name] = col
                        clinical_df[var_name] = df[col]
                        print(f"  [OK] Found {var_name}: {col}")
                    break
    
    # Also check in characteristics_ch1 or similar columns
    char_cols = [c for c in df.columns if 'characteristic' in str(c).lower() or 'ch1' in str(c).lower()]
    
    for col in char_cols:
        if col in df.columns:
            # Try to extract information from characteristics
            for idx, val in df[col].items():
                if pd.notna(val):
                    val_str = str(val).lower()
                    
                    # Extract age
                    age_match = re.search(r'age[:\s]*(\d+)', val_str)
                    if age_match and 'Age' not in found_vars:
                        age = int(age_match.group(1))
                        if 'Age' not in clinical_df.columns:
                            clinical_df['Age'] = np.nan
                        clinical_df.loc[idx, 'Age'] = age
                        if 'Age' not in found_vars:
                            found_vars['Age'] = col
                            print(f"  [OK] Extracted Age from {col}")
                    
                    # Extract sex/gender
                    if ('female' in val_str or 'f' in val_str) and 'Sex' not in found_vars:
                        if 'Sex' not in clinical_df.columns:
                            clinical_df['Sex'] = 'Unknown'
                        clinical_df.loc[idx, 'Sex'] = 'Female'
                        if 'Sex' not in found_vars:
                            found_vars['Sex'] = col
                            print(f"  [OK] Extracted Sex from {col}")
                    elif ('male' in val_str or 'm' in val_str) and 'Sex' not in found_vars:
                        if 'Sex' not in clinical_df.columns:
                            clinical_df['Sex'] = 'Unknown'
                        clinical_df.loc[idx, 'Sex'] = 'Male'
                        if 'Sex' not in found_vars:
                            found_vars['Sex'] = col
                            print(f"  [OK] Extracted Sex from {col}")
                    
                    # Extract hypertension
                    if ('hypertension' in val_str or 'htn' in val_str):
                        if 'Hypertension' not in clinical_df.columns:
                            clinical_df['Hypertension'] = np.nan
                        # Set 1 if present and not explicitly 'no'
                        clinical_df.loc[idx, 'Hypertension'] = 0 if 'no hypertension' in val_str or 'htn: 0' in val_str else 1
                        if 'Hypertension' not in found_vars:
                            found_vars['Hypertension'] = col
                            print(f"  [OK] Extracted Hypertension from {col}")

                    # Extract aneurysm size (mm)
                    size_match = re.search(r'(size|diameter)[:\s]*([0-9]+\.?[0-9]*)\s*mm', val_str)
                    if size_match:
                        size_val = float(size_match.group(2))
                        if 'Size' not in clinical_df.columns:
                            clinical_df['Size'] = np.nan
                        clinical_df.loc[idx, 'Size'] = size_val
                        if 'Size' not in found_vars:
                            found_vars['Size'] = col
                            print(f"  [OK] Extracted Size from {col}")

                    # Extract location/site
                    if ('location' in val_str or 'site' in val_str) and 'Location' not in clinical_df.columns:
                        # take the full value as location descriptor
                        clinical_df.loc[idx, 'Location'] = val
                        if 'Location' not in found_vars:
                            found_vars['Location'] = col
                            print(f"  [OK] Extracted Location from {col}")
    
    # Check source_name_ch1 for additional info
    source_cols = [c for c in df.columns if 'source' in str(c).lower()]
    for col in source_cols:
        if col in df.columns:
            print(f"  [INFO] Found source column: {col}")
            # May contain tissue type or other info
    
    print(f"\n  [OK] Extracted {len(found_vars)} clinical variables")
    
    return clinical_df

def merge_with_existing_metadata(clinical_df, existing_metadata_file):
    """Merge with existing metadata"""
    print(f"\nMerging with existing metadata...")
    
    if not os.path.exists(existing_metadata_file):
        print(f"  [WARNING] Existing metadata file not found: {existing_metadata_file}")
        return clinical_df
    
    existing = pd.read_csv(existing_metadata_file, sep='\t')
    print(f"  [OK] Loaded existing metadata: {existing.shape}")
    
    # Merge
    if 'Sample' in existing.columns and 'Sample' in clinical_df.columns:
        merged = existing.merge(clinical_df, on='Sample', how='outer', suffixes=('_old', '_new'))
        print(f"  [OK] Merged metadata: {merged.shape}")
        return merged
    else:
        print(f"  [WARNING] Cannot merge - Sample column mismatch")
        return clinical_df

def main():
    """Main extraction pipeline"""
    print("=" * 80)
    print("EXTRACTING CLINICAL METADATA FROM GSE122897")
    print("=" * 80)
    
    gse_id = "GSE122897"
    data_dir = f"data/sc_rnaseq/{gse_id}"
    os.makedirs(data_dir, exist_ok=True)
    
    # Step 1: Download series matrix
    matrix_file = download_series_matrix(gse_id)
    
    if matrix_file is None or not os.path.exists(matrix_file):
        print("\n[WARNING] Could not download series matrix")
        print("  Trying to use existing metadata file...")
        existing_file = f"{data_dir}/sample_metadata.tsv"
        if os.path.exists(existing_file):
            df = pd.read_csv(existing_file, sep='\t')
            clinical_df = extract_clinical_variables(df)
        else:
            print("[ERROR] No metadata available")
            return
    else:
        # Step 2: Parse series matrix
        df, metadata = parse_series_matrix(matrix_file)
        
        if df is None:
            print("\n[ERROR] Could not parse series matrix")
            return
        
        # Step 3: Extract clinical variables
        clinical_df = extract_clinical_variables(df)
    
    if clinical_df is None:
        print("\n[ERROR] Could not extract clinical variables")
        return
    
    # Step 4: Merge with existing metadata
    existing_file = f"{data_dir}/sample_metadata.tsv"
    if os.path.exists(existing_file):
        clinical_df = merge_with_existing_metadata(clinical_df, existing_file)
    
    # Step 5: Save results
    output_file = f"{data_dir}/sample_metadata_with_clinical.tsv"
    clinical_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved clinical metadata to {output_file}")
    
    # Summary
    print("\n" + "=" * 80)
    print("CLINICAL METADATA EXTRACTION SUMMARY")
    print("=" * 80)
    print(f"\nSamples: {len(clinical_df)}")
    print(f"\nClinical variables found:")
    clinical_cols = [c for c in clinical_df.columns if c not in ['Sample']]
    for col in clinical_cols:
        non_null = clinical_df[col].notna().sum()
        print(f"  {col}: {non_null}/{len(clinical_df)} samples ({non_null/len(clinical_df)*100:.1f}%)")
        if non_null > 0:
            print(f"    Values: {clinical_df[col].value_counts().head(5).to_dict()}")
    
    print(f"\n[OK] Clinical metadata extraction complete!")
    print(f"\nNext step: Integrate clinical features into multi-omics model")

if __name__ == "__main__":
    main()



