"""
Helper script to parse GEO series matrix and extract metadata for blood transcriptomics
"""
import gzip
import pandas as pd
import os

def parse_geo_series_matrix(matrix_file):
    """Parse GEO series matrix file"""
    
    print(f"Parsing GEO series matrix: {matrix_file}")
    
    metadata_dict = {}
    sample_names = []
    
    try:
        with gzip.open(matrix_file, 'rt', encoding='utf-8', errors='ignore') as f:
            for line in f:
                line = line.strip()
                if not line or not line.startswith('!'):
                    continue
                
                # Parse metadata lines
                if line.startswith('!Sample_'):
                    parts = line.split('\t', 1)
                    if len(parts) >= 2:
                        key = parts[0].replace('!Sample_', '').strip()
                        values = parts[1].split('\t')
                        
                        # Clean up values (remove quotes)
                        values = [v.strip('"').strip("'").strip() for v in values]
                        
                        metadata_dict[key] = values
                        
                        # Get sample names from title or geo_accession
                        if key == 'geo_accession':
                            sample_names = values
                        elif key == 'title' and not sample_names:
                            sample_names = values
        return metadata_dict, sample_names
        
    except Exception as e:
        print(f"Error parsing series matrix: {e}")
        return None, None

def create_blood_metadata(metadata_dict, sample_names):
    """Create metadata DataFrame for blood samples"""
    
    if not metadata_dict or not sample_names:
        return None
    
    # Create metadata DataFrame
    metadata = pd.DataFrame({'Sample': sample_names})
    
    # Try to extract useful metadata
    # Common columns: title, characteristics_ch1, source_name_ch1, etc.
    for key, values in metadata_dict.items():
        if len(values) == len(sample_names):
            metadata[key] = values
    
    # Try to identify status/group column
    # Look for characteristics or title columns that might indicate case/control
    status_col = None
    for col in metadata.columns:
        if 'characteristics' in col.lower() or 'title' in col.lower():
            unique_vals = metadata[col].unique()
            # Check if it might be status/group info
            unique_vals_str = [str(v).lower() for v in unique_vals if pd.notna(v)]
            if any('aneurysm' in v or 'control' in v or 'case' in v or 'normal' in v 
                   for v in unique_vals_str):
                status_col = col
                break
    
    # Create Status column if we found something useful
    if status_col:
        metadata['Status'] = metadata[status_col]
        # Try to clean up status values
        metadata['Status'] = metadata['Status'].str.replace('"', '').str.strip()
        print(f"Found status column: {status_col}")
    else:
        # Default: unknown
        metadata['Status'] = 'Unknown'
        print("Could not identify status column - set to 'Unknown'")
    
    return metadata

def main():
    """Main function"""
    
    matrix_file = "data/blood_rnaseq/GSE159610/GSE159610_series_matrix.txt.gz"
    output_file = "data/blood_rnaseq/GSE159610/sample_metadata.tsv"
    
    if not os.path.exists(matrix_file):
        print(f"Series matrix file not found: {matrix_file}")
        return
    
    # Parse series matrix
    metadata_dict, sample_names = parse_geo_series_matrix(matrix_file)
    
    if metadata_dict is None:
        print("Failed to parse series matrix")
        return
    
    print(f"Found {len(sample_names)} samples")
    print(f"Found {len(metadata_dict)} metadata fields")
    
    # Create metadata DataFrame
    metadata = create_blood_metadata(metadata_dict, sample_names)
    
    if metadata is not None:
        # Save metadata
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        metadata.to_csv(output_file, sep='\t', index=False)
        print(f"\n[OK] Saved metadata to {output_file}")
        print(f"\nMetadata columns: {list(metadata.columns)}")
        print(f"\nStatus values: {metadata['Status'].value_counts().to_dict()}")
        print(f"\nFirst few rows:")
        print(metadata.head())

if __name__ == "__main__":
    main()






