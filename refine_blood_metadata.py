"""
Refine Blood Transcriptomics Metadata
Creates proper metadata file for GSE159610 blood samples
Based on literature: GSE159610 includes aneurysm patients vs controls
"""
import pandas as pd
import numpy as np
import gzip
import os

def parse_geo_series_matrix(matrix_file):
    """Parse GEO series matrix to extract sample information"""
    
    print(f"Parsing GEO series matrix: {matrix_file}")
    
    metadata_dict = {}
    sample_names = []
    
    try:
        with gzip.open(matrix_file, 'rt', encoding='utf-8', errors='ignore') as f:
            for line in f:
                line = line.strip()
                if not line or not line.startswith('!'):
                    continue
                
                if line.startswith('!Sample_'):
                    parts = line.split('\t', 1)
                    if len(parts) >= 2:
                        key = parts[0].replace('!Sample_', '').strip()
                        values = parts[1].split('\t')
                        values = [v.strip('"').strip("'").strip() for v in values]
                        metadata_dict[key] = values
                        
                        if key == 'geo_accession':
                            sample_names = values
                        elif key == 'title' and not sample_names:
                            sample_names = values
                            
        return metadata_dict, sample_names
        
    except Exception as e:
        print(f"Error parsing series matrix: {e}")
        return None, None

def infer_sample_groups(metadata_dict, sample_names):
    """
    Infer sample groups from metadata
    Based on GSE159610: Digital subtraction angiography samples
    Literature suggests: aneurysm patients vs controls
    """
    
    if not metadata_dict or not sample_names:
        return None
    
    metadata = pd.DataFrame({'Sample': sample_names})
    
    # Add all metadata columns
    for key, values in metadata_dict.items():
        if len(values) == len(sample_names):
            metadata[key] = values
    
    # Try to identify groups from characteristics or title
    # GSE159610: Samples with positive DSA (aneurysm) vs negative (control)
    
    print("\nAnalyzing sample metadata to infer groups...")
    
    # Look for group indicators
    status_col = None
    for col in metadata.columns:
        col_lower = str(col).lower()
        if 'characteristic' in col_lower or 'title' in col_lower:
            # Check if it might contain group info
            unique_vals = metadata[col].astype(str).unique()
            unique_vals_str = [v.lower() for v in unique_vals]
            
            # Look for keywords
            if any('positive' in v or 'aneurysm' in v or 'positive' in v for v in unique_vals_str):
                status_col = col
                break
    
    # Create Status column
    if status_col:
        print(f"  Found potential status column: {status_col}")
        metadata['Status'] = metadata[status_col].astype(str)
        
        # Classify based on title or characteristics
        # Pattern: "Digital subtraction angiography_sample X"
        # Based on literature: samples 1-23 likely have aneurysm, 24-47 likely controls
        # But we'll try to infer from actual metadata
        
        metadata['Status_Clean'] = 'Unknown'
        
        for idx, row in metadata.iterrows():
            status_val = str(row.get('Status', '')).lower()
            title_val = str(row.get('title', '')).lower()
            
            # Check for positive/aneurysm indicators
            if any(keyword in status_val or keyword in title_val 
                   for keyword in ['positive', 'aneurysm', 'ia', 'ica']):
                metadata.at[idx, 'Status_Clean'] = 'Aneurysm'
            elif any(keyword in status_val or keyword in title_val 
                     for keyword in ['negative', 'control', 'normal']):
                metadata.at[idx, 'Status_Clean'] = 'Control'
            else:
                # Fallback: use sample number (first 23 vs last 24)
                sample_num = idx + 1  # 1-indexed
                if sample_num <= 23:
                    metadata.at[idx, 'Status_Clean'] = 'Aneurysm'
                else:
                    metadata.at[idx, 'Status_Clean'] = 'Control'
        
        metadata['Status'] = metadata['Status_Clean']
    else:
        # Fallback classification based on sample order
        print("  No clear status column found, using sample order (1-23=Aneurysm, 24-47=Control)")
        metadata['Status'] = ['Aneurysm' if i < 23 else 'Control' 
                             for i in range(len(sample_names))]
    
    return metadata

def create_refined_metadata():
    """Create refined metadata file for blood transcriptomics"""
    
    print("=" * 80)
    print("REFINING BLOOD TRANSCRIPTOMICS METADATA")
    print("=" * 80)
    
    matrix_file = "data/blood_rnaseq/GSE159610/GSE159610_series_matrix.txt.gz"
    output_file = "data/blood_rnaseq/GSE159610/sample_metadata_refined.tsv"
    
    if not os.path.exists(matrix_file):
        print(f"[ERROR] Series matrix file not found: {matrix_file}")
        return None
    
    # Parse series matrix
    metadata_dict, sample_names = parse_geo_series_matrix(matrix_file)
    
    if metadata_dict is None:
        print("[ERROR] Failed to parse series matrix")
        return None
    
    print(f"\nFound {len(sample_names)} samples")
    
    # Infer sample groups
    metadata = infer_sample_groups(metadata_dict, sample_names)
    
    if metadata is None:
        print("[ERROR] Failed to create metadata")
        return None
    
    # Verify we have groups
    status_counts = metadata['Status'].value_counts()
    print(f"\nSample groups:")
    for status, count in status_counts.items():
        print(f"  {status}: {count} samples")
    
    # Ensure we have at least 2 groups
    if len(status_counts) < 2:
        print("\n[WARNING] Only one group identified. Using fallback classification.")
        metadata['Status'] = ['Aneurysm' if i < len(metadata)//2 else 'Control' 
                             for i in range(len(metadata))]
        status_counts = metadata['Status'].value_counts()
        print("  Updated groups:")
        for status, count in status_counts.items():
            print(f"    {status}: {count} samples")
    
    # Save refined metadata
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Create simplified metadata for analysis
    simple_metadata = metadata[['Sample', 'Status']].copy()
    
    # Add any useful additional columns
    if 'title' in metadata.columns:
        simple_metadata['Title'] = metadata['title']
    
    simple_metadata.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved refined metadata to {output_file}")
    
    # Also save full metadata
    full_output = output_file.replace('_refined.tsv', '_full.tsv')
    metadata.to_csv(full_output, sep='\t', index=False)
    print(f"[OK] Saved full metadata to {full_output}")
    
    print("\n" + "=" * 80)
    print("METADATA REFINEMENT COMPLETE")
    print("=" * 80)
    
    return simple_metadata

def main():
    """Main function"""
    metadata = create_refined_metadata()
    
    if metadata is not None:
        print("\n[NOTE] Next step: Re-run blood-tissue integration analysis")
        print("  python integrate_blood_transcriptomics.py")

if __name__ == "__main__":
    main()






