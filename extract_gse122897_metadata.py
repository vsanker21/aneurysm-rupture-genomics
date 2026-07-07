"""
Extract sample metadata from GSE122897 series matrix file
"""
import gzip
import pandas as pd
import re
from pathlib import Path

def extract_metadata_from_series_matrix(series_matrix_file):
    """Extract sample metadata from GEO series matrix file"""
    
    print("=" * 80)
    print("EXTRACTING GSE122897 SAMPLE METADATA")
    print("=" * 80)
    
    print(f"\nReading: {series_matrix_file}")
    
    # Read the file
    with gzip.open(series_matrix_file, 'rt', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    
    print(f"✓ Read {len(lines)} lines")
    
    # Parse series matrix format
    # Format: !Sample_<attribute> = <value1> <value2> ...
    metadata_dict = {}
    sample_names = []
    
    for line in lines:
        line = line.strip()
        
        if line.startswith('!Sample_'):
            # Extract attribute name and values
            # Format: !Sample_<attr>    "value1"   "value2" ...
            parts = line.split('\t', 1)
            if len(parts) == 2:
                attr_name = parts[0].replace('!Sample_', '').strip()
                # Parse tab-separated quoted values
                values_str = parts[1]
                # Split by tab and remove quotes
                values = [v.strip().strip('"') for v in values_str.split('\t')]
                
                metadata_dict[attr_name] = values
                
                # Get sample names from geo_accession or title
                if attr_name == 'geo_accession':
                    sample_names = values
                    print(f"\n✓ Found {len(sample_names)} samples from geo_accession")
                elif attr_name == 'title' and not sample_names:
                    # Use titles as sample identifiers if geo_accession not found
                    sample_names = values
                    print(f"\n✓ Found {len(sample_names)} samples from title")
    
    if not sample_names:
        # Try to get from title
        if 'title' in metadata_dict:
            sample_names = metadata_dict['title']
        else:
            print("⚠️  Could not find sample names")
            return None
    
    # Create metadata DataFrame
    metadata_df = pd.DataFrame({
        'Sample': sample_names
    })
    
    # Add all other attributes as columns
    for attr, values in metadata_dict.items():
        if attr not in ['geo_accession', 'title']:
            if len(values) == len(sample_names):
                metadata_df[attr] = values
    
    print(f"\n✓ Created metadata DataFrame: {metadata_df.shape}")
    print(f"\nColumns: {list(metadata_df.columns)}")
    
    # Display first few rows
    print(f"\nFirst 5 rows:")
    print(metadata_df.head().to_string())
    
    # Look for rupture/unruptured information
    print("\n" + "=" * 80)
    print("SEARCHING FOR RUPTURE STATUS INFORMATION")
    print("=" * 80)
    
    # Check all columns for rupture-related keywords
    rupture_keywords = ['rupture', 'ruptured', 'unrupture', 'unruptured', 'status', 
                       'group', 'condition', 'phenotype', 'characteristics', 'type']
    
    relevant_columns = []
    for col in metadata_df.columns:
        col_lower = col.lower()
        if any(keyword in col_lower for keyword in rupture_keywords):
            relevant_columns.append(col)
    
    if relevant_columns:
        print(f"\n✓ Found potentially relevant columns:")
        for col in relevant_columns:
            print(f"\n  {col}:")
            unique_vals = metadata_df[col].unique()
            print(f"    Unique values: {unique_vals}")
            print(f"    Value counts:")
            print(metadata_df[col].value_counts().to_string())
    else:
        print("\n⚠️  No obvious rupture status column found")
        print("\nAll columns:")
        for col in metadata_df.columns:
            unique_vals = metadata_df[col].unique()[:5]  # First 5 unique values
            print(f"  {col}: {unique_vals}")
    
    # Try to infer rupture status from available data
    print("\n" + "=" * 80)
    print("ATTEMPTING TO INFER RUPTURE STATUS")
    print("=" * 80)
    
    # Check if we can create a Group column
    if 'characteristics_ch1' in metadata_df.columns:
        print("\nChecking 'characteristics_ch1' column...")
        char_values = metadata_df['characteristics_ch1'].unique()
        print(f"Unique values: {char_values}")
        
        # Try to extract rupture status
        metadata_df['Group'] = 'Unknown'
        
        for idx, row in metadata_df.iterrows():
            char = str(row.get('characteristics_ch1', '')).lower()
            if 'rupture' in char or 'ruptured' in char:
                if 'unrupture' in char or 'unrupt' in char:
                    metadata_df.at[idx, 'Group'] = 'Unruptured'
                else:
                    metadata_df.at[idx, 'Group'] = 'Ruptured'
            elif 'control' in char:
                metadata_df.at[idx, 'Group'] = 'Control'
    
    # Check other potential columns
    for col in metadata_df.columns:
        if col != 'Sample' and col != 'Group':
            col_lower = str(metadata_df[col].iloc[0]).lower() if len(metadata_df) > 0 else ''
            if 'rupture' in col_lower or any(x in col_lower for x in ['ruptured', 'unruptured']):
                print(f"\nFound rupture info in column: {col}")
                # Try to extract
                for idx, row in metadata_df.iterrows():
                    val = str(row[col]).lower()
                    if 'unrupture' in val or 'unrupt' in val:
                        metadata_df.at[idx, 'Group'] = 'Unruptured'
                    elif 'rupture' in val and 'unrupture' not in val:
                        metadata_df.at[idx, 'Group'] = 'Ruptured'
    
    # Try to extract from tissue/status column
    if 'tissue/status' in metadata_df.columns or 'tissue_status' in metadata_df.columns:
        col_name = 'tissue/status' if 'tissue/status' in metadata_df.columns else 'tissue_status'
        print(f"\nExtracting groups from '{col_name}' column...")
        
        metadata_df['Group'] = 'Unknown'
        
        for idx, row in metadata_df.iterrows():
            status = str(row[col_name]).lower()
            if 'ruptured' in status and 'unruptured' not in status:
                metadata_df.at[idx, 'Group'] = 'Ruptured'
            elif 'unruptured' in status:
                metadata_df.at[idx, 'Group'] = 'Unruptured'
            elif 'control' in status or 'cortical artery' in status:
                metadata_df.at[idx, 'Group'] = 'Control'
    
    # Also try to infer from sample titles (CO = Control, IA = Intracranial Aneurysm)
    if 'title' in metadata_df.columns:
        print("\nInferring from sample titles (CO = Control, IA = Aneurysm)...")
        for idx, row in metadata_df.iterrows():
            title = str(row.get('title', '')).upper()
            if title.startswith('CO'):
                if 'Group' not in metadata_df.columns or metadata_df.at[idx, 'Group'] == 'Unknown':
                    metadata_df.at[idx, 'Group'] = 'Control'
            elif title.startswith('IA'):
                # Need to check if ruptured or unruptured from other columns
                if 'Group' not in metadata_df.columns or metadata_df.at[idx, 'Group'] == 'Unknown':
                    # Will be set from tissue/status if available
                    pass
    
    # Check Group column
    if 'Group' in metadata_df.columns:
        group_counts = metadata_df['Group'].value_counts()
        print(f"\n✓ Inferred groups:")
        print(group_counts.to_string())
        
        if 'Unknown' in group_counts and group_counts['Unknown'] == len(metadata_df):
            print("\n⚠️  Could not infer groups automatically")
            print("   Please manually review metadata and update Group column")
        else:
            print("\n✓ Successfully inferred sample groups!")
    else:
        print("\n⚠️  Could not create Group column")
        print("   Please manually add Group column with 'Ruptured', 'Unruptured', or 'Control'")
    
    return metadata_df

def main():
    """Main function"""
    
    # Check Downloads folder
    downloads_dir = Path("C:/Users/vsanker/Downloads")
    series_matrix_file = downloads_dir / "GSE122897_series_matrix.txt.gz"
    
    if not series_matrix_file.exists():
        print(f"❌ File not found: {series_matrix_file}")
        return
    
    # Extract metadata
    metadata_df = extract_metadata_from_series_matrix(str(series_matrix_file))
    
    if metadata_df is not None:
        # Save metadata
        output_file = "data/sc_rnaseq/GSE122897/sample_metadata.tsv"
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        metadata_df.to_csv(output_file, sep='\t', index=False)
        print(f"\n✓ Saved metadata to: {output_file}")
        
        # Also save to results
        results_file = "results/gse122897_sample_metadata.tsv"
        os.makedirs("results", exist_ok=True)
        metadata_df.to_csv(results_file, sep='\t', index=False)
        print(f"✓ Saved metadata to: {results_file}")
        
        return metadata_df
    
    return None

if __name__ == "__main__":
    import os
    metadata = main()

