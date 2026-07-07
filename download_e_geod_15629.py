"""
Download E-GEOD-15629 from ArrayExpress
19 samples (8 ruptured + 6 unruptured + 5 controls) for external validation
"""
import requests
import os
import json
from pathlib import Path
import time

print("=" * 80)
print("DOWNLOADING E-GEOD-15629 FROM ARRAYEXPRESS")
print("=" * 80)

# Create output directory
output_dir = Path("data/microarray/E-GEOD-15629")
output_dir.mkdir(parents=True, exist_ok=True)

experiment_id = "E-GEOD-15629"
base_url = "https://www.ebi.ac.uk/arrayexpress/json/v3/experiments"

print(f"\nDataset: {experiment_id}")
print(f"Description: Global expression profiling of ruptured and unruptured")
print(f"             human intracranial aneurysms")
print(f"Samples: 8 ruptured + 6 unruptured + 5 controls = 19 total")

# Try to get experiment metadata
print(f"\n[Step 1] Fetching experiment metadata...")
try:
    # ArrayExpress API v3
    api_url = f"{base_url}/{experiment_id}"
    response = requests.get(api_url, timeout=30)
    
    if response.status_code == 200:
        print(f"  ✓ API accessible")
        try:
            data = response.json()
            print(f"  ✓ Metadata retrieved")
            
            # Save metadata
            metadata_file = output_dir / "experiment_metadata.json"
            with open(metadata_file, 'w') as f:
                json.dump(data, f, indent=2)
            print(f"  ✓ Saved metadata to {metadata_file}")
            
            # Extract file information if available
            if 'experiments' in data and len(data['experiments']) > 0:
                exp = data['experiments'][0]
                if 'files' in exp:
                    files = exp['files']
                    print(f"  ✓ Found {len(files)} files available")
                    
                    # List files
                    for i, file_info in enumerate(files[:10], 1):  # Show first 10
                        filename = file_info.get('name', 'Unknown')
                        size = file_info.get('size', 'Unknown')
                        print(f"    {i}. {filename} ({size})")
                    
        except json.JSONDecodeError:
            print(f"  ⚠ API returned non-JSON response, trying alternative method")
    else:
        print(f"  ⚠ API returned status {response.status_code}")
        
except Exception as e:
    print(f"  ⚠ Error accessing API: {e}")

# Try alternative: Direct ArrayExpress website scraping or FTP
print(f"\n[Step 2] Attempting direct download...")

# ArrayExpress FTP/HTTP patterns
download_urls = [
    # ArrayExpress FTP pattern
    f"https://www.ebi.ac.uk/arrayexpress/files/{experiment_id}/{experiment_id}.raw.1.zip",
    f"https://www.ebi.ac.uk/arrayexpress/files/{experiment_id}/{experiment_id}.processed.1.zip",
    # Alternative patterns
    f"https://www.ebi.ac.uk/biostudies/files/{experiment_id}/{experiment_id}.zip",
]

# Also try accessing via website
website_url = f"https://www.ebi.ac.uk/arrayexpress/experiments/{experiment_id}/"
print(f"  Info: Visit {website_url} for manual download")
print(f"  Note: Some ArrayExpress datasets require manual download from website")

# Try to download using requests
success_count = 0
for url in download_urls:
    filename = url.split('/')[-1]
    output_path = output_dir / filename
    
    if output_path.exists():
        print(f"  ⚠ File already exists: {filename}")
        continue
    
    try:
        print(f"  Trying: {filename}...")
        response = requests.head(url, timeout=10, allow_redirects=True)
        
        if response.status_code == 200:
            print(f"    ✓ File accessible, downloading...")
            file_response = requests.get(url, stream=True, timeout=60)
            
            if file_response.status_code == 200:
                with open(output_path, 'wb') as f:
                    for chunk in file_response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                
                file_size = output_path.stat().st_size
                print(f"    ✓ Downloaded: {filename} ({file_size:,} bytes)")
                success_count += 1
            else:
                print(f"    ✗ Download failed: HTTP {file_response.status_code}")
        else:
            print(f"    ✗ File not accessible: HTTP {response.status_code}")
            
    except Exception as e:
        print(f"    ✗ Error: {e}")

# Create README with download instructions
print(f"\n[Step 3] Creating README with download instructions...")
readme_content = f"""# E-GEOD-15629 Dataset

## Dataset Information
- **Accession**: E-GEOD-15629
- **Platform**: ArrayExpress
- **Description**: Global expression profiling of ruptured and unruptured human intracranial aneurysms
- **Samples**: 19 total
  - 8 ruptured intracranial aneurysms
  - 6 unruptured intracranial aneurysms
  - 5 control intracranial arteries
- **Technology**: Microarray
- **Purpose**: External validation dataset

## Manual Download Instructions

### Option 1: ArrayExpress Website
1. Visit: https://www.ebi.ac.uk/arrayexpress/experiments/{experiment_id}/
2. Click on "Download" or "Files" tab
3. Download all available files (usually .CEL files for microarray)

### Option 2: BioStudies (ArrayExpress migrated to BioStudies)
1. Visit: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{experiment_id}
2. Download files from the Files section

### Option 3: GEO (if cross-listed)
- Some ArrayExpress datasets are also in GEO
- Check: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15629

## Expected Files
- Raw microarray data (.CEL files)
- Processed expression data (if available)
- Sample metadata
- Experiment description

## Usage
This dataset is used for external validation of rupture prediction models.

## Notes
- Microarray technology (different from RNA-seq)
- Useful for cross-platform validation
- Smaller sample size (19 samples) but valuable for validation
"""
readme_path = output_dir / "README.md"
with open(readme_path, 'w') as f:
    f.write(readme_content)
print(f"  ✓ Created README at {readme_path}")

# Summary
print(f"\n" + "=" * 80)
print("DOWNLOAD SUMMARY")
print("=" * 80)
print(f"\nDataset: {experiment_id}")
print(f"Output directory: {output_dir}")
print(f"Files downloaded: {success_count}")
print(f"\nStatus:")
if success_count > 0:
    print(f"  ✓ Successfully downloaded {success_count} file(s)")
else:
    print(f"  ⚠ No files downloaded automatically")
    print(f"  → Please download manually from:")
    print(f"     https://www.ebi.ac.uk/arrayexpress/experiments/{experiment_id}/")
    print(f"  → Or check BioStudies:")
    print(f"     https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{experiment_id}")

print(f"\n✓ README created with download instructions")


