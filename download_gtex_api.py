"""
Download GTEx vessel tissue expression data using the GTEx Portal API
Note: The direct Google Cloud Storage URL is broken (404)
"""

import requests
import json
from pathlib import Path
import time

BASE_DIR = Path("data/gtex")
BASE_DIR.mkdir(parents=True, exist_ok=True)

# GTEx Portal API base URL
API_BASE = "https://gtexportal.org/api/v2"

# Vessel tissues we need
VESSEL_TISSUES = [
    "Artery_Aorta",
    "Artery_Coronary", 
    "Artery_Tibial"
]

def get_tissue_list():
    """Get list of available tissues from GTEx API"""
    print("Fetching available tissues from GTEx API...")
    try:
        url = f"{API_BASE}/dataset/tissueInfo?datasetId=gtex_v8"
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            data = response.json()
            tissues = [t['tissueSiteDetailId'] for t in data.get('data', [])]
            print(f"  ✅ Found {len(tissues)} tissues")
            return tissues
        else:
            print(f"  ❌ Status: {response.status_code}")
            return []
    except Exception as e:
        print(f"  ❌ Error: {e}")
        return []

def get_gene_list():
    """Get list of genes (using a sample query to understand structure)"""
    print("\nFetching gene list...")
    try:
        # Try to get a sample of genes
        url = f"{API_BASE}/reference/gene?pageSize=100"
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            data = response.json()
            genes = data.get('data', [])
            print(f"  ✅ Found {len(genes)} genes (sample)")
            return genes
        else:
            print(f"  ❌ Status: {response.status_code}")
            return []
    except Exception as e:
        print(f"  ❌ Error: {e}")
        return []

def download_tissue_expression(tissue_id, dataset_id="gtex_v8"):
    """Download median gene expression for a specific tissue"""
    print(f"\nDownloading expression data for {tissue_id}...")
    
    # First, get a list of genes to query
    # The API requires gencodeId, so we need to get gene IDs first
    try:
        # Get reference genes
        gene_url = f"{API_BASE}/reference/gene?pageSize=1000"
        gene_response = requests.get(gene_url, timeout=60)
        
        if gene_response.status_code != 200:
            print(f"  ❌ Could not fetch gene list: {gene_response.status_code}")
            return None
        
        genes = gene_response.json().get('data', [])
        if not genes:
            print(f"  ❌ No genes found")
            return None
        
        print(f"  Found {len(genes)} genes, downloading expression for {tissue_id}...")
        
        # For bulk download, we might need to query in batches
        # Let's try getting expression for all genes at once
        # Note: This might be rate-limited, so we'll do it in batches
        
        all_expression = []
        batch_size = 100
        
        for i in range(0, len(genes), batch_size):
            batch = genes[i:i+batch_size]
            gencode_ids = [g.get('gencodeId') for g in batch if g.get('gencodeId')]
            
            if not gencode_ids:
                continue
            
            # Query expression for this batch
            expr_url = f"{API_BASE}/expression/medianGeneExpression"
            params = {
                "gencodeId": ",".join(gencode_ids[:10]),  # API might limit batch size
                "tissueSiteDetailId": tissue_id,
                "datasetId": dataset_id
            }
            
            try:
                expr_response = requests.get(expr_url, params=params, timeout=60)
                if expr_response.status_code == 200:
                    expr_data = expr_response.json().get('data', [])
                    all_expression.extend(expr_data)
                    print(f"    Batch {i//batch_size + 1}: {len(expr_data)} genes")
                else:
                    print(f"    Batch {i//batch_size + 1}: Status {expr_response.status_code}")
                
                time.sleep(0.5)  # Rate limiting
                
            except Exception as e:
                print(f"    Batch {i//batch_size + 1}: Error {e}")
                continue
        
        if all_expression:
            output_file = BASE_DIR / f"GTEx_v8_{tissue_id}_expression.json"
            with open(output_file, 'w', encoding='utf-8') as f:
                json.dump(all_expression, f, indent=2)
            print(f"  ✅ Saved {len(all_expression)} gene expression values to {output_file.name}")
            return output_file
        else:
            print(f"  ❌ No expression data retrieved")
            return None
            
    except Exception as e:
        print(f"  ❌ Error: {e}")
        return None

def main():
    print("="*80)
    print("GTEx Vessel Tissue Expression Download via API")
    print("="*80)
    print("\nNote: Direct Google Cloud Storage URL is broken (404)")
    print("Using GTEx Portal API instead\n")
    
    # Check available tissues
    tissues = get_tissue_list()
    
    # Filter for vessel tissues
    vessel_tissues_available = [t for t in VESSEL_TISSUES if t in tissues]
    
    if not vessel_tissues_available:
        print(f"\n⚠️  Warning: None of the requested vessel tissues found in API")
        print(f"   Requested: {VESSEL_TISSUES}")
        print(f"   Available tissues (first 10): {tissues[:10]}")
        print("\n   Trying to download with exact tissue IDs anyway...")
        vessel_tissues_available = VESSEL_TISSUES
    
    # Download expression for each vessel tissue
    downloaded_files = []
    for tissue in vessel_tissues_available:
        result = download_tissue_expression(tissue)
        if result:
            downloaded_files.append(result)
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    
    if downloaded_files:
        print(f"\n✅ Successfully downloaded {len(downloaded_files)} tissue expression files:")
        for f in downloaded_files:
            size = f.stat().st_size
            print(f"   {f.name} ({size:,} bytes)")
    else:
        print("\n❌ No files downloaded")
        print("\nAlternative options:")
        print("1. Visit GTEx Portal: https://gtexportal.org/home/datasets")
        print("2. Use GTEx Portal web interface (requires JavaScript)")
        print("3. Apply for dbGaP access: https://www.ncbi.nlm.nih.gov/gap/")
        print("4. Check GTEx GitHub: https://github.com/broadinstitute/gtex-pipeline")
    
    print("\n" + "="*80)

if __name__ == "__main__":
    main()






