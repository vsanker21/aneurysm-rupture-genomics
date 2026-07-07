"""
Search GEO for human intracranial aneurysm single-cell datasets
"""

import requests
from bs4 import BeautifulSoup
import time

def search_geo(query, max_results=20):
    """Search GEO database"""
    url = "https://www.ncbi.nlm.nih.gov/gds"
    params = {
        'term': query,
        'report': 'Full',
        'display': '500'
    }
    
    try:
        response = requests.get(url, params=params, timeout=30)
        if response.status_code == 200:
            soup = BeautifulSoup(response.text, 'html.parser')
            # Look for dataset links
            links = soup.find_all('a', href=lambda x: x and 'GSE' in str(x))
            datasets = []
            for link in links[:max_results]:
                href = link.get('href', '')
                if 'GSE' in href and 'acc=' in href:
                    acc = href.split('acc=')[1].split('&')[0] if 'acc=' in href else None
                    if acc and acc.startswith('GSE'):
                        datasets.append(acc)
            return list(set(datasets))
    except Exception as e:
        print(f"Search error: {e}")
    return []

# Search queries
queries = [
    "human intracranial aneurysm single cell",
    "cerebral aneurysm scRNA-seq human",
    "aneurysm endothelial single cell human",
    "intracranial aneurysm tissue single cell RNA-seq"
]

print("="*80)
print("SEARCHING GEO FOR HUMAN INTRACRANIAL ANEURYSM SINGLE-CELL DATASETS")
print("="*80)

all_datasets = set()

for query in queries:
    print(f"\nSearching: '{query}'...")
    results = search_geo(query)
    if results:
        print(f"  Found {len(results)} datasets: {', '.join(results[:10])}")
        all_datasets.update(results)
    else:
        print("  No results found")
    time.sleep(1)  # Be polite to the server

print(f"\n\nTotal unique datasets found: {len(all_datasets)}")
if all_datasets:
    print(f"Datasets to check: {', '.join(sorted(all_datasets)[:20])}")






