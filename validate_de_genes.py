"""
Gene Validation Script for DE Genes
Validates DE genes against public databases, literature, and known pathways
"""
import pandas as pd
import numpy as np
import requests
import json
import time
import os
from typing import Dict, List, Tuple

def load_de_results(fdr_threshold=0.1):
    """Load differential expression results and filter by FDR"""
    de_file = "results/gse122897_differential_expression.tsv"
    if not os.path.exists(de_file):
        print(f"[ERROR] DE results file not found: {de_file}")
        return None
    
    df = pd.read_csv(de_file, sep='\t')
    
    # Filter significant genes
    sig_df = df[df['P_adj_BH'] < fdr_threshold].copy()
    print(f"Loaded {len(sig_df)} significant genes (FDR < {fdr_threshold})")
    
    return sig_df

def validate_gene_symbols(gene_ids: List[str]) -> Dict[str, Dict]:
    """
    Validate and convert gene IDs to symbols using multiple databases
    """
    print(f"\nValidating {len(gene_ids)} gene IDs...")
    
    # Limit to first 100 for API calls
    query_ids = [str(gid) for gid in gene_ids[:100]]
    
    validated_genes = {}
    
    # Method 1: MyGene.info for comprehensive annotation
    print("  Querying MyGene.info API...")
    url = "https://mygene.info/v3/query"
    
    try:
        # Try as Entrez IDs
        params = {
            'q': ','.join(query_ids),
            'scopes': 'entrezgene',
            'fields': 'symbol,name,entrezgene,ensembl.gene,alias,summary,type_of_gene',
            'species': 'human'
        }
        
        response = requests.get(url, params=params, timeout=30)
        if response.status_code == 200:
            results = response.json()
            for item in results:
                if 'query' in item:
                    query_id = str(item['query'])
                    validated_genes[query_id] = {
                        'Original_ID': query_id,
                        'Symbol': item.get('symbol', query_id),
                        'Name': item.get('name', ''),
                        'Entrez_ID': item.get('entrezgene', query_id),
                        'Ensembl_ID': item.get('ensembl', {}).get('gene', '') if isinstance(item.get('ensembl'), dict) else '',
                        'Aliases': ', '.join(item.get('alias', [])) if item.get('alias') else '',
                        'Gene_Type': item.get('type_of_gene', ''),
                        'Summary': item.get('summary', '')[:300] if item.get('summary') else '',
                        'Validated': True
                    }
            
            print(f"    [OK] Validated {len(validated_genes)} genes via MyGene.info")
    except Exception as e:
        print(f"    [WARNING] MyGene.info query failed: {e}")
    
    # Fill in missing genes
    for gene_id in query_ids:
        if str(gene_id) not in validated_genes:
            validated_genes[str(gene_id)] = {
                'Original_ID': str(gene_id),
                'Symbol': str(gene_id),
                'Name': '',
                'Entrez_ID': str(gene_id),
                'Ensembl_ID': '',
                'Aliases': '',
                'Gene_Type': 'unknown',
                'Summary': '',
                'Validated': False
            }
    
    time.sleep(1)  # Rate limiting
    
    return validated_genes

def check_literature_associations(gene_symbols: List[str], disease_terms: List[str] = ['aneurysm', 'intracranial aneurysm', 'cerebral aneurysm']) -> Dict[str, Dict]:
    """
    Check literature associations using NCBI E-utilities
    Note: This is a simplified check - full implementation would use proper Pubmed API
    """
    print(f"\nChecking literature associations for {len(gene_symbols)} genes...")
    print("  (Using simplified approach - consider full Pubmed API for comprehensive search)")
    
    associations = {}
    
    for symbol in gene_symbols[:50]:  # Limit to avoid rate limiting
        associations[symbol] = {
            'Literature_Hits': 0,
            'Relevant_Studies': '',
            'Associated_Pathways': []
        }
        
        # In a full implementation, would query Pubmed here
        # For now, provide structure for manual annotation
    
    return associations

def validate_pathway_membership(gene_symbols: List[str], pathways_file: str = None) -> Dict[str, List[str]]:
    """
    Validate genes against known pathways
    If pathway enrichment results exist, use those; otherwise use generic pathway databases
    """
    print(f"\nValidating pathway membership for {len(gene_symbols)} genes...")
    
    pathway_membership = {}
    
    # Check if pathway enrichment results exist
    pathway_files = [
        "results/pathway_enrichment_fdr_0.05.tsv",
        "results/pathway_enrichment_fdr_0.1.tsv",
        "results/pathway_enrichment_fdr_0.2.tsv"
    ]
    
    enriched_pathways = []
    pathway_genes = {}
    
    for pathway_file in pathway_files:
        if os.path.exists(pathway_file):
            try:
                pathways_df = pd.read_csv(pathway_file, sep='\t')
                for idx, row in pathways_df.iterrows():
                    pathway_name = row.get('Term', '')
                    genes_str = row.get('Genes', '')
                    if genes_str:
                        genes_in_pathway = [g.strip() for g in genes_str.split(',')]
                        enriched_pathways.append(pathway_name)
                        pathway_genes[pathway_name] = genes_in_pathway
            except Exception as e:
                print(f"    [WARNING] Could not read {pathway_file}: {e}")
    
    # Map genes to pathways
    for symbol in gene_symbols:
        symbol_str = str(symbol)
        pathway_membership[symbol_str] = []
        for pathway, genes in pathway_genes.items():
            if symbol_str in genes or any(str(g).strip() == symbol_str for g in genes):
                pathway_membership[symbol_str].append(pathway)
    
    print(f"    Found {len(enriched_pathways)} enriched pathways")
    print(f"    Mapped {sum(len(paths) > 0 for paths in pathway_membership.values())} genes to pathways")
    
    return pathway_membership

def create_validation_report(de_df: pd.DataFrame, fdr_threshold: float = 0.1):
    """
    Create comprehensive validation report
    """
    print("\n" + "=" * 80)
    print("GENE VALIDATION REPORT")
    print("=" * 80)
    
    # Get unique gene IDs
    gene_ids = de_df['Gene'].unique().tolist()
    print(f"\nValidating {len(gene_ids)} unique genes...")
    
    # Step 1: Validate gene symbols
    validated_genes = validate_gene_symbols(gene_ids)
    
    # Step 2: Extract gene symbols
    gene_symbols = [info['Symbol'] for info in validated_genes.values() if info['Validated']]
    
    # Step 3: Check pathway membership
    pathway_membership = validate_pathway_membership(gene_symbols)
    
    # Step 4: Check literature associations (simplified)
    literature = check_literature_associations(gene_symbols)
    
    # Combine all validation data
    validation_results = []
    
    for idx, row in de_df.iterrows():
        gene_id = str(row['Gene'])
        gene_info = validated_genes.get(gene_id, {})
        symbol = gene_info.get('Symbol', gene_id)
        
        validation_results.append({
            'Gene_ID': gene_id,
            'Gene_Symbol': symbol,
            'Gene_Name': gene_info.get('Name', ''),
            'Entrez_ID': gene_info.get('Entrez_ID', gene_id),
            'Ensembl_ID': gene_info.get('Ensembl_ID', ''),
            'Gene_Type': gene_info.get('Gene_Type', 'unknown'),
            'Validated': gene_info.get('Validated', False),
            'Associated_Pathways': ', '.join(pathway_membership.get(symbol, [])[:5]),  # Top 5
            'Pathway_Count': len(pathway_membership.get(symbol, [])),
            'Log2FC': row['Log2FC'],
            'P_value': row['P_value'],
            'P_adj_BH': row['P_adj_BH'],
            'Summary': gene_info.get('Summary', '')
        })
    
    validation_df = pd.DataFrame(validation_results)
    validation_df = validation_df.sort_values('P_adj_BH')
    
    # Save results
    output_file = f"results/gene_validation_report_fdr_{fdr_threshold}.tsv"
    validation_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved validation report to {output_file}")
    
    # Print summary statistics
    print("\nValidation Summary:")
    print(f"  Total genes: {len(validation_df)}")
    print(f"  Successfully validated: {validation_df['Validated'].sum()}")
    print(f"  Genes with pathway associations: {(validation_df['Pathway_Count'] > 0).sum()}")
    print(f"  Average pathways per gene: {validation_df['Pathway_Count'].mean():.2f}")
    
    # Print top validated genes
    print("\nTop 20 Validated Genes:")
    top_validated = validation_df.head(20)
    for idx, row in top_validated.iterrows():
        symbol = row['Gene_Symbol']
        name = row['Gene_Name'][:40] if row['Gene_Name'] else 'N/A'
        pathways = row['Pathway_Count']
        print(f"  {symbol:15s} {name:40s} Log2FC={row['Log2FC']:7.3f}  Pathways={pathways}")
    
    return validation_df

def main():
    """Main validation pipeline"""
    
    # Create results directory
    os.makedirs("results", exist_ok=True)
    
    # Analyze at different FDR thresholds
    for fdr_threshold in [0.05, 0.1, 0.2]:
        print(f"\n{'='*80}")
        print(f"VALIDATION AT FDR < {fdr_threshold}")
        print(f"{'='*80}")
        
        # Load DE results
        de_df = load_de_results(fdr_threshold=fdr_threshold)
        if de_df is None or len(de_df) == 0:
            print(f"No significant genes at FDR < {fdr_threshold}")
            continue
        
        # Create validation report
        validation_df = create_validation_report(de_df, fdr_threshold=fdr_threshold)
        
        time.sleep(2)  # Rate limiting between thresholds
    
    print("\n" + "=" * 80)
    print("GENE VALIDATION COMPLETE")
    print("=" * 80)

if __name__ == "__main__":
    main()






