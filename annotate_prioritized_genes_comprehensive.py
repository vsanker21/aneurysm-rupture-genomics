"""
Comprehensive Gene Annotation and Mapping for Prioritized Genes
Uses multiple methods to get gene symbols, names, and functions
"""
import pandas as pd
import numpy as np
import requests
import time
import os

def annotate_genes_mygene_api(gene_ids, batch_size=50):
    """
    Annotate genes using MyGene.info API with batching
    """
    print(f"\nAnnotating {len(gene_ids)} genes using MyGene.info API...")
    print(f"  Using batch size: {batch_size}")
    
    url = "https://mygene.info/v3/query"
    annotations = {}
    
    # Process in batches
    for batch_start in range(0, len(gene_ids), batch_size):
        batch_end = min(batch_start + batch_size, len(gene_ids))
        batch_ids = [str(gid).rstrip('.0') for gid in gene_ids[batch_start:batch_end]]
        
        print(f"  Processing batch {batch_start//batch_size + 1}: genes {batch_start+1}-{batch_end}")
        
        # Try multiple ID types
        id_types = ['entrezgene', 'ensembl.gene', 'symbol']
        
        for id_type in id_types:
            try:
                params = {
                    'q': ','.join(batch_ids),
                    'scopes': id_type,
                    'fields': 'symbol,name,entrezgene,ensembl.gene,alias,summary,type_of_gene,go,pathway,omim,uniprot',
                    'species': 'human'
                }
                
                response = requests.get(url, params=params, timeout=30)
                if response.status_code == 200:
                    results = response.json()
                    for item in results:
                        if 'query' in item:
                            query_id = str(item['query']).rstrip('.0')
                            
                            # Extract GO terms
                            go_terms = []
                            if 'go' in item and isinstance(item['go'], dict):
                                for go_type in ['BP', 'CC', 'MF']:
                                    if go_type in item['go']:
                                        terms = item['go'][go_type]
                                        if isinstance(terms, list):
                                            go_terms.extend([t.get('term', '') for t in terms[:5]])
                            
                            # Extract pathways
                            pathways = []
                            if 'pathway' in item:
                                if isinstance(item['pathway'], dict):
                                    if 'kegg' in item['pathway']:
                                        pathways.extend([p.get('name', '') for p in item['pathway']['kegg'][:5]])
                                    if 'reactome' in item['pathway']:
                                        pathways.extend([p.get('name', '') for p in item['pathway']['reactome'][:5]])
                            
                            # Extract OMIM
                            omim_diseases = []
                            if 'omim' in item and isinstance(item['omim'], dict):
                                if 'name' in item['omim']:
                                    omim_diseases.append(item['omim']['name'])
                            
                            annotations[query_id] = {
                                'Gene_ID': query_id,
                                'Symbol': item.get('symbol', query_id),
                                'Name': item.get('name', ''),
                                'Entrez_ID': item.get('entrezgene', query_id),
                                'Ensembl_ID': item.get('ensembl', {}).get('gene', '') if isinstance(item.get('ensembl'), dict) else '',
                                'Aliases': ', '.join(item.get('alias', [])) if item.get('alias') else '',
                                'Gene_Type': item.get('type_of_gene', ''),
                                'Summary': item.get('summary', '')[:500] if item.get('summary') else '',
                                'GO_Terms': '; '.join(go_terms[:10]),
                                'Pathways': '; '.join(pathways[:10]),
                                'OMIM_Disease': '; '.join(omim_diseases) if omim_diseases else '',
                                'UniProt_ID': item.get('uniprot', {}).get('Swiss-Prot', '') if isinstance(item.get('uniprot'), dict) else ''
                            }
                    
                    if annotations:
                        print(f"    [OK] Annotated {len([k for k in annotations.keys() if k in batch_ids])} genes using {id_type}")
                        break
            except Exception as e:
                print(f"    [WARNING] Batch annotation failed for {id_type}: {e}")
                continue
        
        time.sleep(2)  # Rate limiting between batches
    
    # Fill in missing genes
    for gene_id in gene_ids:
        gene_id_str = str(gene_id).rstrip('.0')
        if gene_id_str not in annotations:
            annotations[gene_id_str] = {
                'Gene_ID': gene_id_str,
                'Symbol': gene_id_str,
                'Name': '',
                'Entrez_ID': gene_id_str,
                'Ensembl_ID': '',
                'Aliases': '',
                'Gene_Type': 'unknown',
                'Summary': '',
                'GO_Terms': '',
                'Pathways': '',
                'OMIM_Disease': '',
                'UniProt_ID': ''
            }
    
    print(f"\n  [OK] Total annotations: {len(annotations)} genes")
    return annotations

def map_genes_to_known_loci(annotated_df):
    """
    Map genes to known IA GWAS loci based on chromosome coordinates
    Known loci from Bakker et al. 2020
    """
    print("\nMapping genes to known IA GWAS loci...")
    
    # Known IA GWAS loci (Bakker et al. 2020)
    known_loci = {
        '9p21.3': {
            'chr': '9',
            'start': 21900000,
            'end': 22100000,
            'genes': ['CDKN2B-AS1', 'CDKN2A', 'CDKN2B', 'CDKN2A-AS1'],
            'description': 'Cell cycle regulation'
        },
        '8q11.23': {
            'chr': '8',
            'start': 48000000,
            'end': 49000000,
            'genes': ['SOX17'],
            'description': 'Endothelial function'
        },
        '4q31': {
            'chr': '4',
            'start': 147000000,
            'end': 148000000,
            'genes': ['EDNRA'],
            'description': 'Endothelin signaling'
        },
        '10q24': {
            'chr': '10',
            'start': 101000000,
            'end': 102000000,
            'genes': ['CFDP1'],
            'description': 'Unknown function'
        },
        '2q13': {
            'chr': '2',
            'start': 189000000,
            'end': 190000000,
            'genes': ['COL3A1'],
            'description': 'Collagen, ECM'
        },
        '18q11': {
            'chr': '18',
            'start': 64000000,
            'end': 65000000,
            'genes': ['LOXL1'],
            'description': 'Lysyl oxidase, ECM crosslinking'
        }
    }
    
    # Check if genes match known loci genes by symbol
    annotated_df['In_Known_Locus'] = False
    annotated_df['Associated_Locus'] = ''
    annotated_df['Locus_Description'] = ''
    
    for idx, row in annotated_df.iterrows():
        symbol = str(row.get('Symbol', '')).upper()
        gene_id = str(row.get('Gene_ID', ''))
        
        for locus_name, locus_info in known_loci.items():
            locus_genes = [g.upper() for g in locus_info['genes']]
            if symbol in locus_genes or any(lg in symbol for lg in locus_genes):
                annotated_df.at[idx, 'In_Known_Locus'] = True
                annotated_df.at[idx, 'Associated_Locus'] = locus_name
                annotated_df.at[idx, 'Locus_Description'] = locus_info['description']
                print(f"  [OK] {symbol} ({gene_id}) mapped to {locus_name}")
                break
    
    mapped_count = annotated_df['In_Known_Locus'].sum()
    print(f"\n  [OK] Mapped {mapped_count} genes to known IA GWAS loci")
    
    return annotated_df

def main():
    """Main annotation pipeline"""
    
    print("=" * 80)
    print("COMPREHENSIVE GENE ANNOTATION AND MAPPING")
    print("=" * 80)
    
    # Load prioritized genes
    prioritized_file = "results/prioritized_genes_summary.tsv"
    if not os.path.exists(prioritized_file):
        print(f"[ERROR] Prioritized genes file not found: {prioritized_file}")
        return
    
    prioritized_df = pd.read_csv(prioritized_file, sep='\t')
    print(f"\nLoaded {len(prioritized_df)} prioritized genes")
    
    # Get all prioritized genes for annotation
    gene_ids = prioritized_df['Gene_ID'].astype(str).str.rstrip('.0').tolist()
    
    print(f"\nAnnotating {len(gene_ids)} prioritized genes...")
    
    # Annotate genes
    annotations = annotate_genes_mygene_api(gene_ids, batch_size=50)
    
    # Create annotated DataFrame
    annotated_list = []
    for idx, row in prioritized_df.iterrows():
        gene_id_str = str(row['Gene_ID']).rstrip('.0')
        annotation = annotations.get(gene_id_str, {})
        
        annotated_row = {
            'Gene_ID': gene_id_str,
            'Gene_Symbol': annotation.get('Symbol', gene_id_str),
            'Gene_Name': annotation.get('Name', ''),
            'Entrez_ID': annotation.get('Entrez_ID', gene_id_str),
            'Ensembl_ID': annotation.get('Ensembl_ID', ''),
            'Gene_Type': annotation.get('Gene_Type', ''),
            'Aliases': annotation.get('Aliases', ''),
            'Summary': annotation.get('Summary', ''),
            'GO_Terms': annotation.get('GO_Terms', ''),
            'Pathways': annotation.get('Pathways', ''),
            'OMIM_Disease': annotation.get('OMIM_Disease', ''),
            'UniProt_ID': annotation.get('UniProt_ID', ''),
            'Log2FC': row['Log2FC'],
            'P_adj_BH': row['P_adj_BH'],
            'Priority_Score': row['Priority_Score'],
            'In_GWAS_Locus': row.get('In_GWAS_Locus', False),
            'GWAS_Association': row.get('GWAS_Association', 'No')
        }
        annotated_list.append(annotated_row)
    
    annotated_df = pd.DataFrame(annotated_list)
    annotated_df = annotated_df.sort_values('Priority_Score', ascending=False)
    
    # Map to known loci
    annotated_df = map_genes_to_known_loci(annotated_df)
    
    # Save annotated results
    output_file = "results/prioritized_genes_annotated_comprehensive.tsv"
    annotated_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved comprehensive annotations to {output_file}")
    
    # Print top annotated genes
    print(f"\n\nTop 20 Prioritized Genes with Annotations:")
    print(f"{'Rank':<6} {'Symbol':<15} {'Name':<50} {'Log2FC':<10} {'P_adj':<12} {'Priority':<12} {'Locus':<12}")
    print("-" * 130)
    
    for idx, (i, row) in enumerate(annotated_df.head(20).iterrows(), 1):
        symbol = row.get('Gene_Symbol', str(row['Gene_ID']))
        name = row.get('Gene_Name', '')
        if len(name) > 48:
            name = name[:45] + "..."
        log2fc = row['Log2FC']
        p_adj = row['P_adj_BH']
        priority = row['Priority_Score']
        locus = row.get('Associated_Locus', '')
        
        symbol_str = symbol[:13] if isinstance(symbol, str) and len(symbol) > 13 else symbol
        name_str = name[:48] if len(name) > 48 else name
        
        log2fc_str = f"{log2fc:.3f}" if pd.notna(log2fc) else "N/A"
        p_adj_str = f"{p_adj:.4f}" if pd.notna(p_adj) else "N/A"
        priority_str = f"{priority:.2f}" if pd.notna(priority) else "N/A"
        locus_str = locus[:10] if locus else ""
        
        print(f"{idx:<6} {str(symbol_str):<15} {name_str:<50} {log2fc_str:>9}  {p_adj_str:>11}  {priority_str:>11}  {locus_str:<12}")
    
    # Summary statistics
    print(f"\n\nAnnotation Summary:")
    print(f"  Total genes annotated: {len(annotated_df)}")
    print(f"  Genes with symbols: {(annotated_df['Gene_Symbol'] != annotated_df['Gene_ID']).sum()}")
    print(f"  Genes with names: {annotated_df['Gene_Name'].notna().sum()}")
    print(f"  Genes with GO terms: {(annotated_df['GO_Terms'] != '').sum()}")
    print(f"  Genes with pathways: {(annotated_df['Pathways'] != '').sum()}")
    print(f"  Genes in known GWAS loci: {annotated_df['In_Known_Locus'].sum()}")
    
    # Genes in known loci
    if annotated_df['In_Known_Locus'].sum() > 0:
        print(f"\n  Genes in Known IA GWAS Loci:")
        loci_genes = annotated_df[annotated_df['In_Known_Locus'] == True]
        for idx, row in loci_genes.iterrows():
            symbol = row.get('Gene_Symbol', row['Gene_ID'])
            locus = row.get('Associated_Locus', '')
            desc = row.get('Locus_Description', '')
            print(f"    {symbol} -> {locus} ({desc})")
    
    print("\n" + "=" * 80)
    print("COMPREHENSIVE ANNOTATION COMPLETE")
    print("=" * 80)
    
    return annotated_df

if __name__ == "__main__":
    main()






