"""
Annotate top genes from GSE122897 analysis for literature review
This script annotates genes and creates/updates the literature search template
"""
import pandas as pd
import numpy as np
import requests
import time
import os

def annotate_genes_mygene_api(gene_ids, batch_size=50):
    """
    Annotate genes using MyGene.info API - query individually for reliability
    """
    print(f"\nAnnotating {len(gene_ids)} genes using MyGene.info API...")
    
    url = "https://mygene.info/v3/gene"
    annotations = {}
    
    # Query genes individually for better reliability
    for i, gene_id in enumerate(gene_ids, 1):
        gene_id_str = str(gene_id).rstrip('.0')
        print(f"  [{i}/{len(gene_ids)}] Querying gene {gene_id_str}...", end=' ')
        
        # Try querying by Entrez ID first
        try:
            response = requests.get(f"{url}/{gene_id_str}", params={'fields': 'symbol,name,entrezgene,ensembl.gene,alias'}, timeout=30)
            if response.status_code == 200:
                item = response.json()
                if item and 'symbol' in item:
                    annotations[gene_id_str] = {
                        'Gene_ID': gene_id_str,
                        'Symbol': item.get('symbol', gene_id_str),
                        'Name': item.get('name', ''),
                        'Entrez_ID': item.get('entrezgene', gene_id_str),
                        'Ensembl_ID': item.get('ensembl', {}).get('gene', '') if isinstance(item.get('ensembl'), dict) else '',
                        'Aliases': ', '.join(item.get('alias', [])) if item.get('alias') else ''
                    }
                    print(f"[OK] {item.get('symbol', gene_id_str)}")
                else:
                    print(f"[NOT FOUND]")
            else:
                print(f"[ERROR {response.status_code}]")
        except Exception as e:
            print(f"[ERROR: {str(e)[:30]}]")
        
        time.sleep(0.5)  # Rate limiting
    
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
                'Aliases': ''
            }
    
    print(f"\n  [OK] Total annotations: {len(annotations)} genes")
    return annotations

def create_literature_search_template(annotated_top_genes, output_file):
    """
    Create literature search template with PubMed search queries
    """
    print(f"\nCreating literature search template for top {len(annotated_top_genes)} genes...")
    
    template_rows = []
    
    for rank, (idx, row) in enumerate(annotated_top_genes.iterrows(), 1):
        gene_id = str(row['Gene'])
        symbol = row.get('Gene_Symbol', gene_id)
        name = row.get('Gene_Name', '')
        log2fc = row.get('Log2FC', 0)
        p_adj = row.get('P_adj_BH', 1.0)
        
        # Calculate priority score (higher is better)
        priority_score = abs(log2fc) * (-np.log10(max(p_adj, 1e-10)))
        
        # Create search queries
        # Query 1: Direct gene search
        search_query_1 = f'"{symbol}" OR "{gene_id}"'
        if name:
            search_query_1 += f' OR "{name}"'
        
        # Query 2: Gene + aneurysm
        search_query_2 = f'("{symbol}" OR "{gene_id}") AND (aneurysm OR intracranial aneurysm OR cerebral aneurysm OR IA)'
        
        # Query 3: Gene + vascular
        search_query_3 = f'("{symbol}" OR "{gene_id}") AND (vascular OR vessel OR artery OR endothelium OR smooth muscle)'
        
        # Query 4: Gene + stroke
        search_query_4 = f'("{symbol}" OR "{gene_id}") AND (stroke OR cerebrovascular OR subarachnoid hemorrhage OR SAH)'
        
        # Query 5: Gene + biomarker
        search_query_5 = f'("{symbol}" OR "{gene_id}") AND (biomarker OR serum OR blood OR plasma OR CSF)'
        
        template_rows.append({
            'Rank': rank,
            'Gene_ID': gene_id,
            'Gene_Symbol': symbol,
            'Gene_Name': name,
            'Log2FC': log2fc,
            'P_adj_BH': p_adj,
            'Priority_Score': priority_score,
            'Search_Query_1_Direct': search_query_1,
            'Search_Query_2_Aneurysm': search_query_2,
            'Search_Query_3_Vascular': search_query_3,
            'Search_Query_4_Stroke': search_query_4,
            'Search_Query_5_Biomarker': search_query_5
        })
    
    template_df = pd.DataFrame(template_rows)
    template_df = template_df.sort_values('Priority_Score', ascending=False)
    template_df['Rank'] = range(1, len(template_df) + 1)
    
    # Save template
    template_df.to_csv(output_file, sep='\t', index=False)
    print(f"[OK] Saved literature search template to {output_file}")
    
    return template_df

def main():
    """Main annotation pipeline"""
    
    print("=" * 80)
    print("GENE ANNOTATION AND MAPPING FOR LITERATURE REVIEW")
    print("=" * 80)
    
    # Load differential expression results to get top genes
    de_file = "results/gse122897_differential_expression.tsv"
    if not os.path.exists(de_file):
        print(f"[ERROR] Differential expression file not found: {de_file}")
        print("   Please run analyze_gse122897.py first")
        return
    
    de_df = pd.read_csv(de_file, sep='\t')
    print(f"\nLoaded {len(de_df)} genes from differential expression results")
    
    # Get top 10 most significant genes (by P_adj_BH, then by absolute Log2FC)
    de_df['Abs_Log2FC'] = de_df['Log2FC'].abs()
    top_genes_df = de_df.nsmallest(10, 'P_adj_BH').nlargest(10, 'Abs_Log2FC').head(10).copy()
    print(f"Selected top 10 genes for annotation")
    
    # Get unique gene IDs
    gene_ids = top_genes_df['Gene'].astype(str).str.rstrip('.0').unique().tolist()
    print(f"Unique genes to annotate: {len(gene_ids)}")
    print(f"Gene IDs: {gene_ids}")
    
    # Annotate genes
    annotations = annotate_genes_mygene_api(gene_ids, batch_size=50)
    
    # Merge annotations with top genes
    annotated_list = []
    for idx, row in top_genes_df.iterrows():
        gene_id_str = str(row['Gene']).rstrip('.0')
        annotation = annotations.get(gene_id_str, {})
        
        annotated_row = row.to_dict()
        annotated_row['Gene_Symbol'] = annotation.get('Symbol', gene_id_str)
        annotated_row['Gene_Name'] = annotation.get('Name', '')
        annotated_row['Entrez_ID'] = annotation.get('Entrez_ID', gene_id_str)
        annotated_row['Ensembl_ID'] = annotation.get('Ensembl_ID', '')
        annotated_row['Aliases'] = annotation.get('Aliases', '')
        annotated_list.append(annotated_row)
    
    # Print annotation status
    print(f"\nAnnotation status:")
    for row in annotated_list:
        gene_id = str(row['Gene']).rstrip('.0')
        symbol = row.get('Gene_Symbol', '')
        name = row.get('Gene_Name', '')
        print(f"  {gene_id}: {symbol} - {name}")
    
    annotated_df = pd.DataFrame(annotated_list)
    
    # Sort by priority (absolute log2FC and significance)
    annotated_df = annotated_df.sort_values(['P_adj_BH', 'Abs_Log2FC'], ascending=[True, False])
    
    # Save annotated top genes
    annotated_output = "results/gse122897_top_genes_annotated.tsv"
    annotated_df.to_csv(annotated_output, sep='\t', index=False)
    print(f"\n[OK] Saved annotated top genes to {annotated_output}")
    
    # Create literature search template for top 10
    top_10 = annotated_df.head(10).copy()
    template_file = "results/literature_search_template_top10.tsv"
    template_df = create_literature_search_template(top_10, template_file)
    
    # Print summary
    print(f"\n\nTop 10 Genes for Literature Review:")
    print(f"{'Rank':<6} {'Symbol':<15} {'Name':<50} {'Log2FC':<10} {'P_adj':<12}")
    print("-" * 100)
    
    for rank, (idx, row) in enumerate(template_df.iterrows(), 1):
        symbol = row.get('Gene_Symbol', str(row['Gene_ID']))
        name = row.get('Gene_Name', '')
        if len(name) > 48:
            name = name[:45] + "..."
        log2fc = row['Log2FC']
        p_adj = row['P_adj_BH']
        
        print(f"{rank:<6} {str(symbol):<15} {name:<50} {log2fc:>9.3f}  {p_adj:>11.4e}")
    
    print(f"\n\nAnnotation Summary:")
    print(f"  Total genes annotated: {len(annotated_df)}")
    print(f"  Genes with symbols: {(annotated_df['Gene_Symbol'] != annotated_df['Gene'].astype(str)).sum()}")
    print(f"  Genes with names: {annotated_df['Gene_Name'].notna().sum()}")
    
    print("\n" + "=" * 80)
    print("ANNOTATION COMPLETE - Ready for literature review")
    print("=" * 80)
    
    return annotated_df, template_df

if __name__ == "__main__":
    main()

