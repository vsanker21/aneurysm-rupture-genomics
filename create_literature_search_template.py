"""
Create Literature Search Template for Top 10 Prioritized Genes
Generates PubMed search queries and literature review templates
"""
import pandas as pd
import os

def create_literature_search_template():
    """Create literature search template for top 10 prioritized genes"""
    
    print("=" * 80)
    print("CREATING LITERATURE SEARCH TEMPLATE")
    print("=" * 80)
    
    # Load prioritized genes with annotations
    merged_file = "results/prioritized_genes_annotated_merged.tsv"
    
    if not os.path.exists(merged_file):
        print(f"[ERROR] Merged genes file not found: {merged_file}")
        return None
    
    merged_df = pd.read_csv(merged_file, sep='\t')
    print(f"\nLoaded prioritized genes with annotations: {len(merged_df)} genes")
    
    # Get top 10 prioritized genes
    top_10 = merged_df.head(10).copy()
    
    print(f"\nCreating literature search template for top 10 prioritized genes...")
    
    # Create literature search template
    search_template = []
    
    for idx, (i, row) in enumerate(top_10.iterrows(), 1):
        gene_id = str(row.get('Gene_ID', 'N/A')).rstrip('.0')
        symbol = row.get('Gene_Symbol', gene_id)
        name = row.get('Gene_Name', '')
        log2fc = row.get('Log2FC', 'N/A')
        p_adj = row.get('P_adj_BH', 'N/A')
        priority = row.get('Priority_Score', 'N/A')
        
        # Create search queries
        search_queries = []
        
        # Query 1: Direct gene search
        if pd.notna(symbol) and symbol != gene_id:
            search_queries.append(f'"{symbol}" OR "{name}" OR "{gene_id}"')
        else:
            search_queries.append(f'"{gene_id}"')
        
        # Query 2: Aneurysm-specific
        if pd.notna(symbol) and symbol != gene_id:
            search_queries.append(f'("{symbol}" OR "{gene_id}") AND (aneurysm OR intracranial aneurysm OR cerebral aneurysm OR IA)')
        else:
            search_queries.append(f'"{gene_id}" AND (aneurysm OR intracranial aneurysm OR cerebral aneurysm OR IA)')
        
        # Query 3: Vascular biology
        if pd.notna(symbol) and symbol != gene_id:
            search_queries.append(f'("{symbol}" OR "{gene_id}") AND (vascular OR vessel OR artery OR endothelium OR smooth muscle)')
        else:
            search_queries.append(f'"{gene_id}" AND (vascular OR vessel OR artery OR endothelium OR smooth muscle)')
        
        # Query 4: Stroke/cerebrovascular
        if pd.notna(symbol) and symbol != gene_id:
            search_queries.append(f'("{symbol}" OR "{gene_id}") AND (stroke OR cerebrovascular OR subarachnoid hemorrhage OR SAH)')
        else:
            search_queries.append(f'"{gene_id}" AND (stroke OR cerebrovascular OR subarachnoid hemorrhage OR SAH)')
        
        # Query 5: Biomarker potential
        if pd.notna(symbol) and symbol != gene_id:
            search_queries.append(f'("{symbol}" OR "{gene_id}") AND (biomarker OR serum OR blood OR plasma OR CSF)')
        else:
            search_queries.append(f'"{gene_id}" AND (biomarker OR serum OR blood OR plasma OR CSF)')
        
        search_template.append({
            'Rank': idx,
            'Gene_ID': gene_id,
            'Gene_Symbol': symbol if pd.notna(symbol) else gene_id,
            'Gene_Name': name if pd.notna(name) else '',
            'Log2FC': log2fc,
            'P_adj_BH': p_adj,
            'Priority_Score': priority,
            'Search_Query_1_Direct': search_queries[0],
            'Search_Query_2_Aneurysm': search_queries[1],
            'Search_Query_3_Vascular': search_queries[2],
            'Search_Query_4_Stroke': search_queries[3],
            'Search_Query_5_Biomarker': search_queries[4]
        })
    
    # Save search template
    template_df = pd.DataFrame(search_template)
    output_file = "results/literature_search_template_top10.tsv"
    template_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved literature search template to {output_file}")
    
    # Create formatted output for manual use
    print(f"\n\n" + "=" * 80)
    print("LITERATURE SEARCH QUERIES FOR TOP 10 PRIORITIZED GENES")
    print("=" * 80)
    
    for item in search_template:
        print(f"\n{'='*80}")
        print(f"Rank {item['Rank']}: Gene {item['Gene_ID']} ({item['Gene_Symbol']})")
        print(f"{'='*80}")
        print(f"  Priority Score: {item['Priority_Score']}")
        print(f"  Log2FC: {item['Log2FC']}")
        print(f"  P_adj_BH: {item['P_adj_BH']}")
        print(f"  Gene Name: {item['Gene_Name']}")
        print(f"\n  PubMed Search Queries:")
        print(f"    1. Direct: {item['Search_Query_1_Direct']}")
        print(f"    2. Aneurysm: {item['Search_Query_2_Aneurysm']}")
        print(f"    3. Vascular: {item['Search_Query_3_Vascular']}")
        print(f"    4. Stroke: {item['Search_Query_4_Stroke']}")
        print(f"    5. Biomarker: {item['Search_Query_5_Biomarker']}")
        
        print(f"\n  Key Questions:")
        print(f"    - What is this gene's function?")
        print(f"    - Any associations with vascular diseases?")
        print(f"    - Expression in cerebrovascular tissue?")
        print(f"    - Biomarker potential in blood/serum?")
        print(f"    - Functional studies available?")
    
    # Create summary table
    print(f"\n\n" + "=" * 80)
    print("SUMMARY TABLE - TOP 10 PRIORITIZED GENES")
    print("=" * 80)
    print(f"{'Rank':<6} {'Gene_ID':<12} {'Symbol':<15} {'Name':<40} {'Log2FC':<10} {'P_adj':<12} {'Priority':<12}")
    print("-" * 110)
    
    for item in search_template:
        name = item['Gene_Name'][:38] if item['Gene_Name'] and len(item['Gene_Name']) > 38 else (item['Gene_Name'] or '')
        print(f"{item['Rank']:<6} {item['Gene_ID']:<12} {str(item['Gene_Symbol']):<15} {name:<40} {item['Log2FC']:<10} {item['P_adj_BH']:<12} {item['Priority_Score']:<12}")
    
    # Create PubMed URLs for easy access
    print(f"\n\n" + "=" * 80)
    print("PUBMED SEARCH URLs (Copy and paste into browser)")
    print("=" * 80)
    
    for item in search_template:
        symbol = item['Gene_Symbol']
        gene_id = item['Gene_ID']
        
        # URL-encode search query
        query_direct = f"{symbol}+OR+{gene_id}" if symbol != gene_id else gene_id
        query_aneurysm = f"({symbol}+OR+{gene_id})+AND+(aneurysm+OR+intracranial+aneurysm)" if symbol != gene_id else f"{gene_id}+AND+(aneurysm+OR+intracranial+aneurysm)"
        
        print(f"\nRank {item['Rank']}: {symbol} ({gene_id})")
        print(f"  Direct: https://pubmed.ncbi.nlm.nih.gov/?term={query_direct}")
        print(f"  Aneurysm: https://pubmed.ncbi.nlm.nih.gov/?term={query_aneurysm}")
    
    return template_df

def main():
    """Main function"""
    template_df = create_literature_search_template()
    
    print("\n" + "=" * 80)
    print("LITERATURE SEARCH TEMPLATE COMPLETE")
    print("=" * 80)
    
    if template_df is not None:
        print(f"\n[NOTE] Search template saved to results/literature_search_template_top10.tsv")
        print("  Next steps:")
        print("  1. Copy search queries to PubMed")
        print("  2. Review top 10-20 results per query")
        print("  3. Document findings in literature review table")
        print("  4. Identify functional studies for validation")

if __name__ == "__main__":
    main()






