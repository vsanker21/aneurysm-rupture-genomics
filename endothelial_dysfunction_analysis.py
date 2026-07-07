"""
Endothelial Dysfunction Analysis
Implements cell-type deconvolution and calculates EC dysfunction score
"""
import pandas as pd
import numpy as np
import os
from scipy import stats
from scipy.stats import ttest_ind
import warnings
warnings.filterwarnings('ignore')

# Key endothelial dysfunction markers from literature
EC_DYSFUNCTION_GENES = {
    'ANGPT2': 'Angiopoietin-2 - destabilizes EC-PC interaction',
    'ICAM1': 'Intercellular adhesion - inflammatory',
    'VCAM1': 'Vascular cell adhesion - transendothelial migration',
    'THBS1': 'Thrombospondin - angiogenic inhibitor',
    'PERP': 'p53-dependent apoptosis - EC apoptosis',
    'LYVE1': 'Lymphatic vessel marker - abnormal lymphangiogenesis',
    'CDH5': 'VE-cadherin - adherens junction integrity',
    'KDR': 'VEGFR2 - angiogenic signaling',
    'TEK': 'TIE2 - endothelial receptor tyrosine kinase',
    'CXCL8': 'IL-8 - neutrophil recruitment, inflammation',
    'VWF': 'von Willebrand factor - EC marker',
    'PECAM1': 'CD31 - EC marker',
    'ENG': 'Endoglin - EC marker'
}

# EC marker genes for identification
EC_MARKER_GENES = ['PECAM1', 'VWF', 'CDH5', 'ENG', 'KDR', 'TEK']

def load_tissue_expression():
    """Load tissue RNA-seq expression data"""
    # Try GSE122897 first (ruptured vs unruptured)
    expr_file = "data/sc_rnaseq/GSE122897/GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
    
    if not os.path.exists(expr_file):
        # Try alternative files
        alt_files = [
            "data/sc_rnaseq/GSE158558/GSE158558_expression.tsv.gz",
            "data/sc_rnaseq/GSE66238/GSE66238_expression.tsv.gz"
        ]
        for alt_file in alt_files:
            if os.path.exists(alt_file):
                expr_file = alt_file
                break
    
    if not os.path.exists(expr_file):
        print(f"[ERROR] Tissue expression file not found")
        return None, None
    
    print(f"\nLoading tissue expression data: {os.path.basename(expr_file)}")
    
    try:
        expr_df = pd.read_csv(expr_file, sep='\t', compression='gzip', index_col=0)
        print(f"  [OK] Expression matrix: {expr_df.shape[0]} genes, {expr_df.shape[1]} samples")
        
        # Load metadata
        metadata_file = "data/sc_rnaseq/GSE122897/sample_metadata.tsv"
        if os.path.exists(metadata_file):
            metadata = pd.read_csv(metadata_file, sep='\t')
            print(f"  [OK] Metadata: {metadata.shape}")
        else:
            metadata = None
            print(f"  [WARNING] Metadata not found, will infer from sample names")
        
        return expr_df, metadata
    except Exception as e:
        print(f"[ERROR] Failed to load expression data: {e}")
        return None, None

def map_gene_symbols_to_ids(gene_symbols):
    """Map gene symbols to Entrez IDs using MyGene.info API"""
    import requests
    import time
    
    mapping = {}
    url = "https://mygene.info/v3/query"
    
    for symbol in gene_symbols:
        try:
            params = {
                'q': symbol,
                'scopes': 'symbol',
                'fields': 'entrezgene',
                'species': 'human'
            }
            response = requests.get(url, params=params, timeout=10)
            if response.status_code == 200:
                data = response.json()
                if 'hits' in data and len(data['hits']) > 0:
                    entrez_id = str(data['hits'][0].get('entrezgene', ''))
                    if entrez_id:
                        mapping[symbol] = entrez_id
            time.sleep(0.3)
        except:
            continue
    
    return mapping

def simple_ec_deconvolution(expr_df, ec_markers=None):
    """
    Simple EC deconvolution using marker gene expression
    Since we don't have scRNA-seq, we'll use EC marker expression as proxy
    """
    if ec_markers is None:
        ec_markers = EC_MARKER_GENES
    
    print(f"\nEstimating endothelial cell signature using marker genes...")
    print(f"  EC markers: {', '.join(ec_markers)}")
    
    # Find available markers
    available_markers = [m for m in ec_markers if m in expr_df.index]
    
    if len(available_markers) == 0:
        # Try to find by gene ID - map symbols to IDs
        print(f"  [INFO] Marker genes not found by symbol, mapping to gene IDs...")
        symbol_to_id = map_gene_symbols_to_ids(ec_markers)
        
        # Try to find by ID
        expr_index_str = [str(idx) for idx in expr_df.index]
        for symbol, gene_id in symbol_to_id.items():
            if gene_id in expr_index_str:
                idx_pos = expr_index_str.index(gene_id)
                available_markers.append(expr_df.index[idx_pos])
        
        if len(available_markers) == 0:
            print(f"  [WARNING] Could not find EC markers in expression data")
            return None
    
    print(f"  Found {len(available_markers)} out of {len(ec_markers)} markers")
    
    # Calculate EC score as mean of marker expression
    ec_scores = expr_df.loc[available_markers].mean(axis=0)
    
    # Normalize (z-score)
    ec_scores = (ec_scores - ec_scores.mean()) / ec_scores.std()
    
    return ec_scores

def calculate_ec_dysfunction_score(expr_df, metadata=None):
    """
    Calculate EC dysfunction score based on dysfunction marker expression
    """
    print(f"\n" + "=" * 80)
    print("CALCULATING ENDOTHELIAL DYSFUNCTION SCORE")
    print("=" * 80)
    
    # Get dysfunction genes available in data
    dysfunction_genes = list(EC_DYSFUNCTION_GENES.keys())
    available_genes = [g for g in dysfunction_genes if g in expr_df.index]
    
    # If not found by symbol, try mapping to gene IDs
    gene_id_to_symbol = {}
    if len(available_genes) < 3:
        print(f"  [INFO] Only {len(available_genes)} dysfunction genes found by symbol, mapping to IDs...")
        symbol_to_id = map_gene_symbols_to_ids(dysfunction_genes)
        
        # Try to find by ID
        expr_index_str = [str(idx) for idx in expr_df.index]
        for symbol, gene_id in symbol_to_id.items():
            if symbol not in available_genes and gene_id in expr_index_str:
                idx_pos = expr_index_str.index(gene_id)
                gene_id_in_expr = expr_df.index[idx_pos]
                available_genes.append(gene_id_in_expr)
                # Map gene ID to symbol for later use
                gene_id_to_symbol[str(gene_id_in_expr)] = symbol
    
    if len(available_genes) < 3:
        print(f"[WARNING] Only {len(available_genes)} dysfunction genes found")
        print(f"  Available: {available_genes}")
        return None
    
    print(f"\nUsing {len(available_genes)} EC dysfunction markers:")
    for gene in available_genes:
        # Get symbol for display - try direct lookup first, then mapping
        gene_str = str(gene)
        if gene in EC_DYSFUNCTION_GENES:
            symbol = gene
            desc = EC_DYSFUNCTION_GENES[gene]
        elif gene_str in gene_id_to_symbol:
            symbol = gene_id_to_symbol[gene_str]
            desc = EC_DYSFUNCTION_GENES.get(symbol, 'Unknown function')
        else:
            symbol = gene_str
            desc = 'Unknown function'
        print(f"  {symbol}: {desc}")
    
    # Extract expression of dysfunction genes
    dysfunction_expr = expr_df.loc[available_genes].T  # Samples as rows
    
    # Log transform
    dysfunction_expr = np.log2(dysfunction_expr + 1)
    
    # Calculate dysfunction score as weighted sum
    # Upregulated genes (ANGPT2, ICAM1, VCAM1) contribute positively
    # Downregulated genes (CDH5, ENG) contribute negatively
    
    weights = {
        'ANGPT2': 1.0,  # Upregulated = dysfunction
        'ICAM1': 1.0,
        'VCAM1': 1.0,
        'THBS1': 0.5,
        'PERP': 0.5,
        'CXCL8': 1.0,
        'CDH5': -1.0,  # Downregulated = dysfunction
        'ENG': -1.0,
        'KDR': -0.5,
        'TEK': -0.5
    }
    
    dysfunction_score = np.zeros(len(dysfunction_expr))
    
    for gene in available_genes:
        # Get symbol for weight lookup
        gene_str = str(gene)
        if gene in EC_DYSFUNCTION_GENES:
            symbol = gene
        elif gene_str in gene_id_to_symbol:
            symbol = gene_id_to_symbol[gene_str]
        else:
            symbol = gene_str
        
        weight = weights.get(symbol, 0.5)
        if weight > 0:
            # Higher expression = more dysfunction
            dysfunction_score += weight * dysfunction_expr[gene].values
        else:
            # Lower expression = more dysfunction
            dysfunction_score -= abs(weight) * dysfunction_expr[gene].values
    
    # Normalize
    dysfunction_score = (dysfunction_score - dysfunction_score.mean()) / dysfunction_score.std()
    
    # Create results DataFrame
    results = pd.DataFrame({
        'Sample': dysfunction_expr.index,
        'EC_Dysfunction_Score': dysfunction_score
    })
    
    # Add individual gene expression (use symbol for column name)
    for gene in available_genes:
        # Get symbol for column name
        gene_str = str(gene)
        if gene in EC_DYSFUNCTION_GENES:
            symbol = gene
        elif gene_str in gene_id_to_symbol:
            symbol = gene_id_to_symbol[gene_str]
        else:
            symbol = gene_str
        results[f'EC_{symbol}'] = dysfunction_expr[gene].values
    
    return results

def compare_ec_dysfunction_by_group(expr_df, metadata, dysfunction_scores):
    """
    Compare EC dysfunction between ruptured and unruptured groups
    """
    if metadata is None or 'Group' not in metadata.columns:
        print("[WARNING] Cannot compare by group - metadata not available")
        return None
    
    print(f"\nComparing EC dysfunction by group...")
    
    # Merge with metadata
    merged = dysfunction_scores.merge(metadata[['Sample', 'Group']], on='Sample', how='left')
    
    # Get groups
    groups = merged['Group'].unique()
    if 'Ruptured' not in groups or 'Unruptured' not in groups:
        print(f"  [WARNING] Groups not found: {groups}")
        return None
    
    ruptured_scores = merged[merged['Group'] == 'Ruptured']['EC_Dysfunction_Score'].values
    unruptured_scores = merged[merged['Group'] == 'Unruptured']['EC_Dysfunction_Score'].values
    
    # Statistical test
    t_stat, p_value = ttest_ind(ruptured_scores, unruptured_scores)
    
    print(f"\nEC Dysfunction Comparison:")
    print(f"  Ruptured (n={len(ruptured_scores)}): mean={ruptured_scores.mean():.3f}")
    print(f"  Unruptured (n={len(unruptured_scores)}): mean={unruptured_scores.mean():.3f}")
    print(f"  T-statistic: {t_stat:.3f}, P-value: {p_value:.4e}")
    
    return {
        'ruptured_mean': ruptured_scores.mean(),
        'unruptured_mean': unruptured_scores.mean(),
        't_statistic': t_stat,
        'p_value': p_value,
        'effect_size': ruptured_scores.mean() - unruptured_scores.mean()
    }

def main():
    """Main EC dysfunction analysis pipeline"""
    
    print("=" * 80)
    print("ENDOTHELIAL DYSFUNCTION ANALYSIS")
    print("=" * 80)
    
    os.makedirs("results", exist_ok=True)
    
    # Step 1: Load tissue expression
    expr_df, metadata = load_tissue_expression()
    if expr_df is None:
        return
    
    # Step 2: Estimate EC signature
    ec_scores = simple_ec_deconvolution(expr_df)
    
    # Step 3: Calculate EC dysfunction score
    dysfunction_scores = calculate_ec_dysfunction_score(expr_df, metadata)
    if dysfunction_scores is None:
        return
    
    # Save results
    output_file = "results/ec_dysfunction_scores.tsv"
    dysfunction_scores.to_csv(output_file, sep='\t', index=False)
    print(f"\n[OK] Saved EC dysfunction scores to {output_file}")
    
    # Step 4: Compare by group
    comparison = compare_ec_dysfunction_by_group(expr_df, metadata, dysfunction_scores)
    
    # Step 5: Identify top dysregulated genes
    if metadata is not None and 'Group' in metadata.columns:
        print(f"\n" + "=" * 80)
        print("TOP EC DYSFUNCTION MARKERS")
        print("=" * 80)
        
        merged = dysfunction_scores.merge(metadata[['Sample', 'Group']], on='Sample', how='left')
        
        ec_genes = [col for col in dysfunction_scores.columns if col.startswith('EC_') and col != 'EC_Dysfunction_Score']
        
        gene_comparisons = []
        for gene_col in ec_genes:
            gene_name = gene_col.replace('EC_', '')
            ruptured = merged[merged['Group'] == 'Ruptured'][gene_col].values
            unruptured = merged[merged['Group'] == 'Unruptured'][gene_col].values
            
            if len(ruptured) > 0 and len(unruptured) > 0:
                t_stat, p_val = ttest_ind(ruptured, unruptured)
                log2fc = np.log2((ruptured.mean() + 1e-6) / (unruptured.mean() + 1e-6))
                
                # Get description for the gene
                desc = EC_DYSFUNCTION_GENES.get(gene_name, 'Unknown function')
                
                gene_comparisons.append({
                    'Gene': gene_name,
                    'Mean_Ruptured': ruptured.mean(),
                    'Mean_Unruptured': unruptured.mean(),
                    'Log2FC': log2fc,
                    'P_value': p_val,
                    'T_statistic': t_stat,
                    'Description': desc
                })
        
        gene_df = pd.DataFrame(gene_comparisons)
        gene_df = gene_df.sort_values('P_value')
        
        print(f"\nTop EC Dysfunction Markers:")
        print(f"{'Gene':<12} {'Log2FC':<10} {'P_value':<12} {'Description':<50}")
        print("-" * 90)
        for idx, row in gene_df.head(10).iterrows():
            desc = row.get('Description', EC_DYSFUNCTION_GENES.get(row['Gene'], ''))
            if len(desc) > 48:
                desc = desc[:45] + "..."
            print(f"{row['Gene']:<12} {row['Log2FC']:>9.3f}  {row['P_value']:>11.4e}  {desc:<50}")
        
        # Save gene comparisons
        gene_file = "results/ec_dysfunction_markers.tsv"
        gene_df.to_csv(gene_file, sep='\t', index=False)
        print(f"\n[OK] Saved EC dysfunction markers to {gene_file}")
    
    print("\n" + "=" * 80)
    print("EC DYSFUNCTION ANALYSIS COMPLETE")
    print("=" * 80)
    
    return dysfunction_scores, comparison

if __name__ == "__main__":
    main()

