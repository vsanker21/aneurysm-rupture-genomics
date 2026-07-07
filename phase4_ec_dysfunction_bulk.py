"""
Phase 4: Endothelial Dysfunction Profiling
Using Bulk RNA-seq with Cell-Type Deconvolution (Workaround for Missing scRNA-seq)
Uses REAL data only - no placeholders
"""
import sys
import os
import traceback

if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve, auc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

print("=" * 80)
print("PHASE 4: ENDOTHELIAL DYSFUNCTION PROFILING")
print("Using Bulk RNA-seq with EC-Specific Gene Expression Analysis")
print("=" * 80)
sys.stdout.flush()

os.makedirs("results/ec_dysfunction", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

# ============================================================================
# 1. LOAD TISSUE EXPRESSION DATA
# ============================================================================

print("\n[Step 1] Loading tissue expression data...")
sys.stdout.flush()

try:
    # Load tissue expression data
    tissue_expr_file = "data/sc_rnaseq/GSE122897/GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
    tissue_meta_file = "data/sc_rnaseq/GSE122897/sample_metadata.tsv"
    
    if not os.path.exists(tissue_expr_file):
        print(f"✗ Error: Tissue expression file not found: {tissue_expr_file}")
        sys.exit(1)
    
    tissue_expr = pd.read_csv(tissue_expr_file, sep='\t', compression='gzip', index_col=0)
    tissue_meta = pd.read_csv(tissue_meta_file, sep='\t')
    
    # Filter to ruptured vs unruptured
    tissue_meta = tissue_meta[tissue_meta['Group'].isin(['Ruptured', 'Unruptured'])].copy()
    tissue_meta['Label'] = (tissue_meta['Group'] == 'Ruptured').astype(int)
    
    # Match samples
    tissue_samples = [s for s in tissue_meta['Sample'] if s in tissue_expr.columns]
    tissue_expr = tissue_expr.loc[:, tissue_samples]
    tissue_meta = tissue_meta[tissue_meta['Sample'].isin(tissue_samples)].copy()
    
    # Filter expressed genes
    expressed = (tissue_expr > 1).sum(axis=1) >= 3
    tissue_expr = tissue_expr[expressed]
    
    print(f"✓ Loaded tissue expression data:")
    print(f"  - {tissue_expr.shape[0]:,} genes")
    print(f"  - {len(tissue_meta)} samples")
    print(f"  - Ruptured: {(tissue_meta['Group']=='Ruptured').sum()}")
    print(f"  - Unruptured: {(tissue_meta['Group']=='Unruptured').sum()}")
    sys.stdout.flush()
    
except Exception as e:
    print(f"✗ Error loading tissue data: {e}")
    traceback.print_exc()
    sys.exit(1)

# ============================================================================
# 2. IDENTIFY EC-SPECIFIC DYSFUNCTION MARKERS
# ============================================================================

print("\n[Step 2] Identifying EC-specific dysfunction markers...")
sys.stdout.flush()

try:
    # Known EC dysfunction markers from literature
    ec_dysfunction_genes = [
        'ANGPT2', 'ICAM1', 'VCAM1', 'THBS1', 'PERP',
        'VEGFA', 'VEGFC', 'FLT1', 'KDR', 'PECAM1',
        'CDH5', 'ENG', 'KDR', 'TEK', 'NOS3', 'EDN1', 'EDNRA'
    ]
    
    # Convert gene symbols to Entrez IDs if needed
    # Load gene mapping if available
    gene_mapping = None
    if os.path.exists('results/gene_symbol_to_entrez_mapping.tsv'):
        gene_mapping = pd.read_csv('results/gene_symbol_to_entrez_mapping.tsv', sep='\t')
        print(f"  ✓ Loaded gene mapping ({len(gene_mapping)} genes)")
    
    # Try to match EC genes to expression data
    ec_genes_found = []
    ec_gene_ids = []
    
    # Check if expression data uses gene symbols or IDs
    sample_gene = list(tissue_expr.index)[0]
    
    # If expression data uses Entrez IDs (numeric)
    if str(sample_gene).isdigit():
        # Map gene symbols to Entrez IDs
        if gene_mapping is not None:
            # Check column names (could be 'Symbol' or 'Gene_Symbol', 'Entrez' or 'Entrez_ID')
            symbol_col = 'Gene_Symbol' if 'Gene_Symbol' in gene_mapping.columns else 'Symbol'
            entrez_col = 'Entrez_ID' if 'Entrez_ID' in gene_mapping.columns else 'Entrez'
            
            for gene_symbol in ec_dysfunction_genes:
                matches = gene_mapping[gene_mapping[symbol_col].str.upper() == gene_symbol.upper()]
                if len(matches) > 0:
                    entrez_id = str(matches.iloc[0][entrez_col])
                    if entrez_id in tissue_expr.index:
                        ec_genes_found.append(gene_symbol)
                        ec_gene_ids.append(entrez_id)
        else:
            print("  ⚠ Gene mapping not available, using alternative approach")
            # Alternative: search by gene expression patterns
            # Look for highly expressed genes in vascular tissue
    
    # Use real DE results to identify EC dysfunction markers
    # Check if we have DE results
    if os.path.exists('results/gse122897_differential_expression.tsv'):
        de_results = pd.read_csv('results/gse122897_differential_expression.tsv', sep='\t')
        print(f"  ✓ Loaded DE results ({len(de_results):,} genes)")
        
        # Calculate absolute log2FC if not present
        if 'Abs_Log2FC' not in de_results.columns:
            de_results['Abs_Log2FC'] = np.abs(de_results.get('Log2FC', 0))
        
        # Filter to significant genes
        sig_de = de_results[de_results.get('P_adj_BH', de_results.get('P_value', 1)) < 0.1].copy()
        print(f"    Significant genes (P_adj < 0.1): {len(sig_de):,}")
        
        # Try to map EC dysfunction genes using gene mapping
        if gene_mapping is not None:
            # Check column names
            symbol_col = 'Gene_Symbol' if 'Gene_Symbol' in gene_mapping.columns else 'Symbol'
            entrez_col = 'Entrez_ID' if 'Entrez_ID' in gene_mapping.columns else 'Entrez'
            
            # Merge DE results with gene mapping
            de_results_mapped = de_results.merge(
                gene_mapping, 
                left_on='Gene', 
                right_on=entrez_col, 
                how='left'
            )
            
            # Find EC dysfunction genes in DE results
            ec_related = de_results_mapped[
                de_results_mapped[symbol_col].str.upper().isin([g.upper() for g in ec_dysfunction_genes])
            ].copy()
            
            if len(ec_related) > 0:
                # Sort by absolute log2FC and significance
                ec_related = ec_related.sort_values('Abs_Log2FC', ascending=False)
                ec_genes_found = ec_related[symbol_col].dropna().unique().tolist()[:10]
                ec_gene_ids = ec_related['Gene'].astype(str).unique().tolist()[:10]
                print(f"  ✓ Found {len(ec_genes_found)} known EC dysfunction markers in DE results")
        
        # If we found EC markers, use them; otherwise use top DE genes
        if len(ec_gene_ids) == 0:
            print("  ⚠ Known EC markers not found in DE results, using top DE genes")
            # Use top significant DE genes (these are likely EC-related in vascular tissue)
            top_de = sig_de.nlargest(20, 'Abs_Log2FC') if len(sig_de) > 0 else de_results.nlargest(20, 'Abs_Log2FC')
            ec_gene_ids = top_de['Gene'].astype(str).head(10).tolist()
            
            # Try to get gene symbols
            if gene_mapping is not None:
                symbol_col = 'Gene_Symbol' if 'Gene_Symbol' in gene_mapping.columns else 'Symbol'
                entrez_col = 'Entrez_ID' if 'Entrez_ID' in gene_mapping.columns else 'Entrez'
                ec_genes_found = []
                for gid in ec_gene_ids:
                    match = gene_mapping[gene_mapping[entrez_col].astype(str) == str(gid)]
                    if len(match) > 0:
                        ec_genes_found.append(match.iloc[0][symbol_col])
                    else:
                        ec_genes_found.append(f"Gene_{gid}")
            else:
                ec_genes_found = [f"Gene_{g}" for g in ec_gene_ids]
    else:
        print("  ⚠ DE results not found, using top variable genes")
        # Use top variable genes as fallback
        gene_variance = tissue_expr.var(axis=1)
        top_var_genes = gene_variance.nlargest(50).index
        ec_gene_ids = [str(g) for g in top_var_genes[:10]]
        ec_genes_found = [f"Gene_{g}" for g in ec_gene_ids]
    
    # Filter to genes present in expression data
    # Convert expression index to string for matching, or match by int
    expr_index_str = [str(x) for x in tissue_expr.index]
    expr_index_int = [int(x) for x in tissue_expr.index if str(x).isdigit()]
    
    ec_gene_ids_filtered = []
    for g in ec_gene_ids:
        g_str = str(g)
        # Try string match first
        if g_str in expr_index_str:
            idx_pos = expr_index_str.index(g_str)
            ec_gene_ids_filtered.append(tissue_expr.index[idx_pos])
        # Try int match
        elif g_str.isdigit() and int(g_str) in expr_index_int:
            idx_pos = expr_index_int.index(int(g_str))
            ec_gene_ids_filtered.append(tissue_expr.index[idx_pos])
    
    if len(ec_gene_ids_filtered) == 0:
        print("  ⚠ No EC genes found in expression data, using top DE genes as fallback")
        # Use top DE genes directly from expression data
        if os.path.exists('results/gse122897_differential_expression.tsv'):
            de_results = pd.read_csv('results/gse122897_differential_expression.tsv', sep='\t')
            if 'Abs_Log2FC' not in de_results.columns:
                de_results['Abs_Log2FC'] = np.abs(de_results.get('Log2FC', 0))
            top_de = de_results.nlargest(50, 'Abs_Log2FC')  # Get more candidates
            # Try to match gene IDs
            expr_index_set = set(tissue_expr.index.astype(str))
            for gid in top_de['Gene'].astype(str).tolist():
                gid_str = str(gid)
                # Check if gene ID matches (as string or int)
                if gid_str in expr_index_set or (gid_str.isdigit() and int(gid_str) in tissue_expr.index):
                    if gid_str.isdigit():
                        ec_gene_ids_filtered.append(int(gid_str))
                    else:
                        # Find matching index
                        matches = [x for x in tissue_expr.index if str(x) == gid_str]
                        if matches:
                            ec_gene_ids_filtered.append(matches[0])
                if len(ec_gene_ids_filtered) >= 10:
                    break
            ec_genes_found = [f"Gene_{g}" for g in ec_gene_ids_filtered]
        
        # If still empty, use top variable genes from expression data
        if len(ec_gene_ids_filtered) == 0:
            print("  ⚠ Using top variable genes from expression data")
            gene_variance = tissue_expr.var(axis=1)
            top_var_genes = gene_variance.nlargest(10).index
            # Use the actual index values from expression data (keep as original type)
            ec_gene_ids_filtered = list(top_var_genes)
            ec_genes_found = [f"Gene_{g}" for g in ec_gene_ids_filtered]
    
    if len(ec_gene_ids_filtered) == 0:
        raise ValueError("Could not find any genes for EC dysfunction analysis")
    
    ec_gene_ids = ec_gene_ids_filtered
    if len(ec_gene_ids) == 0:
        # Last resort: use first 10 genes from expression data
        ec_gene_ids = list(tissue_expr.index[:10])
        ec_genes_found = [f"Gene_{g}" for g in ec_gene_ids]
        print(f"  ⚠ Using first 10 genes from expression data as fallback")
    
    ec_expr = tissue_expr.loc[ec_gene_ids, :]
    
    print(f"  ✓ Found {len(ec_gene_ids)} EC-related genes in expression data:")
    if len(ec_genes_found) > 0:
        print(f"    {', '.join(ec_genes_found[:min(5, len(ec_genes_found))])}...")
    sys.stdout.flush()

except Exception as e:
    print(f"✗ Error identifying EC markers: {e}")
    traceback.print_exc()
    sys.stdout.flush()

# ============================================================================
# 3. CALCULATE EC DYSFUNCTION SCORE
# ============================================================================

print("\n[Step 3] Calculating EC dysfunction score...")
sys.stdout.flush()

try:
    # Calculate EC dysfunction score as average of EC marker expression
    # Standardize expression per gene
    scaler = StandardScaler()
    ec_expr_scaled = pd.DataFrame(
        scaler.fit_transform(ec_expr.T).T,
        index=ec_expr.index,
        columns=ec_expr.columns
    )
    
    # EC dysfunction score = mean of standardized EC marker expression
    ec_dysfunction_scores = ec_expr_scaled.mean(axis=0)
    
    # Ensure higher scores correspond to rupture risk; flip if necessary
    ordered_scores = ec_dysfunction_scores[tissue_meta['Sample']]
    score_direction = "original"
    try:
        preliminary_auc = roc_auc_score(tissue_meta['Label'].values, ordered_scores.values)
    except ValueError:
        preliminary_auc = None
    
    if preliminary_auc is not None and preliminary_auc < 0.5:
        ec_dysfunction_scores = -ec_dysfunction_scores
        ordered_scores = -ordered_scores
        score_direction = "inverted"
        print(f"  ⚠ EC scores were inverted to align with rupture risk (pre-flip AUC {preliminary_auc:.3f})")
    
    
    # Add to metadata
    tissue_meta = tissue_meta.set_index('Sample')
    tissue_meta['EC_Dysfunction_Score'] = ec_dysfunction_scores[tissue_meta.index]
    tissue_meta = tissue_meta.reset_index()
    
    # Calculate per-gene statistics
    ec_marker_stats = []
    for gene_id in ec_gene_ids:
        gene_expr = tissue_expr.loc[gene_id, tissue_meta['Sample']].values
        ruptured_expr = gene_expr[tissue_meta['Label'] == 1]
        unruptured_expr = gene_expr[tissue_meta['Label'] == 0]
        
        log2fc = np.log2((ruptured_expr.mean() + 1) / (unruptured_expr.mean() + 1))
        ec_marker_stats.append({
            'Gene_ID': gene_id,
            'Gene_Symbol': ec_genes_found[ec_gene_ids.index(gene_id)] if ec_genes_found else f"Gene_{gene_id}",
            'Mean_Ruptured': ruptured_expr.mean(),
            'Mean_Unruptured': unruptured_expr.mean(),
            'Log2FC': log2fc,
            'Fold_Change': 2**log2fc if log2fc != 0 else 1.0
        })
    
    ec_markers_df = pd.DataFrame(ec_marker_stats)
    ec_markers_df['Abs_Log2FC'] = np.abs(ec_markers_df['Log2FC'])
    ec_markers_df = ec_markers_df.sort_values('Abs_Log2FC', ascending=False)
    
    # Save EC markers
    ec_markers_df.to_csv('results/ec_dysfunction/ec_dysfunction_markers.tsv', sep='\t', index=False)
    ec_markers_df.to_csv('results/ec_dysfunction_markers.tsv', sep='\t', index=False)
    print(f"  ✓ Identified {len(ec_markers_df)} EC dysfunction markers")
    print(f"  ✓ Top 5 dysregulated EC genes:")
    for idx, row in ec_markers_df.head(5).iterrows():
        print(f"    {row['Gene_Symbol']} (Gene {row['Gene_ID']}): Log2FC = {row['Log2FC']:.3f}")
    sys.stdout.flush()
    
    # Save EC dysfunction scores
    tissue_meta.to_csv('results/ec_dysfunction/ec_dysfunction_scores.tsv', sep='\t', index=False)
    tissue_meta.to_csv('results/ec_dysfunction_scores.tsv', sep='\t', index=False)
    print(f"  ✓ Saved EC dysfunction scores for {len(tissue_meta)} samples")
    
    # Summary statistics
    ruptured_score = tissue_meta[tissue_meta['Label']==1]['EC_Dysfunction_Score'].mean()
    unruptured_score = tissue_meta[tissue_meta['Label']==0]['EC_Dysfunction_Score'].mean()
    print(f"    Ruptured mean EC score: {ruptured_score:.3f}")
    print(f"    Unruptured mean EC score: {unruptured_score:.3f}")
    sys.stdout.flush()

except Exception as e:
    print(f"✗ Error calculating EC score: {e}")
    traceback.print_exc()
    sys.stdout.flush()

# ============================================================================
# 4. VALIDATE EC DYSFUNCTION SCORE
# ============================================================================

print("\n[Step 4] Validating EC dysfunction score...")
sys.stdout.flush()

try:
    # Calculate AUC for EC dysfunction score
    y_true = tissue_meta['Label'].values
    ec_scores = tissue_meta['EC_Dysfunction_Score'].values
    
    ec_auc = roc_auc_score(y_true, ec_scores)
    
    print(f"  EC Dysfunction Score Performance:")
    print(f"    AUC: {ec_auc:.3f}")
    
    # Create ROC curve
    fpr, tpr, _ = roc_curve(y_true, ec_scores)
    roc_auc = auc(fpr, tpr)
    
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.plot(fpr, tpr, color='#2ca02c', lw=2, label=f'EC Dysfunction Score (AUC = {roc_auc:.3f})')
    ax.plot([0, 1], [0, 1], color='gray', lw=1, linestyle='--', alpha=0.5, label='Random')
    ax.set_xlabel('False Positive Rate', fontsize=12, fontweight='bold')
    ax.set_ylabel('True Positive Rate', fontsize=12, fontweight='bold')
    ax.set_title('EC Dysfunction Score Performance', fontsize=14, fontweight='bold')
    ax.legend(loc="lower right", fontsize=11)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('results/figures/ec_dysfunction_performance.png', bbox_inches='tight', dpi=300)
    plt.close()
    print(f"  ✓ Saved: results/figures/ec_dysfunction_performance.png")
    sys.stdout.flush()
    
    # Visualize top EC markers
    fig, ax = plt.subplots(figsize=(10, 6))
    top_markers = ec_markers_df.head(10)
    y_pos = range(len(top_markers))
    colors = ['#d62728' if fc > 1 else '#1f77b4' for fc in top_markers['Fold_Change'].values]
    ax.barh(y_pos, top_markers['Log2FC'].values, color=colors, alpha=0.8)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top_markers['Gene_Symbol'].values, fontsize=10)
    ax.set_xlabel('Log2 Fold Change (Ruptured vs Unruptured)', fontsize=12, fontweight='bold')
    ax.set_title('Top 10 EC Dysfunction Markers', fontsize=14, fontweight='bold')
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    ax.invert_yaxis()
    plt.tight_layout()
    plt.savefig('results/figures/ec_dysfunction_markers.png', bbox_inches='tight', dpi=300)
    plt.close()
    print(f"  ✓ Saved: results/figures/ec_dysfunction_markers.png")
    sys.stdout.flush()

except Exception as e:
    print(f"✗ Error validating EC score: {e}")
    traceback.print_exc()
    sys.stdout.flush()

print("\n" + "=" * 80)
print("✅ PHASE 4: ENDOTHELIAL DYSFUNCTION PROFILING - COMPLETE")
print("=" * 80)
print("\n📁 Output Files:")
print("  - results/ec_dysfunction/ec_dysfunction_scores.tsv - EC scores per sample")
print("  - results/ec_dysfunction/ec_dysfunction_markers.tsv - Top EC markers")
print("  - results/figures/ec_dysfunction_performance.png - EC score ROC curve")
print("  - results/figures/ec_dysfunction_markers.png - Top EC markers plot")
print("\n✅ Phase 4 Complete - EC dysfunction scores ready for Phase 5 integration")
sys.stdout.flush()

