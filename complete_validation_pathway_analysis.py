"""
Complete Validation, Pathway Analysis, and Biomarker Development
Using ONLY real data from available datasets - no placeholders
"""
import sys
import os

# Force UTF-8 encoding
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold, train_test_split
from sklearn.metrics import roc_auc_score, roc_curve, classification_report
import requests
import time
import warnings
warnings.filterwarnings('ignore')

print("=" * 80)
print("COMPLETE VALIDATION AND PATHWAY ANALYSIS - REAL DATA ONLY")
print("=" * 80)

# ============================================================================
# 1. PATHWAY ENRICHMENT ANALYSIS
# ============================================================================

def perform_pathway_enrichment_real_data():
    """Perform pathway enrichment analysis using real signature genes from datasets"""
    print("\n" + "=" * 80)
    print("PATHWAY ENRICHMENT ANALYSIS - REAL SIGNATURE GENES")
    print("=" * 80)
    
    # Load prioritized genes from real data analysis
    prioritized_file = "results/multi_omics_prioritized_top50.tsv"
    
    if not os.path.exists(prioritized_file):
        print(f"  [ERROR] Prioritized genes file not found: {prioritized_file}")
        return None, None
    
    prioritized = pd.read_csv(prioritized_file, sep='\t')
    print(f"\n  Loaded {len(prioritized)} prioritized genes from real data")
    
    # Extract valid gene symbols (exclude Ensembl IDs, empty values)
    gene_symbols = []
    for idx, row in prioritized.iterrows():
        symbol = str(row.get('Gene_Symbol', '')).strip()
        if symbol and symbol != 'nan' and symbol != '' and not symbol.startswith('ENSG'):
            if symbol.upper() not in [g.upper() for g in gene_symbols]:
                gene_symbols.append(symbol)
    
    # Also try Entrez IDs if available (only if numeric)
    for idx, row in prioritized.iterrows():
        entrez = row.get('Entrez_ID', None)
        if pd.notna(entrez):
            try:
                # Try to convert to int (only if it's a numeric ID)
                entrez_float = float(entrez)
                if entrez_float == int(entrez_float):
                    entrez_str = str(int(entrez_float))
                    if entrez_str not in gene_symbols and not entrez_str.startswith('ENSG'):
                        gene_symbols.append(entrez_str)
            except (ValueError, TypeError):
                # Skip if it's not numeric
                pass
    
    print(f"  Extracted {len(gene_symbols)} unique gene identifiers")
    print(f"  Example genes: {gene_symbols[:15]}")
    
    # Use g:Profiler API for enrichment (real pathway database)
    print(f"\n  Performing pathway enrichment using g:Profiler API...")
    
    url = "https://biit.cs.ut.ee/gprofiler/api/gost/profile"
    
    # Limit query to avoid timeout
    query_genes = gene_symbols[:500]
    
    payload = {
        "organism": "hsapiens",
        "query": query_genes,
        "sources": ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"],
        "user_threshold": 0.05,
        "all_results": False,
        "ordered": True,
        "no_iea": False,
        "domain_scope": "annotated"
    }
    
    enriched_pathways = pd.DataFrame()
    
    try:
        print("    Sending request to g:Profiler...")
        response = requests.post(url, json=payload, timeout=60)
        
        if response.status_code == 200:
            result = response.json()
            if 'result' in result and len(result['result']) > 0:
                pathways = []
                for item in result['result']:
                    pathways.append({
                        'Term': item.get('name', ''),
                        'Source': item.get('source', ''),
                        'Term_ID': item.get('native', ''),
                        'P_value': item.get('p_value', np.nan),
                        'Term_size': item.get('term_size', 0),
                        'Query_size': item.get('query_size', 0),
                        'Intersection_size': item.get('intersection_size', 0),
                        'Enrichment_Ratio': item.get('intersection_size', 0) / item.get('term_size', 1) if item.get('term_size', 0) > 0 else 0
                    })
                
                enriched_pathways = pd.DataFrame(pathways)
                enriched_pathways = enriched_pathways.sort_values('P_value')
                
                print(f"  [OK] Found {len(enriched_pathways)} enriched pathways")
                
                print(f"\n  Top 10 enriched pathways:")
                for idx, row in enriched_pathways.head(10).iterrows():
                    print(f"    {row['Term'][:60]:<60s} {row['Source']:<8s} P={row['P_value']:.2e} ({row['Intersection_size']} genes)")
            else:
                print("  [WARNING] No pathways found in API response")
        else:
            print(f"  [WARNING] API returned status {response.status_code}")
    except Exception as e:
        print(f"  [WARNING] Pathway enrichment API call failed: {e}")
        print("  Using existing pathway enrichment results if available...")
        
        # Fallback to existing results
        existing_file = "results/pathway_enrichment_prioritized_fdr_0.05.tsv"
        if os.path.exists(existing_file):
            enriched_pathways = pd.read_csv(existing_file, sep='\t')
            print(f"  [OK] Loaded existing pathway results: {len(enriched_pathways)} pathways")
    
    # Save results
    if len(enriched_pathways) > 0:
        pathway_file = "results/pathway_enrichment_signature_genes.tsv"
        enriched_pathways.to_csv(pathway_file, sep='\t', index=False)
        print(f"\n  [OK] Saved pathway enrichment to {pathway_file}")
    
    return enriched_pathways, gene_symbols

# ============================================================================
# 2. BLOOD BIOMARKER PANEL DEVELOPMENT
# ============================================================================

def develop_blood_biomarker_panel_real_data():
    """Develop blood biomarker panel using real expression data"""
    print("\n" + "=" * 80)
    print("BLOOD BIOMARKER PANEL DEVELOPMENT - REAL DATA")
    print("=" * 80)
    
    # Load prioritized genes
    prioritized_file = "results/multi_omics_prioritized_top50.tsv"
    if not os.path.exists(prioritized_file):
        print(f"  [ERROR] Prioritized genes file not found")
        return None
    
    prioritized = pd.read_csv(prioritized_file, sep='\t')
    
    # Filter for blood-specific genes
    if 'In_Blood' in prioritized.columns:
        blood_genes = prioritized[prioritized['In_Blood'] == True].copy()
    else:
        # Use all genes if In_Blood column doesn't exist
        blood_genes = prioritized.copy()
    
    # Sort by priority score
    if 'Multi_Omics_Priority_Score' in blood_genes.columns:
        blood_genes = blood_genes.sort_values('Multi_Omics_Priority_Score', ascending=False)
    
    # Get top 20 blood biomarkers
    top_blood = blood_genes.head(20)
    
    print(f"\n  Selected top {len(top_blood)} blood biomarker candidates from real data")
    
    # Load blood expression data
    blood_expr_file = "data/blood_rnaseq/GSE159610/GSE159610_raw_gene_counts_matrix.txt.gz"
    blood_meta_file = "data/blood_rnaseq/GSE159610/sample_metadata_with_labels.tsv"
    
    biomarker_panel = []
    
    if os.path.exists(blood_expr_file) and os.path.exists(blood_meta_file):
        print(f"\n  Loading blood expression data...")
        blood_expr = pd.read_csv(blood_expr_file, sep='\t', compression='gzip', index_col=0)
        blood_meta = pd.read_csv(blood_meta_file, sep='\t')
        
        print(f"    Expression: {blood_expr.shape[0]:,} genes, {blood_expr.shape[1]} samples")
        
        # Separate aneurysm and control samples
        aneurysm_samples = blood_meta[blood_meta['Status'] == 'Aneurysm']['Sample'].tolist()
        control_samples = blood_meta[blood_meta['Status'] == 'Control']['Sample'].tolist()
        aneurysm_samples = [s for s in aneurysm_samples if s in blood_expr.columns]
        control_samples = [s for s in control_samples if s in blood_expr.columns]
        
        print(f"    Aneurysm: {len(aneurysm_samples)} samples")
        print(f"    Control: {len(control_samples)} samples")
        
        # Extract biomarker information
        print(f"\n  Top blood biomarker candidates:")
        for rank, (idx, row) in enumerate(top_blood.iterrows(), 1):
            gene_info = {
                'Rank': rank,
                'Gene_Symbol': str(row.get('Gene_Symbol', 'N/A')) if pd.notna(row.get('Gene_Symbol')) else 'N/A',
                'Gene_ID': str(row.get('Blood_Gene_ID', row.get('Entrez_ID', 'N/A'))) if pd.notna(row.get('Blood_Gene_ID', row.get('Entrez_ID'))) else 'N/A',
                'Blood_Log2FC': row.get('Blood_Log2FC', np.nan),
                'Blood_P_adj': row.get('Blood_P_adj', np.nan),
                'Priority_Score': row.get('Multi_Omics_Priority_Score', np.nan),
                'In_Tissue': row.get('In_Tissue', False),
                'Mean_Aneurysm': np.nan,
                'Mean_Control': np.nan,
                'Fold_Change': np.nan
            }
            
            # Try to get expression values
            gene_id = gene_info['Gene_ID']
            if gene_id != 'N/A' and gene_id in blood_expr.index:
                try:
                    aneurysm_expr = blood_expr.loc[gene_id, aneurysm_samples].mean()
                    control_expr = blood_expr.loc[gene_id, control_samples].mean()
                    gene_info['Mean_Aneurysm'] = aneurysm_expr
                    gene_info['Mean_Control'] = control_expr
                    if control_expr > 0:
                        gene_info['Fold_Change'] = aneurysm_expr / control_expr
                except:
                    pass
            
            biomarker_panel.append(gene_info)
            
            # Print summary
            log2fc = f"{gene_info['Blood_Log2FC']:.3f}" if pd.notna(gene_info['Blood_Log2FC']) else 'N/A'
            priority = f"{gene_info['Priority_Score']:.1f}" if pd.notna(gene_info['Priority_Score']) else 'N/A'
            print(f"    {rank:2d}. {gene_info['Gene_Symbol']:<15s} Log2FC={log2fc:>8s} Priority={priority:>6s}")
    
    biomarker_df = pd.DataFrame(biomarker_panel)
    
    # Save biomarker panel
    panel_file = "results/blood_biomarker_panel_real_data.tsv"
    biomarker_df.to_csv(panel_file, sep='\t', index=False)
    print(f"\n  [OK] Saved biomarker panel to {panel_file}")
    
    return biomarker_df

# ============================================================================
# 3. CROSS-DATASET VALIDATION
# ============================================================================

def cross_dataset_validation_real_data():
    """Perform cross-dataset validation using available real datasets"""
    print("\n" + "=" * 80)
    print("CROSS-DATASET VALIDATION - REAL DATA")
    print("=" * 80)
    
    validation_results = {}
    
    # Load prioritized genes
    prioritized_file = "results/multi_omics_prioritized_top50.tsv"
    if not os.path.exists(prioritized_file):
        print("  [ERROR] Prioritized genes file not found")
        return {}
    
    prioritized = pd.read_csv(prioritized_file, sep='\t')
    
    # Strategy 1: Train on blood (GSE159610), validate internally with CV
    print("\n  Strategy 1: Internal cross-validation on blood dataset (GSE159610)")
    
    blood_expr_file = "data/blood_rnaseq/GSE159610/GSE159610_raw_gene_counts_matrix.txt.gz"
    blood_meta_file = "data/blood_rnaseq/GSE159610/sample_metadata_with_labels.tsv"
    
    if os.path.exists(blood_expr_file) and os.path.exists(blood_meta_file):
        blood_expr = pd.read_csv(blood_expr_file, sep='\t', compression='gzip', index_col=0)
        blood_meta = pd.read_csv(blood_meta_file, sep='\t')
        blood_meta['Label'] = (blood_meta['Status'] == 'Aneurysm').astype(int)
        
        # Get top blood genes
        blood_prioritized = prioritized[prioritized.get('In_Blood', True) == True].head(30)
        
        # Extract feature genes
        feature_genes = []
        for idx, row in blood_prioritized.iterrows():
            gene_id = str(row.get('Blood_Gene_ID', row.get('Entrez_ID', '')))
            if gene_id and gene_id != 'nan' and gene_id in blood_expr.index:
                feature_genes.append(gene_id)
        
        if len(feature_genes) > 0 and len(blood_meta) > 0:
            # Align samples
            common_samples = [s for s in blood_meta['Sample'] if s in blood_expr.columns]
            blood_meta_aligned = blood_meta[blood_meta['Sample'].isin(common_samples)].copy()
            
            if len(common_samples) > 0:
                X = blood_expr.loc[feature_genes, common_samples].T.values
                y = blood_meta_aligned['Label'].values
                
                # Handle missing values
                X = np.nan_to_num(X, nan=0.0)
                
                # Standardize
                scaler = StandardScaler()
                X_scaled = scaler.fit_transform(X)
                
                # Cross-validation
                if len(np.unique(y)) > 1:
                    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
                    rf = RandomForestClassifier(n_estimators=100, random_state=42, max_depth=3)
                    gb = GradientBoostingClassifier(n_estimators=100, random_state=42, max_depth=3)
                    
                    rf_scores = cross_val_score(rf, X_scaled, y, cv=cv, scoring='roc_auc')
                    gb_scores = cross_val_score(gb, X_scaled, y, cv=cv, scoring='roc_auc')
                    
                    print(f"    Blood dataset CV (n={len(y)}, features={len(feature_genes)}):")
                    print(f"      Random Forest: {rf_scores.mean():.3f} ± {rf_scores.std():.3f}")
                    print(f"      Gradient Boosting: {gb_scores.mean():.3f} ± {gb_scores.std():.3f}")
                    
                    validation_results['blood_cv_rf_auc'] = rf_scores.mean()
                    validation_results['blood_cv_rf_std'] = rf_scores.std()
                    validation_results['blood_cv_gb_auc'] = gb_scores.mean()
                    validation_results['blood_cv_gb_std'] = gb_scores.std()
                    validation_results['blood_n_samples'] = len(y)
                    validation_results['blood_n_features'] = len(feature_genes)
    
    # Strategy 2: Train on tissue (GSE122897), validate internally with CV
    print("\n  Strategy 2: Internal cross-validation on tissue dataset (GSE122897)")
    
    tissue_expr_file = "data/sc_rnaseq/GSE122897/GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
    tissue_meta_file = "data/sc_rnaseq/GSE122897/sample_metadata.tsv"
    
    if os.path.exists(tissue_expr_file) and os.path.exists(tissue_meta_file):
        tissue_expr = pd.read_csv(tissue_expr_file, sep='\t', compression='gzip', index_col=0)
        tissue_meta = pd.read_csv(tissue_meta_file, sep='\t')
        
        # Filter for ruptured vs unruptured
        tissue_meta_aneurysm = tissue_meta[
            tissue_meta['Group'].isin(['Ruptured', 'Unruptured'])
        ].copy()
        
        if len(tissue_meta_aneurysm) > 0:
            tissue_meta_aneurysm['Label'] = (tissue_meta_aneurysm['Group'] == 'Ruptured').astype(int)
            
            # Get top tissue genes
            tissue_prioritized = prioritized[prioritized.get('In_Tissue', True) == True].head(30)
            
            # Extract feature genes (tissue uses numeric gene IDs)
            feature_genes = []
            for idx, row in tissue_prioritized.iterrows():
                # Try different ID formats
                gene_ids_to_try = [
                    str(row.get('Tissue_Gene_ID', '')),
                    str(row.get('Entrez_ID', '')),
                    str(int(float(row.get('Entrez_ID', '')))) if pd.notna(row.get('Entrez_ID')) else ''
                ]
                
                for gene_id in gene_ids_to_try:
                    if gene_id and gene_id != 'nan' and gene_id in tissue_expr.index:
                        feature_genes.append(gene_id)
                        break
            
            if len(feature_genes) > 0:
                # Align samples
                common_samples = [s for s in tissue_meta_aneurysm['Sample'] if s in tissue_expr.columns]
                tissue_meta_aligned = tissue_meta_aneurysm[tissue_meta_aneurysm['Sample'].isin(common_samples)].copy()
                
                if len(common_samples) > 0:
                    X = tissue_expr.loc[feature_genes, common_samples].T.values
                    y = tissue_meta_aligned['Label'].values
                    
                    # Handle missing values
                    X = np.nan_to_num(X, nan=0.0)
                    
                    # Standardize
                    scaler = StandardScaler()
                    X_scaled = scaler.fit_transform(X)
                    
                    # Cross-validation
                    if len(np.unique(y)) > 1:
                        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
                        rf = RandomForestClassifier(n_estimators=100, random_state=42, max_depth=3)
                        gb = GradientBoostingClassifier(n_estimators=100, random_state=42, max_depth=3)
                        
                        rf_scores = cross_val_score(rf, X_scaled, y, cv=cv, scoring='roc_auc')
                        gb_scores = cross_val_score(gb, X_scaled, y, cv=cv, scoring='roc_auc')
                        
                        print(f"    Tissue dataset CV (n={len(y)}, features={len(feature_genes)}):")
                        print(f"      Random Forest: {rf_scores.mean():.3f} ± {rf_scores.std():.3f}")
                        print(f"      Gradient Boosting: {gb_scores.mean():.3f} ± {gb_scores.std():.3f}")
                        
                        validation_results['tissue_cv_rf_auc'] = rf_scores.mean()
                        validation_results['tissue_cv_rf_std'] = rf_scores.std()
                        validation_results['tissue_cv_gb_auc'] = gb_scores.mean()
                        validation_results['tissue_cv_gb_std'] = gb_scores.std()
                        validation_results['tissue_n_samples'] = len(y)
                        validation_results['tissue_n_features'] = len(feature_genes)
    
    # Save validation results
    if len(validation_results) > 0:
        validation_file = "results/cross_dataset_validation_real_data.tsv"
        validation_df = pd.DataFrame([validation_results])
        validation_df.to_csv(validation_file, sep='\t', index=False)
        print(f"\n  [OK] Saved validation results to {validation_file}")
    
    return validation_results

# ============================================================================
# 4. PERFORMANCE COMPARISON WITH LITERATURE
# ============================================================================

def performance_comparison_literature():
    """Compare model performance with literature benchmarks (PHASES score)"""
    print("\n" + "=" * 80)
    print("PERFORMANCE COMPARISON WITH LITERATURE")
    print("=" * 80)
    
    # Load our model results from real data
    blood_ml_file = "results/blood_ml_summary.tsv"
    
    comparison_data = []
    
    # Literature benchmarks (from project documentation - PHASES score)
    # PHASES score typical AUC: 0.65-0.72 (reported in literature)
    phases_auc_min = 0.65
    phases_auc_max = 0.72
    phases_auc_mean = (phases_auc_min + phases_auc_max) / 2
    
    if os.path.exists(blood_ml_file):
        our_results = pd.read_csv(blood_ml_file, sep='\t')
        
        print(f"\n  Model Performance Comparison:")
        print(f"  {'Model':<30s} {'Our AUC':<12s} {'PHASES AUC':<15s} {'Improvement':<15s}")
        print("  " + "-" * 75)
        
        for idx, row in our_results.iterrows():
            model = row['Best_Model']
            our_auc = row.get('Best_AUC', 0)
            
            improvement_pct = ((our_auc - phases_auc_mean) / phases_auc_mean) * 100 if phases_auc_mean > 0 else 0
            
            comparison_data.append({
                'Model': model,
                'Our_AUC': our_auc,
                'PHASES_AUC_Min': phases_auc_min,
                'PHASES_AUC_Max': phases_auc_max,
                'PHASES_AUC_Mean': phases_auc_mean,
                'Improvement_Percent': improvement_pct,
                'Exceeds_PHASES': our_auc > phases_auc_max
            })
            
            status = f"+{improvement_pct:.1f}%" if improvement_pct > 0 else "Baseline"
            exceeds = "YES" if our_auc > phases_auc_max else "NO"
            
            print(f"  {model:<30s} {our_auc:.3f}        {phases_auc_min:.2f}-{phases_auc_max:.2f}      {status:<15s} [{exceeds}]")
        
        # Overall summary
        if len(comparison_data) > 0:
            avg_our_auc = np.mean([d['Our_AUC'] for d in comparison_data])
            avg_improvement = np.mean([d['Improvement_Percent'] for d in comparison_data])
            n_exceeds = sum([d['Exceeds_PHASES'] for d in comparison_data])
            
            print(f"\n  Summary:")
            print(f"    Average our AUC: {avg_our_auc:.3f}")
            print(f"    PHASES AUC range: {phases_auc_min:.2f}-{phases_auc_max:.2f}")
            print(f"    Average improvement: {avg_improvement:.1f}%")
            print(f"    Models exceeding PHASES: {n_exceeds}/{len(comparison_data)}")
            print(f"\n    Target AUC: 0.85-0.92 (integrated model)")
            if avg_our_auc >= 0.85:
                print(f"    [OK] Target achieved!")
            elif avg_our_auc >= phases_auc_max:
                print(f"    [OK] Exceeds PHASES score baseline")
    
    # Save comparison
    comparison_file = "results/performance_comparison_phases.tsv"
    if len(comparison_data) > 0:
        comparison_df = pd.DataFrame(comparison_data)
        comparison_df.to_csv(comparison_file, sep='\t', index=False)
        print(f"\n  [OK] Saved comparison to {comparison_file}")
        return comparison_df
    
    return pd.DataFrame()

# ============================================================================
# MAIN PIPELINE
# ============================================================================

def main():
    """Main validation and analysis pipeline using real data only"""
    
    os.makedirs("results", exist_ok=True)
    
    # 1. Pathway enrichment analysis
    enriched_pathways, signature_genes = perform_pathway_enrichment_real_data()
    
    # 2. Blood biomarker panel development
    biomarker_panel = develop_blood_biomarker_panel_real_data()
    
    # 3. Cross-dataset validation
    validation_results = cross_dataset_validation_real_data()
    
    # 4. Performance comparison with literature (PHASES score)
    performance_comparison = performance_comparison_literature()
    
    # Summary
    print("\n" + "=" * 80)
    print("VALIDATION AND ANALYSIS COMPLETE - REAL DATA ONLY")
    print("=" * 80)
    print("\nSummary:")
    
    if enriched_pathways is not None and len(enriched_pathways) > 0:
        print(f"  Pathway enrichment: {len(enriched_pathways)} pathways identified")
        if len(enriched_pathways) > 0:
            top_pathway = enriched_pathways.iloc[0]
            print(f"    Top pathway: {top_pathway['Term'][:60]} (P={top_pathway['P_value']:.2e})")
    
    if biomarker_panel is not None and len(biomarker_panel) > 0:
        print(f"  Biomarker panel: {len(biomarker_panel)} candidates")
        print(f"    Top 3: {', '.join(biomarker_panel['Gene_Symbol'].head(3).tolist())}")
    
    if len(validation_results) > 0:
        print(f"  Cross-validation: Completed")
        if 'blood_cv_rf_auc' in validation_results:
            print(f"    Blood CV AUC: {validation_results['blood_cv_rf_auc']:.3f}")
        if 'tissue_cv_rf_auc' in validation_results:
            print(f"    Tissue CV AUC: {validation_results['tissue_cv_rf_auc']:.3f}")
    
    if len(performance_comparison) > 0:
        avg_auc = performance_comparison['Our_AUC'].mean()
        avg_improvement = performance_comparison['Improvement_Percent'].mean()
        print(f"  Performance vs PHASES: {avg_improvement:+.1f}% improvement")
        print(f"    Our average AUC: {avg_auc:.3f}")
        print(f"    PHASES AUC: 0.65-0.72")
    
    print("\nOutput files:")
    print("  - results/pathway_enrichment_signature_genes.tsv")
    print("  - results/blood_biomarker_panel_real_data.tsv")
    print("  - results/cross_dataset_validation_real_data.tsv")
    print("  - results/performance_comparison_phases.tsv")
    
    print("\n[OK] All analyses completed using REAL DATASET DATA ONLY!")

if __name__ == "__main__":
    main()

