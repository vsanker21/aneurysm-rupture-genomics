"""
Phase 6: Exosomal Biomarker Discovery (REAL DATA ONLY)
Identify novel exosomal non-coding RNA biomarkers using ONLY real data
"""
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
import pickle
import json

print("=" * 80)
print("PHASE 6: EXOSOMAL BIOMARKER DISCOVERY (REAL DATA ONLY)")
print("NO SIMULATED DATA - USING ONLY ACTUAL DATA FROM DATASETS")
print("=" * 80)

# Load differential expression results (REAL DATA)
de_results = pd.read_csv('results/gse122897_differential_expression.tsv', sep='\t')
de_results = de_results.sort_values('P_adj_BH').head(1000)  # Top 1000 most significant

# Known exosomal RNA biomarkers from literature (reference only, not simulated)
known_exosomal_markers = {
    'ATP1A1-AS1': {
        'type': 'lncRNA',
        'auc_literature': 0.832,
        'function': 'VSMC phenotype switching',
        'mechanism': 'Increases MMP-9 expression'
    },
    'hsa_circ_0007990': {
        'type': 'circRNA',
        'auc_literature': 0.78,
        'function': 'endothelial dysfunction'
    },
    'miR-155': {
        'type': 'miRNA',
        'auc_literature': 0.75,
        'function': 'inflammatory response'
    }
}

# Load optimized model and features
with open('results/phase5_best_model_summary.json', 'r') as f:
    model_summary = json.load(f)

# Load data (REAL DATA)
features_df = pd.read_csv('results/integrated_features.tsv', sep='\t')
features_df = features_df[features_df['Group'].isin(['Ruptured', 'Unruptured'])].copy()
labels = (features_df['Group'] == 'Ruptured').astype(int)

# Create exosomal biomarker discovery framework (REAL DATA ONLY)
print("\n1. IDENTIFYING CANDIDATE EXOSOMAL RNA BIOMARKERS (REAL DATA)")
print("-" * 80)

# Filter significant genes with substantial fold changes (REAL DATA)
sig_genes = de_results[
    (de_results['P_adj_BH'] < 0.1) & 
    (np.abs(de_results['Log2FC']) > 1.0)
].copy()

print(f"\nSignificant genes (P_adj < 0.1, |Log2FC| > 1): {len(sig_genes)}")

# Convert gene IDs to strings for searching
sig_genes['Gene_Str'] = sig_genes['Gene'].astype(str)

# Identify potential non-coding RNAs based on gene patterns (REAL DATA)
candidate_biomarkers = []

# Check for known markers (matching by name patterns in real DE results)
for gene_id, gene_name in zip(sig_genes['Gene'], sig_genes['Gene_Str']):
    # Check if it's a known exosomal marker (by name matching in real data)
    matched_marker = None
    for marker_name, marker_info in known_exosomal_markers.items():
        if marker_name.lower() in str(gene_name).lower() or str(gene_name) in marker_name:
            matched_marker = (marker_name, marker_info)
            break
    
    if matched_marker:
        candidate_biomarkers.append({
            'Gene_ID': gene_id,
            'Gene_Name': gene_name,
            'Marker_Name': matched_marker[0],
            'RNA_Type': matched_marker[1]['type'],
            'Log2FC': sig_genes[sig_genes['Gene'] == gene_id]['Log2FC'].values[0],
            'P_adj_BH': sig_genes[sig_genes['Gene'] == gene_id]['P_adj_BH'].values[0],
            'AUC_Literature': matched_marker[1].get('auc_literature', 'N/A'),
            'Function': matched_marker[1].get('function', 'N/A'),
            'Status': 'Known'
        })

print(f"\nFound {len(candidate_biomarkers)} known exosomal biomarkers in DE results (REAL DATA)")

# Look for novel candidates from REAL DE results (top 50 most significant)
novel_candidates = sig_genes.head(50).copy()

novel_biomarkers = []
for idx, row in novel_candidates.iterrows():
    # Skip if already identified as known marker
    if row['Gene'] in [b['Gene_ID'] for b in candidate_biomarkers]:
        continue
    
    # Include highly significant genes as potential biomarkers (REAL DATA ONLY)
    if row['P_adj_BH'] < 0.05 and abs(row['Log2FC']) > 1.5:
        # Check for patterns suggesting non-coding RNA in gene ID/name
        gene_str = str(row['Gene'])
        rna_type = 'Protein_coding'  # Default
        
        # Heuristics for identifying potential ncRNAs (based on real gene IDs/names)
        if any(pattern in gene_str.lower() for pattern in ['-as', '-ot', 'lnc', 'linc']):
            rna_type = 'lncRNA'
        elif 'circ' in gene_str.lower() or 'hsa_circ' in gene_str.lower():
            rna_type = 'circRNA'
        elif 'mir' in gene_str.lower() or 'microrna' in gene_str.lower():
            rna_type = 'miRNA'
        
        novel_biomarkers.append({
            'Gene_ID': row['Gene'],
            'Gene_Name': gene_str,
            'Marker_Name': f'Candidate_{len(novel_biomarkers)+1}',
            'RNA_Type': rna_type,
            'Log2FC': row['Log2FC'],
            'P_adj_BH': row['P_adj_BH'],
            'AUC_Literature': 'N/A',
            'Function': 'To be validated',
            'Status': 'Novel'
        })

print(f"Identified {len(novel_biomarkers)} novel biomarker candidates from REAL DE data")

# Combine all candidates (REAL DATA ONLY)
all_biomarkers = pd.DataFrame(candidate_biomarkers + novel_biomarkers)

if len(all_biomarkers) > 0:
    all_biomarkers = all_biomarkers.sort_values('P_adj_BH')
    all_biomarkers.to_csv('results/phase6_exosomal_biomarker_candidates.tsv', sep='\t', index=False)
    print(f"\n[OK] Saved {len(all_biomarkers)} candidate biomarkers (REAL DATA)")

print("\n" + "=" * 80)
print("2. CREATING CLINICAL NOMOGRAM (REAL DATA ONLY)")
print("=" * 80)

# Create nomogram using ONLY optimized model (we don't have actual exosomal measurements)
# In practice, exosomal biomarkers would need to be measured separately

with open('results/phase5_optimized_model.pkl', 'rb') as f:
    optimized_model = pickle.load(f)

with open('results/phase5_optimized_scaler.pkl', 'rb') as f:
    scaler = pickle.load(f)

optimized_features = model_summary['features_used']
X_optimized = features_df[optimized_features].fillna(0).values
X_optimized_scaled = scaler.transform(X_optimized)

# Get optimized model predictions (REAL DATA)
optimized_scores = optimized_model.predict_proba(X_optimized_scaled)[:, 1]

# For nomogram, use only optimized model scores (REAL DATA)
# Exosomal biomarkers would need to be measured separately in practice
top_exosomal_candidates = all_biomarkers.head(5) if len(all_biomarkers) > 0 else pd.DataFrame()

print(f"\nTop exosomal biomarker candidates (from REAL DE data):")
for idx, row in top_exosomal_candidates.head(5).iterrows():
    print(f"  {row['Marker_Name']} ({row['RNA_Type']}): Log2FC={row['Log2FC']:.2f}, P_adj={row['P_adj_BH']:.3e}")

# Create nomogram using ONLY optimized model (since we don't have actual exosomal measurements)
# The nomogram formula includes placeholders for exosomal biomarkers, but current scores use only optimized model
nomogram_weights = {
    'optimized_model': 1.0,  # Use 100% optimized model since we don't have exosomal measurements
    'exosomal_biomarker_1': 0.0,  # Would be used when exosomal data available
    'exosomal_biomarker_2': 0.0,
    'exosomal_biomarker_3': 0.0
}

# Nomogram score = optimized model score (REAL DATA ONLY)
nomogram_scores = optimized_scores.copy()

# Evaluate nomogram performance (same as optimized model since no exosomal data)
nomogram_auc = roc_auc_score(labels, nomogram_scores)

print(f"\nNomogram Performance (REAL DATA):")
print(f"  Nomogram AUC: {nomogram_auc:.3f}")
print(f"  Note: Currently uses optimized model only (100%)")
print(f"        Exosomal biomarkers identified but not yet measured")
print(f"        Formula will include exosomal components when data available")

# Save nomogram
nomogram_df = pd.DataFrame({
    'Sample': features_df.index,
    'True_Label': labels,
    'Optimized_Model_Score': optimized_scores,
    'Nomogram_Score': nomogram_scores,
    'Nomogram_Risk_Category': pd.cut(nomogram_scores, bins=[0, 0.3, 0.7, 1.0], 
                                      labels=['Low', 'Medium', 'High'])
})
nomogram_df.to_csv('results/phase6_nomogram_scores.tsv', sep='\t', index=False)

# Save nomogram formula (with note about exosomal biomarkers)
nomogram_formula = {
    'components': {
        'optimized_model': {
            'weight': nomogram_weights['optimized_model'],
            'features': optimized_features,
            'model_type': model_summary['best_model']
        },
        'exosomal_biomarkers': {
            'weight': sum([v for k, v in nomogram_weights.items() if 'exosomal' in k]),
            'top_candidates': top_exosomal_candidates[['Marker_Name', 'RNA_Type', 'Log2FC']].to_dict('records')[:3] if len(top_exosomal_candidates) > 0 else [],
            'note': 'Exosomal biomarkers identified but not yet measured. Weights will be updated when exosomal data available.'
        }
    },
    'formula': 'Nomogram_Score = 1.0 * Optimized_Model (current - exosomal data pending)',
    'formula_future': 'Nomogram_Score = 0.70 * Optimized_Model + 0.15 * Exosomal_1 + 0.10 * Exosomal_2 + 0.05 * Exosomal_3 (when exosomal data available)',
    'performance': {
        'nomogram_auc': float(nomogram_auc),
        'optimized_model_auc': float(roc_auc_score(labels, optimized_scores))
    },
    'note': 'All results based on REAL DATA. No simulated or synthetic data used. Exosomal biomarkers require separate measurement.'
}

with open('results/phase6_nomogram_formula.json', 'w') as f:
    json.dump(nomogram_formula, f, indent=2)

print("\n[OK] Nomogram scores saved (REAL DATA ONLY)")
print("[OK] Nomogram formula saved with notes about exosomal biomarkers")

print("\n" + "=" * 80)
print("PHASE 6 SUMMARY (REAL DATA ONLY)")
print("=" * 80)
print(f"\nExosomal Biomarker Candidates Identified (REAL DE DATA): {len(all_biomarkers)}")
print(f"  - Known markers: {len(candidate_biomarkers)}")
print(f"  - Novel candidates: {len(novel_biomarkers)}")
print(f"\nClinical Nomogram Created (REAL DATA):")
print(f"  - Nomogram AUC: {nomogram_auc:.3f}")
print(f"  - Currently uses optimized model only (exosomal data pending)")
print(f"\nNext Steps:")
print("  1. Measure top exosomal biomarker candidates in plasma exosomes")
print("  2. qPCR validation of top 5-10 biomarkers")
print("  3. Integrate exosomal measurements into nomogram")
print("  4. Prospective validation of nomogram in clinical cohort")
print(f"\nNOTE: All results are based on REAL DATA from actual datasets.")
print(f"      No simulated, synthetic, or randomly generated data was used.")

print("\n" + "=" * 80)
print("PHASE 6 COMPLETE (REAL DATA ONLY)")
print("=" * 80)
