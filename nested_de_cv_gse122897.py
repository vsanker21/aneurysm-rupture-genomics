"""
Nested DE within Cross-Validation for GSE122897.
Performs differential expression on the training set only in each CV fold,
then computes the DE signature score from those fold-specific genes.
This avoids circularity (DE genes from full data used in model evaluation).
"""
import pandas as pd
import numpy as np
import os
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib.pyplot as plt

# EC and ECM are predefined; DE is data-derived within each fold
EC_GENES = ['ANGPT2', 'ICAM1', 'VCAM1', 'THBS1', 'PERP']
ECM_GENES = ['MMP9', 'MMP2', 'TIMP1', 'COL3A1', 'COL1A1']

def load_data():
    """Load TPM and metadata"""
    expr_path = "data/sc_rnaseq/GSE122897/GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
    meta_paths = [
        "data/sc_rnaseq/GSE122897/sample_metadata_with_clinical.tsv",
        "data/sc_rnaseq/GSE122897/sample_metadata.tsv",
    ]
    if not os.path.exists(expr_path):
        raise FileNotFoundError(expr_path)
    expr = pd.read_csv(expr_path, sep='\t', compression='gzip', index_col=0)
    meta = None
    for p in meta_paths:
        if os.path.exists(p):
            meta = pd.read_csv(p, sep='\t')
            break
    if meta is None or 'Group' not in meta.columns:
        raise ValueError("Metadata with Group column not found")
    meta = meta[meta['Group'].isin(['Ruptured', 'Unruptured'])].copy()
    meta = meta[meta['Sample'].isin(expr.columns)]
    common = [c for c in meta['Sample'] if c in expr.columns]
    expr = expr[common]
    meta = meta[meta['Sample'].isin(common)].reset_index(drop=True)
    if 'Sample' not in meta.columns:
        meta['Sample'] = common
    return expr, meta

def _map_genes(expr_df, symbols):
    """Map symbols to expr index (NCBI/Entrez IDs)."""
    # Prefer gene_symbol_to_entrez_mapping (Symbol -> Entrez)
    for mapping_file in ["results/gene_symbol_to_entrez_mapping.tsv", "results/gene_id_mapping_unified.tsv"]:
        if not os.path.exists(mapping_file):
            continue
        m = pd.read_csv(mapping_file, sep='\t')
        id_col = 'Entrez_ID' if 'Entrez_ID' in m.columns else 'Original_ID'
        sym_col = 'Gene_Symbol' if 'Gene_Symbol' in m.columns else None
        if not sym_col or id_col not in m.columns:
            continue
        m = m.dropna(subset=[sym_col, id_col])
        sym_to_id = dict(zip(m[sym_col].astype(str).str.upper(), m[id_col].astype(int).astype(str)))
        idx_str = set(expr_df.index.astype(str))
        hit_ids = []
        for s in symbols:
            uid = str(s).strip().upper()
            eid = sym_to_id.get(uid)
            if eid and eid in idx_str:
                hit_ids.append(eid)
        if hit_ids:
            return [i for i in expr_df.index if str(i) in hit_ids]
    return []

def run_de_on_subset(expr, train_samples, meta, top_n=20):
    """
    Run Welch t-test DE on train samples only. Return top_n gene IDs by FDR and |Log2FC|.
    """
    groups = meta['Group'].unique()
    group1 = meta[meta['Group'] == groups[0]]['Sample'].tolist()
    group2 = meta[meta['Group'] == groups[1]]['Sample'].tolist()
    g1 = [s for s in group1 if s in train_samples]
    g2 = [s for s in group2 if s in train_samples]
    if len(g1) < 2 or len(g2) < 2:
        return []
    ex = expr[[c for c in g1 + g2 if c in expr.columns]]
    ex1 = ex[g1].values
    ex2 = ex[g2].values
    ex1_log = np.log2(ex1 + 1)
    ex2_log = np.log2(ex2 + 1)
    t, p = ttest_ind(ex1_log, ex2_log, axis=1)
    m1 = np.mean(ex1, axis=1)
    m2 = np.mean(ex2, axis=1)
    log2fc = np.log2((m1 + 1e-6) / (m2 + 1e-6))
    p_for_corr = np.nan_to_num(p, nan=1.0)
    _, padj, _, _ = multipletests(p_for_corr, method='fdr_bh')
    padj = np.where(np.isnan(p), np.nan, padj)
    df = pd.DataFrame({
        'gene_idx': expr.index,
        'Log2FC': log2fc,
        'P_adj': padj,
        'Abs_Log2FC': np.abs(log2fc),
    })
    df = df.dropna(subset=['P_adj'])
    df = df.sort_values(['P_adj', 'Abs_Log2FC'], ascending=[True, False])
    top = df.head(top_n)['gene_idx'].tolist()
    return top

def compute_de_score(expr, gene_ids, samples):
    """Mean log2(TPM+1) over gene_ids for given samples."""
    ex = expr.loc[[g for g in gene_ids if g in expr.index], [s for s in samples if s in expr.columns]]
    if ex.empty:
        return np.zeros(len(samples))
    return np.log2(ex + 1).mean(axis=0).values

def compute_ec_ecm_scores(expr, meta):
    """Compute EC and ECM scores (predefined gene sets)."""
    ec_ids = _map_genes(expr, EC_GENES)
    ecm_ids = _map_genes(expr, ECM_GENES)
    ec_score = np.log2(expr.loc[ec_ids] + 1).mean(axis=0) if ec_ids else pd.Series(0, index=expr.columns)
    ecm_score = np.log2(expr.loc[ecm_ids] + 1).mean(axis=0) if ecm_ids else pd.Series(0, index=expr.columns)
    return ec_score, ecm_score

def main():
    print("Nested DE within Cross-Validation (GSE122897)")
    print("=" * 55)
    expr, meta = load_data()
    y = (meta['Group'] == 'Ruptured').astype(int).values
    samples = meta['Sample'].tolist()
    n_splits = 5
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)

    ec_score, ecm_score = compute_ec_ecm_scores(expr, meta)
    ec_arr = ec_score.loc[samples].values
    ecm_arr = ecm_score.loc[samples].values

    fold_aucs = []
    oof_pred = np.full(len(y), np.nan)
    for fold, (train_idx, test_idx) in enumerate(cv.split(np.zeros(len(y)), y)):
        train_samples = [samples[i] for i in train_idx]
        test_samples = [samples[i] for i in test_idx]
        de_genes = run_de_on_subset(expr, train_samples, meta, top_n=25)
        if len(de_genes) < 2:
            de_train = np.zeros(len(train_samples))
            de_test = np.zeros(len(test_samples))
        else:
            de_train = compute_de_score(expr, de_genes, train_samples)
            de_test = compute_de_score(expr, de_genes, test_samples)

        X_train = np.column_stack([
            ec_arr[train_idx], ecm_arr[train_idx], de_train
        ])
        X_test = np.column_stack([
            ec_arr[test_idx], ecm_arr[test_idx], de_test
        ])
        scaler = StandardScaler()
        X_train_s = scaler.fit_transform(X_train)
        X_test_s = scaler.transform(X_test)
        clf = LogisticRegression(C=1.0, max_iter=1000, random_state=42)
        clf.fit(X_train_s, y[train_idx])
        pred = clf.predict_proba(X_test_s)[:, 1]
        oof_pred[test_idx] = pred
        auc = roc_auc_score(y[test_idx], pred)
        fold_aucs.append(auc)
        print(f"  Fold {fold+1}: AUC = {auc:.3f} (DE genes in fold: {len(de_genes)})")

    mean_auc = np.mean(fold_aucs)
    std_auc = np.std(fold_aucs)
    print(f"\nNested CV AUC: {mean_auc:.3f} ± {std_auc:.3f}")

    # Save summary
    out_path = "results/gse122897_nested_de_cv_summary.tsv"
    os.makedirs("results", exist_ok=True)
    pd.DataFrame({
        'Metric': ['Nested_CV_AUC_Mean', 'Nested_CV_AUC_Std', 'N_Folds'],
        'Value': [mean_auc, std_auc, n_splits],
    }).to_csv(out_path, sep='\t', index=False)
    print(f"Saved: {out_path}")

    valid = ~np.isnan(oof_pred)
    if valid.sum() == len(y):
        pred_df = pd.DataFrame({
            'Sample': samples,
            'True_Label': y,
            'Nested_OOF_Prediction': oof_pred,
        })
        pred_path = "results/gse122897_nested_de_cv_oof_predictions.tsv"
        pred_df.to_csv(pred_path, sep='\t', index=False)
        print(f"Saved: {pred_path}")

        fpr, tpr, thresholds = roc_curve(y, oof_pred)
        roc_df = pd.DataFrame({
            'FPR': fpr,
            'TPR': tpr,
            'Threshold': thresholds,
        })
        roc_path = "results/gse122897_nested_de_cv_roc_curve.tsv"
        roc_df.to_csv(roc_path, sep='\t', index=False)
        print(f"Saved: {roc_path}")

        fig_dir = "results/figures/ultra_hd"
        os.makedirs(fig_dir, exist_ok=True)
        fig_path = os.path.join(fig_dir, "nested_cv_roc_curve_ultra_hd.png")
        plt.figure(figsize=(7, 6), dpi=300)
        plt.plot(fpr, tpr, color="#1f77b4", linewidth=2.5, label=f"Nested CV OOF ROC (AUC = {mean_auc:.3f})")
        plt.plot([0, 1], [0, 1], linestyle="--", color="gray", linewidth=1.5, label="Chance")
        plt.xlabel("False positive rate")
        plt.ylabel("True positive rate")
        plt.title("Leakage-controlled nested CV ROC")
        plt.legend(loc="lower right")
        plt.tight_layout()
        plt.savefig(fig_path, bbox_inches="tight")
        plt.close()
        print(f"Saved: {fig_path}")
    print("\nNote: Nested CV typically yields lower AUC than non-nested (avoids circularity).")

if __name__ == "__main__":
    main()
