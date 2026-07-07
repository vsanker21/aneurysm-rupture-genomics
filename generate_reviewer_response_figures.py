"""
Generate figures for reviewer-addressed analyses.
Creates:
  1. Pathway enrichment bar plot from gse122897_pathway_enrichment_filtered.tsv (for S6)
  2. DESeq2 volcano plot from gse122897_deseq2_results.tsv
  3. Nested vs non-nested CV AUC comparison bar plot
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

# Publication-quality style
try:
    plt.style.use('seaborn-v0_8-paper')
except Exception:
    try:
        plt.style.use('seaborn-paper')
    except Exception:
        plt.style.use('default')

plt.rcParams.update({
    'font.size': 14,
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 11,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

OUT_DIR = Path("results/figures/ultra_hd")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def plot_pathway_enrichment_gse122897():
    """Pathway enrichment bar plot from GSE122897 g:Profiler results (matches Table S5)."""
    path = Path("results/gse122897_pathway_enrichment_filtered.tsv")
    if not path.exists():
        print(f"  [SKIP] {path} not found")
        return None

    df = pd.read_csv(path, sep='\t')
    if len(df) == 0:
        print("  [SKIP] No pathways in file")
        return None

    # Use Term_Name and P_value (g:Profiler format)
    term_col = 'Term_Name' if 'Term_Name' in df.columns else 'Term'
    df = df.dropna(subset=[term_col, 'P_value'])
    df = df.sort_values('P_value', ascending=True).head(15)

    df['neg_log10_p'] = -np.log10(df['P_value'].clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(10, 8))
    y_pos = np.arange(len(df))
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(df)))
    bars = ax.barh(y_pos, df['neg_log10_p'], color=colors, alpha=0.85, edgecolor='black', linewidth=0.8)

    labels = [str(x)[:55] + ('...' if len(str(x)) > 55 else '') for x in df[term_col]]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=11)
    ax.set_xlabel('-log10(P-value)', fontweight='bold')
    ax.set_title('Pathway Enrichment (GSE122897 DE Genes)\nRanked gene list, g:Profiler GO/KEGG/Reactome, FDR < 0.05', fontweight='bold')
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3, axis='x')

    out = OUT_DIR / "pathway_enrichment_gse122897_ultra_hd.png"
    plt.tight_layout()
    plt.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  [OK] Saved: {out}")
    return str(out)


def plot_deseq2_volcano():
    """Volcano plot from DESeq2 results (raw counts, padj)."""
    path = Path("results/gse122897_deseq2_results.tsv")
    if not path.exists():
        print(f"  [SKIP] {path} not found")
        return None

    df = pd.read_csv(path, sep='\t')
    # Handle quoted headers
    df.columns = [c.strip('"') for c in df.columns]

    fc_col = 'log2FoldChange' if 'log2FoldChange' in df.columns else 'Log2FC'
    padj_col = 'padj' if 'padj' in df.columns else 'P_adj'
    if fc_col not in df.columns or padj_col not in df.columns:
        print("  [SKIP] Missing log2FoldChange or padj columns")
        return None

    df = df.dropna(subset=[fc_col, padj_col])
    df['neg_log10_padj'] = -np.log10(df[padj_col].clip(lower=1e-300))

    # Significance categories
    sig = (df[padj_col] < 0.05) & (np.abs(df[fc_col]) >= 0.5)
    up = sig & (df[fc_col] > 0)
    down = sig & (df[fc_col] < 0)

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(df.loc[~sig, fc_col], df.loc[~sig, 'neg_log10_padj'],
               c='#888888', s=15, alpha=0.6, label='Not significant (padj ≥ 0.05)')
    ax.scatter(df.loc[up, fc_col], df.loc[up, 'neg_log10_padj'],
               c='#b71c1c', s=25, alpha=0.8, label='Upregulated (padj < 0.05)')
    ax.scatter(df.loc[down, fc_col], df.loc[down, 'neg_log10_padj'],
               c='#1565c0', s=25, alpha=0.8, label='Downregulated (padj < 0.05)')

    # Label top genes from annotated file
    annot_path = Path("results/gse122897_deseq2_top_annotated.tsv")
    if annot_path.exists():
        top = pd.read_csv(annot_path, sep='\t').head(12)
        if 'Gene_Symbol' in top.columns and 'Log2FC' in top.columns:
            df['Gene_ID_str'] = df['Gene_ID'].astype(str).str.strip('"')
            for _, row in top.iterrows():
                sym = str(row.get('Gene_Symbol', '')).strip()
                if sym and sym != 'nan':
                    fc = float(row.get('Log2FC', row.get('log2FoldChange', 0)))
                    gid = str(row.get('Gene_ID', '')).strip()
                    match = df[df['Gene_ID_str'] == gid]
                    if len(match) > 0:
                        yval = match['neg_log10_padj'].iloc[0]
                    else:
                        yval = -np.log10(float(row.get('P_adj', 1e-10)) + 1e-300)
                    ax.annotate(sym, (fc, yval), fontsize=9, alpha=0.9,
                                xytext=(5, 5), textcoords='offset points')

    ax.axhline(-np.log10(0.05), color='gray', linestyle='--', linewidth=1)
    ax.axvline(0, color='black', linewidth=0.5)
    ax.set_xlabel('log2 fold change (ruptured vs unruptured)', fontweight='bold')
    ax.set_ylabel('-log10(adjusted p-value)', fontweight='bold')
    ax.set_title('DESeq2 Differential Expression (GSE122897, raw counts)\nn=43: 22 ruptured, 21 unruptured', fontweight='bold')
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)

    out = OUT_DIR / "deseq2_volcano_ultra_hd.png"
    plt.tight_layout()
    plt.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  [OK] Saved: {out}")
    return str(out)


def plot_nested_vs_non_nested_cv():
    """Bar plot comparing nested DE-within-CV vs non-nested CV AUC."""
    nested_path = Path("results/gse122897_nested_de_cv_summary.tsv")
    # Non-nested: from manuscript (0.798 ± 0.156)
    non_nested_auc = 0.798
    non_nested_std = 0.156

    nested_auc = 0.635
    nested_std = 0.173
    if nested_path.exists():
        ndf = pd.read_csv(nested_path, sep='\t')
        if 'Metric' in ndf.columns and 'Value' in ndf.columns:
            d = dict(zip(ndf['Metric'], ndf['Value']))
            nested_auc = float(d.get('Nested_CV_AUC_Mean', nested_auc))
            nested_std = float(d.get('Nested_CV_AUC_Std', nested_std))

    fig, ax = plt.subplots(figsize=(6, 5))
    x = [0, 1]
    aucs = [non_nested_auc, nested_auc]
    stds = [non_nested_std, nested_std]
    labels = ['Non-nested CV\n(DE on full cohort)', 'Nested CV\n(DE within each fold)']
    colors = ['#2e7d32', '#1565c0']

    bars = ax.bar(x, aucs, yerr=stds, color=colors, capsize=6, edgecolor='black', linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel('Cross-validation AUC (mean ± SD)', fontweight='bold')
    ax.set_title('Effect of Circularity: Nested vs Non-nested Design\n(GSE122897, 5-fold stratified CV)', fontweight='bold')
    ax.set_ylim(0, 1)
    ax.axhline(0.5, color='gray', linestyle='--', linewidth=1)
    ax.grid(True, alpha=0.3, axis='y')

    for i, (a, s) in enumerate(zip(aucs, stds)):
        ax.text(i, a + s + 0.03, f'{a:.2f}±{s:.2f}', ha='center', fontsize=11, fontweight='bold')

    out = OUT_DIR / "nested_vs_non_nested_cv_ultra_hd.png"
    plt.tight_layout()
    plt.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  [OK] Saved: {out}")
    return str(out)


def main():
    print("=" * 60)
    print("GENERATING REVIEWER RESPONSE FIGURES")
    print("=" * 60)

    plot_pathway_enrichment_gse122897()
    plot_deseq2_volcano()
    plot_nested_vs_non_nested_cv()

    print("\nDone. Figures saved to results/figures/ultra_hd/")


if __name__ == "__main__":
    main()
