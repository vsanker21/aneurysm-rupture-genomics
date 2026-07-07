from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm, ttest_ind
from statsmodels.stats.multitest import multipletests

from advanced_public_signature_transfer import _load_gse13353_gene_level
from manuscript_results_detail import annotate_deseq2_with_symbols, load_gene_symbol_map

BASE = Path(r"e:\Aneurysms Genomics Project")
OUT_DIR = BASE / "results" / "rupture_meta_analysis"
FIG_DIR = BASE / "results" / "figures" / "ultra_hd"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)


def _effective_n(n_case: int, n_control: int) -> float:
    return float(np.sqrt((n_case * n_control) / max(n_case + n_control, 1)))


def _load_gse122897() -> tuple[pd.DataFrame, int, int]:
    de = pd.read_csv(BASE / "results" / "gse122897_deseq2_results.tsv", sep="\t")
    sym = load_gene_symbol_map(BASE)
    de = annotate_deseq2_with_symbols(de, sym)
    out = de[["display_symbol", "log2FoldChange", "pvalue", "padj"]].rename(
        columns={
            "display_symbol": "Gene",
            "log2FoldChange": "Log2FC",
            "pvalue": "P_value",
            "padj": "FDR_BH",
        }
    )
    meta = pd.read_csv(BASE / "data" / "sc_rnaseq" / "GSE122897" / "sample_metadata.tsv", sep="\t")
    meta = meta[meta["Group"].isin(["Ruptured", "Unruptured"])].copy()
    n_r = int((meta["Group"] == "Ruptured").sum())
    n_u = int((meta["Group"] == "Unruptured").sum())
    return out.drop_duplicates("Gene"), n_r, n_u


def _load_gse54083() -> tuple[pd.DataFrame, int, int]:
    out = pd.read_csv(BASE / "results" / "external_validation_gse54083" / "gse54083_gene_level_validation.tsv", sep="\t")
    out = out[["Gene", "Log2FC_RvU", "P_value", "FDR_BH"]].rename(columns={"Log2FC_RvU": "Log2FC"})
    meta = pd.read_csv(BASE / "results" / "external_validation_gse54083" / "gse54083_sample_metadata.tsv", sep="\t")
    meta = meta[meta["Group"].isin(["Ruptured", "Unruptured"])].copy()
    n_r = int((meta["Group"] == "Ruptured").sum())
    n_u = int((meta["Group"] == "Unruptured").sum())
    return out.drop_duplicates("Gene"), n_r, n_u


def _compute_gse13353_gene_level() -> tuple[pd.DataFrame, int, int]:
    expr, meta = _load_gse13353_gene_level()
    transform_needed = float(np.nanpercentile(expr.values, 95)) > 100.0
    expr_use = np.log2(expr + 1.0) if transform_needed else expr.copy()
    r_cols = meta.loc[meta["Group"] == "Ruptured", "Sample"].tolist()
    u_cols = meta.loc[meta["Group"] == "Unruptured", "Sample"].tolist()

    rows = []
    for gene, row in expr_use.iterrows():
        rv = pd.to_numeric(row[r_cols], errors="coerce").dropna().values.astype(float)
        uv = pd.to_numeric(row[u_cols], errors="coerce").dropna().values.astype(float)
        if len(rv) < 2 or len(uv) < 2:
            continue
        t_stat, p_val = ttest_ind(rv, uv, equal_var=False, nan_policy="omit")
        rows.append(
            {
                "Gene": gene,
                "Mean_Ruptured": float(np.nanmean(rv)),
                "Mean_Unruptured": float(np.nanmean(uv)),
                "Log2FC": float(np.nanmean(rv) - np.nanmean(uv)),
                "T_statistic": float(t_stat) if pd.notna(t_stat) else np.nan,
                "P_value": float(p_val) if pd.notna(p_val) else np.nan,
            }
        )
    out = pd.DataFrame(rows).sort_values("P_value", ascending=True).reset_index(drop=True)
    if out.empty:
        raise RuntimeError("GSE13353 gene-level analysis produced no rows.")
    out["FDR_BH"] = multipletests(out["P_value"].fillna(1.0).values, method="fdr_bh")[1]
    out.to_csv(OUT_DIR / "gse13353_gene_level_meta_input.tsv", sep="\t", index=False)
    return out, len(r_cols), len(u_cols)


def _stouffer_meta(df: pd.DataFrame, weight_cols: list[str], z_cols: list[str]) -> pd.Series:
    weights = df[weight_cols].to_numpy(dtype=float)
    zs = df[z_cols].to_numpy(dtype=float)
    numerator = np.nansum(weights * zs, axis=1)
    denominator = np.sqrt(np.nansum(weights**2, axis=1))
    meta_z = numerator / denominator
    meta_p = 2.0 * norm.sf(np.abs(meta_z))
    return pd.Series(meta_z), pd.Series(meta_p)


def _plot_top_consensus(df: pd.DataFrame, out_path: Path) -> None:
    top = df.head(20).copy().sort_values("Meta_Z", ascending=True)
    colors = np.where(top["Consensus_direction"] == "Up in ruptured", "#b54a4a", "#4b7bb8")

    fig, ax = plt.subplots(figsize=(10, 9))
    ax.barh(top["Gene"], top["Meta_Z"], color=colors, alpha=0.9)
    ax.axvline(0.0, color="black", lw=1)
    ax.set_xlabel("Weighted Stouffer meta-Z")
    ax.set_title("Consensus rupture genes across direct public rupture cohorts")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    g122, n122_r, n122_u = _load_gse122897()
    g540, n540_r, n540_u = _load_gse54083()
    g133, n133_r, n133_u = _compute_gse13353_gene_level()

    merged = (
        g122.rename(columns={"Log2FC": "Log2FC_GSE122897", "P_value": "P_GSE122897", "FDR_BH": "FDR_GSE122897"})
        .merge(
            g540.rename(columns={"Log2FC": "Log2FC_GSE54083", "P_value": "P_GSE54083", "FDR_BH": "FDR_GSE54083"}),
            on="Gene",
            how="inner",
        )
        .merge(
            g133.rename(columns={"Log2FC": "Log2FC_GSE13353", "P_value": "P_GSE13353", "FDR_BH": "FDR_GSE13353"}),
            on="Gene",
            how="inner",
        )
    )
    if merged.empty:
        raise RuntimeError("No overlapping genes across GSE122897, GSE54083, and GSE13353.")

    eps = 1e-300
    for dataset in ["GSE122897", "GSE54083", "GSE13353"]:
        pcol = f"P_{dataset}"
        lcol = f"Log2FC_{dataset}"
        zcol = f"Z_{dataset}"
        merged[zcol] = np.sign(merged[lcol]).replace(0, np.nan) * norm.isf(merged[pcol].clip(lower=eps) / 2.0)

    merged["Weight_GSE122897"] = _effective_n(n122_r, n122_u)
    merged["Weight_GSE54083"] = _effective_n(n540_r, n540_u)
    merged["Weight_GSE13353"] = _effective_n(n133_r, n133_u)

    merged["All_three_directionally_concordant"] = (
        (np.sign(merged["Log2FC_GSE122897"]) == np.sign(merged["Log2FC_GSE54083"]))
        & (np.sign(merged["Log2FC_GSE122897"]) == np.sign(merged["Log2FC_GSE13353"]))
    )
    merged["Consensus_direction"] = np.where(
        merged["All_three_directionally_concordant"] & (merged["Log2FC_GSE122897"] > 0),
        "Up in ruptured",
        np.where(
            merged["All_three_directionally_concordant"] & (merged["Log2FC_GSE122897"] < 0),
            "Down in ruptured",
            "Discordant",
        ),
    )

    merged["Meta_Z"], merged["Meta_P"] = _stouffer_meta(
        merged,
        ["Weight_GSE122897", "Weight_GSE54083", "Weight_GSE13353"],
        ["Z_GSE122897", "Z_GSE54083", "Z_GSE13353"],
    )
    merged = merged.sort_values(["Meta_P", "All_three_directionally_concordant"], ascending=[True, False]).reset_index(drop=True)
    merged["Meta_FDR_BH"] = multipletests(merged["Meta_P"].fillna(1.0).values, method="fdr_bh")[1]

    concordant = merged[merged["All_three_directionally_concordant"]].copy()
    concordant = concordant.sort_values("Meta_P", ascending=True).reset_index(drop=True)
    top = concordant.head(25).copy()

    merged.to_csv(OUT_DIR / "rupture_consensus_meta_all.tsv", sep="\t", index=False)
    concordant.to_csv(OUT_DIR / "rupture_consensus_meta_concordant.tsv", sep="\t", index=False)
    top.to_csv(OUT_DIR / "rupture_consensus_meta_top25.tsv", sep="\t", index=False)

    figure_path = FIG_DIR / "rupture_consensus_meta_ultra_hd.png"
    _plot_top_consensus(top, figure_path)

    summary = {
        "datasets": {
            "GSE122897": {"n_ruptured": n122_r, "n_unruptured": n122_u},
            "GSE54083": {"n_ruptured": n540_r, "n_unruptured": n540_u},
            "GSE13353": {"n_ruptured": n133_r, "n_unruptured": n133_u},
        },
        "n_overlap_all_three": int(len(merged)),
        "n_directionally_concordant_all_three": int(len(concordant)),
        "n_meta_fdr_lt_0_05": int((merged["Meta_FDR_BH"] < 0.05).sum()),
        "n_concordant_meta_fdr_lt_0_05": int((concordant["Meta_FDR_BH"] < 0.05).sum()),
        "top_concordant_genes": top[["Gene", "Consensus_direction", "Meta_Z", "Meta_P", "Meta_FDR_BH"]]
        .head(10)
        .to_dict(orient="records"),
    }
    (OUT_DIR / "rupture_consensus_meta_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    print("\nTop concordant genes:")
    print(
        top[
            [
                "Gene",
                "Consensus_direction",
                "Meta_Z",
                "Meta_P",
                "Meta_FDR_BH",
                "Log2FC_GSE122897",
                "Log2FC_GSE54083",
                "Log2FC_GSE13353",
            ]
        ]
        .head(15)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
