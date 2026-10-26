from __future__ import annotations

import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, norm

matplotlib.use("Agg")

import matplotlib.pyplot as plt

try:
    from project_paths import get_project_root
    BASE = get_project_root()
except ImportError:
    BASE = Path(__file__).resolve().parent
OUT_DIR = BASE / "results" / "blood_bridge_panel"
OUT_GENE_TSV = OUT_DIR / "blood_bridge_panel_gene_summary.tsv"
OUT_PERF_TSV = OUT_DIR / "blood_bridge_panel_performance.tsv"
OUT_SCORE_TSV = OUT_DIR / "blood_bridge_panel_sample_scores.tsv"
OUT_SUMMARY_JSON = OUT_DIR / "blood_bridge_panel_summary.json"
OUT_FIG = BASE / "results" / "figures" / "ultra_hd" / "blood_bridge_panel_transport_ultra_hd.png"

PANEL_SIZE = 8


def _auc_from_scores(labels: np.ndarray, scores: np.ndarray) -> float:
    labels = np.asarray(labels, dtype=bool)
    scores = np.asarray(scores, dtype=float)
    pos = labels.sum()
    neg = (~labels).sum()
    if pos == 0 or neg == 0:
        return float("nan")
    order = scores.argsort(kind="mergesort")
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(scores) + 1, dtype=float)
    pos_ranks = ranks[labels].sum()
    return float((pos_ranks - pos * (pos + 1) / 2.0) / (pos * neg))


def _load_bridge_candidates(base: Path) -> pd.DataFrame:
    shortlist = pd.read_csv(base / "results" / "phases_g" / "shortlist_evidence.tsv", sep="\t")
    g106520 = pd.read_csv(base / "results" / "external_validation_gse106520" / "gse106520_shortlist_overlap.tsv", sep="\t")

    keep = shortlist[
        [
            "symbol",
            "DESeq2_log2FC",
            "DESeq2_padj",
            "GSE159610_Log2FC_aneurysm_vs_control",
            "GSE159610_P",
            "direction_match_rupture_and_aneurysm",
        ]
    ].rename(columns={"symbol": "Gene"})
    keep = keep.dropna(subset=["Gene", "GSE159610_Log2FC_aneurysm_vs_control"]).copy()

    g106_keep = g106520[
        [
            "Gene",
            "Log2FC_AvC",
            "P_value",
            "DirectionMatch_GSE122897",
        ]
    ].rename(columns={"Log2FC_AvC": "GSE106520_Log2FC_aneurysm_vs_control", "P_value": "GSE106520_P"})

    df = keep.merge(g106_keep, on="Gene", how="inner")
    df["Bridge_direction_match_GSE159610"] = df["direction_match_rupture_and_aneurysm"] == True  # noqa: E712
    df["Bridge_direction_match_GSE106520"] = df["DirectionMatch_GSE122897"] == True  # noqa: E712
    df = df[df["Bridge_direction_match_GSE159610"] & df["Bridge_direction_match_GSE106520"]].copy()

    p1 = np.clip(pd.to_numeric(df["GSE159610_P"], errors="coerce").fillna(1.0).to_numpy(dtype=float), 1e-300, 1.0)
    p2 = np.clip(pd.to_numeric(df["GSE106520_P"], errors="coerce").fillna(1.0).to_numpy(dtype=float), 1e-300, 1.0)
    df["Bridge_Stouffer_Z"] = (norm.isf(p1 / 2.0) + norm.isf(p2 / 2.0)) / np.sqrt(2.0)
    df["Bridge_meta_P"] = 2.0 * norm.sf(np.abs(df["Bridge_Stouffer_Z"]))
    df["abs_discovery_log2FC"] = pd.to_numeric(df["DESeq2_log2FC"], errors="coerce").abs()
    df = df.sort_values(["Bridge_meta_P", "abs_discovery_log2FC", "Gene"], ascending=[True, False, True]).reset_index(drop=True)
    df["Panel_rank"] = np.arange(1, len(df) + 1)
    return df


def _load_gse159610_expression(base: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = pd.read_csv(
        base / "data" / "blood_rnaseq" / "GSE159610" / "GSE159610_raw_gene_counts_matrix.txt.gz",
        sep="\t",
        compression="gzip",
    )
    sample_cols = [c for c in counts.columns if c.startswith("DSA")]
    expr = counts[["gene_name", *sample_cols]].copy()
    expr["gene_name"] = expr["gene_name"].astype(str).str.strip()
    expr = expr[expr["gene_name"].notna() & (expr["gene_name"] != "") & (expr["gene_name"] != "nan")]
    expr = expr.groupby("gene_name", as_index=True)[sample_cols].sum(min_count=1)
    lib_sizes = expr.sum(axis=0).replace(0, np.nan)
    expr = np.log2(expr.div(lib_sizes, axis=1) * 1_000_000.0 + 1.0)
    expr.index.name = "Gene"

    meta = pd.read_csv(base / "data" / "blood_rnaseq" / "GSE159610" / "sample_metadata_with_labels.tsv", sep="\t")
    meta["sample_number"] = meta["Title"].astype(str).str.extract(r"sample\s+(\d+)$").astype(int)
    meta["Sample"] = "DSA" + meta["sample_number"].astype(str)
    meta = meta[["Sample", "Group"]].copy()
    return expr, meta


def _load_simple_expression(
    matrix_path: Path,
    metadata_path: Path,
    *,
    sample_col: str,
    group_col: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = pd.read_csv(matrix_path, sep="\t").set_index("Gene")
    meta = pd.read_csv(metadata_path, sep="\t")
    meta = meta[[sample_col, group_col]].rename(columns={sample_col: "Sample", group_col: "Group"})
    meta["Sample"] = meta["Sample"].astype(str).str.extract(r"(\S+)$", expand=False).fillna(meta["Sample"].astype(str))
    meta["Group"] = meta["Group"].astype(str)
    return expr, meta


def _zscore_rows(frame: pd.DataFrame) -> pd.DataFrame:
    means = frame.mean(axis=1)
    sds = frame.std(axis=1, ddof=0).replace(0, np.nan)
    return frame.sub(means, axis=0).div(sds, axis=0)


def _score_cohort(
    expr: pd.DataFrame,
    meta: pd.DataFrame,
    *,
    panel: pd.DataFrame,
    dataset: str,
    role: str,
    positive_label: str,
    negative_label: str,
) -> tuple[dict[str, object], pd.DataFrame]:
    available_samples = [c for c in expr.columns if c in set(meta["Sample"])]
    expr = expr.loc[:, available_samples].copy()
    meta = meta[meta["Sample"].isin(available_samples)].copy()
    meta = meta[meta["Group"].isin([positive_label, negative_label])].copy()
    meta = meta.drop_duplicates(subset=["Sample"]).set_index("Sample")
    expr = expr.loc[:, expr.columns.isin(meta.index)].copy()
    meta = meta.loc[expr.columns]

    panel = panel.copy()
    panel["sign"] = np.sign(pd.to_numeric(panel["DESeq2_log2FC"], errors="coerce"))
    panel = panel[panel["sign"] != 0]
    genes = [g for g in panel["Gene"] if g in expr.index]
    signs = panel.set_index("Gene").loc[genes, "sign"].astype(float)

    panel_expr = expr.loc[genes].apply(pd.to_numeric, errors="coerce")
    panel_z = _zscore_rows(panel_expr).dropna(how="all")
    genes = [g for g in genes if g in panel_z.index]
    signed = panel_z.loc[genes].mul(signs.loc[genes], axis=0)
    scores = signed.mean(axis=0, skipna=True)

    labels = meta["Group"] == positive_label
    auc = _auc_from_scores(labels.to_numpy(dtype=bool), scores.to_numpy(dtype=float))
    pos_scores = scores[labels]
    neg_scores = scores[~labels]
    mw_p = float(mannwhitneyu(pos_scores, neg_scores, alternative="two-sided").pvalue) if len(pos_scores) and len(neg_scores) else float("nan")

    summary = {
        "Dataset": dataset,
        "Role": role,
        "Positive_group": positive_label,
        "Negative_group": negative_label,
        "Samples": f"{int(labels.sum())} positive / {int((~labels).sum())} comparator",
        "Panel_genes_measured": int(len(genes)),
        "AUC": auc,
        "Mann_Whitney_P": mw_p,
        "Mean_score_positive": float(pos_scores.mean()),
        "Mean_score_negative": float(neg_scores.mean()),
        "Score_difference": float(pos_scores.mean() - neg_scores.mean()),
    }

    sample_scores = pd.DataFrame(
        {
            "Dataset": dataset,
            "Role": role,
            "Sample": scores.index,
            "Group": meta.loc[scores.index, "Group"].values,
            "Panel_score": scores.values,
        }
    )
    return summary, sample_scores


def _load_context_gene_effects(base: Path, genes: list[str]) -> pd.DataFrame:
    g367 = pd.read_csv(base / "results" / "external_validation_gse36791" / "gse36791_shortlist_overlap.tsv", sep="\t")
    g733 = pd.read_csv(base / "results" / "external_validation_gse73378" / "gse73378_shortlist_overlap.tsv", sep="\t")

    g367 = g367[["Gene", "Log2FC_RvC", "P_value", "DirectionMatch_GSE122897"]].rename(
        columns={
            "Log2FC_RvC": "GSE36791_Log2FC_postrupture",
            "P_value": "GSE36791_P",
            "DirectionMatch_GSE122897": "GSE36791_DirectionMatch",
        }
    )
    g733 = g733[["Gene", "Log2FC_CasevControl", "P_value", "DirectionMatch_GSE122897"]].rename(
        columns={
            "Log2FC_CasevControl": "GSE73378_Log2FC_pastSAH",
            "P_value": "GSE73378_P",
            "DirectionMatch_GSE122897": "GSE73378_DirectionMatch",
        }
    )
    merged = pd.DataFrame({"Gene": genes}).merge(g367, on="Gene", how="left").merge(g733, on="Gene", how="left")
    return merged


def _save_figure(panel: pd.DataFrame, perf: pd.DataFrame, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    heat = panel[
        [
            "Gene",
            "DESeq2_log2FC",
            "GSE159610_Log2FC_aneurysm_vs_control",
            "GSE106520_Log2FC_aneurysm_vs_control",
            "GSE36791_Log2FC_postrupture",
            "GSE73378_Log2FC_pastSAH",
        ]
    ].copy()
    heat = heat.rename(
        columns={
            "DESeq2_log2FC": "GSE122897 tissue rupture",
            "GSE159610_Log2FC_aneurysm_vs_control": "GSE159610 blood bridge",
            "GSE106520_Log2FC_aneurysm_vs_control": "GSE106520 blood bridge",
            "GSE36791_Log2FC_postrupture": "GSE36791 post-rupture",
            "GSE73378_Log2FC_pastSAH": "GSE73378 past SAH",
        }
    ).set_index("Gene")

    perf_plot = perf.copy()
    perf_plot["label"] = perf_plot["Dataset"] + "\n" + perf_plot["Samples"]
    colors = perf_plot["Role"].map(
        {
            "Blood bridge": "#4c72b0",
            "Blood context only": "#55a868",
        }
    ).fillna("#808080")

    fig = plt.figure(figsize=(14, 8))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1.0])

    ax0 = fig.add_subplot(gs[0, 0])
    vmax = float(np.nanmax(np.abs(heat.to_numpy(dtype=float)))) if not heat.empty else 1.0
    vmax = max(vmax, 0.25)
    im = ax0.imshow(heat.to_numpy(dtype=float), aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax0.set_xticks(np.arange(heat.shape[1]))
    ax0.set_xticklabels(list(heat.columns), rotation=30, ha="right")
    ax0.set_yticks(np.arange(heat.shape[0]))
    ax0.set_yticklabels(list(heat.index))
    ax0.set_title("A. Discovery-to-blood effect directions in the reduced panel", loc="left", fontsize=12, fontweight="bold")
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            val = heat.iloc[i, j]
            if pd.notna(val):
                ax0.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8)
    cbar = fig.colorbar(im, ax=ax0, fraction=0.046, pad=0.04)
    cbar.set_label("Log2 fold change", fontsize=9)

    ax1 = fig.add_subplot(gs[0, 1])
    y = np.arange(len(perf_plot))[::-1]
    ax1.barh(y, perf_plot["AUC"], color=list(colors))
    ax1.axvline(0.5, color="black", linestyle="--", linewidth=1)
    ax1.set_yticks(y)
    ax1.set_yticklabels(list(perf_plot["label"]), fontsize=9)
    ax1.set_xlim(0, 1)
    ax1.set_xlabel("Panel-score AUROC")
    ax1.set_title("B. Sample-level transportability of the reduced blood panel", loc="left", fontsize=12, fontweight="bold")
    for yi, auc, delta in zip(y, perf_plot["AUC"], perf_plot["Score_difference"]):
        ax1.text(min(float(auc) + 0.02, 0.98), yi, f"AUC {float(auc):.2f}; delta {float(delta):.2f}", va="center", fontsize=8)

    from matplotlib.patches import Patch

    ax1.legend(
        handles=[
            Patch(facecolor="#4c72b0", label="Bridge cohorts used for panel derivation"),
            Patch(facecolor="#55a868", label="Post-rupture / SAH context only"),
        ],
        frameon=False,
        loc="lower right",
        fontsize=8,
    )

    fig.suptitle("Reduced blood-testable panel anchored on GSE159610 and GSE106520", fontsize=14, fontweight="bold")
    fig.text(
        0.5,
        0.01,
        "The panel is selected only from genes directionally concordant in discovery tissue and both blood bridge cohorts. "
        "GSE36791 and GSE73378 are displayed as post-rupture/SAH context and are not used for panel selection.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.04, 1, 0.95])
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_blood_bridge_panel_assets(base: Path = BASE, panel_size: int = PANEL_SIZE) -> dict[str, object]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    candidates = _load_bridge_candidates(base)
    panel = candidates.head(panel_size).copy()
    context = _load_context_gene_effects(base, panel["Gene"].tolist())
    panel = panel.merge(context, on="Gene", how="left")

    g159610_expr, g159610_meta = _load_gse159610_expression(base)
    g106520_expr, g106520_meta = _load_simple_expression(
        base / "results" / "external_validation_gse106520" / "gse106520_gene_level_log2fpkm.tsv",
        base / "results" / "external_validation_gse106520" / "gse106520_sample_metadata.tsv",
        sample_col="Sample",
        group_col="Group",
    )
    g36791_expr, g36791_meta = _load_simple_expression(
        base / "results" / "external_validation_gse36791" / "gse36791_gene_level_log2_expression.tsv",
        base / "results" / "external_validation_gse36791" / "gse36791_sample_metadata.tsv",
        sample_col="Sample",
        group_col="Group",
    )
    g73378_expr, g73378_meta = _load_simple_expression(
        base / "results" / "external_validation_gse73378" / "gse73378_gene_level_expression.tsv",
        base / "results" / "external_validation_gse73378" / "gse73378_sample_metadata.tsv",
        sample_col="Title",
        group_col="Group",
    )

    perf_rows: list[dict[str, object]] = []
    score_frames: list[pd.DataFrame] = []
    for args in [
        (g159610_expr, g159610_meta, "GSE159610", "Blood bridge", "Aneurysm", "Control"),
        (g106520_expr, g106520_meta, "GSE106520", "Blood bridge", "Aneurysm", "Control"),
        (g36791_expr, g36791_meta, "GSE36791", "Blood context only", "Ruptured", "Control"),
        (g73378_expr, g73378_meta, "GSE73378", "Blood context only", "Case", "Control"),
    ]:
        perf_row, score_df = _score_cohort(
            args[0],
            args[1],
            panel=panel,
            dataset=args[2],
            role=args[3],
            positive_label=args[4],
            negative_label=args[5],
        )
        perf_rows.append(perf_row)
        score_frames.append(score_df)

    perf = pd.DataFrame(perf_rows)
    perf["AUC"] = perf["AUC"].astype(float)
    perf["Mann_Whitney_P"] = perf["Mann_Whitney_P"].astype(float)
    perf["Mean_score_positive"] = perf["Mean_score_positive"].astype(float)
    perf["Mean_score_negative"] = perf["Mean_score_negative"].astype(float)
    perf["Score_difference"] = perf["Score_difference"].astype(float)
    perf = perf.sort_values(["Role", "Dataset"]).reset_index(drop=True)

    scores = pd.concat(score_frames, ignore_index=True)

    gene_cols = [
        "Panel_rank",
        "Gene",
        "DESeq2_log2FC",
        "DESeq2_padj",
        "GSE159610_Log2FC_aneurysm_vs_control",
        "GSE159610_P",
        "GSE106520_Log2FC_aneurysm_vs_control",
        "GSE106520_P",
        "Bridge_Stouffer_Z",
        "Bridge_meta_P",
        "GSE36791_Log2FC_postrupture",
        "GSE36791_P",
        "GSE73378_Log2FC_pastSAH",
        "GSE73378_P",
    ]
    panel[gene_cols].to_csv(OUT_GENE_TSV, sep="\t", index=False)
    perf.to_csv(OUT_PERF_TSV, sep="\t", index=False)
    scores.to_csv(OUT_SCORE_TSV, sep="\t", index=False)
    _save_figure(panel, perf, OUT_FIG)

    summary = {
        "selection_rule": (
            "Top genes ranked by combined Stouffer Z across GSE159610 and GSE106520 after requiring directional concordance "
            "with the GSE122897 tissue shortlist in both blood bridge cohorts."
        ),
        "panel_size": int(len(panel)),
        "panel_genes": panel["Gene"].tolist(),
        "bridge_datasets": ["GSE159610", "GSE106520"],
        "context_only_datasets": ["GSE36791", "GSE73378"],
        "panel_figure": str(OUT_FIG),
        "gene_table": str(OUT_GENE_TSV),
        "performance_table": str(OUT_PERF_TSV),
        "sample_score_table": str(OUT_SCORE_TSV),
        "cohort_performance": perf.to_dict(orient="records"),
    }
    OUT_SUMMARY_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


if __name__ == "__main__":
    result = build_blood_bridge_panel_assets()
    print(json.dumps(result, indent=2))
