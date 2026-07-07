from __future__ import annotations

import csv
import gzip
import io
import json
import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import average_precision_score, roc_auc_score

BASE = Path(r"e:\Aneurysms Genomics Project")
OUT_DIR = BASE / "results" / "advanced_validation"
FIG_DIR = BASE / "results" / "figures" / "ultra_hd"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

RANDOM_SEED = 42
BOOTSTRAP_N = 2000
PERMUTATION_N = 2000

MODEL_FEATURES = ["DE_Signature_Score", "EC_Dysfunction_Score", "ECM_Remodeling_Score"]
ECM_GENES = ["MMP9", "MMP2", "TIMP1", "COL3A1", "COL1A1", "ELN", "FBN1", "LAMA3"]
EC_MARKERS = {
    "ANGPT2": 1.2,
    "ICAM1": 1.1,
    "VCAM1": 1.1,
    "THBS1": 0.5,
    "PERP": 0.5,
    "CXCL8": 1.0,
    "CDH5": -1.0,
    "ENG": -1.0,
    "KDR": -0.5,
    "TEK": -0.5,
    "VWF": -0.8,
    "PECAM1": 0.3,
    "LYVE1": 0.3,
    "SELE": 0.8,
    "SELP": 0.8,
    "NOS3": -0.6,
    "NOS2": 0.7,
}
TOP_DE_GENE_IDS = ["6611", "124220", "119587", "8685", "3909", "7045"]


def _load_de_signature_symbols() -> list[str]:
    shortlist = pd.read_csv(BASE / "results" / "phases_g" / "shortlist_evidence.tsv", sep="\t")
    if "Gene_ID" not in shortlist.columns or "symbol" not in shortlist.columns:
        raise RuntimeError("shortlist_evidence.tsv is missing Gene_ID or symbol columns.")
    shortlist["Gene_ID"] = shortlist["Gene_ID"].astype(str).str.replace(r"\.0$", "", regex=True)
    subset = shortlist[shortlist["Gene_ID"].isin(TOP_DE_GENE_IDS)].copy()
    symbols = [str(x).strip() for x in subset["symbol"].dropna().tolist() if str(x).strip()]
    symbols = [x for x in symbols if x not in {"nan", "NA", "None"}]
    ordered = []
    for symbol in symbols:
        if symbol not in ordered:
            ordered.append(symbol)
    return ordered


DE_SIGNATURE_GENES = _load_de_signature_symbols()


def _zscore(series: pd.Series) -> pd.Series:
    std = float(series.std(ddof=0))
    if std < 1e-8:
        return pd.Series(np.zeros(len(series)), index=series.index, dtype=float)
    return (series - float(series.mean())) / std


def _as_gene_expression(df: pd.DataFrame) -> pd.DataFrame:
    expr = df.copy()
    expr.index = expr.index.astype(str).str.strip()
    expr = expr[~expr.index.duplicated(keep="first")]
    expr.columns = [str(c).strip() for c in expr.columns]
    return expr


def _score_mean_signature(expr: pd.DataFrame, genes: list[str], log2_input: bool) -> tuple[pd.Series, list[str]]:
    hits = [g for g in genes if g in expr.index]
    if len(hits) < 2:
        raise RuntimeError(f"Too few genes present for mean signature: {hits}")
    values = expr.loc[hits]
    if not log2_input:
        values = np.log2(values + 1.0)
    score = values.mean(axis=0)
    return _zscore(score), hits


def _score_weighted_signature(expr: pd.DataFrame, weights: dict[str, float], log2_input: bool) -> tuple[pd.Series, list[str]]:
    accum = None
    hits: list[str] = []
    for gene, weight in weights.items():
        if gene not in expr.index:
            continue
        vals = expr.loc[gene]
        vals = vals if log2_input else np.log2(vals + 1.0)
        term = vals * float(weight)
        accum = term if accum is None else accum + term
        hits.append(gene)
    if accum is None or len(hits) < 2:
        raise RuntimeError(f"Too few genes present for weighted signature: {hits}")
    return _zscore(accum), hits


def _binary_metrics(y_true: np.ndarray, y_score: np.ndarray, rng: np.random.Generator) -> dict[str, float]:
    auc = float(roc_auc_score(y_true, y_score))
    auprc = float(average_precision_score(y_true, y_score))
    p_u = float(mannwhitneyu(y_score[y_true == 1], y_score[y_true == 0], alternative="two-sided").pvalue)

    boot = []
    n = len(y_true)
    idx = np.arange(n)
    for _ in range(BOOTSTRAP_N):
        sample_idx = rng.choice(idx, size=n, replace=True)
        y_b = y_true[sample_idx]
        if len(np.unique(y_b)) < 2:
            continue
        boot.append(roc_auc_score(y_b, y_score[sample_idx]))
    if boot:
        ci_low, ci_high = np.percentile(np.array(boot, dtype=float), [2.5, 97.5])
    else:
        ci_low, ci_high = np.nan, np.nan

    perm = []
    for _ in range(PERMUTATION_N):
        y_perm = rng.permutation(y_true)
        if len(np.unique(y_perm)) < 2:
            continue
        perm.append(roc_auc_score(y_perm, y_score))
    perm = np.array(perm, dtype=float)
    perm_p = float((1.0 + np.sum(perm >= auc)) / (1.0 + len(perm))) if len(perm) else np.nan

    return {
        "AUC": auc,
        "AUPRC": auprc,
        "MW_pvalue": p_u,
        "AUC_bootstrap_ci_low": float(ci_low),
        "AUC_bootstrap_ci_high": float(ci_high),
        "AUC_permutation_pvalue": perm_p,
    }


def _parse_geo_series_matrix(path: Path) -> tuple[pd.DataFrame, dict[str, list[list[str]]]]:
    meta: dict[str, list[list[str]]] = {}
    in_table = False
    header = None
    rows = []
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                in_table = False
                continue
            if in_table:
                parsed = next(csv.reader([line], delimiter="\t", quotechar='"'))
                if header is None:
                    header = parsed
                else:
                    rows.append(parsed)
                continue
            if line.startswith("!Sample_"):
                parts = next(csv.reader([line], delimiter="\t", quotechar='"'))
                key = parts[0].replace("!Sample_", "")
                meta.setdefault(key, []).append(parts[1:])
    if header is None:
        raise RuntimeError(f"Could not parse series matrix table from {path}")
    df = pd.DataFrame(rows, columns=header).rename(columns={header[0]: "ProbeID"})
    samples = header[1:]
    for sample in samples:
        df[sample] = pd.to_numeric(df[sample], errors="coerce")
    return df.set_index("ProbeID"), meta


def _load_gpl570_annotation() -> pd.DataFrame:
    path = BASE / "data" / "microarray" / "GSE26969" / "GPL570.annot.gz"
    lines = gzip.open(path, "rt", encoding="utf-8", errors="ignore").read().splitlines()
    begin = next(i for i, line in enumerate(lines) if line.startswith("!platform_table_begin"))
    end = next(i for i, line in enumerate(lines) if line.startswith("!platform_table_end"))
    annot = pd.read_csv(io.StringIO("\n".join(lines[begin + 1 : end])), sep="\t", low_memory=False)
    annot = annot.rename(columns={"ID": "ProbeID"})
    gene_col = "Gene symbol"
    annot["Gene"] = (
        annot[gene_col]
        .astype(str)
        .str.split("///")
        .str[0]
        .str.strip()
        .replace({"---": np.nan, "nan": np.nan, "": np.nan})
    )
    annot = annot.dropna(subset=["Gene"])[["ProbeID", "Gene"]].drop_duplicates("ProbeID")
    return annot


def _load_gse13353_gene_level() -> tuple[pd.DataFrame, pd.DataFrame]:
    series = BASE / "External Validation" / "GSE13353_series_matrix.txt.gz"
    probe_expr, meta_raw = _parse_geo_series_matrix(series)
    samples = probe_expr.columns.tolist()

    def _col(name: str, idx: int = 0) -> list[str]:
        rows = meta_raw.get(name, [])
        return rows[idx] if len(rows) > idx else [""] * len(samples)

    titles = _col("title", 0)
    chars = _col("characteristics_ch1", 0)
    phenotype = []
    for title, char in zip(titles, chars):
        text = f"{title} {char}".lower()
        if "unruptured" in text:
            phenotype.append("Unruptured")
        elif "ruptured" in text:
            phenotype.append("Ruptured")
        else:
            phenotype.append("Unknown")
    meta = pd.DataFrame({"Sample": samples, "Title": titles, "Group": phenotype})
    meta = meta[meta["Group"].isin(["Ruptured", "Unruptured"])].copy()
    annot = _load_gpl570_annotation()
    merged = probe_expr.reset_index().merge(annot, on="ProbeID", how="inner")
    gene_expr = merged.groupby("Gene")[meta["Sample"].tolist()].mean().sort_index()
    return _as_gene_expression(gene_expr), meta


def _load_gse54083() -> tuple[pd.DataFrame, pd.DataFrame, bool]:
    expr = pd.read_csv(
        BASE / "results" / "external_validation_gse54083" / "gse54083_gene_level_log2_expression.tsv",
        sep="\t",
        index_col=0,
    )
    meta = pd.read_csv(BASE / "results" / "external_validation_gse54083" / "gse54083_sample_metadata.tsv", sep="\t")
    meta = meta[meta["Group"].isin(["Ruptured", "Unruptured"])].copy()
    return _as_gene_expression(expr[meta["Sample"].tolist()]), meta, True


def _load_gse13353() -> tuple[pd.DataFrame, pd.DataFrame, bool]:
    expr, meta = _load_gse13353_gene_level()
    transform_needed = float(np.nanpercentile(expr.values, 95)) > 100.0
    return expr[meta["Sample"].tolist()], meta, not transform_needed


def _load_gse26969() -> tuple[pd.DataFrame, pd.DataFrame, bool]:
    expr = pd.read_csv(
        BASE / "results" / "external_validation_gse26969" / "gse26969_gene_level_expression.tsv",
        sep="\t",
        index_col=0,
    )
    meta = pd.read_csv(BASE / "results" / "external_validation_gse26969" / "gse26969_sample_metadata.tsv", sep="\t")
    return _as_gene_expression(expr[meta["Sample"].tolist()]), meta, True


def _load_gse36791() -> tuple[pd.DataFrame, pd.DataFrame, bool]:
    expr = pd.read_csv(
        BASE / "results" / "external_validation_gse36791" / "gse36791_gene_level_log2_expression.tsv",
        sep="\t",
        index_col=0,
    )
    meta = pd.read_csv(BASE / "results" / "external_validation_gse36791" / "gse36791_sample_metadata.tsv", sep="\t")
    return _as_gene_expression(expr[meta["Sample"].tolist()]), meta, True


def _load_gse159610() -> tuple[pd.DataFrame, pd.DataFrame, bool]:
    expr = pd.read_csv(
        BASE / "data" / "blood_rnaseq" / "GSE159610" / "GSE159610_raw_gene_counts_matrix.txt.gz",
        sep="\t",
        compression="gzip",
    )
    expr = expr[["gene_name"] + [c for c in expr.columns if c.startswith("DSA")]].copy()
    expr = expr.rename(columns={"gene_name": "Gene"})
    expr["Gene"] = expr["Gene"].astype(str).str.strip()
    expr = expr[(expr["Gene"] != "") & (expr["Gene"] != "nan")].copy()
    expr = expr.groupby("Gene", as_index=False)[[c for c in expr.columns if c.startswith("DSA")]].sum()
    expr = expr.set_index("Gene").sort_index()

    meta = pd.read_csv(BASE / "data" / "blood_rnaseq" / "GSE159610" / "sample_metadata_with_labels.tsv", sep="\t")
    meta = meta.copy().reset_index(drop=True)
    sample_cols = sorted([c for c in expr.columns if c.startswith("DSA")], key=lambda x: int(x.replace("DSA", "")))
    if len(sample_cols) < len(meta):
        raise RuntimeError("Not enough DSA columns to align GSE159610 metadata.")
    meta = meta.iloc[: len(sample_cols)].copy()
    meta["Sample"] = sample_cols[: len(meta)]
    meta["Group"] = meta["Status"].replace({"Aneurysm": "Aneurysm", "Control": "Control"})
    return _as_gene_expression(expr[meta["Sample"].tolist()]), meta, False


def _load_model_and_scaler():
    with open(BASE / "results" / "phase5_optimized_model.pkl", "rb") as handle:
        model = pickle.load(handle)
    with open(BASE / "results" / "phase5_optimized_scaler.pkl", "rb") as handle:
        scaler = pickle.load(handle)
    return model, scaler


def _evaluate_dataset(
    dataset: str,
    role: str,
    expr: pd.DataFrame,
    meta: pd.DataFrame,
    log2_input: bool,
    positive_group: str,
    negative_group: str,
    model,
    scaler,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, dict[str, object]]:
    meta = meta[meta["Group"].isin([positive_group, negative_group])].copy()
    meta["Sample"] = meta["Sample"].astype(str)
    expr = expr.loc[:, meta["Sample"].tolist()].copy()

    de_score, de_hits = _score_mean_signature(expr, DE_SIGNATURE_GENES, log2_input=log2_input)
    ec_score, ec_hits = _score_weighted_signature(expr, EC_MARKERS, log2_input=log2_input)
    ecm_score, ecm_hits = _score_mean_signature(expr, ECM_GENES, log2_input=log2_input)

    scores = pd.DataFrame(
        {
            "Sample": expr.columns,
            "DE_Signature_Score": de_score.reindex(expr.columns).values,
            "EC_Dysfunction_Score": ec_score.reindex(expr.columns).values,
            "ECM_Remodeling_Score": ecm_score.reindex(expr.columns).values,
        }
    )
    scores = meta.merge(scores, on="Sample", how="left")
    X = scores[MODEL_FEATURES].astype(float).values
    X_scaled = scaler.transform(X)
    scores["Transferred_Rupture_Probability"] = model.predict_proba(X_scaled)[:, 1]
    scores["Dataset"] = dataset
    scores["Role"] = role
    scores["Is_positive_group"] = (scores["Group"] == positive_group).astype(int)

    y_true = scores["Is_positive_group"].to_numpy(dtype=int)
    y_score = scores["Transferred_Rupture_Probability"].to_numpy(dtype=float)
    metrics = _binary_metrics(y_true, y_score, rng)

    summary = {
        "Dataset": dataset,
        "Role": role,
        "Positive_group": positive_group,
        "Negative_group": negative_group,
        "Samples": int(len(scores)),
        "N_positive_group": int((scores["Group"] == positive_group).sum()),
        "N_negative_group": int((scores["Group"] == negative_group).sum()),
        "DE_genes_measured": int(len(de_hits)),
        "EC_genes_measured": int(len(ec_hits)),
        "ECM_genes_measured": int(len(ecm_hits)),
        "DE_genes_used": ", ".join(de_hits),
        "EC_genes_used": ", ".join(ec_hits),
        "ECM_genes_used": ", ".join(ecm_hits),
        "Mean_probability_positive": float(scores.loc[scores["Group"] == positive_group, "Transferred_Rupture_Probability"].mean()),
        "Mean_probability_negative": float(scores.loc[scores["Group"] == negative_group, "Transferred_Rupture_Probability"].mean()),
        **metrics,
    }
    return scores, summary


def _save_transfer_figure(summary_df: pd.DataFrame, out_path: Path) -> None:
    plot_df = summary_df.copy()
    plot_df["Label"] = plot_df["Dataset"] + " | " + plot_df["Role"]
    plot_df = plot_df.sort_values("AUC", ascending=True)

    fig, axes = plt.subplots(1, 2, figsize=(16, 9))
    ax = axes[0]
    ax.barh(plot_df["Label"], plot_df["AUC"], color="#2f6db2", alpha=0.9)
    xerr = np.vstack(
        [
            plot_df["AUC"] - plot_df["AUC_bootstrap_ci_low"],
            plot_df["AUC_bootstrap_ci_high"] - plot_df["AUC"],
        ]
    )
    ax.errorbar(plot_df["AUC"], plot_df["Label"], xerr=xerr, fmt="none", ecolor="black", capsize=4, lw=1.2)
    ax.axvline(0.5, color="gray", linestyle="--", lw=1)
    ax.set_xlabel("Transferred model AUROC")
    ax.set_title("A. External transfer of the three-score rupture model")

    ax2 = axes[1]
    width = 0.38
    y = np.arange(len(plot_df))
    ax2.barh(y - width / 2, plot_df["Mean_probability_positive"], height=width, color="#b54a4a", label="Positive group")
    ax2.barh(y + width / 2, plot_df["Mean_probability_negative"], height=width, color="#5a7f3d", label="Comparator group")
    ax2.set_yticks(y)
    ax2.set_yticklabels(plot_df["Label"])
    ax2.set_xlabel("Mean transferred rupture probability")
    ax2.set_title("B. Mean transferred probabilities by group")
    ax2.legend(frameon=False)

    fig.suptitle("Advanced public-data signature transfer across tissue and blood cohorts", fontsize=16)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    model, scaler = _load_model_and_scaler()
    rng = np.random.default_rng(RANDOM_SEED)

    dataset_specs = [
        ("GSE54083", "direct tissue rupture transfer", _load_gse54083, "Ruptured", "Unruptured"),
        ("GSE13353", "direct tissue rupture transfer", _load_gse13353, "Ruptured", "Unruptured"),
        ("GSE26969", "orthogonal tissue context", _load_gse26969, "Aneurysm", "Control"),
        ("GSE36791", "orthogonal blood rupture/SAH context", _load_gse36791, "Ruptured", "Control"),
        ("GSE159610", "orthogonal blood aneurysm context", _load_gse159610, "Aneurysm", "Control"),
    ]

    sample_frames = []
    summary_rows = []

    for dataset, role, loader, positive_group, negative_group in dataset_specs:
        expr, meta, log2_input = loader()
        scores, summary = _evaluate_dataset(
            dataset=dataset,
            role=role,
            expr=expr,
            meta=meta,
            log2_input=log2_input,
            positive_group=positive_group,
            negative_group=negative_group,
            model=model,
            scaler=scaler,
            rng=rng,
        )
        sample_frames.append(scores)
        summary_rows.append(summary)

    sample_df = pd.concat(sample_frames, ignore_index=True)
    summary_df = pd.DataFrame(summary_rows).sort_values("AUC", ascending=False).reset_index(drop=True)

    sample_df.to_csv(OUT_DIR / "external_signature_transfer_sample_scores.tsv", sep="\t", index=False)
    summary_df.to_csv(OUT_DIR / "external_signature_transfer_summary.tsv", sep="\t", index=False)
    (OUT_DIR / "external_signature_transfer_summary.json").write_text(
        json.dumps(summary_rows, indent=2),
        encoding="utf-8",
    )

    figure_path = FIG_DIR / "external_signature_transfer_ultra_hd.png"
    _save_transfer_figure(summary_df, figure_path)

    print("Saved:")
    print(f"  - {OUT_DIR / 'external_signature_transfer_sample_scores.tsv'}")
    print(f"  - {OUT_DIR / 'external_signature_transfer_summary.tsv'}")
    print(f"  - {figure_path}")
    print("\nSummary:")
    print(
        summary_df[
            [
                "Dataset",
                "Role",
                "Samples",
                "AUC",
                "AUC_bootstrap_ci_low",
                "AUC_bootstrap_ci_high",
                "AUC_permutation_pvalue",
                "Mean_probability_positive",
                "Mean_probability_negative",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
