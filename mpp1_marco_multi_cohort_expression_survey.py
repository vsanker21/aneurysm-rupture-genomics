"""
MPP1 and MARCO — multi-cohort public expression survey.

Programmatically loads expression (with GEO download fallbacks where needed),
computes group means, Welch log2FC, Cohen's d, and direction concordance vs
GSE122897 tissue discovery (ruptured vs unruptured).

Comparisons are labelled explicitly; orthogonal blood cohorts (post-rupture SAH,
prior SAH history) are reported separately and not pooled with prospective
rupture tissue meta-analysis.
"""

from __future__ import annotations

import gzip
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm, ttest_ind

from advanced_public_signature_transfer import (
    _as_gene_expression,
    _load_gse13353_gene_level,
    _load_gse54083,
    _load_gse26969,
    _load_gse36791,
)
from blood_bridge_panel import _load_gse159610_expression
from blood_bridge_panel_geoparse_check import ensure_gse106520_expression, ensure_gse159610_counts
from public_cohort_context_utils import add_gse122897_reference, compute_gene_level_welch

try:
    from project_paths import get_project_root
    BASE = get_project_root()
except ImportError:
    BASE = Path(__file__).resolve().parent
OUT_DIR = BASE / "results" / "mpp1_marco_survey"
FIG_DIR = BASE / "results" / "figures" / "ultra_hd"
GENES = ["MPP1", "MARCO"]


@dataclass(frozen=True)
class ComparisonSpec:
    cohort_id: str
    comparison_id: str
    biospecimen: str
    matrix_type: str
    evidence_tier: str
    positive_label: str
    negative_label: str
    case_definition: str
    control_definition: str
    interpretation_note: str
    loader: Callable[[Path], tuple[pd.DataFrame, pd.DataFrame]]


def _cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    pooled = np.sqrt(((len(a) - 1) * np.var(a, ddof=1) + (len(b) - 1) * np.var(b, ddof=1)) / (len(a) + len(b) - 2))
    if pooled < 1e-12:
        return float("nan")
    return float((np.mean(a) - np.mean(b)) / pooled)


def _discovery_signs(base: Path) -> pd.Series:
    panel = pd.read_csv(base / "results" / "blood_bridge_panel" / "blood_bridge_panel_gene_summary.tsv", sep="\t")
    signs = np.sign(pd.to_numeric(panel.set_index("Gene")["DESeq2_log2FC"], errors="coerce"))
    de = pd.read_csv(base / "results" / "gse122897_deseq2_results.tsv", sep="\t")
    mapping = pd.read_csv(base / "results" / "gene_id_mapping_unified.tsv", sep="\t")
    mapping["Gene_ID"] = mapping["Original_ID"].astype(str).str.replace(r"\.0$", "", regex=True)
    de["Gene_ID"] = de["Gene_ID"].astype(str).str.replace(r"\.0$", "", regex=True)
    de = de.merge(mapping[["Gene_ID", "Gene_Symbol"]].drop_duplicates("Gene_ID"), on="Gene_ID", how="left")
    for gene in GENES:
        row = de[de["Gene_Symbol"] == gene]
        if not row.empty and gene not in signs.index:
            signs[gene] = np.sign(float(row.iloc[0]["log2FoldChange"]))
    return signs.reindex(GENES).astype(float)


def _align(expr: pd.DataFrame, meta: pd.DataFrame, positive: str, negative: str) -> tuple[pd.DataFrame, list[str], list[str]]:
    meta = meta[meta["Group"].isin([positive, negative])].copy()
    meta = meta.drop_duplicates(subset=["Sample"])
    cols = [c for c in expr.columns if c in set(meta["Sample"])]
    expr = expr.loc[:, cols].copy()
    meta = meta.set_index("Sample").loc[cols]
    pos_cols = meta.index[meta["Group"] == positive].tolist()
    neg_cols = meta.index[meta["Group"] == negative].tolist()
    return expr, pos_cols, neg_cols


def _analyze_genes(
    expr: pd.DataFrame,
    pos_cols: list[str],
    neg_cols: list[str],
    spec: ComparisonSpec,
    discovery: pd.Series,
    *,
    log2_input: bool,
) -> pd.DataFrame:
    rows = []
    for gene in GENES:
        if gene not in expr.index:
            rows.append(
                {
                    "Gene": gene,
                    "Cohort": spec.cohort_id,
                    "Comparison": spec.comparison_id,
                    "Status": "gene_absent",
                }
            )
            continue
        vals = pd.to_numeric(expr.loc[gene, pos_cols + neg_cols], errors="coerce")
        pv = vals[pos_cols].values.astype(float)
        nv = vals[neg_cols].values.astype(float)
        if not log2_input:
            pv = np.log2(pv + 1.0)
            nv = np.log2(nv + 1.0)
        if np.sum(~np.isnan(pv)) < 2 or np.sum(~np.isnan(nv)) < 2:
            rows.append(
                {
                    "Gene": gene,
                    "Cohort": spec.cohort_id,
                    "Comparison": spec.comparison_id,
                    "Status": "insufficient_samples",
                }
            )
            continue
        t_stat, p_val = ttest_ind(pv[~np.isnan(pv)], nv[~np.isnan(nv)], equal_var=False, nan_policy="omit")
        log2fc = float(np.nanmean(pv) - np.nanmean(nv))
        disc = float(discovery.get(gene, np.nan))
        rows.append(
            {
                "Gene": gene,
                "Cohort": spec.cohort_id,
                "Comparison": spec.comparison_id,
                "Biospecimen": spec.biospecimen,
                "Matrix": spec.matrix_type,
                "Evidence_tier": spec.evidence_tier,
                "Status": "ok",
                "N_positive": len(pos_cols),
                "N_negative": len(neg_cols),
                "Mean_positive": float(np.nanmean(pv)),
                "Mean_negative": float(np.nanmean(nv)),
                "Log2FC_positive_minus_negative": log2fc,
                "Cohens_d": _cohens_d(pv, nv),
                "T_statistic": float(t_stat) if pd.notna(t_stat) else np.nan,
                "P_value": float(p_val) if pd.notna(p_val) else np.nan,
                "Discovery_GSE122897_sign": disc,
                "Direction_match_discovery": bool(np.sign(log2fc) == np.sign(disc)) if disc not in (0, np.nan) else False,
                "Signed_log2FC_vs_discovery": log2fc * disc if disc not in (0, np.nan) else np.nan,
                "Case_definition": spec.case_definition,
                "Control_definition": spec.control_definition,
                "Interpretation_note": spec.interpretation_note,
            }
        )
    return pd.DataFrame(rows)


def _gse122897_entrez_to_symbol(base: Path) -> dict[str, str]:
    mapping = pd.read_csv(base / "results" / "gene_id_mapping_unified.tsv", sep="\t")
    mapping["Gene_ID"] = mapping["Original_ID"].astype(str).str.replace(r"\.0$", "", regex=True)
    entrez_map = (
        mapping.dropna(subset=["Gene_Symbol"])
        .drop_duplicates("Gene_ID")
        .set_index("Gene_ID")["Gene_Symbol"]
        .astype(str)
        .to_dict()
    )
    shortlist = pd.read_csv(base / "results" / "phases_g" / "shortlist_evidence.tsv", sep="\t")
    if {"Gene_ID", "symbol"}.issubset(shortlist.columns):
        for _, row in shortlist.dropna(subset=["Gene_ID", "symbol"]).iterrows():
            gid = str(row["Gene_ID"]).replace(".0", "")
            entrez_map[gid] = str(row["symbol"]).strip()
    entrez_map.setdefault("4354", "MPP1")
    entrez_map.setdefault("8685", "MARCO")
    return entrez_map


def _maybe_log2_transform(expr: pd.DataFrame) -> pd.DataFrame:
    if float(np.nanpercentile(expr.values, 95)) > 50.0:
        return np.log2(expr + 1.0)
    return expr


def _load_gse122897(base: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    tpm = pd.read_csv(
        base / "data" / "sc_rnaseq" / "GSE122897" / "GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz",
        sep="\t",
        compression="gzip",
        index_col=0,
    )
    entrez_map = _gse122897_entrez_to_symbol(base)
    tpm.index = tpm.index.astype(str).str.replace(r"\.0$", "", regex=True)
    tpm.index = tpm.index.map(lambda x: entrez_map.get(x, x))
    tpm = tpm.groupby(level=0).mean()
    expr = np.log2(tpm + 1.0)
    meta = pd.read_csv(base / "data" / "sc_rnaseq" / "GSE122897" / "sample_metadata.tsv", sep="\t")
    return _as_gene_expression(expr), meta


def _load_gse13353(base: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr, meta = _load_gse13353_gene_level()
    return _as_gene_expression(_maybe_log2_transform(expr)), meta


def _load_gse15629(base: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = pd.read_csv(
        base / "results" / "external_validation_gse15629" / "gse15629_gene_level_expression.tsv",
        sep="\t",
        index_col=0,
    )
    meta = pd.read_csv(
        base / "results" / "external_validation_gse15629" / "gse15629_sample_metadata.tsv",
        sep="\t",
    )
    return _as_gene_expression(_maybe_log2_transform(expr)), meta


def _load_gse106520(base: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr, meta = ensure_gse106520_expression()
    return _as_gene_expression(expr), meta


def _load_gse159610(base: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    ensure_gse159610_counts()
    return _load_gse159610_expression(base)


def _parse_series_matrix_metadata(path: Path) -> pd.DataFrame:
    import gzip as _gzip

    lines = _gzip.open(path, "rt", encoding="utf-8", errors="ignore").read().splitlines()
    sample_ids: list[str] = []
    titles: list[str] = []
    for line in lines:
        if line.startswith("!Sample_geo_accession"):
            sample_ids = [x.strip().strip('"') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_title"):
            titles = [x.strip().strip('"') for x in line.split("\t")[1:]]
    return pd.DataFrame({"Sample": sample_ids, "Title": titles})


def _parse_series_matrix_sample_table(path: Path) -> pd.DataFrame:
    """Parse GEO series matrix sample metadata including stacked characteristics rows."""
    import gzip as _gzip

    lines = _gzip.open(path, "rt", encoding="utf-8", errors="ignore").read().splitlines()
    sample_ids: list[str] = []
    titles: list[str] = []
    char_rows: list[list[str]] = []
    for line in lines:
        if line.startswith("!Sample_geo_accession"):
            sample_ids = [x.strip().strip('"') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_title"):
            titles = [x.strip().strip('"') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1"):
            char_rows.append([x.strip().strip('"') for x in line.split("\t")[1:]])

    rows = []
    for i, sample in enumerate(sample_ids):
        rec: dict[str, str] = {"Sample": sample, "Title": titles[i] if i < len(titles) else ""}
        for vals in char_rows:
            if i >= len(vals):
                continue
            text = vals[i]
            if ":" in text:
                key, value = text.split(":", 1)
                rec[key.strip().lower().replace(" ", "_")] = value.strip()
        rows.append(rec)
    return pd.DataFrame(rows)


def _load_gse6551_gpl570(base: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    from advanced_public_signature_transfer import _load_gpl570_annotation, _parse_geo_series_matrix

    matrix = base / "data" / "microarray" / "GSE6551" / "GSE6551-GPL570_series_matrix.txt.gz"
    probe_expr, _ = _parse_geo_series_matrix(matrix)
    annot = _load_gpl570_annotation()
    merged = probe_expr.reset_index().merge(annot, on="ProbeID", how="inner")
    gene_expr = merged.groupby("Gene")[probe_expr.columns.tolist()].mean().sort_index()
    meta = _parse_series_matrix_sample_table(matrix)
    tissue = meta.get("tissue", pd.Series([""] * len(meta))).astype(str).str.lower()
    meta["Group_aneurysm_vs_artery"] = np.where(
        tissue.str.contains("non-aneurysm"),
        "ArteryControl",
        np.where(tissue.str.contains("aneurysm"), "Aneurysm", "Unknown"),
    )
    disease = meta.get("source_disease_status", pd.Series([""] * len(meta))).astype(str).str.lower()
    meta["Group_rupture"] = np.where(
        disease.str.contains("unruptured"),
        "Unruptured",
        np.where(disease.str.contains("ruptured"), "Ruptured", "Unknown"),
    )
    return _as_gene_expression(_maybe_log2_transform(gene_expr)), meta


def _load_gse50867(base: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    from advanced_public_signature_transfer import _parse_geo_series_matrix

    matrix = base / "data" / "microarray" / "GSE50867" / "GSE50867_series_matrix.txt.gz"
    probe_expr, _ = _parse_geo_series_matrix(matrix)
    meta = _parse_series_matrix_metadata(matrix)
    title = meta["Title"].astype(str).str.lower()
    meta["Group"] = np.where(
        title.str.startswith("ruptured"),
        "Ruptured",
        np.where(title.str.startswith("unruptured"), "Unruptured", np.where(title.str.startswith("control"), "Control", "Unknown")),
    )
    meta = meta[meta["Group"].isin(["Ruptured", "Unruptured", "Control"])].copy()
    return _as_gene_expression(probe_expr), meta


def _load_gse66238(base: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    tpm = pd.read_csv(
        base / "data" / "sc_rnaseq" / "GSE66238" / "GSE66238_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz",
        sep="\t",
        compression="gzip",
        index_col=0,
    )
    entrez_map = _gse122897_entrez_to_symbol(base)
    tpm.index = tpm.index.astype(str).str.replace(r"\.0$", "", regex=True)
    tpm.index = tpm.index.map(lambda x: entrez_map.get(x, x))
    expr = np.log2(tpm.groupby(level=0).mean() + 1.0)
    meta = _parse_series_matrix_metadata(
        base / "data" / "sc_rnaseq" / "GSE66238" / "GSE66238_series_matrix.txt.gz"
    )
    meta["Group"] = np.where(
        meta["Title"].str.contains("BrainAneurysm", case=False, na=False),
        "Aneurysm",
        np.where(meta["Title"].str.contains("Control", case=False, na=False), "Control", "Unknown"),
    )
    meta = meta[meta["Group"].isin(["Aneurysm", "Control"])].copy()
    return _as_gene_expression(expr), meta


def _load_gse158558(base: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    tpm = pd.read_csv(
        base / "data" / "sc_rnaseq" / "GSE158558" / "GSE158558_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz",
        sep="\t",
        compression="gzip",
        index_col=0,
    )
    entrez_map = _gse122897_entrez_to_symbol(base)
    tpm.index = tpm.index.astype(str).str.replace(r"\.0$", "", regex=True)
    tpm.index = tpm.index.map(lambda x: entrez_map.get(x, x))
    expr = np.log2(tpm.groupby(level=0).mean() + 1.0)
    meta = _parse_series_matrix_metadata(
        base / "data" / "sc_rnaseq" / "GSE158558" / "GSE158558_series_matrix.txt.gz"
    )
    meta["Group"] = "Aneurysm"
    return _as_gene_expression(expr), meta


def _load_gse75436(base: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    from advanced_public_signature_transfer import _load_gpl570_annotation, _parse_geo_series_matrix

    matrix = base / "data" / "microarray" / "GSE75436" / "GSE75436_series_matrix.txt.gz"
    probe_expr, meta_raw = _parse_geo_series_matrix(matrix)
    annot = _load_gpl570_annotation()
    merged = probe_expr.reset_index().merge(annot, on="ProbeID", how="inner")
    gene_expr = merged.groupby("Gene")[probe_expr.columns.tolist()].mean().sort_index()
    samples = probe_expr.columns.tolist()

    def _col(name: str, idx: int = 0) -> list[str]:
        rows = meta_raw.get(name, [])
        return rows[idx] if len(rows) > idx else [""] * len(samples)

    titles = _col("title", 0)
    tissue_rows = meta_raw.get("characteristics_ch1", [])
    tissue_vals = tissue_rows[0] if tissue_rows else [""] * len(samples)
    groups = []
    for tissue in tissue_vals:
        text = str(tissue).lower()
        if "intracranial aneurysm" in text:
            groups.append("Aneurysm")
        elif "superficial temporal" in text or "temporal artery" in text:
            groups.append("Control")
        else:
            groups.append("Unknown")
    meta = pd.DataFrame({"Sample": samples, "Title": titles, "Group": groups})
    meta = meta[meta["Group"].isin(["Aneurysm", "Control"])].copy()
    return _as_gene_expression(_maybe_log2_transform(gene_expr)), meta


def _load_gse73378(base: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = pd.read_csv(
        base / "results" / "external_validation_gse73378" / "gse73378_gene_level_expression.tsv",
        sep="\t",
        index_col=0,
    )
    meta = pd.read_csv(
        base / "results" / "external_validation_gse73378" / "gse73378_sample_metadata.tsv",
        sep="\t",
    )
    if "Sample" not in meta.columns:
        meta["Sample"] = meta["Title"].astype(str).str.split().str[-1].str.strip()
    meta["Group"] = meta["Group"].replace({"Case": "Prior_SAH", "Control": "Control"})
    return _as_gene_expression(expr), meta


def _stouffer_by_gene(df: pd.DataFrame, gene: str, subset: pd.DataFrame) -> dict[str, float]:
    sub = subset[(subset["Gene"] == gene) & (subset["Status"] == "ok")].copy()
    if sub.empty:
        return {"Meta_Z": np.nan, "Meta_P": np.nan, "N_cohorts": 0}
    zs = []
    for _, row in sub.iterrows():
        p = float(np.clip(row["P_value"], 1e-300, 1.0))
        z = norm.isf(p / 2.0) * np.sign(row["Log2FC_positive_minus_negative"])
        zs.append(z)
    meta_z = float(np.sum(zs) / np.sqrt(len(zs)))
    meta_p = float(2.0 * norm.sf(abs(meta_z)))
    n_concordant = int((sub["Direction_match_discovery"] == True).sum())  # noqa: E712
    return {
        "Meta_Z": meta_z,
        "Meta_P": meta_p,
        "N_cohorts": len(sub),
        "N_direction_concordant": n_concordant,
        "Fraction_concordant": n_concordant / len(sub) if len(sub) else np.nan,
        "Mean_signed_log2FC": float(sub["Signed_log2FC_vs_discovery"].mean()),
    }


def _plot_forest(results: pd.DataFrame, out_path: Path) -> None:
    ok = results[results["Status"] == "ok"].copy()
    if ok.empty:
        return
    ok["Label"] = ok["Cohort"] + "\n" + ok["Comparison"].str.replace("_", " ")
    ok = ok.sort_values(["Gene", "Evidence_tier", "Cohort"])

    fig, axes = plt.subplots(1, 2, figsize=(14, max(6, 0.45 * len(ok) / 2)), sharey=False)
    colors = {"tier1_direct_rupture_tissue": "#b54a4a", "tier2_tissue_context": "#e09b4d",
              "tier4_blood_aneurysm_control": "#4b7bb8", "tier3_blood_rupture_sah_context": "#7a5cad"}

    for ax, gene in zip(axes, GENES):
        sub = ok[ok["Gene"] == gene]
        y = np.arange(len(sub))
        for idx, (_, row) in enumerate(sub.iterrows()):
            fc = row["Log2FC_positive_minus_negative"]
            color = colors.get(row["Evidence_tier"], "#666666")
            ax.barh(idx, fc, color=color, alpha=0.85, height=0.6)
        ax.axvline(0, color="black", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(sub["Label"], fontsize=8)
        ax.set_xlabel("Log2FC (positive group − negative group)")
        ax.set_title(f"{gene} across public cohorts")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _build_specs(base: Path) -> list[ComparisonSpec]:
    def wrap(loader: Callable[[Path], tuple[pd.DataFrame, pd.DataFrame]]) -> Callable[[Path], tuple[pd.DataFrame, pd.DataFrame]]:
        return loader

    return [
        ComparisonSpec(
            "GSE122897",
            "ruptured_vs_unruptured",
            "aneurysm_wall_tissue",
            "RNA-seq log2(TPM+1)",
            "tier1_direct_rupture_tissue",
            "Ruptured",
            "Unruptured",
            "Ruptured aneurysm wall",
            "Unruptured aneurysm wall",
            "Primary discovery cohort for tissue rupture signature.",
            wrap(_load_gse122897),
        ),
        ComparisonSpec(
            "GSE54083",
            "ruptured_vs_unruptured",
            "aneurysm_wall_tissue",
            "microarray log2",
            "tier1_direct_rupture_tissue",
            "Ruptured",
            "Unruptured",
            "Ruptured aneurysm wall",
            "Unruptured aneurysm wall",
            "Independent microarray replication of direct rupture contrast.",
            lambda b: _load_gse54083()[:2],
        ),
        ComparisonSpec(
            "GSE13353",
            "ruptured_vs_unruptured",
            "aneurysm_wall_tissue",
            "microarray log2",
            "tier1_direct_rupture_tissue",
            "Ruptured",
            "Unruptured",
            "Ruptured aneurysm wall",
            "Unruptured aneurysm wall",
            "Small direct rupture tissue support cohort.",
            wrap(_load_gse13353),
        ),
        ComparisonSpec(
            "GSE15629",
            "ruptured_vs_unruptured",
            "aneurysm_wall_tissue",
            "microarray log2",
            "tier2_tissue_context",
            "Ruptured",
            "Unruptured",
            "Ruptured vessel wall",
            "Unruptured vessel wall",
            "Mixed-design cohort; rupture subset only (excludes controls).",
            wrap(_load_gse15629),
        ),
        ComparisonSpec(
            "GSE15629",
            "aneurysm_any_vs_control",
            "aneurysm_wall_tissue",
            "microarray log2",
            "tier2_tissue_context",
            "AneurysmAny",
            "Control",
            "All ruptured + unruptured aneurysm samples",
            "Non-aneurysm vessel controls",
            "Orthogonal aneurysm presence context; not a pure rupture contrast.",
            wrap(_load_gse15629),
        ),
        ComparisonSpec(
            "GSE26969",
            "aneurysm_vs_control_vessel",
            "aneurysm_wall_tissue",
            "microarray log2",
            "tier2_tissue_context",
            "Aneurysm",
            "Control",
            "Unruptured aneurysm tissue",
            "Superficial temporal artery (normal vessel)",
            "Unruptured aneurysm vs normal artery; rupture status not tested.",
            lambda b: _load_gse26969()[:2],
        ),
        ComparisonSpec(
            "GSE159610",
            "aneurysm_vs_control",
            "whole_blood",
            "RNA-seq log2(CPM+1)",
            "tier4_blood_aneurysm_control",
            "Aneurysm",
            "Control",
            "Intracranial aneurysm (blood)",
            "Healthy control (blood)",
            "Blood bridge: aneurysm detection context, not rupture status.",
            wrap(_load_gse159610),
        ),
        ComparisonSpec(
            "GSE106520",
            "aneurysm_vs_control",
            "neutrophil_blood",
            "RNA-seq log2(FPKM+1)",
            "tier4_blood_aneurysm_control",
            "Aneurysm",
            "Control",
            "Unruptured intracranial aneurysm (neutrophils)",
            "Healthy control (neutrophils)",
            "Cell-type-specific blood bridge; not ruptured vs unruptured.",
            wrap(_load_gse106520),
        ),
        ComparisonSpec(
            "GSE36791",
            "ruptured_sah_vs_control",
            "whole_blood",
            "microarray log2",
            "tier3_blood_rupture_sah_context",
            "Ruptured",
            "Control",
            "Ruptured aneurysm with SAH",
            "Headache controls without aneurysm rupture",
            "Post-rupture systemic blood response; do not treat as prospective rupture validation.",
            lambda b: _load_gse36791()[:2],
        ),
        ComparisonSpec(
            "GSE73378",
            "prior_sah_vs_control",
            "whole_blood",
            "microarray log2",
            "tier3_blood_rupture_sah_context",
            "Prior_SAH",
            "Control",
            "History of aneurysmal SAH",
            "Controls without SAH history",
            "Long-horizon post-SAH blood context; orthogonal to unruptured risk.",
            wrap(_load_gse73378),
        ),
        ComparisonSpec(
            "GSE66238",
            "aneurysm_vs_control_vessel",
            "aneurysm_wall_tissue",
            "RNA-seq log2(TPM+1)",
            "tier2_tissue_context",
            "Aneurysm",
            "Control",
            "Cerebral aneurysm dome",
            "Superficial temporal artery (control vessel)",
            "Not ruptured vs unruptured; aneurysm vs STA control from cuffdiff design.",
            wrap(_load_gse66238),
        ),
        ComparisonSpec(
            "GSE158558",
            "aneurysm_wall_case_series",
            "aneurysm_wall_tissue",
            "RNA-seq log2(TPM+1)",
            "tier2_tissue_context",
            "Aneurysm",
            "Aneurysm",
            "Intracranial aneurysm wall (n=7)",
            "Same group (descriptive levels only)",
            "No rupture/control contrast in GEO metadata; reports wall expression levels only.",
            wrap(_load_gse158558),
        ),
        ComparisonSpec(
            "GSE6551",
            "aneurysm_vs_artery",
            "aneurysm_wall_tissue",
            "microarray log2 (GPL570)",
            "tier2_tissue_context",
            "Aneurysm",
            "ArteryControl",
            "Pooled aneurysm wall RNA",
            "Contralateral non-aneurysm artery (paired vessel)",
            "Legacy 2007 tissue array; aneurysm vs contralateral artery (not healthy control).",
            wrap(_load_gse6551_gpl570),
        ),
        ComparisonSpec(
            "GSE6551",
            "ruptured_vs_unruptured",
            "aneurysm_wall_tissue",
            "microarray log2 (GPL570)",
            "tier1_direct_rupture_tissue",
            "Ruptured",
            "Unruptured",
            "Aneurysm tissue pools with ruptured status",
            "Aneurysm tissue pools with unruptured status",
            "GSE6551 rupture labels on aneurysm-tissue pools only (small n, pooled RNA).",
            wrap(_load_gse6551_gpl570),
        ),
        ComparisonSpec(
            "GSE50867",
            "ruptured_vs_unruptured",
            "blood_plasma",
            "miRNA array (GPL17725)",
            "tier3_blood_rupture_sah_context",
            "Ruptured",
            "Unruptured",
            "Plasma pool — ruptured IA",
            "Plasma pool — unruptured IA",
            "Circulating miRNA pools (GPL17725); MPP1/MARCO mRNA not on platform.",
            wrap(_load_gse50867),
        ),
        ComparisonSpec(
            "GSE50867",
            "aneurysm_vs_control",
            "blood_plasma",
            "miRNA array (GPL17725)",
            "tier4_blood_aneurysm_control",
            "AneurysmAny",
            "Control",
            "Pooled ruptured + unruptured plasma",
            "Healthy control plasma",
            "Combined IA pools vs control; miRNA-only platform.",
            wrap(_load_gse50867),
        ),
        ComparisonSpec(
            "GSE75436",
            "aneurysm_vs_control_vessel",
            "aneurysm_wall_tissue",
            "microarray log2",
            "tier2_tissue_context",
            "Aneurysm",
            "Control",
            "Intracranial aneurysm tissue",
            "Matched superficial temporal artery",
            "Tissue microarray; paired aneurysm vs STA (not blood, not rupture-stratified).",
            wrap(_load_gse75436),
        ),
    ]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    discovery = _discovery_signs(BASE)

    all_rows: list[pd.DataFrame] = []
    sample_rows: list[dict[str, object]] = []

    for spec in _build_specs(BASE):
        key = f"{spec.cohort_id}|{spec.comparison_id}"
        print(f"Processing {key}...")
        try:
            expr, meta = spec.loader(BASE)
        except Exception as exc:  # noqa: BLE001
            all_rows.append(
                pd.DataFrame(
                    [
                        {
                            "Gene": gene,
                            "Cohort": spec.cohort_id,
                            "Comparison": spec.comparison_id,
                            "Status": f"load_failed: {exc}",
                        }
                        for gene in GENES
                    ]
                )
            )
            continue

        if spec.comparison_id == "aneurysm_any_vs_control" and spec.cohort_id == "GSE15629":
            meta = meta.copy()
            meta.loc[meta["Group"].isin(["Ruptured", "Unruptured"]), "Group"] = "AneurysmAny"

        if spec.cohort_id == "GSE6551":
            if spec.comparison_id == "aneurysm_vs_artery":
                meta = meta.copy()
                meta["Group"] = meta["Group_aneurysm_vs_artery"]
            elif spec.comparison_id == "ruptured_vs_unruptured":
                meta = meta[meta["Group_aneurysm_vs_artery"] == "Aneurysm"].copy()
                meta["Group"] = meta["Group_rupture"]

        if spec.cohort_id == "GSE50867" and spec.comparison_id == "aneurysm_vs_control":
            meta = meta.copy()
            meta.loc[meta["Group"].isin(["Ruptured", "Unruptured"]), "Group"] = "AneurysmAny"

        if spec.comparison_id == "aneurysm_wall_case_series":
            result = pd.DataFrame(
                [
                    {
                        "Gene": gene,
                        "Cohort": spec.cohort_id,
                        "Comparison": spec.comparison_id,
                        "Status": "no_contrast",
                        "Interpretation_note": spec.interpretation_note,
                    }
                    for gene in GENES
                ]
            )
            all_rows.append(result)
            continue

        expr, pos_cols, neg_cols = _align(expr, meta, spec.positive_label, spec.negative_label)
        log2_input = True

        result = _analyze_genes(expr, pos_cols, neg_cols, spec, discovery, log2_input=log2_input)
        all_rows.append(result)

        for gene in GENES:
            if gene not in expr.index:
                continue
            for sample in pos_cols + neg_cols:
                val = float(pd.to_numeric(expr.loc[gene, sample], errors="coerce"))
                if not log2_input:
                    val = float(np.log2(val + 1.0))
                grp = spec.positive_label if sample in pos_cols else spec.negative_label
                sample_rows.append(
                    {
                        "Gene": gene,
                        "Cohort": spec.cohort_id,
                        "Comparison": spec.comparison_id,
                        "Sample": sample,
                        "Group": grp,
                        "Expression_log2_scale": val,
                    }
                )

    results = pd.concat(all_rows, ignore_index=True)
    ok_only = pd.concat(all_rows, ignore_index=True)
    ok_part = ok_only[ok_only["Status"] == "ok"].copy()
    failed = ok_only[ok_only["Status"] != "ok"].copy()
    if not ok_part.empty:
        ok_part = add_gse122897_reference(
            ok_part,
            base=BASE,
            effect_col="Log2FC_positive_minus_negative",
        )
    full = pd.concat([ok_part, failed], ignore_index=True, sort=False)

    full.to_csv(OUT_DIR / "mpp1_marco_all_cohort_results.tsv", sep="\t", index=False)
    pd.DataFrame(sample_rows).to_csv(OUT_DIR / "mpp1_marco_sample_level_expression.tsv", sep="\t", index=False)

    meta_blocks = []
    families = {
        "tissue_ruptured_vs_unruptured": (
            "Direct rupture tissue meta-analysis (Tier 1 + GSE15629 rupture subset)",
            ["GSE122897", "GSE54083", "GSE13353"],
            "ruptured_vs_unruptured",
            "tier1_direct_rupture_tissue",
        ),
        "tissue_rupture_including_gse15629": (
            "Tissue rupture meta including GSE15629 rupture subset",
            ["GSE122897", "GSE54083", "GSE13353", "GSE15629"],
            "ruptured_vs_unruptured",
            None,
        ),
        "blood_aneurysm_vs_control": (
            "Blood aneurysm vs control bridge meta",
            ["GSE159610", "GSE106520"],
            "aneurysm_vs_control",
            "tier4_blood_aneurysm_control",
        ),
        "blood_post_rupture_context": (
            "Blood post-rupture / prior SAH context (orthogonal; not prospective rupture)",
            ["GSE36791", "GSE73378"],
            None,
            "tier3_blood_rupture_sah_context",
        ),
    }
    ok = full[full["Status"] == "ok"]
    for gene in GENES:
        for fam_id, (label, cohorts, comp, tier) in families.items():
            sub = ok[ok["Gene"] == gene]
            if comp:
                sub = sub[sub["Comparison"] == comp]
            if tier:
                sub = sub[sub["Evidence_tier"] == tier]
            if cohorts:
                sub = sub[sub["Cohort"].isin(cohorts)]
            stats = _stouffer_by_gene(ok, gene, sub)
            meta_blocks.append({"Gene": gene, "Meta_family": fam_id, "Meta_label": label, **stats})

    meta_df = pd.DataFrame(meta_blocks)
    meta_df.to_csv(OUT_DIR / "mpp1_marco_stouffer_meta.tsv", sep="\t", index=False)

    _plot_forest(ok, FIG_DIR / "mpp1_marco_multi_cohort_forest_ultra_hd.png")

    summary = {
        "genes": GENES,
        "discovery_GSE122897_log2FC": {
            g: float(discovery[g]) for g in GENES if g in discovery.index and pd.notna(discovery[g])
        },
        "n_comparisons_ok": int((ok.groupby(["Cohort", "Comparison"]).size() > 0).sum()),
        "meta": meta_df.to_dict(orient="records"),
    }
    (OUT_DIR / "mpp1_marco_survey_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    md = [
        "# MPP1 and MARCO — multi-cohort public expression survey",
        "",
        "Discovery reference: GSE122897 ruptured vs unruptured aneurysm wall (DESeq2 sign).",
        "All cohorts use Welch t-test on harmonized log2-scale expression.",
        "",
        "## Discovery direction (tissue)",
        "",
    ]
    for g in GENES:
        d = discovery.get(g, np.nan)
        md.append(f"- **{g}**: sign = {d:+.0f} (higher in ruptured in GSE122897)")
    md.extend(
        [
            "",
            "## Interpretation tiers",
            "",
            "1. **Tier 1 tissue (ruptured vs unruptured)** — direct rupture biology; suitable for convergence with discovery.",
            "2. **Tier 2 tissue context** — aneurysm presence or mixed designs; supportive only.",
            "3. **Tier 4 blood (aneurysm vs control)** — bridge to accessible biospecimens; not rupture status.",
            "4. **Tier 3 blood (post-rupture / prior SAH)** — systemic aftermath; do not treat as unruptured-risk validation.",
            "",
            "## Per-cohort results (ok only)",
            "",
            ok[
                [
                    "Gene",
                    "Cohort",
                    "Comparison",
                    "N_positive",
                    "N_negative",
                    "Mean_positive",
                    "Mean_negative",
                    "Log2FC_positive_minus_negative",
                    "P_value",
                    "Direction_match_discovery",
                ]
            ].to_markdown(index=False),
            "",
            "## Stouffer meta-analysis",
            "",
            meta_df.to_markdown(index=False),
        ]
    )
    (OUT_DIR / "mpp1_marco_survey_summary.md").write_text("\n".join(md), encoding="utf-8")

    print("\n=== Summary ===")
    print(ok.pivot_table(index=["Gene", "Cohort", "Comparison"], values="Log2FC_positive_minus_negative"))
    print("\nMeta:")
    print(meta_df.to_string(index=False))
    print(f"\nWrote {OUT_DIR}")


if __name__ == "__main__":
    main()
