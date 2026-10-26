"""
External blood-context analysis in GSE73378:
whole-blood expression profiles from patients with a past aneurysmal SAH history
versus controls.

This cohort is not a prospective rupture-risk dataset. We use it as an
orthogonal blood context that may contain longer-horizon post-aSAH or aneurysm-
related systemic biology.
"""

from __future__ import annotations

import gzip
import io
import json
from pathlib import Path

import numpy as np
import pandas as pd

from public_cohort_context_utils import add_gse122897_reference, clean_symbol, compute_gene_level_welch


BASE = Path(r"e:\Aneurysms Genomics Project")
DATA_DIR = BASE / "data" / "microarray" / "GSE73378"
SERIES_MATRIX = DATA_DIR / "GSE73378-GPL10558_series_matrix.txt.gz"
NORMALIZED = DATA_DIR / "GSE73378_normalized.txt.gz"
OUT_DIR = BASE / "results" / "external_validation_gse73378"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PLATFORM_CANDIDATES = [
    DATA_DIR / "GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt.gz",
    BASE / "data" / "microarray" / "GSE36791" / "raw_extracted" / "GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt.gz",
]


def _parse_series_matrix(path: Path) -> pd.DataFrame:
    lines = gzip.open(path, "rt", encoding="utf-8", errors="ignore").read().splitlines()
    begin = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_begin"))
    sample_lines = [line for line in lines[:begin] if line.startswith("!Sample_")]

    titles: list[str] = []
    sample_ids: list[str] = []
    sources: list[str] = []
    characteristics_rows: list[list[str]] = []

    for line in sample_lines:
        parts = [x.strip().strip('"') for x in line.split("\t")[1:]]
        if line.startswith("!Sample_title"):
            titles = parts
        elif line.startswith("!Sample_geo_accession"):
            sample_ids = parts
        elif line.startswith("!Sample_source_name_ch1"):
            sources = parts
        elif line.startswith("!Sample_characteristics_ch1"):
            characteristics_rows.append(parts)

    meta = pd.DataFrame({"Title": titles, "GSM": sample_ids, "Source": sources})
    meta["Sample"] = meta["Title"].str.split().str[-1].str.strip()

    for row in characteristics_rows:
        if not row:
            continue
        key = row[0].split(":", 1)[0].strip().lower().replace(" ", "_")
        values = []
        for value in row:
            if ":" in value:
                values.append(value.split(":", 1)[1].strip())
            else:
                values.append(value.strip())
        meta[key] = values

    meta["Group"] = meta["disease_status"].astype(str).str.title().replace({"Case": "Case", "Control": "Control"})
    meta["Sex"] = meta["sex"].astype(str).str.strip()
    meta["Tissue"] = meta["tissue"].astype(str).str.strip()
    return meta


def _load_platform_annotation(path: Path) -> pd.DataFrame:
    lines = gzip.open(path, "rt", encoding="utf-8", errors="ignore").read().splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip() == "[Probes]")
    annot = pd.read_csv(io.StringIO("\n".join(lines[start + 1 :])), sep="\t")
    annot = annot.rename(columns={"Array_Address_Id": "ProbeID"})
    annot["Gene"] = annot["Symbol"].map(clean_symbol)
    annot["Gene_title_clean"] = annot["Definition"].astype(str).str.strip()
    annot = annot[annot["Gene"].notna()].copy()
    annot = annot[["ProbeID", "Gene", "Gene_title_clean"]].drop_duplicates(subset=["ProbeID"])
    annot["ProbeID"] = annot["ProbeID"].astype(str).str.replace(r"\.0$", "", regex=True)
    return annot


def _gene_level_matrix(expr_path: Path, annot: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    with gzip.open(expr_path, "rt", encoding="utf-8", errors="ignore") as handle:
        probe_expr = pd.read_csv(handle, sep="\t")
    probe_expr = probe_expr.rename(columns={probe_expr.columns[0]: "ProbeID"})
    probe_expr["ProbeID"] = probe_expr["ProbeID"].astype(str).str.replace(r"\.0$", "", regex=True)
    numeric_cols = [c for c in probe_expr.columns if c != "ProbeID"]
    probe_expr[numeric_cols] = probe_expr[numeric_cols].apply(pd.to_numeric, errors="coerce")
    merged = probe_expr.merge(annot, on="ProbeID", how="inner")
    stats = {
        "n_probes_in_matrix": int(len(probe_expr)),
        "n_annotated_probes_used": int(merged["ProbeID"].nunique()),
        "n_unique_genes_after_annotation": int(merged["Gene"].nunique()),
    }
    gene_expr = merged.groupby("Gene", as_index=True)[numeric_cols].mean().sort_index()
    return gene_expr, stats


def main() -> None:
    for path in (SERIES_MATRIX, NORMALIZED):
        if not path.is_file():
            raise FileNotFoundError(path)
    platform_path = next((path for path in PLATFORM_CANDIDATES if path.is_file()), None)
    if platform_path is None:
        raise FileNotFoundError("No GPL10558 platform annotation file was found in expected locations.")

    meta = _parse_series_matrix(SERIES_MATRIX)
    annot = _load_platform_annotation(platform_path)
    gene_expr, parse_stats = _gene_level_matrix(NORMALIZED, annot)

    available_cols = [col for col in meta["Sample"].tolist() if col in gene_expr.columns]
    meta = meta[meta["Sample"].isin(available_cols)].copy()
    gene_expr = gene_expr[available_cols]

    meta.to_csv(OUT_DIR / "gse73378_sample_metadata.tsv", sep="\t", index=False)
    gene_expr.to_csv(OUT_DIR / "gse73378_gene_level_expression.tsv", sep="\t", index=True)

    case_cols = meta.loc[meta["Group"] == "Case", "Sample"].tolist()
    control_cols = meta.loc[meta["Group"] == "Control", "Sample"].tolist()
    res = compute_gene_level_welch(
        gene_expr,
        group_a_cols=case_cols,
        group_b_cols=control_cols,
        group_a_label="Case",
        group_b_label="Control",
        effect_name="Log2FC_CasevControl",
    )
    merged = add_gse122897_reference(res, base=BASE, effect_col="Log2FC_CasevControl")

    res.to_csv(OUT_DIR / "gse73378_gene_level_context.tsv", sep="\t", index=False)
    merged.to_csv(OUT_DIR / "gse73378_gene_level_context_with_gse122897.tsv", sep="\t", index=False)
    merged.sort_values("P_value", ascending=True).head(30).to_csv(
        OUT_DIR / "gse73378_top30_genes.tsv", sep="\t", index=False
    )

    shortlist_hits = merged[merged["In_PHAGESG_shortlist"] == True].copy()  # noqa: E712
    if not shortlist_hits.empty:
        shortlist_hits.sort_values("P_value", ascending=True).to_csv(
            OUT_DIR / "gse73378_shortlist_overlap.tsv", sep="\t", index=False
        )

    summary = {
        "dataset": "GSE73378",
        "comparison": "whole blood, past aneurysmal SAH history vs control",
        "platform": "Illumina GPL10558 normalized whole-blood array, probe-mapped to genes",
        "n_total_samples_loaded": int(len(meta)),
        "n_case": int(len(case_cols)),
        "n_control": int(len(control_cols)),
        "n_genes_tested": int(len(res)),
        "n_genes_p_lt_0_05": int((res["P_value"] < 0.05).sum()),
        "n_genes_fdr_lt_0_05": int((res["FDR_BH"] < 0.05).sum()),
        "n_overlap_with_gse122897_symbols": int(merged["GSE122897_log2FC_compare"].notna().sum()),
        "n_direction_match_with_gse122897": int(
            merged.loc[merged["GSE122897_log2FC_compare"].notna(), "DirectionMatch_GSE122897"].fillna(False).sum()
        ),
        "n_shortlist_genes_present": int(merged["In_PHAGESG_shortlist"].sum()),
        "n_shortlist_direction_match": int(
            merged.loc[merged["In_PHAGESG_shortlist"] == True, "DirectionMatch_GSE122897"].fillna(False).sum()  # noqa: E712
        ),
        "n_shortlist_p_lt_0_05": int(
            ((merged["In_PHAGESG_shortlist"] == True) & (merged["P_value"] < 0.05)).sum()  # noqa: E712
        ),
        "n_shortlist_fdr_lt_0_05": int(
            ((merged["In_PHAGESG_shortlist"] == True) & (merged["FDR_BH"] < 0.05)).sum()  # noqa: E712
        ),
        "orthogonal_context_note": (
            "This cohort captures blood expression in people with a previous aneurysmal SAH history and should be interpreted as "
            "long-horizon blood context rather than direct validation of prospective rupture biology."
        ),
        **parse_stats,
    }
    (OUT_DIR / "gse73378_context_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
