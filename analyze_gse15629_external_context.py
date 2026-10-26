"""
External tissue-context analysis for GSE15629:
ruptured aneurysm, unruptured aneurysm, and control vessel wall samples.

This cohort is useful in two ways:
1. a secondary direct rupture subset (ruptured vs unruptured), and
2. an orthogonal aneurysm-any vs control tissue context.
"""

from __future__ import annotations

import gzip
import io
import json
import time
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd

from public_cohort_context_utils import (
    add_gse122897_reference,
    compute_gene_level_welch,
    extract_symbol_from_assignment,
)


BASE = Path(r"e:\Aneurysms Genomics Project")
DATA_DIR = BASE / "data" / "microarray" / "GSE15629"
SERIES_MATRIX = DATA_DIR / "GSE15629_series_matrix.txt.gz"
PLATFORM_TEXT = DATA_DIR / "GPL6244_platform_full.txt"
OUT_DIR = BASE / "results" / "external_validation_gse15629"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SERIES_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE15nnn/GSE15629/matrix/GSE15629_series_matrix.txt.gz"
PLATFORM_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6244&targ=self&form=text&view=full"
REQUEST_HEADERS = {"User-Agent": "AneurysmContextGSE15629/1.0"}


def _download_if_missing(url: str, destination: Path, *, binary: bool = True) -> None:
    if destination.is_file() and destination.stat().st_size > 0:
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    last_error: Exception | None = None
    for attempt in range(3):
        try:
            request = Request(url, headers=REQUEST_HEADERS)
            with urlopen(request, timeout=300) as response:
                payload = response.read()
            break
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            if attempt == 2:
                raise
            time.sleep(2.0 * (attempt + 1))
    else:
        raise RuntimeError(f"Unable to download {url}: {last_error}")
    if binary:
        destination.write_bytes(payload)
    else:
        destination.write_text(payload.decode("utf-8", errors="ignore"), encoding="utf-8")


def _parse_series_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    lines = gzip.open(path, "rt", encoding="utf-8", errors="ignore").read().splitlines()
    begin = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_begin"))
    end = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_end"))

    sample_ids: list[str] = []
    titles: list[str] = []
    sources: list[str] = []
    group_labels: list[str] = []

    for line in lines[:begin]:
        if line.startswith("!Sample_geo_accession"):
            sample_ids = [x.strip().strip('"') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_title"):
            titles = [x.strip().strip('"') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_source_name_ch1"):
            sources = [x.strip().strip('"') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1") and "type:" in line.lower():
            group_labels = [x.strip().strip('"') for x in line.split("\t")[1:]]

    table = pd.read_csv(io.StringIO("\n".join(lines[begin + 1 : end])), sep="\t")
    table = table.rename(columns={"ID_REF": "ProbeID"})
    table["ProbeID"] = table["ProbeID"].astype(str)
    table = table.set_index("ProbeID")
    table.columns = [str(col).strip().strip('"') for col in table.columns]
    table = table.apply(pd.to_numeric, errors="coerce")

    meta = pd.DataFrame({"Sample": sample_ids, "Title": titles, "Source": sources, "Group_label": group_labels})
    meta["Group"] = meta["Group_label"].astype(str).str.replace("type:", "", regex=False).str.strip().str.title()
    return table, meta


def _load_platform_annotation(path: Path) -> pd.DataFrame:
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    begin = next(i for i, line in enumerate(lines) if line.strip() == "!platform_table_begin")
    end = next(i for i, line in enumerate(lines) if line.strip() == "!platform_table_end")
    annot = pd.read_csv(io.StringIO("\n".join(lines[begin + 1 : end])), sep="\t")
    annot = annot.rename(columns={"ID": "ProbeID"})
    annot["ProbeID"] = annot["ProbeID"].astype(str)
    annot["Gene"] = annot["gene_assignment"].map(extract_symbol_from_assignment)
    annot = annot[annot["Gene"].notna()].copy()
    annot = annot[["ProbeID", "Gene"]].drop_duplicates(subset=["ProbeID"])
    return annot


def _gene_level_matrix(probe_expr: pd.DataFrame, annot: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    transform_needed = float(np.nanpercentile(probe_expr.values, 95)) > 100.0
    expr_use = np.log2(probe_expr + 1.0) if transform_needed else probe_expr.copy()
    df = expr_use.reset_index().rename(columns={"index": "ProbeID"})
    merged = df.merge(annot, on="ProbeID", how="inner")
    stats = {
        "n_probes_in_matrix": int(probe_expr.shape[0]),
        "n_annotated_probes_used": int(merged["ProbeID"].nunique()),
        "n_unique_genes_after_annotation": int(merged["Gene"].nunique()),
        "log2_transform_applied": bool(transform_needed),
    }
    numeric_cols = [c for c in merged.columns if c not in {"ProbeID", "Gene"}]
    gene_expr = merged.groupby("Gene", as_index=True)[numeric_cols].mean().sort_index()
    return gene_expr, stats


def main() -> None:
    _download_if_missing(SERIES_URL, SERIES_MATRIX, binary=True)
    _download_if_missing(PLATFORM_URL, PLATFORM_TEXT, binary=False)

    probe_expr, meta = _parse_series_matrix(SERIES_MATRIX)
    annot = _load_platform_annotation(PLATFORM_TEXT)
    gene_expr, parse_stats = _gene_level_matrix(probe_expr, annot)

    meta.to_csv(OUT_DIR / "gse15629_sample_metadata.tsv", sep="\t", index=False)
    gene_expr.to_csv(OUT_DIR / "gse15629_gene_level_expression.tsv", sep="\t", index=True)

    ruptured_cols = meta.loc[meta["Group"] == "Ruptured", "Sample"].tolist()
    unruptured_cols = meta.loc[meta["Group"] == "Unruptured", "Sample"].tolist()
    control_cols = meta.loc[meta["Group"] == "Control", "Sample"].tolist()
    aneurysm_cols = meta.loc[meta["Group"].isin(["Ruptured", "Unruptured"]), "Sample"].tolist()

    rupture_res = compute_gene_level_welch(
        gene_expr,
        group_a_cols=ruptured_cols,
        group_b_cols=unruptured_cols,
        group_a_label="Ruptured",
        group_b_label="Unruptured",
        effect_name="Log2FC_RvU",
    )
    rupture_merged = add_gse122897_reference(rupture_res, base=BASE, effect_col="Log2FC_RvU")

    context_res = compute_gene_level_welch(
        gene_expr,
        group_a_cols=aneurysm_cols,
        group_b_cols=control_cols,
        group_a_label="AneurysmAny",
        group_b_label="Control",
        effect_name="Log2FC_AneurysmAnyvControl",
    )
    context_merged = add_gse122897_reference(context_res, base=BASE, effect_col="Log2FC_AneurysmAnyvControl")

    rupture_res.to_csv(OUT_DIR / "gse15629_rupture_subset.tsv", sep="\t", index=False)
    rupture_merged.to_csv(OUT_DIR / "gse15629_rupture_subset_with_gse122897.tsv", sep="\t", index=False)
    context_res.to_csv(OUT_DIR / "gse15629_aneurysm_vs_control.tsv", sep="\t", index=False)
    context_merged.to_csv(OUT_DIR / "gse15629_aneurysm_vs_control_with_gse122897.tsv", sep="\t", index=False)

    rupture_shortlist = rupture_merged[rupture_merged["In_PHAGESG_shortlist"] == True].copy()  # noqa: E712
    if not rupture_shortlist.empty:
        rupture_shortlist.sort_values("P_value", ascending=True).to_csv(
            OUT_DIR / "gse15629_shortlist_overlap_rupture.tsv", sep="\t", index=False
        )

    context_shortlist = context_merged[context_merged["In_PHAGESG_shortlist"] == True].copy()  # noqa: E712
    if not context_shortlist.empty:
        context_shortlist.sort_values("P_value", ascending=True).to_csv(
            OUT_DIR / "gse15629_shortlist_overlap_aneurysm_vs_control.tsv", sep="\t", index=False
        )

    summary = {
        "dataset": "GSE15629",
        "platform": "Affymetrix GPL6244 series-matrix values mapped to transcript-cluster genes",
        "n_total_samples_loaded": int(len(meta)),
        "n_ruptured": int(len(ruptured_cols)),
        "n_unruptured": int(len(unruptured_cols)),
        "n_control": int(len(control_cols)),
        "n_genes_tested_rupture_subset": int(len(rupture_res)),
        "n_genes_p_lt_0_05_rupture_subset": int((rupture_res["P_value"] < 0.05).sum()),
        "n_genes_fdr_lt_0_05_rupture_subset": int((rupture_res["FDR_BH"] < 0.05).sum()),
        "n_overlap_with_gse122897_symbols_rupture_subset": int(rupture_merged["GSE122897_log2FC_compare"].notna().sum()),
        "n_direction_match_with_gse122897_rupture_subset": int(
            rupture_merged.loc[
                rupture_merged["GSE122897_log2FC_compare"].notna(), "DirectionMatch_GSE122897"
            ].fillna(False).sum()
        ),
        "n_shortlist_genes_present_rupture_subset": int(rupture_merged["In_PHAGESG_shortlist"].sum()),
        "n_shortlist_direction_match_rupture_subset": int(
            rupture_merged.loc[rupture_merged["In_PHAGESG_shortlist"] == True, "DirectionMatch_GSE122897"].fillna(False).sum()  # noqa: E712
        ),
        "n_genes_tested_aneurysm_vs_control": int(len(context_res)),
        "n_genes_p_lt_0_05_aneurysm_vs_control": int((context_res["P_value"] < 0.05).sum()),
        "n_genes_fdr_lt_0_05_aneurysm_vs_control": int((context_res["FDR_BH"] < 0.05).sum()),
        "n_shortlist_genes_present_aneurysm_vs_control": int(context_merged["In_PHAGESG_shortlist"].sum()),
        "n_shortlist_direction_match_aneurysm_vs_control": int(
            context_merged.loc[context_merged["In_PHAGESG_shortlist"] == True, "DirectionMatch_GSE122897"].fillna(False).sum()  # noqa: E712
        ),
        "context_note": (
            "This cohort supports both a secondary mixed-design rupture subset and an orthogonal aneurysm-vessel context; "
            "the latter should not be presented as direct rupture replication."
        ),
        **parse_stats,
    }
    (OUT_DIR / "gse15629_context_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
