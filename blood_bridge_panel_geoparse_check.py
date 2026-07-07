"""
Blood panel gene-level verification (GEO supplementary + local fallbacks).

Downloads expression data programmatically when missing, then runs:
  1. Within-cohort Welch DE (aneurysm vs control): log2FC, t, P
  2. Single-gene AUROC with fixed discovery tissue sign (GSE122897 DESeq2)
  3. Leave-one-out AUROC drop using the same discovery-signed panel score

GSE106520 labels come from results/external_validation_gse106520/gse106520_sample_metadata.tsv
(not keyword parsing on GEO soft metadata).

GSE159610 uses supplementary raw gene counts; GSE106520 uses per-sample FPKM tables.
"""

from __future__ import annotations

import gzip
import json
import re
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urljoin
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import norm

from blood_bridge_panel import (
    BASE,
    _auc_from_scores,
    _load_gse159610_expression,
)
from blood_bridge_panel_gene_attribution import (
    PANEL_GENES,
    _gene_score,
    _load_panel_signs,
    _panel_score,
    _prepare_cohort,
)
from public_cohort_context_utils import compute_gene_level_welch

OUT_DIR = BASE / "results" / "blood_bridge_panel"
REQUEST_HEADERS = {"User-Agent": "BloodPanelVerification/1.0"}

GSE159610_DIR = BASE / "data" / "blood_rnaseq" / "GSE159610"
GSE106520_DIR = BASE / "data" / "blood_rnaseq" / "GSE106520"
GSE159610_COUNTS = GSE159610_DIR / "GSE159610_raw_gene_counts_matrix.txt.gz"
GSE106520_MATRIX = GSE106520_DIR / "GSE106520_series_matrix.txt.gz"
GSE106520_META = BASE / "results" / "external_validation_gse106520" / "gse106520_sample_metadata.tsv"
GSE106520_EXPR = BASE / "results" / "external_validation_gse106520" / "gse106520_gene_level_log2fpkm.tsv"

COHORT_CONFIG = {
    "GSE159610": {
        "description": "Whole blood RNA-seq: aneurysm vs control",
        "positive": "Aneurysm",
        "negative": "Control",
    },
    "GSE106520": {
        "description": "Neutrophil RNA-seq: aneurysm vs control (project metadata)",
        "positive": "Aneurysm",
        "negative": "Control",
    },
}


class LinkCollector(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag != "a":
            return
        for key, value in attrs:
            if key == "href" and value:
                self.links.append(value)


def _series_prefix(accession: str) -> str:
    match = re.match(r"^(GSE)(\d+)$", accession)
    if not match:
        raise ValueError(accession)
    return f"GSE{match.group(2)[:-3]}nnn"


def _download(url: str, destination: Path) -> None:
    if destination.is_file() and destination.stat().st_size > 0:
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = Request(url, headers=REQUEST_HEADERS)
    with urlopen(request, timeout=300) as response, destination.open("wb") as handle:
        handle.write(response.read())


def _list_ftp_files(url: str) -> list[str]:
    request = Request(url, headers=REQUEST_HEADERS)
    with urlopen(request, timeout=60) as response:
        html = response.read().decode("utf-8", errors="ignore")
    parser = LinkCollector()
    parser.feed(html)
    out = []
    for link in parser.links:
        if link in {"../", "./", ""}:
            continue
        out.append(urljoin(url, link))
    return sorted(set(out))


def ensure_gse159610_counts() -> Path:
    """Download raw gene count matrix from GEO supplementary if missing."""
    if GSE159610_COUNTS.is_file():
        return GSE159610_COUNTS

    prefix = _series_prefix("GSE159610")
    suppl_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/GSE159610/suppl/"
    urls = _list_ftp_files(suppl_url)
    targets = [u for u in urls if "raw_gene_counts_matrix" in u.lower() and u.endswith(".gz")]
    if not targets:
        raise FileNotFoundError(f"No raw count matrix found at {suppl_url}")
    _download(targets[0], GSE159610_COUNTS)
    return GSE159610_COUNTS


def _parse_gse106520_series_matrix(path: Path) -> pd.DataFrame:
    lines = gzip.open(path, "rt", encoding="utf-8", errors="ignore").read().splitlines()
    begin = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_begin"))
    sample_lines = [line for line in lines[:begin] if line.startswith("!Sample_")]

    sample_ids: list[str] = []
    titles: list[str] = []
    supp_files: list[str] = []

    for line in sample_lines:
        parts = [x.strip().strip('"') for x in line.split("\t")[1:]]
        if line.startswith("!Sample_geo_accession"):
            sample_ids = parts
        elif line.startswith("!Sample_title"):
            titles = parts
        elif line.startswith("!Sample_supplementary_file_1"):
            supp_files = parts

    return pd.DataFrame({"Sample": sample_ids, "Title": titles, "Supplementary_file": supp_files})


def ensure_gse106520_expression() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load or rebuild GSE106520 expression from supplementary per-sample FPKM files."""
    if GSE106520_EXPR.is_file() and GSE106520_META.is_file():
        expr = pd.read_csv(GSE106520_EXPR, sep="\t", index_col=0)
        meta = pd.read_csv(GSE106520_META, sep="\t")
        meta = meta[meta["Group"].isin(["Aneurysm", "Control"])].copy()
        return expr, meta

    matrix_url = (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE106nnn/GSE106520/matrix/GSE106520_series_matrix.txt.gz"
    )
    _download(matrix_url, GSE106520_MATRIX)

    series_meta = _parse_gse106520_series_matrix(GSE106520_MATRIX)
    matrices: list[pd.DataFrame] = []
    rows_meta: list[dict[str, str]] = []

    for _, row in series_meta.iterrows():
        title = str(row["Title"])
        if not (title.lower().startswith("aneurysm_") or title.lower().startswith("control_")):
            continue
        group = "Aneurysm" if title.lower().startswith("aneurysm_") else "Control"
        url = str(row["Supplementary_file"]).replace("ftp://", "https://")
        filename = url.rsplit("/", 1)[-1]
        local_path = GSE106520_DIR / filename
        _download(url, local_path)

        with gzip.open(local_path, "rt", encoding="utf-8", errors="ignore") as handle:
            df = pd.read_csv(handle, sep="\t")
        value_col = df.columns[-1]
        sample_col = str(row["Sample"])
        df = df[["gene_short_name", value_col]].rename(columns={"gene_short_name": "Gene", sample_col: value_col})
        df["Gene"] = df["Gene"].astype(str).str.strip()
        df = df[(df["Gene"] != "") & (df["Gene"] != "-")].copy()
        df[sample_col] = pd.to_numeric(df[sample_col], errors="coerce")
        df = df.groupby("Gene", as_index=False)[sample_col].mean()
        matrices.append(df)
        rows_meta.append({"Sample": sample_col, "Title": title, "Group": group, "Local_sample_file": filename})

    merged = matrices[0]
    for df in matrices[1:]:
        merged = merged.merge(df, on="Gene", how="outer")
    expr = np.log2(merged.set_index("Gene").sort_index() + 1.0)
    meta = pd.DataFrame(rows_meta)
    GSE106520_META.parent.mkdir(parents=True, exist_ok=True)
    meta.to_csv(GSE106520_META, sep="\t", index=False)
    expr.to_csv(GSE106520_EXPR, sep="\t", index=True)
    return expr, meta


def load_cohort(cohort_id: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    cfg = COHORT_CONFIG[cohort_id]
    if cohort_id == "GSE159610":
        ensure_gse159610_counts()
        return _load_gse159610_expression(BASE)
    expr, meta = ensure_gse106520_expression()
    return expr, meta


def within_cohort_de(
    expr: pd.DataFrame,
    labels: pd.Series,
    genes: list[str],
    cohort_id: str,
) -> pd.DataFrame:
    pos_cols = labels[labels].index.tolist()
    neg_cols = labels[~labels].index.tolist()
    sub = expr.loc[[g for g in genes if g in expr.index], pos_cols + neg_cols]
    res = compute_gene_level_welch(
        sub,
        group_a_cols=pos_cols,
        group_b_cols=neg_cols,
        group_a_label=COHORT_CONFIG[cohort_id]["positive"],
        group_b_label=COHORT_CONFIG[cohort_id]["negative"],
        effect_name="Log2FC_aneurysm_vs_control",
    )
    return res.rename(columns={"T_statistic": "T_statistic", "P_value": "P_value"})[
        ["Gene", "Log2FC_aneurysm_vs_control", "T_statistic", "P_value"]
    ]


def bridge_stouffer(p159: float, p106: float) -> tuple[float, float]:
    p1 = float(np.clip(p159, 1e-300, 1.0))
    p2 = float(np.clip(p106, 1e-300, 1.0))
    z = (norm.isf(p1 / 2.0) + norm.isf(p2 / 2.0)) / np.sqrt(2.0)
    meta_p = float(2.0 * norm.sf(abs(z)))
    return z, meta_p


def analyze_cohort(
    cohort_id: str,
    expr: pd.DataFrame,
    meta: pd.DataFrame,
    genes: list[str],
    signs: pd.Series,
) -> dict[str, pd.DataFrame]:
    cfg = COHORT_CONFIG[cohort_id]
    expr, labels = _prepare_cohort(expr, meta, cfg["positive"], cfg["negative"])

    de = within_cohort_de(expr, labels, genes, cohort_id)
    de["Cohort"] = cohort_id

    auc_rows = []
    for gene in genes:
        if gene not in expr.index:
            continue
        scores = _gene_score(expr, gene, float(signs[gene]))
        auc = _auc_from_scores(labels.to_numpy(dtype=bool), scores.to_numpy(dtype=float))
        pos = scores[labels]
        neg = scores[~labels]
        mw_p = float(stats.mannwhitneyu(pos, neg, alternative="two-sided").pvalue)
        auc_rows.append(
            {
                "Gene": gene,
                "Cohort": cohort_id,
                "Single_gene_AUROC": auc,
                "Mann_Whitney_P": mw_p,
                "Discovery_sign": float(signs[gene]),
            }
        )
    auc_df = pd.DataFrame(auc_rows)

    full_scores = _panel_score(expr, genes, signs)
    full_auc = _auc_from_scores(labels.to_numpy(dtype=bool), full_scores.to_numpy(dtype=float))

    loo_rows = []
    for gene in genes:
        if gene not in expr.index:
            continue
        reduced = [g for g in genes if g != gene]
        loo_scores = _panel_score(expr, reduced, signs)
        loo_auc = _auc_from_scores(labels.to_numpy(dtype=bool), loo_scores.to_numpy(dtype=float))
        loo_rows.append(
            {
                "Gene": gene,
                "Cohort": cohort_id,
                "Full_panel_AUROC": full_auc,
                "Leave_one_out_AUROC": loo_auc,
                "AUROC_drop": full_auc - loo_auc,
            }
        )
    loo_df = pd.DataFrame(loo_rows)

    return {"de": de, "auroc": auc_df, "loo": loo_df, "full_panel_auroc": full_auc}


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    signs = _load_panel_signs(BASE)
    genes = [g for g in PANEL_GENES if g in signs.index]

    de_parts: list[pd.DataFrame] = []
    auc_parts: list[pd.DataFrame] = []
    loo_parts: list[pd.DataFrame] = []
    summary: dict[str, object] = {}

    for cohort_id in COHORT_CONFIG:
        print(f"\n{'=' * 60}\n{cohort_id}: {COHORT_CONFIG[cohort_id]['description']}\n{'=' * 60}")
        expr, meta = load_cohort(cohort_id)
        result = analyze_cohort(cohort_id, expr, meta, genes, signs)
        de_parts.append(result["de"])
        auc_parts.append(result["auroc"])
        loo_parts.append(result["loo"])
        summary[cohort_id] = {"full_panel_auroc": result["full_panel_auroc"]}

        print(result["de"].to_string(index=False))
        print(f"\nFull panel AUROC: {result['full_panel_auroc']:.3f}")
        print("\nMPP1 / MARCO:")
        focus = result["auroc"].merge(result["loo"], on=["Gene", "Cohort"], how="left")
        print(focus[focus["Gene"].isin(["MPP1", "MARCO"])].to_string(index=False))

    de_all = pd.concat(de_parts, ignore_index=True)
    auc_all = pd.concat(auc_parts, ignore_index=True)
    loo_all = pd.concat(loo_parts, ignore_index=True)

    p159 = de_all[de_all["Cohort"] == "GSE159610"].set_index("Gene")["P_value"]
    p106 = de_all[de_all["Cohort"] == "GSE106520"].set_index("Gene")["P_value"]
    stouffer_rows = []
    for gene in genes:
        if gene in p159.index and gene in p106.index:
            z, mp = bridge_stouffer(float(p159[gene]), float(p106[gene]))
            stouffer_rows.append({"Gene": gene, "Bridge_Stouffer_Z": z, "Bridge_meta_P": mp})
    stouffer_df = pd.DataFrame(stouffer_rows)

    de_path = OUT_DIR / "blood_bridge_panel_verification_de.tsv"
    auc_path = OUT_DIR / "blood_bridge_panel_verification_auroc.tsv"
    loo_path = OUT_DIR / "blood_bridge_panel_verification_loo.tsv"
    de_all.to_csv(de_path, sep="\t", index=False)
    auc_all.to_csv(auc_path, sep="\t", index=False)
    loo_all.to_csv(loo_path, sep="\t", index=False)
    stouffer_df.to_csv(OUT_DIR / "blood_bridge_panel_verification_stouffer.tsv", sep="\t", index=False)

    # Compare to prior attribution outputs
    prior_de = OUT_DIR / "blood_bridge_panel_within_cohort_de.tsv"
    prior_auc = OUT_DIR / "blood_bridge_panel_single_gene_auroc.tsv"
    prior_loo = OUT_DIR / "blood_bridge_panel_leave_one_out_auroc.tsv"
    comparison = {"de_max_abs_diff": None, "auroc_max_abs_diff": None, "loo_max_abs_diff": None}
    if prior_de.is_file():
        old = pd.read_csv(prior_de, sep="\t")
        merged = de_all.merge(
            old,
            on=["Gene", "Cohort"],
            suffixes=("_new", "_prior"),
        )
        comparison["de_max_abs_diff"] = float(
            (merged["Log2FC_aneurysm_vs_control_new"] - merged["Log2FC_aneurysm_vs_control_prior"]).abs().max()
        )
    if prior_auc.is_file():
        old = pd.read_csv(prior_auc, sep="\t")
        merged = auc_all.merge(old, on=["Gene", "Cohort"], suffixes=("_new", "_prior"))
        comparison["auroc_max_abs_diff"] = float(
            (merged["Single_gene_AUROC_new"] - merged["Single_gene_AUROC_prior"]).abs().max()
        )
    if prior_loo.is_file():
        old = pd.read_csv(prior_loo, sep="\t")
        merged = loo_all.merge(old, on=["Gene", "Cohort"], suffixes=("_new", "_prior"))
        comparison["loo_max_abs_diff"] = float((merged["AUROC_drop_new"] - merged["AUROC_drop_prior"]).abs().max())

    (OUT_DIR / "blood_bridge_panel_verification_summary.json").write_text(
        json.dumps({"cohorts": summary, "comparison_to_prior": comparison}, indent=2),
        encoding="utf-8",
    )

    md_lines = [
        "# Blood bridge panel verification (discovery-signed)",
        "",
        "Data: programmatic GEO supplementary download with local fallbacks.",
        "GSE106520 labels: `gse106520_sample_metadata.tsv` (explicit Aneurysm/Control).",
        "AUROC: discovery tissue sign from GSE122897; **no** `max(AUC, 1-AUC)` flipping.",
        "LOO: discovery-signed panel score (same as `blood_bridge_panel.py`).",
        "",
        f"Full panel AUROC — GSE159610: **{summary['GSE159610']['full_panel_auroc']:.3f}**, "
        f"GSE106520: **{summary['GSE106520']['full_panel_auroc']:.3f}**.",
        "",
        "## Comparison to prior attribution outputs",
        "",
        f"- Max |Δlog2FC|: {comparison.get('de_max_abs_diff')}",
        f"- Max |ΔAUROC|: {comparison.get('auroc_max_abs_diff')}",
        f"- Max |ΔLOO drop|: {comparison.get('loo_max_abs_diff')}",
        "",
        "## MPP1 and MARCO",
        "",
    ]
    for gene in ["MPP1", "MARCO"]:
        md_lines.append(f"### {gene}")
        md_lines.append(de_all[de_all["Gene"] == gene].to_markdown(index=False))
        md_lines.append(auc_all[auc_all["Gene"] == gene].to_markdown(index=False))
        md_lines.append(loo_all[loo_all["Gene"] == gene].to_markdown(index=False))
        md_lines.append("")

    (OUT_DIR / "blood_bridge_panel_verification_summary.md").write_text("\n".join(md_lines), encoding="utf-8")

    print(f"\nWrote:\n  {de_path}\n  {auc_path}\n  {loo_path}")
    print(f"Comparison to prior pipeline: {comparison}")


if __name__ == "__main__":
    main()
