"""
External tissue-context analysis in GSE26969 (Affymetrix GPL570):
unruptured intracranial aneurysm tissue versus normal superficial temporal artery.

The dataset is not a rupture validation cohort. We use it as orthogonal tissue
context to ask whether genes implicated in GSE122897 rupture biology are already
dysregulated in aneurysm wall relative to non-aneurysmal artery.
"""

from __future__ import annotations

import gzip
import io
import json
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

BASE = Path(r"e:\Aneurysms Genomics Project")
SERIES_MATRIX = BASE / "data" / "microarray" / "GSE26969" / "GSE26969_series_matrix.txt.gz"
GPL570_ANNOT = BASE / "data" / "microarray" / "GSE26969" / "GPL570.annot.gz"
OUT_DIR = BASE / "results" / "external_validation_gse26969"
OUT_DIR.mkdir(parents=True, exist_ok=True)

GPL570_URL = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"


def _clean_symbol(raw: object) -> str | None:
    if pd.isna(raw):
        return None
    text = str(raw).strip()
    if not text or text in {"---", "nan", "NA", "None"}:
        return None
    parts = [part.strip() for part in text.split("///") if part.strip() and part.strip() != "---"]
    unique: list[str] = []
    for part in parts:
        if part not in unique:
            unique.append(part)
    if len(unique) != 1:
        return None
    return unique[0]


def _download_if_missing(url: str, dest: Path) -> None:
    if dest.is_file() and dest.stat().st_size > 0:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    resp = requests.get(url, timeout=300)
    resp.raise_for_status()
    dest.write_bytes(resp.content)


def _parse_series_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    lines = gzip.open(path, "rt", encoding="utf-8", errors="ignore").read().splitlines()
    begin = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_begin"))
    end = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_end"))

    sample_ids: list[str] = []
    titles: list[str] = []
    sources: list[str] = []
    statuses: list[str] = []

    for line in lines[:begin]:
        if line.startswith("!Sample_geo_accession"):
            sample_ids = [x.strip().strip('"') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_title"):
            titles = [x.strip().strip('"') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_source_name_ch1"):
            sources = [x.strip().strip('"') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1") and "status:" in line.lower():
            statuses = [x.strip().strip('"') for x in line.split("\t")[1:]]

    table = pd.read_csv(io.StringIO("\n".join(lines[begin + 1 : end])), sep="\t")
    table = table.rename(columns={"ID_REF": "ProbeID"}).set_index("ProbeID")
    table.columns = [str(col).strip().strip('"') for col in table.columns]
    table = table.apply(pd.to_numeric, errors="coerce")

    meta = pd.DataFrame(
        {
            "Sample": sample_ids,
            "Title": titles,
            "Source": sources,
            "Status": statuses,
        }
    )

    def infer_group(row: pd.Series) -> str:
        text = " ".join([str(row.get("Title", "")), str(row.get("Source", "")), str(row.get("Status", ""))]).lower()
        if "aneurysmal tissue" in text or "unruptured aneurysm" in text:
            return "Aneurysm"
        return "Control"

    meta["Group"] = meta.apply(infer_group, axis=1)
    return table, meta


def _load_gpl570_annotation(path: Path) -> pd.DataFrame:
    _download_if_missing(GPL570_URL, path)
    lines = gzip.open(path, "rt", encoding="utf-8", errors="ignore").read().splitlines()
    begin = next(i for i, line in enumerate(lines) if line.startswith("!platform_table_begin"))
    end = next(i for i, line in enumerate(lines) if line.startswith("!platform_table_end"))
    annot = pd.read_csv(io.StringIO("\n".join(lines[begin + 1 : end])), sep="\t")
    annot = annot.rename(columns={"ID": "ProbeID"})
    annot["Gene"] = annot["Gene symbol"].map(_clean_symbol)
    annot["Gene_title_clean"] = annot["Gene title"].astype(str).str.split("///").str[0].str.strip()
    annot = annot[annot["Gene"].notna()].copy()
    annot = annot[["ProbeID", "Gene", "Gene_title_clean"]].drop_duplicates(subset=["ProbeID"])
    return annot


def _gene_level_matrix(probe_expr: pd.DataFrame, annot: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    df = probe_expr.reset_index().rename(columns={"index": "ProbeID"})
    merged = df.merge(annot, on="ProbeID", how="inner")
    stats = {
        "n_probes_in_matrix": int(probe_expr.shape[0]),
        "n_annotated_probes_used": int(merged["ProbeID"].nunique()),
        "n_unique_genes_after_annotation": int(merged["Gene"].nunique()),
    }
    numeric_cols = [c for c in merged.columns if c not in {"ProbeID", "Gene", "Gene_title_clean"}]
    gene_expr = merged.groupby("Gene", as_index=True)[numeric_cols].mean().sort_index()
    return gene_expr, stats


def _load_gse122897_deseq2() -> pd.DataFrame:
    de = pd.read_csv(BASE / "results" / "gse122897_deseq2_results.tsv", sep="\t")
    mapping = pd.read_csv(BASE / "results" / "gene_id_mapping_unified.tsv", sep="\t")
    mapping = mapping[["Original_ID", "Gene_Symbol"]].dropna(subset=["Gene_Symbol"]).copy()
    mapping["Gene_ID"] = mapping["Original_ID"].astype(str).str.replace(r"\.0$", "", regex=True)
    de["Gene_ID"] = de["Gene_ID"].astype(str).str.replace(r"\.0$", "", regex=True)
    de = de.merge(mapping[["Gene_ID", "Gene_Symbol"]].drop_duplicates("Gene_ID"), on="Gene_ID", how="left")
    de["Gene"] = de["Gene_Symbol"].where(de["Gene_Symbol"].notna(), de["Gene_ID"])
    return de[["Gene", "log2FoldChange", "padj", "baseMean"]].rename(
        columns={
            "log2FoldChange": "GSE122897_log2FC",
            "padj": "GSE122897_padj",
            "baseMean": "GSE122897_baseMean",
        }
    )


def _load_shortlist_effects() -> pd.DataFrame:
    shortlist_path = BASE / "results" / "phases_g" / "shortlist_evidence.tsv"
    if not shortlist_path.is_file():
        return pd.DataFrame(columns=["Gene", "Shortlist_GSE122897_log2FC", "Shortlist_GSE122897_padj"])
    sl = pd.read_csv(shortlist_path, sep="\t")
    if "symbol" not in sl.columns:
        return pd.DataFrame(columns=["Gene", "Shortlist_GSE122897_log2FC", "Shortlist_GSE122897_padj"])
    keep = ["symbol"]
    if "DESeq2_log2FC" in sl.columns:
        keep.append("DESeq2_log2FC")
    if "DESeq2_padj" in sl.columns:
        keep.append("DESeq2_padj")
    sl = sl[keep].dropna(subset=["symbol"]).drop_duplicates(subset=["symbol"]).rename(columns={"symbol": "Gene"})
    rename = {}
    if "DESeq2_log2FC" in sl.columns:
        rename["DESeq2_log2FC"] = "Shortlist_GSE122897_log2FC"
    if "DESeq2_padj" in sl.columns:
        rename["DESeq2_padj"] = "Shortlist_GSE122897_padj"
    return sl.rename(columns=rename)


def main() -> None:
    if not SERIES_MATRIX.is_file():
        raise FileNotFoundError(SERIES_MATRIX)

    probe_expr, meta = _parse_series_matrix(SERIES_MATRIX)
    annot = _load_gpl570_annotation(GPL570_ANNOT)
    gene_expr, parse_stats = _gene_level_matrix(probe_expr, annot)

    meta.to_csv(OUT_DIR / "gse26969_sample_metadata.tsv", sep="\t", index=False)
    gene_expr.to_csv(OUT_DIR / "gse26969_gene_level_expression.tsv", sep="\t", index=True)

    aneurysm_cols = meta.loc[meta["Group"] == "Aneurysm", "Sample"].tolist()
    control_cols = meta.loc[meta["Group"] == "Control", "Sample"].tolist()

    rows = []
    for gene, row in gene_expr.iterrows():
        av = pd.to_numeric(row[aneurysm_cols], errors="coerce").dropna().values.astype(float)
        cv = pd.to_numeric(row[control_cols], errors="coerce").dropna().values.astype(float)
        if len(av) < 2 or len(cv) < 2:
            continue
        t_stat, p_val = ttest_ind(av, cv, equal_var=False, nan_policy="omit")
        log2fc = float(np.nanmean(av) - np.nanmean(cv))
        rows.append(
            {
                "Gene": gene,
                "Mean_Aneurysm": float(np.nanmean(av)),
                "Mean_Control": float(np.nanmean(cv)),
                "Log2FC_AvC": log2fc,
                "T_statistic": float(t_stat) if pd.notna(t_stat) else np.nan,
                "P_value": float(p_val) if pd.notna(p_val) else np.nan,
            }
        )

    res = pd.DataFrame(rows).sort_values("P_value", ascending=True).reset_index(drop=True)
    res["FDR_BH"] = multipletests(res["P_value"].fillna(1.0).values, method="fdr_bh")[1]

    de122 = _load_gse122897_deseq2()
    sl_effects = _load_shortlist_effects()
    merged = res.merge(de122, on="Gene", how="left")
    merged = merged.merge(sl_effects, on="Gene", how="left")
    merged["GSE122897_log2FC_compare"] = merged["Shortlist_GSE122897_log2FC"].where(
        merged["Shortlist_GSE122897_log2FC"].notna(), merged["GSE122897_log2FC"]
    )
    merged["GSE122897_padj_compare"] = merged["Shortlist_GSE122897_padj"].where(
        merged["Shortlist_GSE122897_padj"].notna(), merged["GSE122897_padj"]
    )
    merged["DirectionMatch_GSE122897"] = (
        np.sign(merged["Log2FC_AvC"]).replace(0, np.nan)
        == np.sign(merged["GSE122897_log2FC_compare"]).replace(0, np.nan)
    )
    merged["In_PHAGESG_shortlist"] = merged["Shortlist_GSE122897_log2FC"].notna()

    res.to_csv(OUT_DIR / "gse26969_gene_level_context.tsv", sep="\t", index=False)
    merged.to_csv(OUT_DIR / "gse26969_gene_level_context_with_gse122897.tsv", sep="\t", index=False)
    merged.sort_values("P_value", ascending=True).head(30).to_csv(
        OUT_DIR / "gse26969_top30_genes.tsv", sep="\t", index=False
    )

    shortlist_hits = merged[merged["In_PHAGESG_shortlist"] == True].copy()  # noqa: E712
    if not shortlist_hits.empty:
        shortlist_hits = shortlist_hits.sort_values("P_value", ascending=True)
        shortlist_hits.to_csv(OUT_DIR / "gse26969_shortlist_overlap.tsv", sep="\t", index=False)

    summary = {
        "dataset": "GSE26969",
        "comparison": "unruptured aneurysm tissue vs superficial temporal artery control",
        "platform": "Affymetrix GPL570 series-matrix RMA values",
        "n_total_samples_loaded": int(len(meta)),
        "n_aneurysm": int(len(aneurysm_cols)),
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
            "This cohort benchmarks aneurysm wall against non-aneurysmal artery and is not a rupture validation dataset."
        ),
        **parse_stats,
    }
    (OUT_DIR / "gse26969_context_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print("\nTop 15 genes:")
    print(
        merged[
            ["Gene", "Log2FC_AvC", "P_value", "FDR_BH", "GSE122897_log2FC_compare", "DirectionMatch_GSE122897"]
        ]
        .head(15)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
