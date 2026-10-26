"""
External validation in GSE54083 (Agilent microarray):
ruptured intracranial aneurysm dome (RIA) vs unruptured dome (UIA).

This script parses the raw sample text files archived in GSE54083_RAW.tar,
builds a gene-level matrix from non-control probes, performs Welch t-tests on
log2(processed_signal + 1), adjusts p-values by Benjamini-Hochberg, and
compares the resulting direction of effect with GSE122897 DESeq2.
"""

from __future__ import annotations

import gzip
import io
import json
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

BASE = Path(r"e:\Aneurysms Genomics Project")
RAW_TAR = BASE / "data" / "microarray" / "GSE54083" / "GSE54083_RAW.tar"
OUT_DIR = BASE / "results" / "external_validation_gse54083"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def _sample_group(name: str) -> str | None:
    tag = name.upper()
    if "_RIA" in tag:
        return "Ruptured"
    if "_UIA" in tag:
        return "Unruptured"
    if "_STA" in tag:
        return "Control"
    return None


def _read_feature_block(raw_bytes: bytes) -> pd.DataFrame:
    lines = gzip.decompress(raw_bytes).decode("utf-8", errors="ignore").splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("FEATURES\t"))
    df = pd.read_csv(io.StringIO("\n".join(lines[start:])), sep="\t")
    return df


def _load_gse54083_matrix() -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_series: dict[str, pd.Series] = {}
    sample_meta: list[dict[str, str]] = []

    with tarfile.open(RAW_TAR, "r") as tar:
        for member in tar.getmembers():
            if not member.name.lower().endswith(".txt.gz"):
                continue
            group = _sample_group(member.name)
            if group is None:
                continue

            handle = tar.extractfile(member)
            if handle is None:
                continue
            df = _read_feature_block(handle.read())

            # Non-control probes only; aggregate probes to gene-level mean on log scale.
            df = df[df["ControlType"] == 0].copy()
            df["GeneName"] = df["GeneName"].astype(str).str.strip()
            df = df[
                df["GeneName"].notna()
                & (df["GeneName"] != "")
                & (~df["GeneName"].str.startswith("DarkCorner", na=False))
                & (~df["GeneName"].str.startswith("GE_", na=False))
            ].copy()
            df["gProcessedSignal"] = pd.to_numeric(df["gProcessedSignal"], errors="coerce")
            df = df[df["gProcessedSignal"].notna()].copy()
            df["log2_signal"] = np.log2(df["gProcessedSignal"].clip(lower=0) + 1.0)

            sample_id = Path(member.name).stem.replace(".txt", "")
            gene_signal = df.groupby("GeneName", as_index=True)["log2_signal"].mean()
            sample_series[sample_id] = gene_signal
            sample_meta.append({"Sample": sample_id, "Group": group, "ArchiveMember": member.name})

    expr = pd.DataFrame(sample_series).sort_index()
    meta = pd.DataFrame(sample_meta).sort_values(["Group", "Sample"]).reset_index(drop=True)
    return expr, meta


def _load_gse122897_deseq2() -> pd.DataFrame:
    de = pd.read_csv(BASE / "results" / "gse122897_deseq2_results.tsv", sep="\t")
    mapping = pd.read_csv(BASE / "results" / "gene_id_mapping_unified.tsv", sep="\t")
    if "Original_ID" in mapping.columns and "Gene_Symbol" in mapping.columns:
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
    if not RAW_TAR.is_file():
        raise FileNotFoundError(RAW_TAR)

    expr, meta = _load_gse54083_matrix()
    expr.to_csv(OUT_DIR / "gse54083_gene_level_log2_expression.tsv", sep="\t", index=True)
    meta.to_csv(OUT_DIR / "gse54083_sample_metadata.tsv", sep="\t", index=False)

    r_cols = meta.loc[meta["Group"] == "Ruptured", "Sample"].tolist()
    u_cols = meta.loc[meta["Group"] == "Unruptured", "Sample"].tolist()

    rows = []
    for gene, row in expr.iterrows():
        rv = pd.to_numeric(row[r_cols], errors="coerce").dropna().values.astype(float)
        uv = pd.to_numeric(row[u_cols], errors="coerce").dropna().values.astype(float)
        if len(rv) < 2 or len(uv) < 2:
            continue
        t_stat, p_val = ttest_ind(rv, uv, equal_var=False, nan_policy="omit")
        log2fc = float(np.nanmean(rv) - np.nanmean(uv))
        rows.append(
            {
                "Gene": gene,
                "Mean_Ruptured_log2_signal": float(np.nanmean(rv)),
                "Mean_Unruptured_log2_signal": float(np.nanmean(uv)),
                "Log2FC_RvU": log2fc,
                "T_statistic": float(t_stat) if pd.notna(t_stat) else np.nan,
                "P_value": float(p_val) if pd.notna(p_val) else np.nan,
            }
        )

    res = pd.DataFrame(rows).sort_values("P_value", ascending=True).reset_index(drop=True)
    if not res.empty:
        res["FDR_BH"] = multipletests(res["P_value"].fillna(1.0).values, method="fdr_bh")[1]
    else:
        res["FDR_BH"] = []

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
        np.sign(merged["Log2FC_RvU"]).replace(0, np.nan)
        == np.sign(merged["GSE122897_log2FC_compare"]).replace(0, np.nan)
    )
    merged["In_PHAGESG_shortlist"] = merged["Shortlist_GSE122897_log2FC"].notna()

    merged.to_csv(OUT_DIR / "gse54083_gene_level_validation_with_gse122897.tsv", sep="\t", index=False)
    res.to_csv(OUT_DIR / "gse54083_gene_level_validation.tsv", sep="\t", index=False)

    top = merged.sort_values("P_value", ascending=True).head(30)
    top.to_csv(OUT_DIR / "gse54083_top30_genes.tsv", sep="\t", index=False)

    shortlist_hits = merged[merged["In_PHAGESG_shortlist"] == True].copy()  # noqa: E712
    if not shortlist_hits.empty:
        shortlist_hits = shortlist_hits.sort_values("P_value", ascending=True)
        shortlist_hits.to_csv(OUT_DIR / "gse54083_shortlist_overlap.tsv", sep="\t", index=False)

    summary = {
        "dataset": "GSE54083",
        "platform": "Agilent microarray",
        "n_total_samples_loaded": int(len(meta)),
        "n_ruptured": int(len(r_cols)),
        "n_unruptured": int(len(u_cols)),
        "n_controls_sta": int((meta["Group"] == "Control").sum()),
        "n_genes_tested": int(len(res)),
        "n_genes_fdr_lt_0_05": int((res["FDR_BH"] < 0.05).sum()),
        "n_genes_p_lt_0_05": int((res["P_value"] < 0.05).sum()),
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
    }
    (OUT_DIR / "gse54083_validation_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print("\nTop 15 genes:")
    print(top[["Gene", "Log2FC_RvU", "P_value", "FDR_BH", "GSE122897_log2FC", "DirectionMatch_GSE122897"]].head(15).to_string(index=False))


if __name__ == "__main__":
    main()
