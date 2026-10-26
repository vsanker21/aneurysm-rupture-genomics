from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests


def clean_symbol(raw: object) -> str | None:
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


def extract_symbol_from_assignment(raw: object) -> str | None:
    if pd.isna(raw):
        return None
    text = str(raw).strip()
    if not text or text in {"---", "nan", "NA", "None"}:
        return None

    symbols: list[str] = []
    for block in text.split("///"):
        parts = [part.strip() for part in block.split("//")]
        if len(parts) < 2:
            continue
        candidate = parts[1].strip()
        if not candidate or candidate in {"---", "NA", "None"}:
            continue
        if candidate not in symbols:
            symbols.append(candidate)
    if len(symbols) != 1:
        return None
    return symbols[0]


def load_gse122897_deseq2(base: Path) -> pd.DataFrame:
    de = pd.read_csv(base / "results" / "gse122897_deseq2_results.tsv", sep="\t")
    mapping = pd.read_csv(base / "results" / "gene_id_mapping_unified.tsv", sep="\t")
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


def load_shortlist_effects(base: Path) -> pd.DataFrame:
    shortlist_path = base / "results" / "phases_g" / "shortlist_evidence.tsv"
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


def add_gse122897_reference(result_df: pd.DataFrame, *, base: Path, effect_col: str) -> pd.DataFrame:
    de122 = load_gse122897_deseq2(base)
    sl_effects = load_shortlist_effects(base)
    merged = result_df.merge(de122, on="Gene", how="left")
    merged = merged.merge(sl_effects, on="Gene", how="left")
    merged["GSE122897_log2FC_compare"] = merged["Shortlist_GSE122897_log2FC"].where(
        merged["Shortlist_GSE122897_log2FC"].notna(), merged["GSE122897_log2FC"]
    )
    merged["GSE122897_padj_compare"] = merged["Shortlist_GSE122897_padj"].where(
        merged["Shortlist_GSE122897_padj"].notna(), merged["GSE122897_padj"]
    )
    merged["DirectionMatch_GSE122897"] = (
        np.sign(merged[effect_col]).replace(0, np.nan)
        == np.sign(merged["GSE122897_log2FC_compare"]).replace(0, np.nan)
    )
    merged["In_PHAGESG_shortlist"] = merged["Shortlist_GSE122897_log2FC"].notna()
    return merged


def compute_gene_level_welch(
    gene_expr: pd.DataFrame,
    *,
    group_a_cols: list[str],
    group_b_cols: list[str],
    group_a_label: str,
    group_b_label: str,
    effect_name: str,
) -> pd.DataFrame:
    rows: list[dict[str, float | str]] = []
    for gene, row in gene_expr.iterrows():
        av = pd.to_numeric(row[group_a_cols], errors="coerce").dropna().values.astype(float)
        bv = pd.to_numeric(row[group_b_cols], errors="coerce").dropna().values.astype(float)
        if len(av) < 2 or len(bv) < 2:
            continue
        t_stat, p_val = ttest_ind(av, bv, equal_var=False, nan_policy="omit")
        rows.append(
            {
                "Gene": gene,
                f"Mean_{group_a_label}": float(np.nanmean(av)),
                f"Mean_{group_b_label}": float(np.nanmean(bv)),
                effect_name: float(np.nanmean(av) - np.nanmean(bv)),
                "T_statistic": float(t_stat) if pd.notna(t_stat) else np.nan,
                "P_value": float(p_val) if pd.notna(p_val) else np.nan,
            }
        )
    if not rows:
        raise RuntimeError("No genes passed minimum sample requirements for Welch testing.")
    result = pd.DataFrame(rows).sort_values("P_value", ascending=True).reset_index(drop=True)
    result["FDR_BH"] = multipletests(result["P_value"].fillna(1.0).values, method="fdr_bh")[1]
    return result
