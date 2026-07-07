"""
Load tabular snippets for Results-section numeric detail in the revised manuscript.
"""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

UNRESOLVED_SYMBOL_REPLACEMENTS = {
    "Ccd": "Entrez 124220 (annotation unresolved)",
}


def _fmt_sci(x: float | None, nd: int = 2) -> str:
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return ""
    try:
        xf = float(x)
    except (TypeError, ValueError):
        return str(x)
    ax = abs(xf)
    if ax == 0:
        return "0"
    if ax < 1e-300:
        return "<1e-300"
    if ax >= 1e4 or ax < 10 ** (-nd):
        return f"{xf:.{nd}e}"
    return f"{xf:.{nd}f}".rstrip("0").rstrip(".")


def load_gene_symbol_map(base: Path) -> pd.DataFrame:
    frames = []

    mp = base / "results/gene_id_mapping_unified.tsv"
    if mp.is_file():
        m = pd.read_csv(mp, sep="\t")
        if "Original_ID" in m.columns and "Gene_Symbol" in m.columns:
            m = m[["Original_ID", "Gene_Symbol"]].dropna(subset=["Gene_Symbol"])
            m["Gene_key"] = m["Original_ID"].astype(str).str.replace(r"\.0$", "", regex=True)
            frames.append(m[["Gene_key", "Gene_Symbol"]])

    for p in [
        base / "results/gse122897_deseq2_top_annotated.tsv",
        base / "results/gse122897_top_genes_annotated.tsv",
    ]:
        if not p.is_file():
            continue
        df = pd.read_csv(p, sep="\t")
        if "Gene_ID" in df.columns and "Gene_Symbol" in df.columns:
            df = df[["Gene_ID", "Gene_Symbol"]].dropna(subset=["Gene_Symbol"])
            df["Gene_key"] = df["Gene_ID"].astype(str).str.replace(r"\.0$", "", regex=True)
            frames.append(df[["Gene_key", "Gene_Symbol"]])
        elif "Gene" in df.columns and "Gene_Symbol" in df.columns:
            df = df[["Gene", "Gene_Symbol"]].dropna(subset=["Gene_Symbol"])
            df["Gene_key"] = df["Gene"].astype(str).str.replace(r"\.0$", "", regex=True)
            frames.append(df[["Gene_key", "Gene_Symbol"]])

    if not frames:
        return pd.DataFrame(columns=["Gene_key", "Gene_Symbol"])

    m = pd.concat(frames, ignore_index=True)
    m["Gene_Symbol"] = m["Gene_Symbol"].replace(UNRESOLVED_SYMBOL_REPLACEMENTS)
    return m.drop_duplicates(subset=["Gene_key"], keep="first")


def annotate_deseq2_with_symbols(de: pd.DataFrame, sym: pd.DataFrame) -> pd.DataFrame:
    out = de.copy()
    out["Gene_key"] = out["Gene_ID"].astype(str).str.replace(r"\.0$", "", regex=True)
    if not sym.empty:
        out = out.merge(sym, on="Gene_key", how="left")
    out["Gene_Symbol"] = out["Gene_Symbol"].replace(UNRESOLVED_SYMBOL_REPLACEMENTS)
    out["display_symbol"] = out["Gene_Symbol"].where(out["Gene_Symbol"].notna(), out["Gene_key"])
    return out


def load_deseq2_top_table(base: Path, n: int = 12) -> pd.DataFrame:
    p = base / "results/gse122897_deseq2_results.tsv"
    if not p.is_file():
        return pd.DataFrame()
    de = pd.read_csv(p, sep="\t")
    sym = load_gene_symbol_map(base)
    de = annotate_deseq2_with_symbols(de, sym)
    de = de.sort_values("padj", ascending=True).head(n)
    return pd.DataFrame(
        {
            "Gene": de["display_symbol"].astype(str),
            "log2FC": de["log2FoldChange"].round(3),
            "padj": de["padj"].map(lambda x: _fmt_sci(x, 2)),
            "baseMean": de["baseMean"].round(2),
        }
    )


def load_lm22_contrast_table(base: Path) -> pd.DataFrame:
    p = base / "results/cibersortx/lm22_ruptured_vs_unruptured.tsv"
    if not p.is_file():
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t")
    df = df[df["P_value"].notna() & (df["P_value"].astype(str) != "")]
    df = df.sort_values("P_value", ascending=True)
    return pd.DataFrame(
        {
            "Cell type": df["Cell_type"],
            "Mean R": df["Mean_Ruptured"].round(3),
            "Mean U": df["Mean_Unruptured"].round(3),
            "Δ(R−U)": df["Diff_R_minus_U"].round(3),
            "P": df["P_value"].map(lambda x: _fmt_sci(x, 3)),
            "FDR": df["FDR_BH"].map(lambda x: _fmt_sci(x, 3)),
        }
    )


def load_immune_adjusted_exemplars(base: Path, n: int = 10) -> pd.DataFrame:
    p = base / "results/cibersortx/deg_immune_covariate_adjusted.tsv"
    if not p.is_file():
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t")
    sym = load_gene_symbol_map(base)
    df = df.copy()
    df["Gene_key"] = df["Gene_ID"].astype(str).str.replace(r"\.0$", "", regex=True)
    if not sym.empty:
        df = df.merge(sym, on="Gene_key", how="left")
    sl = base / "results/phases_g/shortlist_evidence.tsv"
    if sl.is_file():
        s = pd.read_csv(sl, sep="\t")
        if "Gene_ID" in s.columns and "symbol" in s.columns:
            s = s[["Gene_ID", "symbol"]].astype(str)
            s["Gene_key"] = s["Gene_ID"].str.replace(r"\.0$", "", regex=True)
            s = s.drop_duplicates(subset=["Gene_key"])
            df = df.merge(s[["Gene_key", "symbol"]], on="Gene_key", how="left")
            df["Gene"] = df["Gene_Symbol"].where(df["Gene_Symbol"].notna(), df["symbol"])
            df["Gene"] = df["Gene"].where(df["Gene"].notna(), df["Gene_key"])
        else:
            df["Gene"] = df["Gene_Symbol"].where(df["Gene_Symbol"].notna(), df["Gene_key"])
    else:
        df["Gene"] = df["Gene_Symbol"].where(df["Gene_Symbol"].notna(), df["Gene_key"])
    df = df[df["DESeq2_padj"] < 0.05].sort_values("LM_immune_adj_P_Ruptured", ascending=True).head(n)
    return pd.DataFrame(
        {
            "Gene": df["Gene"].astype(str),
            "DESeq2 log2FC": df["DESeq2_log2FC"].round(3),
            "DESeq2 padj": df["DESeq2_padj"].map(lambda x: _fmt_sci(x, 2)),
            "LM coef (rupture)": df["LM_immune_adj_log2_coef_Ruptured"].round(4),
            "LM P (rupture)": df["LM_immune_adj_P_Ruptured"].map(lambda x: _fmt_sci(x, 3)),
        }
    )


def gse159610_shortlist_stats(base: Path) -> tuple[int, int, int]:
    p = base / "results/phases_g/shortlist_evidence.tsv"
    if not p.is_file():
        return 0, 0, 0
    sl = pd.read_csv(p, sep="\t")
    n_total = len(sl)
    col = "GSE159610_Log2FC_aneurysm_vs_control"
    if col not in sl.columns:
        return n_total, 0, 0
    has = int(sl[col].notna().sum())
    dm = sl["direction_match_rupture_and_aneurysm"] == True  # noqa: E712
    n_dir = int((sl[col].notna() & dm).sum())
    return n_total, has, n_dir


def load_gse159610_example_table(base: Path, n: int = 12) -> pd.DataFrame:
    p = base / "results/phases_g/shortlist_evidence.tsv"
    if not p.is_file():
        return pd.DataFrame()
    sl = pd.read_csv(p, sep="\t")
    col_fc = "GSE159610_Log2FC_aneurysm_vs_control"
    col_p = "GSE159610_P"
    if col_fc not in sl.columns:
        return pd.DataFrame()
    m = sl[sl[col_fc].notna()].copy()
    m["abs_fc"] = m[col_fc].abs()
    if col_p in m.columns:
        m["_p_sort"] = pd.to_numeric(m[col_p], errors="coerce")
        m = m.sort_values(["_p_sort", "abs_fc"], ascending=[True, False]).head(n)
    else:
        m = m.sort_values("abs_fc", ascending=False).head(n)
    dm = m["direction_match_rupture_and_aneurysm"].map({True: "Yes", False: "No"}).fillna("")
    return pd.DataFrame(
        {
            "Gene": m["symbol"].astype(str),
            "GSE122897 DESeq2 log2FC": m["DESeq2_log2FC"].round(3),
            "GSE159610 log2FC (A vs C)": m[col_fc].round(3),
            "GSE159610 P": m[col_p].map(lambda x: _fmt_sci(x, 2)) if col_p in m.columns else "",
            "Dir. match": dm,
        }
    )


def load_gse26969_example_table(base: Path, n: int = 12) -> pd.DataFrame:
    p = base / "results/external_validation_gse26969/gse26969_shortlist_overlap.tsv"
    if not p.is_file():
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t")
    if "P_value" not in df.columns:
        return pd.DataFrame()
    df = df.sort_values(["P_value", "DirectionMatch_GSE122897"], ascending=[True, False]).head(n)
    dm = df["DirectionMatch_GSE122897"].map({True: "Yes", False: "No"}).fillna("")
    return pd.DataFrame(
        {
            "Gene": df["Gene"].astype(str),
            "GSE26969 log2FC (A vs C)": df["Log2FC_AvC"].round(3),
            "GSE26969 P": df["P_value"].map(lambda x: _fmt_sci(x, 2)),
            "GSE122897 log2FC": df["GSE122897_log2FC_compare"].round(3),
            "Dir. match": dm,
        }
    )


def load_gse36791_example_table(base: Path, n: int = 12) -> pd.DataFrame:
    p = base / "results/external_validation_gse36791/gse36791_shortlist_overlap.tsv"
    if not p.is_file():
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t")
    if "P_value" not in df.columns:
        return pd.DataFrame()
    df = df.sort_values(["P_value", "DirectionMatch_GSE122897"], ascending=[True, False]).head(n)
    dm = df["DirectionMatch_GSE122897"].map({True: "Yes", False: "No"}).fillna("")
    return pd.DataFrame(
        {
            "Gene": df["Gene"].astype(str),
            "GSE36791 log2FC (R vs C)": df["Log2FC_RvC"].round(3),
            "GSE36791 P": df["P_value"].map(lambda x: _fmt_sci(x, 2)),
            "GSE122897 log2FC": df["GSE122897_log2FC_compare"].round(3),
            "Dir. match": dm,
        }
    )


def load_gse54083_example_table(base: Path, n: int = 12) -> pd.DataFrame:
    p = base / "results/external_validation_gse54083/gse54083_shortlist_overlap.tsv"
    if not p.is_file():
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t")
    if "P_value" not in df.columns:
        return pd.DataFrame()
    df = df.sort_values(["P_value", "DirectionMatch_GSE122897"], ascending=[True, False]).head(n)
    dm = df["DirectionMatch_GSE122897"].map({True: "Yes", False: "No"}).fillna("")
    return pd.DataFrame(
        {
            "Gene": df["Gene"].astype(str),
            "GSE54083 log2FC (R vs U)": df["Log2FC_RvU"].round(3),
            "GSE54083 P": df["P_value"].map(lambda x: _fmt_sci(x, 2)),
            "GSE122897 log2FC": df["GSE122897_log2FC_compare"].round(3),
            "Dir. match": dm,
        }
    )


def load_gse13353_example_table(base: Path) -> pd.DataFrame:
    p = base / "results/external_validation_gse13353/gse13353_validation_table_with_gse122897.tsv"
    if not p.is_file():
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t")
    keep = [
        ("Gene", "Gene"),
        ("Log2FC_RvU", "GSE13353 log2FC (R vs U)"),
        ("P_value", "GSE13353 P"),
        ("GSE122897_log2FC", "GSE122897 log2FC"),
        ("Match", "Dir. match"),
    ]
    out = pd.DataFrame()
    for src, dst in keep:
        if src in df.columns:
            out[dst] = df[src]
    if "GSE13353 P" in out.columns:
        out["GSE13353 P"] = out["GSE13353 P"].map(lambda x: _fmt_sci(x, 2))
    if "Dir. match" in out.columns:
        out["Dir. match"] = out["Dir. match"].map({True: "Yes", False: "No"}).fillna(out["Dir. match"])
    return out


def load_cibersortx_merged_exemplars(base: Path, n: int = 8) -> pd.DataFrame:
    """Representative CIBERSORTx Job 7 mixtures after GSM merge (wide table subset)."""
    p = base / "results/cibersortx/cibersortx_merged_with_group.tsv"
    if not p.is_file():
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t")
    if "Correlation" not in df.columns or "Sample" not in df.columns:
        return pd.DataFrame()
    df = df.sort_values("Correlation", ascending=True)
    k = max(2, n // 2)
    low = df.head(k)
    high = df.tail(k)
    m = pd.concat([low, high], ignore_index=True).drop_duplicates(subset=["Sample"])
    m = m.head(n)
    cols = [
        "Sample",
        "Group",
        "Monocytes",
        "Macrophages M0",
        "Macrophages M1",
        "Macrophages M2",
        "T cells CD8",
        "P-value",
        "Correlation",
        "RMSE",
    ]
    use = [c for c in cols if c in m.columns]
    out = m[use].copy()
    for c in out.columns:
        if c in ("Sample", "Group"):
            continue
        out[c] = pd.to_numeric(out[c], errors="coerce").round(3)
    return out


def load_catalogue_eqtl_examples(base: Path, n: int = 8) -> pd.DataFrame:
    """eQTL Catalogue hits on the shortlist (where rsID was returned)."""
    p = base / "results/phases_g/shortlist_evidence.tsv"
    if not p.is_file():
        return pd.DataFrame()
    sl = pd.read_csv(p, sep="\t")
    if "blood_eqtl_top_rsid" not in sl.columns:
        return pd.DataFrame()
    sub = sl[sl["blood_eqtl_top_rsid"].notna() & (sl["blood_eqtl_top_rsid"].astype(str) != "")].copy()
    sub = sub.sort_values("DESeq2_padj", ascending=True).head(n)
    beta = sub["blood_eqtl_beta"] if "blood_eqtl_beta" in sub.columns else None
    return pd.DataFrame(
        {
            "Gene": sub["symbol"].astype(str),
            "Catalogue rsID": sub["blood_eqtl_top_rsid"].astype(str),
            "Catalogue P": sub["blood_eqtl_top_pvalue"].map(lambda x: _fmt_sci(x, 2))
            if "blood_eqtl_top_pvalue" in sub.columns
            else "",
            "β (blood)": beta.map(lambda x: _fmt_sci(x, 3)) if beta is not None else "",
            "Dataset": sub["blood_eqtl_dataset_id"].astype(str)
            if "blood_eqtl_dataset_id" in sub.columns
            else "",
        }
    )


def load_phases_g_genetics_table(base: Path, n: int = 14) -> pd.DataFrame:
    p = base / "results/phases_g/shortlist_evidence.tsv"
    if not p.is_file():
        return pd.DataFrame()
    sl = pd.read_csv(p, sep="\t")
    sub = sl[sl["eQTLGen_top_SNP"].notna() & (sl["eQTLGen_top_SNP"].astype(str) != "")].copy()
    sub = sub.sort_values("DESeq2_padj", ascending=True).head(n)
    gwas_p = sub["IA_GWAS_P_at_eQTLGen_SNP"] if "IA_GWAS_P_at_eQTLGen_SNP" in sub.columns else None
    gwas_b = sub["IA_GWAS_BETA_at_eQTLGen_SNP"] if "IA_GWAS_BETA_at_eQTLGen_SNP" in sub.columns else None
    return pd.DataFrame(
        {
            "Gene": sub["symbol"].astype(str),
            "DESeq2 log2FC": sub["DESeq2_log2FC"].round(3),
            "eQTLGen SNP": sub["eQTLGen_top_SNP"].astype(str),
            "eQTLGen P": sub["eQTLGen_top_Pvalue"].map(lambda x: _fmt_sci(x, 2)),
            "IA GWAS P@SNP": gwas_p.map(lambda x: _fmt_sci(x, 2)) if gwas_p is not None else "",
            "IA GWAS β@SNP": gwas_b.map(lambda x: _fmt_sci(x, 3)) if gwas_b is not None else "",
        }
    )


def deseq2_top_prose_sentence(base: Path, n: int = 6) -> str:
    p = base / "results/gse122897_deseq2_results.tsv"
    if not p.is_file():
        return ""
    de = pd.read_csv(p, sep="\t")
    sym = load_gene_symbol_map(base)
    de = annotate_deseq2_with_symbols(de, sym)
    de = de.sort_values("padj", ascending=True).head(n)
    parts = []
    for _, r in de.iterrows():
        g = str(r["display_symbol"])
        parts.append(
            f"{g} (log₂FC={float(r['log2FoldChange']):.3f}, padj={_fmt_sci(r['padj'], 2)}, baseMean={float(r['baseMean']):.2f})"
        )
    return "Strongest DESeq2 signals (lowest padj) included " + "; ".join(parts) + "."
