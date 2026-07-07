"""
Complete Extended Data Table 3 / 300-gene shortlist annotations.

- Resolve HGNC-approved gene symbols via local maps, MyGene.info, and NCBI Gene.
- Backfill GSE159610 / GSE106520 bridge fields from project results.
- Mark genes not on a blood platform as NR (not represented).
- Regenerate submission TSV and Tables_Revised.docx Extended Data Table 3.
"""
from __future__ import annotations

import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy.stats import norm, ttest_ind

try:
    from project_paths import get_project_root

    BASE = get_project_root()
except ImportError:
    BASE = Path(__file__).resolve().parent
OUT_TSV = BASE / "results" / "supplementary" / "tables" / "Extended_Data_Table_3_300gene_shortlist.tsv"
SUPP_TSV = BASE / "results/supplementary/tables/Supplementary_Table_S2_300gene_shortlist.tsv"
EVIDENCE = BASE / "results/phases_g/shortlist_evidence.tsv"
G106_OVERLAP = BASE / "results/external_validation_gse106520/gse106520_shortlist_overlap.tsv"
G106_CONTEXT = BASE / "results/external_validation_gse106520/gse106520_gene_level_context_with_gse122897.tsv"
G106_FPKM = BASE / "results/external_validation_gse106520/gse106520_gene_level_log2fpkm.tsv"
G106_META = BASE / "results/external_validation_gse106520/gse106520_sample_metadata.tsv"
G159610_DE = BASE / "results/gse159610_full_de_results_annotated.tsv"

PANEL_GENES = {"DYNLT1", "OASL", "MPP1", "TMEM130", "MARCO", "MCOLN1", "RNASET2", "ANK3"}
WALL_CORE = {"TGFBI", "ANPEP", "MPP1", "MARCO", "CLEC5A", "CPXM2"}

DISPLAY_COLUMNS = [
    "Shortlist_rank",
    "Gene_symbol",
    "DESeq2_log2FC",
    "DESeq2_padj",
    "LM_immune_adj_P_Ruptured",
    "GSE159610_Log2FC_aneurysm_vs_control",
    "GSE159610_P",
    "GSE106520_log2FC",
    "GSE106520_P",
    "Blood_both_direction_concordant",
    "Bridge_Stouffer_Z",
    "Bridge_meta_P",
    "In_8gene_panel",
    "In_6gene_wall_core",
]


def _ranked_shortlist_with_ids() -> pd.DataFrame:
    from generate_supplementary_materials_complete import build_table_s2

    df, _ = build_table_s2()
    ev = pd.read_csv(EVIDENCE, sep="\t")
    g106 = pd.read_csv(G106_OVERLAP, sep="\t")
    sl = ev.rename(columns={"symbol": "Gene_symbol"}) if "symbol" in ev.columns else ev.copy()
    if "Gene_symbol" not in sl.columns:
        sl["Gene_symbol"] = sl.get("symbol")
    merged = sl.merge(g106[["Gene", "P_value"]].rename(columns={"Gene": "Gene_symbol"}), on="Gene_symbol", how="left")
    p1 = pd.to_numeric(merged.get("GSE159610_P"), errors="coerce").fillna(1.0).clip(1e-300, 1.0)
    p2 = pd.to_numeric(merged["P_value"], errors="coerce").fillna(1.0).clip(1e-300, 1.0)
    merged["Bridge_Stouffer_Z"] = (norm.isf(p1 / 2) + norm.isf(p2 / 2)) / np.sqrt(2)
    merged = merged.sort_values(["Bridge_Stouffer_Z", "DESeq2_padj"], ascending=[False, True]).reset_index(drop=True)
    merged["Shortlist_rank"] = np.arange(1, len(merged) + 1)
    id_map = merged.set_index("Shortlist_rank")[["Gene_ID", "ensembl_gene"]]
    out = df.merge(id_map, on="Shortlist_rank", how="left")
    return out


def _load_local_symbol_maps() -> dict[str, str]:
    mapping: dict[str, str] = {}
    for path in [
        BASE / "results/gene_symbol_to_entrez_mapping.tsv",
        BASE / "results/gene_id_mapping_unified.tsv",
        BASE / "results/gse122897_deseq2_top_annotated.tsv",
        BASE / "results/gse122897_top_genes_annotated.tsv",
    ]:
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        if "Entrez_ID" in df.columns and "Gene_Symbol" in df.columns:
            for _, r in df.dropna(subset=["Entrez_ID", "Gene_Symbol"]).iterrows():
                try:
                    eid = str(int(float(r["Entrez_ID"])))
                except (TypeError, ValueError):
                    continue
                sym = str(r["Gene_Symbol"]).strip()
                if sym and sym.lower() != "nan":
                    mapping[eid] = sym
        if "Gene_ID" in df.columns and "Gene_Symbol" in df.columns:
            for _, r in df.dropna(subset=["Gene_ID", "Gene_Symbol"]).iterrows():
                try:
                    eid = str(int(float(r["Gene_ID"])))
                except (TypeError, ValueError):
                    continue
                sym = str(r["Gene_Symbol"]).strip()
                if sym and sym.lower() != "nan":
                    mapping.setdefault(eid, sym)
    ev = pd.read_csv(EVIDENCE, sep="\t")
    for _, r in ev.dropna(subset=["Gene_ID", "symbol"]).iterrows():
        mapping[str(int(r["Gene_ID"]))] = str(r["symbol"]).strip()
    return mapping


def _mygene_batch(entrez_ids: list[str]) -> dict[str, str]:
    url = "https://mygene.info/v3/gene"
    out: dict[str, str] = {}
    for i in range(0, len(entrez_ids), 500):
        chunk = entrez_ids[i : i + 500]
        try:
            r = requests.post(
                url,
                json={"ids": chunk, "fields": "symbol,alias,name", "species": "human"},
                timeout=90,
            )
            r.raise_for_status()
            for item in r.json():
                if not isinstance(item, dict) or item.get("notfound"):
                    continue
                q = str(item.get("query") or item.get("_id") or "").strip()
                sym = str(item.get("symbol") or "").strip()
                if q and sym:
                    out[q] = sym
        except Exception:
            pass
        time.sleep(0.2)
    return out


def _ncbi_gene_symbols(entrez_ids: list[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    for i in range(0, len(entrez_ids), 100):
        chunk = entrez_ids[i : i + 100]
        try:
            r = requests.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                params={"db": "gene", "id": ",".join(chunk), "retmode": "json"},
                timeout=60,
            )
            r.raise_for_status()
            result = r.json().get("result", {})
            for eid in chunk:
                rec = result.get(str(eid), {})
                if not isinstance(rec, dict):
                    continue
                sym = str(rec.get("nomenclaturesymbol") or "").strip()
                if sym:
                    out[str(eid)] = sym
                    continue
                desc = str(rec.get("description") or "")
                m = re.search(r"\b(LOC\d+)\b", desc)
                if m:
                    out[str(eid)] = m.group(1)
        except Exception:
            pass
        time.sleep(0.2)
    return out


def resolve_hgnc_symbols(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["Gene_ID"] = out["Gene_ID"].astype(str).str.replace(r"\.0$", "", regex=True)
    local = _load_local_symbol_maps()
    entrez = out["Gene_ID"].dropna().unique().tolist()
    mygene = _mygene_batch(entrez)
    ncbi = _ncbi_gene_symbols([e for e in entrez if e not in mygene or not mygene.get(e)])

    symbols: list[str] = []
    sources: list[str] = []
    for _, row in out.iterrows():
        eid = str(row["Gene_ID"])
        candidates = [
            ("existing", str(row.get("Gene_symbol", "")).replace("ENTREZ:", "")),
            ("local", local.get(eid, "")),
            ("mygene", mygene.get(eid, "")),
            ("ncbi", ncbi.get(eid, "")),
        ]
        sym = ""
        src = ""
        for name, val in candidates:
            val = str(val).strip()
            if val and val.lower() not in {"nan", "none", eid}:
                sym = val
                src = name
                break
        if not sym:
            sym = f"LOC{eid}" if eid.isdigit() else eid
            src = "fallback"
        symbols.append(sym)
        sources.append(src)

    out["Gene_symbol"] = symbols
    out["symbol_source"] = sources
    out["In_8gene_panel"] = out["Gene_symbol"].isin(PANEL_GENES)
    out["In_6gene_wall_core"] = out["Gene_symbol"].isin(WALL_CORE)
    return out


def _load_g159610_lookup() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    de = pd.read_csv(G159610_DE, sep="\t")
    de["Entrez_ID"] = de["Entrez_ID"].astype(str).str.replace(r"\.0$", "", regex=True)
    de["ens_base"] = de["Gene_ID"].astype(str).str.split(".").str[0]
    by_sym = (
        de.dropna(subset=["Gene_Symbol"])
        .drop_duplicates("Gene_Symbol", keep="first")
        .set_index("Gene_Symbol", drop=False)
    )
    by_eid = (
        de.dropna(subset=["Entrez_ID"])
        .drop_duplicates("Entrez_ID", keep="first")
        .set_index("Entrez_ID", drop=False)
    )
    by_ens = (
        de.dropna(subset=["ens_base"])
        .drop_duplicates("ens_base", keep="first")
        .set_index("ens_base", drop=False)
    )
    return by_sym, by_eid, by_ens


def _load_g106_lookup() -> pd.DataFrame:
    ctx = pd.read_csv(G106_CONTEXT, sep="\t")
    return ctx.drop_duplicates("Gene", keep="first").set_index("Gene")


def _welch_g106_for_gene(gene: str) -> tuple[float, float] | None:
    if not G106_FPKM.exists() or not G106_META.exists():
        return None
    fpkm = pd.read_csv(G106_FPKM, sep="\t", index_col=0)
    if gene not in fpkm.index:
        return None
    meta = pd.read_csv(G106_META, sep="\t")
    aneurysm = meta.loc[meta["Group"] == "Aneurysm", "Sample"].astype(str).tolist()
    control = meta.loc[meta["Group"] == "Control", "Sample"].astype(str).tolist()
    aneurysm = [c for c in aneurysm if c in fpkm.columns]
    control = [c for c in control if c in fpkm.columns]
    if len(aneurysm) < 2 or len(control) < 2:
        return None
    av = pd.to_numeric(fpkm.loc[gene, aneurysm], errors="coerce").dropna().values.astype(float)
    bv = pd.to_numeric(fpkm.loc[gene, control], errors="coerce").dropna().values.astype(float)
    if len(av) < 2 or len(bv) < 2:
        return None
    _, p = ttest_ind(av, bv, equal_var=False)
    return float(np.nanmean(av) - np.nanmean(bv)), float(p)


def backfill_bridge_cohorts(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    by_sym, by_eid, by_ens = _load_g159610_lookup()
    g106 = _load_g106_lookup()

    fc159 = []
    p159 = []
    st159 = []
    fc106 = []
    p106 = []
    st106 = []

    for _, row in out.iterrows():
        sym = row["Gene_symbol"]
        eid = str(row["Gene_ID"])
        ens = str(row.get("ensembl_gene", "")).split(".")[0] if pd.notna(row.get("ensembl_gene")) else ""

        # GSE159610
        l2_159 = row.get("GSE159610_Log2FC_aneurysm_vs_control")
        pv_159 = row.get("GSE159610_P")
        status_159 = "measured"
        if pd.isna(pv_159):
            hit = None
            if sym in by_sym.index:
                hit = by_sym.loc[sym]
            elif eid in by_eid.index:
                hit = by_eid.loc[eid]
            elif ens and ens in by_ens.index:
                hit = by_ens.loc[ens]
            if hit is not None:
                l2_159 = float(hit["Log2FC"])
                pv_159 = float(hit["P_value"])
            else:
                status_159 = "NR"
                l2_159 = np.nan
                pv_159 = np.nan

        # GSE106520
        l2_106 = row.get("GSE106520_log2FC")
        pv_106 = row.get("GSE106520_P")
        status_106 = "measured"
        if pd.isna(pv_106):
            if sym in g106.index:
                l2_106 = float(g106.loc[sym, "Log2FC_AvC"])
                pv_106 = float(g106.loc[sym, "P_value"])
            else:
                welch = _welch_g106_for_gene(sym)
                if welch is not None:
                    l2_106, pv_106 = welch
                else:
                    status_106 = "NR"
                    l2_106 = np.nan
                    pv_106 = np.nan

        fc159.append(l2_159)
        p159.append(pv_159)
        st159.append(status_159)
        fc106.append(l2_106)
        p106.append(pv_106)
        st106.append(status_106)

    out["GSE159610_Log2FC_aneurysm_vs_control"] = fc159
    out["GSE159610_P"] = p159
    out["GSE159610_platform"] = st159
    out["GSE106520_log2FC"] = fc106
    out["GSE106520_P"] = p106
    out["GSE106520_platform"] = st106

    rdir = np.sign(pd.to_numeric(out["DESeq2_log2FC"], errors="coerce"))
    d159 = np.sign(pd.to_numeric(out["GSE159610_Log2FC_aneurysm_vs_control"], errors="coerce"))
    d106 = np.sign(pd.to_numeric(out["GSE106520_log2FC"], errors="coerce"))
    out["Blood_both_direction_concordant"] = (
        (out["GSE159610_platform"] == "measured")
        & (out["GSE106520_platform"] == "measured")
        & (d159 == rdir)
        & (d106 == rdir)
        & (rdir != 0)
    )

    p1 = pd.to_numeric(out["GSE159610_P"], errors="coerce")
    p2 = pd.to_numeric(out["GSE106520_P"], errors="coerce")
    p1_fill = p1.where(out["GSE159610_platform"] == "measured", 1.0).fillna(1.0).clip(1e-300, 1.0)
    p2_fill = p2.where(out["GSE106520_platform"] == "measured", 1.0).fillna(1.0).clip(1e-300, 1.0)
    out["Bridge_Stouffer_Z"] = (norm.isf(p1_fill / 2) + norm.isf(p2_fill / 2)) / np.sqrt(2)
    out["Bridge_meta_P"] = 2 * norm.sf(np.abs(out["Bridge_Stouffer_Z"]))
    return out


def export_table(df: pd.DataFrame) -> pd.DataFrame:
    export = df[DISPLAY_COLUMNS].copy()
    bridge_cols = [
        "GSE159610_Log2FC_aneurysm_vs_control",
        "GSE159610_P",
        "GSE106520_log2FC",
        "GSE106520_P",
    ]
    for col in bridge_cols:
        export[col] = export[col].apply(lambda x: "NR" if pd.isna(x) else x)
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    SUPP_TSV.parent.mkdir(parents=True, exist_ok=True)
    export.to_csv(OUT_TSV, sep="\t", index=False)
    export.to_csv(SUPP_TSV, sep="\t", index=False)
    return export


def audit_bottom_33(df: pd.DataFrame) -> None:
    bottom = df[df["Shortlist_rank"] >= 268]
    bridge = [
        "GSE159610_Log2FC_aneurysm_vs_control",
        "GSE159610_P",
        "GSE106520_log2FC",
        "GSE106520_P",
    ]
    entrez_like = bottom["Gene_symbol"].astype(str).str.startswith("ENTREZ:")
    print(f"Bottom 33 genes: entrez-like symbols={entrez_like.sum()}")
    for col in bridge:
        nr = bottom[col].isna().sum()
        print(f"  {col}: {nr} NaN (NR will display in Word)")
    if "GSE159610_platform" in df.columns:
        print(f"  GSE159610 NR: {(bottom['GSE159610_platform']=='NR').sum()}")
        print(f"  GSE106520 NR: {(bottom['GSE106520_platform']=='NR').sum()}")


def update_shortlist_evidence_symbols(df: pd.DataFrame) -> None:
    """Write resolved HGNC symbols back into shortlist_evidence.tsv."""
    if not EVIDENCE.exists():
        return
    ev = pd.read_csv(EVIDENCE, sep="\t")
    sym_map = (
        df.dropna(subset=["Gene_ID", "Gene_symbol"])
        .assign(Gene_ID=lambda x: x["Gene_ID"].astype(str).str.replace(r"\.0$", "", regex=True))
        .set_index("Gene_ID")["Gene_symbol"]
        .to_dict()
    )
    updated = 0
    for idx, row in ev.iterrows():
        eid = str(int(row["Gene_ID"]))
        if (pd.isna(row.get("symbol")) or str(row.get("symbol")).strip() == "") and eid in sym_map:
            ev.at[idx, "symbol"] = sym_map[eid]
            updated += 1
    if updated:
        ev.to_csv(EVIDENCE, sep="\t", index=False)
        print(f"Updated {updated} symbols in {EVIDENCE.name}")


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser(description="Complete 300-gene shortlist annotations.")
    ap.add_argument(
        "--tsv-only",
        action="store_true",
        help="Write TSV outputs only; do not update submission Word documents.",
    )
    args = ap.parse_args()

    print("Building ranked 300-gene shortlist...")
    df = _ranked_shortlist_with_ids()
    print("Resolving HGNC / official gene symbols...")
    df = resolve_hgnc_symbols(df)
    print("Backfilling bridge cohort statistics...")
    df = backfill_bridge_cohorts(df)
    update_shortlist_evidence_symbols(df)
    export = export_table(df)
    audit_bottom_33(df)

    unresolved = df["Gene_symbol"].astype(str).str.startswith("ENTREZ:").sum()
    print(f"Exported {len(export)} genes to {OUT_TSV.name}")
    print(f"ENTREZ-prefixed symbols remaining: {unresolved}")

    if not args.tsv_only:
        from add_extended_table_3_to_tables_revised import patch_tables_revised

        patch_tables_revised(df[DISPLAY_COLUMNS])
    print("Done.")


if __name__ == "__main__":
    main()
