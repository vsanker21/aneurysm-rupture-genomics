"""Compare frozen 22 vs 21 DESeq2 with explicit-label 21 vs 21.

Does not overwrite as-submitted discovery files.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

BASE = Path(r"E:\Aneurysms Genomics Project")
FROZEN = BASE / "results/deseq2_gse122897_as_submitted_22v21/gse122897_deseq2_results.tsv"
NEW = BASE / "results/deseq2_gse122897_21v21_explicit/gse122897_deseq2_results_21v21.tsv"
OUT = BASE / "results/deseq2_gse122897_21v21_explicit"
IMMUNE = BASE / "results/cibersortx/deg_immune_covariate_adjusted.tsv"
SHORTLIST = BASE / "results/phases_g/shortlist_evidence.tsv"
PANEL_TSV = BASE / "results/blood_bridge_panel/blood_bridge_panel_gene_summary.tsv"

CORE = {
    "TGFBI": "7045",
    "ANPEP": "290",
    "MPP1": "4354",
    "MARCO": "8685",
    "CLEC5A": "23601",
    "CPXM2": "119587",
}
PANEL = ["DYNLT1", "OASL", "MPP1", "TMEM130", "MARCO", "MCOLN1", "RNASET2", "ANK3"]
CORE_IN_SHORTLIST = ["TGFBI", "MPP1", "MARCO", "CPXM2"]


def _gid(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip('"').str.replace(r"\.0$", "", regex=True)


def _load_deseq(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df["Gene_ID"] = _gid(df["Gene_ID"])
    for c in ("log2FoldChange", "padj", "pvalue", "baseMean"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def _counts(df: pd.DataFrame) -> dict:
    sig = df["padj"] < 0.05
    return {
        "n_genes_tested": int(len(df)),
        "n_padj_lt_0.05": int(sig.sum()),
        "n_padj_lt_0.05_and_abs_lfc_ge_0.5": int((sig & (df["log2FoldChange"].abs() >= 0.5)).sum()),
    }


def _rebuild_shortlist(des: pd.DataFrame, immune: pd.DataFrame) -> pd.DataFrame:
    adj = immune.copy()
    adj["Gene_ID"] = _gid(adj["Gene_ID"])
    des_s = des[["Gene_ID", "log2FoldChange", "padj", "baseMean"]].rename(
        columns={"log2FoldChange": "DESeq2_log2FC", "padj": "DESeq2_padj"}
    )
    adj = adj.drop(columns=[c for c in ("DESeq2_log2FC", "DESeq2_padj", "baseMean") if c in adj.columns])
    m = adj.merge(des_s, on="Gene_ID", how="inner")
    m["DESeq2_padj"] = pd.to_numeric(m["DESeq2_padj"], errors="coerce")
    m["LM_immune_adj_P_Ruptured"] = pd.to_numeric(m["LM_immune_adj_P_Ruptured"], errors="coerce")
    short = m[(m["DESeq2_padj"] <= 0.05) & (m["LM_immune_adj_P_Ruptured"] <= 0.05)].copy()
    short = short.sort_values(["DESeq2_padj", "Gene_ID"]).head(300).reset_index(drop=True)
    short["new_rank"] = np.arange(1, len(short) + 1)
    return short


def _stouffer_two(p1, p2) -> tuple[float, float]:
    p1 = float(np.clip(p1, 1e-300, 1.0))
    p2 = float(np.clip(p2, 1e-300, 1.0))
    z = (norm.isf(p1 / 2.0) + norm.isf(p2 / 2.0)) / np.sqrt(2.0)
    return float(z), float(2.0 * norm.sf(abs(z)))


def main() -> None:
    old = _load_deseq(FROZEN)
    new = _load_deseq(NEW)
    old_c, new_c = _counts(old), _counts(new)

    # --- 1. six-gene core ---
    rows = []
    for gene, eid in CORE.items():
        o = old.loc[old["Gene_ID"] == eid]
        n = new.loc[new["Gene_ID"] == eid]
        if o.empty or n.empty:
            rows.append({"Gene": gene, "Entrez": eid, "present": False})
            continue
        o, n = o.iloc[0], n.iloc[0]
        rows.append(
            {
                "Gene": gene,
                "Entrez": eid,
                "present": True,
                "log2FC_22v21": float(o["log2FoldChange"]),
                "padj_22v21": float(o["padj"]),
                "log2FC_21v21": float(n["log2FoldChange"]),
                "padj_21v21": float(n["padj"]),
                "delta_log2FC": float(n["log2FoldChange"] - o["log2FoldChange"]),
                "sign_unchanged": bool(np.sign(o["log2FoldChange"]) == np.sign(n["log2FoldChange"])),
                "still_padj_lt_0.05": bool(n["padj"] < 0.05),
            }
        )
    core = pd.DataFrame(rows)
    core.to_csv(OUT / "core6_22v21_vs_21v21.tsv", sep="\t", index=False)

    # --- 3. shortlist movement (immune P frozen; DESeq2 replaced) ---
    immune = pd.read_csv(IMMUNE, sep="\t")
    sl_old = pd.read_csv(SHORTLIST, sep="\t")
    sl_old["Gene_ID"] = _gid(sl_old["Gene_ID"])
    sl_new = _rebuild_shortlist(new, immune)
    # attach symbols from old shortlist / immune is entrez only — use old shortlist + deseq top map
    sym_map = sl_old.dropna(subset=["symbol"]).drop_duplicates("Gene_ID").set_index("Gene_ID")["symbol"]
    sl_new["symbol"] = sl_new["Gene_ID"].map(sym_map)
    old_ids = set(sl_old["Gene_ID"])
    new_ids = set(sl_new["Gene_ID"])
    gained = sorted(new_ids - old_ids)
    lost = sorted(old_ids - new_ids)
    sl_move = pd.DataFrame(
        {
            "Gene_ID": gained + lost,
            "change": ["gained"] * len(gained) + ["lost"] * len(lost),
        }
    )
    sl_move["symbol"] = sl_move["Gene_ID"].map(sym_map)
    sl_new.to_csv(OUT / "shortlist_300_rebuilt_on_21v21_deseq2.tsv", sep="\t", index=False)
    sl_move.to_csv(OUT / "shortlist_membership_changes.tsv", sep="\t", index=False)

    # --- 2. eight-gene panel ---
    panel_pub = pd.read_csv(PANEL_TSV, sep="\t")
    g159 = pd.read_csv(BASE / "results/gse159610_full_de_results_annotated.tsv", sep="\t")
    g159 = g159.dropna(subset=["Gene_Symbol"]).copy()
    g159["Gene_Symbol"] = g159["Gene_Symbol"].astype(str)
    g159 = g159.sort_values("P_value").drop_duplicates("Gene_Symbol")
    g106 = pd.read_csv(
        BASE / "results/external_validation_gse106520/gse106520_gene_level_context.tsv",
        sep="\t",
    )
    # Prefer original shortlist blood stats when present (same estimator as published panel).
    sl_blood = sl_old.set_index("symbol") if "symbol" in sl_old.columns else None

    cand_rows = []
    for _, r in sl_new.iterrows():
        gene = r["symbol"]
        if not isinstance(gene, str) or not gene:
            continue
        disc_lfc = float(r["DESeq2_log2FC"])
        disc_sign = np.sign(disc_lfc)
        if disc_sign == 0:
            continue
        if sl_blood is not None and gene in sl_blood.index:
            b = sl_blood.loc[gene]
            if isinstance(b, pd.DataFrame):
                b = b.iloc[0]
            lfc159 = pd.to_numeric(b.get("GSE159610_Log2FC_aneurysm_vs_control"), errors="coerce")
            p159 = pd.to_numeric(b.get("GSE159610_P"), errors="coerce")
        else:
            hit = g159[g159["Gene_Symbol"] == gene]
            lfc159 = pd.to_numeric(hit["Log2FC"].iloc[0], errors="coerce") if len(hit) else np.nan
            p159 = pd.to_numeric(hit["P_value"].iloc[0], errors="coerce") if len(hit) else np.nan
        hit106 = g106[g106["Gene"] == gene]
        lfc106 = pd.to_numeric(hit106["Log2FC_AvC"].iloc[0], errors="coerce") if len(hit106) else np.nan
        p106 = pd.to_numeric(hit106["P_value"].iloc[0], errors="coerce") if len(hit106) else np.nan
        if pd.isna(lfc159) or pd.isna(lfc106) or pd.isna(p159) or pd.isna(p106):
            continue
        match159 = bool(np.sign(lfc159) == disc_sign)
        match106 = bool(np.sign(lfc106) == disc_sign)
        z, mp = _stouffer_two(p159, p106)
        cand_rows.append(
            {
                "Gene": gene,
                "Gene_ID": r["Gene_ID"],
                "DESeq2_log2FC": disc_lfc,
                "DESeq2_padj": float(r["DESeq2_padj"]),
                "match_GSE159610": match159,
                "match_GSE106520": match106,
                "both_concordant": match159 and match106,
                "Bridge_Stouffer_Z": z,
                "Bridge_meta_P": mp,
                "abs_discovery_log2FC": abs(disc_lfc),
                "in_published_panel": gene in PANEL,
            }
        )
    cand = pd.DataFrame(cand_rows)
    both = cand[cand["both_concordant"]].sort_values(
        ["Bridge_meta_P", "abs_discovery_log2FC", "Gene"], ascending=[True, False, True]
    ).reset_index(drop=True)
    both["new_panel_rank"] = np.arange(1, len(both) + 1)
    both.to_csv(OUT / "bridge_concordant_genes_on_21v21_shortlist.tsv", sep="\t", index=False)
    new_top8 = both.head(8)["Gene"].tolist()

    panel_status = []
    for gene in PANEL:
        in_new_sl = gene in set(sl_new["symbol"].dropna().astype(str))
        in_new_top8 = gene in new_top8
        o = sl_old.loc[sl_old["symbol"] == gene]
        n = sl_new.loc[sl_new["symbol"] == gene]
        sign_ok = True
        if len(o) and len(n):
            sign_ok = bool(np.sign(o.iloc[0]["DESeq2_log2FC"]) == np.sign(n.iloc[0]["DESeq2_log2FC"]))
        panel_status.append(
            {
                "Gene": gene,
                "in_new_shortlist": in_new_sl,
                "in_new_top8": in_new_top8,
                "discovery_sign_unchanged": sign_ok,
            }
        )
    panel_df = pd.DataFrame(panel_status)
    panel_df.to_csv(OUT / "panel8_status_on_21v21.tsv", sep="\t", index=False)

    core_ok = bool(core["sign_unchanged"].all() and core["still_padj_lt_0.05"].all())
    panel_ok = bool(
        panel_df["in_new_shortlist"].all()
        and panel_df["discovery_sign_unchanged"].all()
        and new_top8 == PANEL
    )
    k4_ok = all(g in set(sl_new["symbol"].dropna().astype(str)) for g in CORE_IN_SHORTLIST)

    decision = "switch_primary_to_21v21" if core_ok and panel_ok else "update_affected_tables_panel_or_core_moved"

    summary = {
        "counts_22v21": old_c,
        "counts_21v21": new_c,
        "core6": core.to_dict(orient="records"),
        "core_signs_and_padj_hold": core_ok,
        "shortlist": {
            "n_old": int(len(sl_old)),
            "n_new": int(len(sl_new)),
            "n_overlap": int(len(old_ids & new_ids)),
            "n_gained": int(len(gained)),
            "n_lost": int(len(lost)),
            "jaccard": float(len(old_ids & new_ids) / len(old_ids | new_ids)) if old_ids | new_ids else 0.0,
            "k4_still_in_shortlist": k4_ok,
            "gained_symbols": sl_move.loc[sl_move["change"] == "gained", "symbol"].fillna(
                sl_move.loc[sl_move["change"] == "gained", "Gene_ID"]
            ).tolist(),
            "lost_symbols": sl_move.loc[sl_move["change"] == "lost", "symbol"].fillna(
                sl_move.loc[sl_move["change"] == "lost", "Gene_ID"]
            ).tolist(),
        },
        "panel": {
            "published_top8": PANEL,
            "new_top8": new_top8,
            "identical_top8": new_top8 == PANEL,
            "all_eight_remain_eligible": bool(panel_df["in_new_shortlist"].all()),
            "all_eight_sign_unchanged": bool(panel_df["discovery_sign_unchanged"].all()),
            "panel_unchanged": panel_ok,
        },
        "decision": decision,
        "rationale": (
            "Primary phenotype is explicit GEO source_name (21 ruptured vs 21 unruptured); "
            "GSM3487892 unknown rupture status excluded. Immune-adjusted P values frozen. "
            "Panel derivation not rebuilt unless the eight-gene set changed under the same selection rule."
        ),
    }
    (OUT / "sensitivity_decision.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
