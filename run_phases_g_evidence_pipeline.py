"""
PHASES-G evidence integration (hypothesis generation; not a clinical score).

Combines, with explicit limitations in the output JSON:
  (1) GSE122897 rupture DE + immune-fraction-adjusted linear model (existing outputs)
  (2) Optional replication bulk cohort GSE159610 (aneurysm vs control; different contrast)
  (3) Blood cis-eQTL evidence via EMBL-EBI eQTL Catalogue (GTEx/Lepik/TwinsUK blood GE)
  (4) Overlap of eQTL sentinel SNPs with Bakker 2020 IA GWAS (rsID lookup; not colocalization)
  (5) Coarse scRNA-style cell-type prioritization from curated marker list (EC/SMC/immune)

Prerequisites (run in order):
  python download_phases_g_public_data.py
  python prepare_cibersortx_mixture_gse122897.py   # if needed
  python immune_deconvolution_lm22_gse122897.py
  python deg_immune_adjusted_gse122897.py --max-genes 15000

Then:
  python run_phases_g_evidence_pipeline.py

Outputs:
  results/phases_g/shortlist_evidence.tsv
  results/phases_g/pipeline_summary.json

Optional tissue-cohort validation (after shortlist exists):
  python run_phases_g_bridge_refresh.py
    -> cohort DE export + tissue_signature_phases_g_bridge.tsv + *_summary.json
  python link_tissue_de_to_phases_g.py --de-tsv your_rupture_vs_unruptured.tsv
    -> bridge TSV + summary JSON alongside it

Completion check (all five objectives + public data + tier-1 artifacts):
  python verify_phases_g_objectives.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

try:
    from project_paths import get_project_root
    BASE = get_project_root()
except ImportError:
    BASE = Path(__file__).resolve().parent
IMMUNE_ADJ = BASE / "results/cibersortx/deg_immune_covariate_adjusted.tsv"
DESEQ = BASE / "results/gse122897_deseq2_results.tsv"
REPL = BASE / "results/gse159610_full_de_results.tsv"
EQTLGEN_CIS_FDR005_DEFAULT = (
    BASE
    / "data/external/eqtlgen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
)
GWAS_DEFAULT = (
    BASE / "data/gwas/bakker_2020/downloaded/IA.GWAS.BakkerMK.2020.sumstats.Stage_2.txt.gz"
)
EQTL_IDS = BASE / "data/external/eqtl_catalogue/blood_ge_dataset_ids.txt"
MARKERS = BASE / "data/reference/scrna_aneurysm_celltype_markers.tsv"
OUT_DIR = BASE / "results/phases_g"
OUT_TSV = OUT_DIR / "shortlist_evidence.tsv"
OUT_JSON = OUT_DIR / "pipeline_summary.json"
OUT_TSV_MOUSE = OUT_DIR / "shortlist_evidence_with_mouse_gse193533_stopgap.tsv"
MOUSE_EXPR = OUT_DIR / "gse193533_mouse_celltype_gene_expression.tsv"

EQTL_DATASET_ASSOCS = (
    "https://www.ebi.ac.uk/eqtl/api/v3/datasets/{dataset_id}/associations"
)
MYGENE_BATCH = "https://mygene.info/v3/gene"


def load_blood_dataset_ids() -> set[str]:
    if not EQTL_IDS.is_file():
        raise FileNotFoundError(
            f"{EQTL_IDS} missing. Run: python download_phases_g_public_data.py"
        )
    return {ln.strip() for ln in EQTL_IDS.read_text().splitlines() if ln.strip()}


def entrez_to_ensembl_symbol(
    entrez_ids: list[str], batch_size: int = 400
) -> pd.DataFrame:
    rows = []
    for i in range(0, len(entrez_ids), batch_size):
        batch = entrez_ids[i : i + batch_size]
        r = requests.post(
            MYGENE_BATCH,
            data={
                "ids": ",".join(batch),
                "fields": "symbol,ensembl.gene,genomic_pos,entrezgene",
                "species": "human",
            },
            timeout=120,
        )
        r.raise_for_status()
        data = r.json()
        if isinstance(data, dict):
            data = [data]
        for row in data:
            q = str(row.get("query", "")).strip()
            sym = row.get("symbol")
            eg = row.get("ensembl")
            ens = None
            if isinstance(eg, dict):
                ens = eg.get("gene")
            elif isinstance(eg, list) and eg:
                ens = eg[0].get("gene") if isinstance(eg[0], dict) else None
            if q and ens:
                rows.append({"Gene_ID": q, "symbol": sym, "ensembl_gene": ens})
        time.sleep(0.15)
    return pd.DataFrame(rows)


def _eqtl_get(url: str, params: dict, timeout: int, retries: int = 6) -> requests.Response:
    """Retries on transient Catalogue failures (500/502/503/504, timeouts)."""
    last_err: BaseException | None = None
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=timeout)
            r.raise_for_status()
            return r
        except (requests.Timeout, requests.ConnectionError) as e:
            last_err = e
            time.sleep(2.5 * (attempt + 1))
        except requests.HTTPError as e:
            code = e.response.status_code if e.response is not None else 0
            if code in (500, 502, 503, 504):
                last_err = e
                time.sleep(3.0 * (attempt + 1))
            else:
                raise
    raise last_err  # type: ignore[misc]


def fetch_blood_eqtl_best(
    ensembl_gene: str,
    blood_ids: set[str],
    max_pages_per_dataset: int = 4,
    request_timeout: int = 180,
) -> dict | None:
    """
    Minimum p-value association across indexed whole-blood GE datasets.

    Uses v3 /datasets/{dataset_id}/associations?gene_id=… per cohort. Global
    /v3/associations?gene_id=&dataset_id= returned rows from other datasets in
    testing. Result is not guaranteed to be the cohort-wide minimum-p cis hit
    unless all pages are read — use FTP/tabix for publication-grade leads.
    """
    best = None
    for ds_id in sorted(blood_ids):
        start = 0
        for _ in range(max(1, max_pages_per_dataset)):
            try:
                r = _eqtl_get(
                    EQTL_DATASET_ASSOCS.format(dataset_id=ds_id),
                    {
                        "gene_id": ensembl_gene,
                        "size": 2500,
                        "start": start,
                    },
                    timeout=request_timeout,
                )
            except requests.RequestException as e:
                print(
                    f"[warn] eQTL skip dataset={ds_id} gene={ensembl_gene}: {e}",
                    flush=True,
                )
                break
            chunk = r.json() or []
            for a in chunk:
                p = float(a.get("pvalue", 1.0))
                if best is None or p < best["pvalue"]:
                    best = {
                        "pvalue": p,
                        "beta": a.get("beta"),
                        "rsid": a.get("rsid"),
                        "variant": a.get("variant"),
                        "dataset_id": a.get("dataset_id"),
                        "study_id": a.get("study_id"),
                    }
            time.sleep(0.1)
            if len(chunk) < 2500:
                break
            start += 2500
    return best


def load_gwas_rows_for_rsids(path: Path, rsids: set[str]) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"GWAS file missing: {path}\nRun download_phases_g_public_data.py")
    need = {r for r in rsids if r and str(r).startswith("rs")}
    hits = []
    for chunk in pd.read_csv(path, sep="\t", compression="gzip", chunksize=400_000, low_memory=False):
        if "SNP" not in chunk.columns:
            raise ValueError(f"Unexpected GWAS columns: {chunk.columns.tolist()}")
        sub = chunk[chunk["SNP"].astype(str).isin(need)]
        if len(sub):
            hits.append(sub)
    if not hits:
        return pd.DataFrame()
    return pd.concat(hits, ignore_index=True).drop_duplicates(subset=["SNP"])


def eqtlgen_top_hit_per_symbol(eqtlgen_gz: Path, symbols_upper: set[str]) -> dict[str, dict]:
    """
    Stream the eQTLGen cis-eQTL file and keep the minimum-Pvalue hit per GeneSymbol
    for the provided symbol set (uppercase).
    """
    import gzip

    best: dict[str, dict] = {}
    with gzip.open(eqtlgen_gz, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(header)}
        need = ["Pvalue", "SNP", "Zscore", "AssessedAllele", "OtherAllele", "Gene", "GeneSymbol", "FDR"]
        miss = [c for c in need if c not in col]
        if miss:
            raise ValueError(f"eQTLGen header missing columns: {miss}")

        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                continue
            sym = parts[col["GeneSymbol"]].strip()
            sym_u = sym.upper()
            if sym_u not in symbols_upper:
                continue
            try:
                p = float(parts[col["Pvalue"]])
            except Exception:
                continue
            row = {
                "eQTLGen_top_SNP": parts[col["SNP"]].strip(),
                "eQTLGen_top_Pvalue": p,
                "eQTLGen_top_Zscore": float(parts[col["Zscore"]]),
                "eQTLGen_top_AssessedAllele": parts[col["AssessedAllele"]].strip(),
                "eQTLGen_top_OtherAllele": parts[col["OtherAllele"]].strip(),
                "eQTLGen_top_FDR": float(parts[col["FDR"]]),
                "eQTLGen_Gene": parts[col["Gene"]].strip(),
                "eQTLGen_GeneSymbol": sym,
            }
            prev = best.get(sym_u)
            if prev is None or p < float(prev["eQTLGen_top_Pvalue"]):
                best[sym_u] = row
    return best


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--max-genes", type=int, default=300, help="Cap eQTL API queries")
    ap.add_argument(
        "--gwas",
        type=Path,
        default=GWAS_DEFAULT,
        help="Bakker Stage_2 summary stats .txt.gz",
    )
    ap.add_argument(
        "--immune-adj",
        type=Path,
        default=IMMUNE_ADJ,
    )
    ap.add_argument("--deseq2-padj-max", type=float, default=0.05)
    ap.add_argument("--lm-p-max", type=float, default=0.05)
    ap.add_argument(
        "--eqtl-pages-per-dataset",
        type=int,
        default=4,
        help="Max 2500-row pages per blood GE dataset (per-dataset v3 associations URL).",
    )
    ap.add_argument(
        "--allow-missing-gwas",
        action="store_true",
        help="If Bakker Stage_2 .txt.gz is absent, skip IA rsID lookup (still writes eQTL/replication).",
    )
    ap.add_argument(
        "--eqtl-timeout",
        type=int,
        default=240,
        help="Seconds per eQTL Catalogue HTTP request (retries on timeout/connection errors).",
    )
    ap.add_argument(
        "--skip-eqtl",
        action="store_true",
        help="Skip eQTL Catalogue + IA rsID overlap (replication + markers still run).",
    )
    ap.add_argument(
        "--eqtlgen-cis-fdr005",
        type=Path,
        default=EQTLGEN_CIS_FDR005_DEFAULT,
        help="Optional eQTLGen cis-eQTL file (FDR<0.05) to add per-gene top blood cis-eQTL + GWAS overlap columns.",
    )
    ap.add_argument(
        "--append-mouse-stopgap",
        action="store_true",
        help="If mouse GSE193533 stopgap summary exists, write an additional TSV with mouse compartment columns appended.",
    )
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if not args.immune_adj.is_file():
        raise FileNotFoundError(
            f"{args.immune_adj} missing. Run deg_immune_adjusted_gse122897.py first."
        )

    adj = pd.read_csv(args.immune_adj, sep="\t")
    adj["Gene_ID"] = adj["Gene_ID"].astype(str).str.replace(".0", "", regex=False)
    adj["DESeq2_padj"] = pd.to_numeric(adj["DESeq2_padj"], errors="coerce")
    adj["LM_immune_adj_P_Ruptured"] = pd.to_numeric(
        adj["LM_immune_adj_P_Ruptured"], errors="coerce"
    )
    short = adj[
        (adj["DESeq2_padj"] <= args.deseq2_padj_max)
        & (adj["LM_immune_adj_P_Ruptured"] <= args.lm_p_max)
    ].copy()
    short = short.sort_values("DESeq2_padj").head(args.max_genes)
    if len(short) == 0:
        raise SystemExit(
            "No genes passed filters. Relax --deseq2-padj-max / --lm-p-max or extend immune_adj run."
        )

    blood_ids = load_blood_dataset_ids()
    idmap = entrez_to_ensembl_symbol(short["Gene_ID"].astype(str).tolist())
    m = short.merge(idmap, on="Gene_ID", how="left")

    markers = pd.read_csv(MARKERS, sep="\t")
    sym_set = {str(s).upper() for s in markers["gene_symbol"]}

    repl_df = None
    if REPL.is_file():
        repl_df = pd.read_csv(REPL, sep="\t")
        repl_df["ensembl_base"] = (
            repl_df["Gene_ID"].astype(str).str.split(".").str[0].str.strip()
        )

    eqtl_rows = []
    rsids = set()
    n_short = len(m)
    for idx, (_, row) in enumerate(m.iterrows()):
        ens = row.get("ensembl_gene")
        if not isinstance(ens, str) or not ens.startswith("ENSG"):
            eqtl_rows.append(
                {
                    "Gene_ID": row["Gene_ID"],
                    "blood_eqtl_top_rsid": np.nan,
                    "blood_eqtl_top_pvalue": np.nan,
                    "blood_eqtl_dataset_id": np.nan,
                }
            )
            continue
        if args.skip_eqtl:
            hit = None
        else:
            print(
                f"[eQTL] gene {idx + 1}/{n_short} Entrez={row['Gene_ID']} {ens}",
                flush=True,
            )
            hit = fetch_blood_eqtl_best(
                ens,
                blood_ids,
                max_pages_per_dataset=args.eqtl_pages_per_dataset,
                request_timeout=args.eqtl_timeout,
            )
        if hit and hit.get("rsid"):
            rsids.add(str(hit["rsid"]))
        eqtl_rows.append(
            {
                "Gene_ID": row["Gene_ID"],
                "blood_eqtl_top_rsid": hit["rsid"] if hit else np.nan,
                "blood_eqtl_top_pvalue": hit["pvalue"] if hit else np.nan,
                "blood_eqtl_beta": hit["beta"] if hit else np.nan,
                "blood_eqtl_dataset_id": hit["dataset_id"] if hit else np.nan,
            }
        )

    eqtl_df = pd.DataFrame(eqtl_rows)
    out = m.merge(eqtl_df, on="Gene_ID", how="left")

    gwas_hits = pd.DataFrame()
    if rsids:
        if args.gwas.is_file():
            gwas_hits = load_gwas_rows_for_rsids(args.gwas, rsids)
        elif args.allow_missing_gwas:
            print(
                f"[warn] GWAS file missing ({args.gwas}); IA overlap skipped.",
                flush=True,
            )
        else:
            raise FileNotFoundError(
                f"GWAS file missing: {args.gwas}\n"
                "Run: python download_phases_g_public_data.py\n"
                "Or pass --allow-missing-gwas to skip IA rsID lookup."
            )
    gmap = {}
    if len(gwas_hits) and "SNP" in gwas_hits.columns:
        for _, gr in gwas_hits.iterrows():
            gmap[str(gr["SNP"])] = {
                "IA_GWAS_P": gr.get("P", gr.get("PVAL", np.nan)),
                "IA_GWAS_BETA": gr.get("BETA", np.nan),
                "IA_GWAS_SE": gr.get("SE", np.nan),
            }

    out["IA_GWAS_P_at_eqtl_snp"] = out["blood_eqtl_top_rsid"].map(
        lambda x: gmap.get(str(x), {}).get("IA_GWAS_P", np.nan) if pd.notna(x) else np.nan
    )
    out["IA_GWAS_BETA_at_eqtl_snp"] = out["blood_eqtl_top_rsid"].map(
        lambda x: gmap.get(str(x), {}).get("IA_GWAS_BETA", np.nan) if pd.notna(x) else np.nan
    )

    # Replication (orthogonal contrast)
    repl_cols = []
    if repl_df is not None:
        for _, row in out.iterrows():
            ensb = row.get("ensembl_gene")
            base = str(ensb).split(".")[0] if pd.notna(ensb) else ""
            sub = repl_df[repl_df["ensembl_base"] == base]
            if len(sub) == 0:
                repl_cols.append(
                    {
                        "GSE159610_Log2FC_aneurysm_vs_control": np.nan,
                        "GSE159610_P": np.nan,
                        "direction_match_rupture_and_aneurysm": np.nan,
                    }
                )
            else:
                r0 = sub.iloc[0]
                l2 = float(r0.get("Log2FC", np.nan))
                pv = float(r0.get("P_value", np.nan))
                rdir = np.sign(float(row.get("DESeq2_log2FC", np.nan)))
                repl_cols.append(
                    {
                        "GSE159610_Log2FC_aneurysm_vs_control": l2,
                        "GSE159610_P": pv,
                        "direction_match_rupture_and_aneurysm": bool(np.sign(l2) == rdir)
                        if np.isfinite(l2) and np.isfinite(rdir) and rdir != 0
                        else np.nan,
                    }
                )
        out = pd.concat([out, pd.DataFrame(repl_cols)], axis=1)

    # Marker-based cell-type hint
    def prio(sym):
        if not isinstance(sym, str) or not sym:
            return "unclassified"
        u = sym.upper()
        hit = markers[markers["gene_symbol"].str.upper() == u]
        if len(hit):
            return hit.iloc[0]["cell_type_priority"]
        return "unclassified"

    out["marker_celltype_hint"] = out["symbol"].map(prio)
    out["is_known_lineage_marker"] = (
        out["symbol"].astype(str).str.upper().isin(sym_set)
    )

    # eQTLGen blood cis-eQTL layer (independent of eQTL Catalogue)
    eqtlgen_added = False
    n_eqtlgen_hit = 0
    n_eqtlgen_snp_in_gwas = 0
    if args.eqtlgen_cis_fdr005 and Path(args.eqtlgen_cis_fdr005).is_file():
        sym_u = set(out["symbol"].dropna().astype(str).str.upper().str.strip())
        best = eqtlgen_top_hit_per_symbol(Path(args.eqtlgen_cis_fdr005), sym_u)
        out["_SYM_U"] = out["symbol"].astype(str).str.upper().str.strip()
        best_df = pd.DataFrame([{"_SYM_U": k, **v} for k, v in best.items()])
        out = out.merge(best_df, on="_SYM_U", how="left").drop(columns=["_SYM_U"])
        eqtlgen_added = True
        n_eqtlgen_hit = int(pd.to_numeric(out["eQTLGen_top_Pvalue"], errors="coerce").notna().sum())

        # GWAS overlap at eQTLGen SNP (rsID intersection)
        rs2 = set(out["eQTLGen_top_SNP"].dropna().astype(str).unique().tolist())
        rs2 = {r for r in rs2 if r.startswith("rs")}
        gwas_hits2 = load_gwas_rows_for_rsids(args.gwas, rs2) if (rs2 and args.gwas.is_file()) else pd.DataFrame()
        gmap2 = {}
        if len(gwas_hits2) and "SNP" in gwas_hits2.columns:
            for _, gr in gwas_hits2.iterrows():
                gmap2[str(gr["SNP"])] = {
                    "IA_GWAS_P": gr.get("P", gr.get("PVAL", np.nan)),
                    "IA_GWAS_BETA": gr.get("BETA", np.nan),
                    "IA_GWAS_SE": gr.get("SE", np.nan),
                    "IA_GWAS_OR": gr.get("OR", np.nan),
                }
        out["IA_GWAS_P_at_eQTLGen_SNP"] = out["eQTLGen_top_SNP"].map(
            lambda x: gmap2.get(str(x), {}).get("IA_GWAS_P", np.nan) if pd.notna(x) else np.nan
        )
        out["IA_GWAS_BETA_at_eQTLGen_SNP"] = out["eQTLGen_top_SNP"].map(
            lambda x: gmap2.get(str(x), {}).get("IA_GWAS_BETA", np.nan) if pd.notna(x) else np.nan
        )
        out["IA_GWAS_SE_at_eQTLGen_SNP"] = out["eQTLGen_top_SNP"].map(
            lambda x: gmap2.get(str(x), {}).get("IA_GWAS_SE", np.nan) if pd.notna(x) else np.nan
        )
        out["IA_GWAS_OR_at_eQTLGen_SNP"] = out["eQTLGen_top_SNP"].map(
            lambda x: gmap2.get(str(x), {}).get("IA_GWAS_OR", np.nan) if pd.notna(x) else np.nan
        )
        n_eqtlgen_snp_in_gwas = int(pd.to_numeric(out["IA_GWAS_P_at_eQTLGen_SNP"], errors="coerce").notna().sum())

    out.to_csv(OUT_TSV, sep="\t", index=False)

    # Optional: append stopgap mouse scRNA compartment evidence if available
    mouse_appended = False
    if args.append_mouse_stopgap and MOUSE_EXPR.is_file():
        try:
            from append_mouse_scrna_stopgap_to_shortlist import build_mouse_gene_summary  # type: ignore

            mm = build_mouse_gene_summary(MOUSE_EXPR)
            out2 = out.merge(mm, on="symbol", how="left")
            out2.to_csv(OUT_TSV_MOUSE, sep="\t", index=False)
            mouse_appended = True
        except Exception as e:
            print(f"[warn] mouse stopgap append failed: {e}", flush=True)

    limitations = [
        "Blood eQTL source = eQTL Catalogue whole-blood GE datasets (not full eQTLGen 3.6GB file).",
        "Per gene and dataset, up to --eqtl-pages-per-dataset windows (2500 rows each) are scanned; true global minimum-p cis hit may lie beyond that cap.",
        "IA GWAS lookup is rsID intersection with eQTL sentinel SNP; not formal colocalization (no COLOC/SMR).",
        "GSE159610 contrast is aneurysm vs control, not ruptured vs unruptured; directional match is exploratory.",
        "Cell-type column uses a small curated marker panel; replace with aneurysm scRNA-derived gene sets for publication.",
    ]
    if args.skip_eqtl:
        limitations.insert(
            0,
            "eQTL Catalogue queries and IA GWAS rsID overlap were not run (--skip-eqtl).",
        )
    if eqtlgen_added:
        limitations.insert(
            0,
            "Added eQTLGen cis-eQTL (FDR<0.05) per-gene top hit + Bakker GWAS rsID overlap at that SNP (rsID intersection only; not colocalization).",
        )
    if rsids and not args.gwas.is_file() and args.allow_missing_gwas:
        limitations.append(
            "IA GWAS overlap was skipped because the Bakker summary-stats file was not on disk."
        )
    if args.append_mouse_stopgap:
        if mouse_appended:
            limitations.append(
                "Mouse scRNA stopgap columns were appended (GSE193533; supportive only; species/model differences; coarse EC/SMC/fibroblast/immune annotation)."
            )
        else:
            limitations.append(
                "Mouse scRNA stopgap append was requested but not written (missing summary TSV or append failed)."
            )

    summary = {
        "n_shortlist_genes": int(len(out)),
        "filters": {
            "DESeq2_padj_max": args.deseq2_padj_max,
            "LM_immune_adj_P_Ruptured_max": args.lm_p_max,
            "max_genes_capped": args.max_genes,
        },
        "gwas_file_used": str(args.gwas.resolve()) if args.gwas.is_file() else None,
        "ia_gwas_overlap_skipped": bool(
            args.skip_eqtl
            or (rsids and not args.gwas.is_file() and args.allow_missing_gwas)
        ),
        "eqtl_request_timeout_s": args.eqtl_timeout,
        "skip_eqtl": bool(args.skip_eqtl),
        "n_with_blood_eqtl_hit": int(out["blood_eqtl_top_pvalue"].notna().sum()),
        "n_blood_eqtl_with_rsid": int(out["blood_eqtl_top_rsid"].notna().sum()),
        "n_eqtl_snps_in_IA_GWAS_file": int(len(gmap)),
        "eqtlgen_added": bool(eqtlgen_added),
        "n_with_eqtlgen_cis_hit": int(n_eqtlgen_hit),
        "n_eqtlgen_snps_in_IA_GWAS_file": int(n_eqtlgen_snp_in_gwas),
        "n_replication_gse159610_matched": int(
            out["GSE159610_Log2FC_aneurysm_vs_control"].notna().sum()
        )
        if "GSE159610_Log2FC_aneurysm_vs_control" in out.columns
        else 0,
        "limitations": limitations,
        "mouse_stopgap_appended": bool(mouse_appended),
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"Wrote {OUT_TSV}")
    if mouse_appended:
        print(f"Wrote {OUT_TSV_MOUSE}")
    print(f"Wrote {OUT_JSON}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    try:
        main()
    except requests.RequestException as e:
        print(f"[error] {e}", file=sys.stderr)
        sys.exit(1)
