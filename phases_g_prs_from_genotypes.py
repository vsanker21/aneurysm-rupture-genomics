"""
Real section-5 PRS construction from **genotype/imputed dosage** × eQTL-aligned weights.

GEO RNA samples do NOT include genotype. This script defines the reproducible
path when you supply dosages aligned to PHASES-labelled individuals:

  PRS_raw_i = sum_j w_j * g_ij

where g_ij ∈ [0,2] counts the rupture-aligned **risk allele** for SNP j, and w_j is
truncated effect magnitude from eQTLGen (same spirit as PHASES-G shortlist docs).

USAGE (after you obtain dosages)
  python phases_g_prs_from_genotypes.py --dosage-long path/to/long.tsv \\
      --clinical path/to/phases_covariates.tsv --out results/phases_g/prs_merged.tsv

  python phases_g_prs_from_genotypes.py --dosage-wide path/to/wide.tsv

  sample_id in long/wide/clinical tables must match the identifier used for PHASES
  on each individual (study subject key)—not GEO GSM labels unless those are how
  you key your registry.

Real cohorts: Prefer PHASES keys in sample_id when exporting from VCF/PLINK/imputation.
If genotype sample names differ, build a two-column map from your pipeline’s ID →
PHASES key before format conversion, then pass --clinical with matching sample_ids.

WEIGHT TABLE ONLY (no genotype yet)
  python phases_g_prs_from_genotypes.py --build-weights-only

Use a frozen weights TSV (skip panel rebuild; no overwrite of --weights):
  python phases_g_prs_from_genotypes.py --from-weights-tsv --dosage-long your_long.tsv

Inputs for --build-weights-only:
  results/phases_g/eqtl_snp_candidate_panel.tsv (from phases_g_practical_analyses_5_8.py)
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent
PANEL_DEFAULT = BASE / "results/phases_g/eqtl_snp_candidate_panel.tsv"
OUT_WEIGHTS = BASE / "results/phases_g/prs_rsid_weights_from_panel.tsv"
OUT_DIR = BASE / "results/phases_g"


def _parse_scientific(s) -> float | None:
    if s is None or (isinstance(s, float) and not math.isfinite(s)):
        return None
    t = str(s).strip()
    if t == "" or t.lower() in ("nan", "none"):
        return None
    try:
        return float(t.replace("e", "E"))
    except ValueError:
        return None


def _tissue_rupture_sign(row: pd.Series) -> float:
    lm = _parse_scientific(row.get("LM_immune_adj_log2_coef_Ruptured"))
    if lm is not None and lm != 0:
        return 1.0 if lm > 0 else -1.0
    d2 = _parse_scientific(row.get("DESeq2_log2FC"))
    if d2 is not None and d2 != 0:
        return 1.0 if d2 > 0 else -1.0
    return 0.0


def _risk_allele(row: pd.Series) -> tuple[str | None, str | None] | tuple[None, None]:
    """
    eQTL Z is conventionally tied to assessed allele vs other allele in cis-eQTL summary.
    If tissue rupture-associated expression is upward, we align risk allele toward
    the allele that cis-eQTL associates with increased expression based on Z sign.
    """
    ass = row.get("eQTLGen_top_AssessedAllele")
    oth = row.get("eQTLGen_top_OtherAllele")
    if pd.isna(ass) or pd.isna(oth) or not str(ass).strip() or not str(oth).strip():
        return None, None
    assessed = str(ass).strip().upper()
    other = str(oth).strip().upper()
    z = _parse_scientific(row.get("eQTLGen_top_Zscore"))
    if z is None:
        return None, None
    rsid = str(row.get("eQTLGen_top_SNP", "")).strip()
    tis = _tissue_rupture_sign(row)
    if tis == 0:
        return None, None

    cis_up = z > 0

    if tis > 0:
        want_expr_up = True
    else:
        want_expr_up = False

    if want_expr_up:
        risk = assessed if cis_up else other
    else:
        risk = other if cis_up else assessed
    return risk.upper(), rsid


def build_rsid_weights(panel_tsv: Path) -> pd.DataFrame:
    df = pd.read_csv(panel_tsv, sep="\t", low_memory=False)
    need = [
        "eQTLGen_top_SNP",
        "eQTLGen_top_Zscore",
        "eQTLGen_top_AssessedAllele",
        "eQTLGen_top_OtherAllele",
    ]
    for c in need:
        if c not in df.columns:
            raise ValueError(f"Panel missing {c}; regenerate eqtl_snp_candidate_panel.tsv")

    rows = []
    for _, r in df.iterrows():
        rsid = str(r.get("eQTLGen_top_SNP", "")).strip()
        if not rsid.startswith("rs"):
            continue
        z = _parse_scientific(r.get("eQTLGen_top_Zscore"))
        if z is None:
            continue
        risk, _ = _risk_allele(r)
        if risk is None:
            continue
        ass = str(r["eQTLGen_top_AssessedAllele"]).strip().upper()
        oth = str(r["eQTLGen_top_OtherAllele"]).strip().upper()
        w = min(abs(float(z)), 10.0)
        rows.append(
            {
                "rsid": rsid,
                "weight": w,
                "risk_allele": risk,
                "assessed_allele": ass,
                "other_allele": oth,
                "eQTLGen_Z": float(z),
                "tissue_rupture_sign": _tissue_rupture_sign(r),
                "Gene_ID": str(r.get("Gene_ID", "")).replace(".0", ""),
                "symbol": r.get("symbol", ""),
            }
        )

    wdf = pd.DataFrame(rows)
    if wdf.empty:
        raise RuntimeError("No rsID rows with alignable alleles (check panel).")

    wdf = wdf.loc[wdf.groupby("rsid")["weight"].idxmax()].reset_index(drop=True)
    return wdf.sort_values("rsid")


def compute_prs_from_long_dosage(
    weights: pd.DataFrame, long_dosage: pd.DataFrame
) -> pd.DataFrame:
    """
    long_dosage columns: sample_id, rsid, dosage
    dosage = count (0–2) of **risk_allele** as defined in weights (after harmonization).
    """
    for c in ("sample_id", "rsid", "dosage"):
        if c not in long_dosage.columns:
            raise ValueError(f"long dosage must have column {c}")
    m = long_dosage.merge(weights[["rsid", "weight"]], on="rsid", how="inner")
    m["contrib"] = m["dosage"].astype(float) * m["weight"].astype(float)
    out = m.groupby("sample_id", as_index=False).agg(
        PRS_raw=("contrib", "sum"),
        n_snps_used=("rsid", "count"),
    )
    mu, sd = out["PRS_raw"].mean(), out["PRS_raw"].std()
    out["PRS_z"] = (out["PRS_raw"] - mu) / sd if sd and sd > 1e-9 else 0.0
    return out


def compute_prs_from_wide_dosage(weights: pd.DataFrame, wide: pd.DataFrame) -> pd.DataFrame:
    if "sample_id" not in wide.columns:
        raise ValueError("wide matrix must have sample_id as first column or named sample_id")
    rsids = [c for c in wide.columns if c != "sample_id"]
    wsub = weights[weights["rsid"].isin(rsids)].set_index("rsid")
    rows = []
    for _, r in wide.iterrows():
        sid = r["sample_id"]
        s = 0.0
        n = 0
        for rs in wsub.index:
            if rs not in wide.columns:
                continue
            g = r[rs]
            if pd.isna(g):
                continue
            s += float(g) * float(wsub.loc[rs, "weight"])
            n += 1
        rows.append({"sample_id": sid, "PRS_raw": s, "n_snps_used": n})
    out = pd.DataFrame(rows)
    mu, sd = out["PRS_raw"].mean(), out["PRS_raw"].std()
    out["PRS_z"] = (out["PRS_raw"] - mu) / sd if sd and sd > 1e-9 else 0.0
    return out


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", type=Path, default=PANEL_DEFAULT)
    ap.add_argument("--build-weights-only", action="store_true")
    ap.add_argument(
        "--from-weights-tsv",
        action="store_true",
        help="Do not rebuild weights from --panel; read existing --weights TSV for PRS only",
    )
    ap.add_argument("--dosage-long", type=Path, default=None)
    ap.add_argument("--dosage-wide", type=Path, default=None)
    ap.add_argument("--weights", type=Path, default=OUT_WEIGHTS)
    ap.add_argument(
        "--clinical",
        type=Path,
        default=None,
        help="TSV merged on sample_id (must match dosage sample_id = PHASES subject key)",
    )
    ap.add_argument("--out", type=Path, default=OUT_DIR / "prs_sample_scores.tsv")
    args = ap.parse_args()

    if args.from_weights_tsv and args.build_weights_only:
        raise SystemExit("Cannot combine --from-weights-tsv with --build-weights-only")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if args.from_weights_tsv:
        if not args.weights.is_file():
            raise SystemExit(f"--from-weights-tsv requires existing file: {args.weights}")
        wread = pd.read_csv(args.weights, sep="\t")
        if "rsid" not in wread.columns or "weight" not in wread.columns:
            raise SystemExit("--weights TSV must contain rsid and weight columns")
    else:
        wdf = build_rsid_weights(args.panel)
        args.weights.parent.mkdir(parents=True, exist_ok=True)
        wdf.to_csv(args.weights, sep="\t", index=False)

        meta = {
            "n_snps_in_weights": int(len(wdf)),
            "weights_tsv": str(args.weights.resolve()),
            "note": (
                "Allele alignment uses eQTLGen assessed/other and Z sign vs tissue rupture "
                "expression sign; independent colocalization still required for causal claims."
            ),
        }
        (OUT_DIR / "prs_weight_build_meta.json").write_text(
            json.dumps(meta, indent=2), encoding="utf-8"
        )

        if args.build_weights_only:
            print(f"Wrote {args.weights}")
            print(json.dumps(meta, indent=2))
            return

        if args.dosage_long is None and args.dosage_wide is None:
            print(
                "Weights written. No --dosage-long / --dosage-wide provided.\n"
                "See data/genotypes/README_GENOTYPE_INPUT.txt for required formats."
            )
            return

        wread = wdf

    if args.dosage_long is None and args.dosage_wide is None:
        if args.from_weights_tsv:
            raise SystemExit("With --from-weights-tsv, provide --dosage-long or --dosage-wide")
        return

    if args.dosage_long is not None:
        d = pd.read_csv(args.dosage_long, sep="\t")
        out = compute_prs_from_long_dosage(wread, d)
    else:
        d = pd.read_csv(args.dosage_wide, sep="\t")
        out = compute_prs_from_wide_dosage(wread, d)

    if args.clinical is not None and args.clinical.is_file():
        clin = pd.read_csv(args.clinical, sep="\t")
        out = out.merge(clin, on="sample_id", how="left")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
