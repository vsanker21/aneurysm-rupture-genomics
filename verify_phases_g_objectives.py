"""
Verify that PHASES-G-style objectives (hypothesis-generation layer) are met.

Maps to run_phases_g_evidence_pipeline.py design:
  (1) GSE122897 rupture DE + immune-fraction-adjusted linear model
  (2) Replication bulk cohort GSE159610 (orthogonal contrast; file + merged columns)
  (3) Blood cis-eQTL via eQTL Catalogue (whole-blood GE index + per-dataset API)
  (4) Bakker 2020 IA GWAS rsID overlap at eQTL sentinel (not colocalization)
  (5) scRNA-style cell-type prioritization (curated marker panel)

Also checks: public data verification, tier-1 supplement artifacts, presence of germline PRS
scripts/weights/README (cohort DNA still required for real patient scores).

Achieved vs pending matrix: results/phases_g/PHASES_G_REQUIREMENTS_STATUS.md

Does NOT require: COLOC/SMR (explicit follow-up; see results/phases_g/stronger_genetics/).

Usage:
  python verify_phases_g_objectives.py
  python verify_phases_g_objectives.py --json-only   # print JSON to stdout only
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent
PY = sys.executable

REQUIRED_SHORTLIST_COLS = [
    "Gene_ID",
    "DESeq2_padj",
    "DESeq2_log2FC",
    "LM_immune_adj_P_Ruptured",
    "symbol",
    "ensembl_gene",
    "blood_eqtl_top_pvalue",
    "blood_eqtl_top_rsid",
    "IA_GWAS_P_at_eqtl_snp",
    "GSE159610_Log2FC_aneurysm_vs_control",
    "GSE159610_P",
    "direction_match_rupture_and_aneurysm",
    "marker_celltype_hint",
    "is_known_lineage_marker",
]


def main() -> int:
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--json-only", action="store_true")
    ap.add_argument(
        "--min-shortlist",
        type=int,
        default=50,
        help="Minimum genes in shortlist_evidence.tsv (default 50).",
    )
    ap.add_argument(
        "--min-eqtl-hits",
        type=int,
        default=10,
        help="Minimum genes with a blood eQTL p-value (default 10).",
    )
    args = ap.parse_args()

    report: dict = {"objectives": [], "errors": [], "warnings": []}

    def ok(name: str, detail: str) -> None:
        report["objectives"].append({"id": name, "status": "OK", "detail": detail})

    def fail(msg: str) -> None:
        report["errors"].append(msg)

    def warn(msg: str) -> None:
        report["warnings"].append(msg)

    # --- Public inputs (programmatic download + verify) ---
    r = subprocess.run(
        [PY, str(BASE / "download_phases_g_public_data.py"), "--verify-only"],
        cwd=BASE,
        capture_output=True,
        text=True,
    )
    if r.returncode != 0:
        fail(
            "Public data verification failed (run: python download_phases_g_public_data.py). "
            f"stderr: {r.stderr[:500]}"
        )
    else:
        ok(
            "Public data",
            "Bakker GWAS + eQTL blood GE index + GSE122897 TPM/raw + markers path; gzip readable",
        )

    # --- Prerequisite: DESeq2 + immune adjustment chain ---
    immune_adj = BASE / "results/cibersortx/deg_immune_covariate_adjusted.tsv"
    if not immune_adj.is_file():
        fail(f"Missing immune-adjusted DE: {immune_adj}")
    else:
        ia = pd.read_csv(immune_adj, sep="\t", nrows=5)
        n_full = sum(1 for _ in open(immune_adj, encoding="utf-8")) - 1
        need = {"Gene_ID", "DESeq2_padj", "LM_immune_adj_P_Ruptured"}
        if not need.issubset(set(ia.columns)):
            fail(f"immune_adj missing columns; need {need}")
        elif n_full < 100:
            fail(f"immune_adj too few rows ({n_full})")
        else:
            ok(
                "(1) Immune-adjusted rupture signal",
                f"deg_immune_covariate_adjusted.tsv n={n_full} genes (LM22 fractions + OLS)",
            )

    deseq = BASE / "results/gse122897_deseq2_results.tsv"
    if not deseq.is_file():
        warn(f"DESeq2 results not found at {deseq} (expected upstream of immune_adj)")
    else:
        ok("(1b) DESeq2 input", "gse122897_deseq2_results.tsv present")

    frac = BASE / "results/cibersortx/lm22_fractions_gse122897.tsv"
    mix = BASE / "results/cibersortx/GSE122897_CIBERSORTx_mixture_TPM.txt"
    if frac.is_file():
        ok("(1c) Deconvolution fractions", "lm22_fractions_gse122897.tsv present")
    elif mix.is_file():
        ok("(1c) Mixture matrix", "CIBERSORTx mixture present; LM22 optional")
    else:
        fail("Missing LM22 fractions and mixture matrix for deconvolution path")

    meta = BASE / "results/gse122897_sample_metadata_supplementary.tsv"
    if not meta.is_file():
        warn(f"Sample meta missing: {meta}")

    # --- Evidence pipeline outputs ---
    shortlist = BASE / "results/phases_g/shortlist_evidence.tsv"
    summary_path = BASE / "results/phases_g/pipeline_summary.json"
    if not shortlist.is_file():
        fail(f"Missing {shortlist}")
    else:
        sh = pd.read_csv(shortlist, sep="\t")
        miss = [c for c in REQUIRED_SHORTLIST_COLS if c not in sh.columns]
        if miss:
            fail(f"shortlist_evidence.tsv missing columns: {miss}")
        elif len(sh) < args.min_shortlist:
            fail(f"shortlist too small: {len(sh)} < {args.min_shortlist}")
        else:
            ok(
                "Shortlist table",
                f"shortlist_evidence.tsv rows={len(sh)}, all required columns present",
            )

    if not summary_path.is_file():
        fail(f"Missing {summary_path}")
    else:
        summ = json.loads(summary_path.read_text(encoding="utf-8"))
        # Genetics layer: accept either (A) eQTL Catalogue blood hits OR (B) eQTLGen cis-eQTL hits.
        # Catalogue can be slow/unreliable in some environments; eQTLGen provides a robust blood cis layer.
        skip_eqtl = bool(summ.get("skip_eqtl"))
        eqtlgen_added = bool(summ.get("eqtlgen_added"))
        n_eq = int(summ.get("n_with_blood_eqtl_hit", 0))
        n_eqtlgen = int(summ.get("n_with_eqtlgen_cis_hit", 0))

        if skip_eqtl and not eqtlgen_added:
            fail(
                "pipeline_summary.json has skip_eqtl=true and eqtlgen_added=false; no blood genetics layer was run. "
                "Run: python run_phases_g_evidence_pipeline.py --max-genes 300 "
                "(or provide the eQTLGen cis-eQTL file and re-run)."
            )

        # Minimum blood genetics evidence: prefer Catalogue count when available, else eQTLGen count.
        if not skip_eqtl:
            if n_eq < args.min_eqtl_hits:
                fail(
                    f"Too few blood eQTL Catalogue hits in summary ({n_eq} < {args.min_eqtl_hits}); "
                    "check eQTL Catalogue connectivity or filters."
                )
        elif eqtlgen_added:
            if n_eqtlgen < args.min_eqtl_hits:
                fail(
                    f"Too few eQTLGen cis-eQTL hits in summary ({n_eqtlgen} < {args.min_eqtl_hits}); "
                    "check eQTLGen file path / symbol mapping."
                )
        if not summ.get("limitations") or len(summ["limitations"]) < 3:
            fail("pipeline_summary.json limitations missing or too short")
        ok(
            "(3)+(4) Blood eQTL + IA GWAS rsID layer",
            f"skip_eqtl={skip_eqtl}; n_blood_eqtl_catalogue={n_eq}; n_eqtlgen={n_eqtlgen}; "
            f"n_eqtl_snps_in_IA_GWAS={summ.get('n_eqtl_snps_in_IA_GWAS_file')}; "
            f"gwas_file_used set",
        )
        n_rep = int(summ.get("n_replication_gse159610_matched", 0))
        repl_file = BASE / "results/gse159610_full_de_results.tsv"
        if not repl_file.is_file():
            warn(f"{repl_file} missing; replication counts may be zero")
        if n_rep < 1:
            warn("No GSE159610 replication matches in summary (optional but recommended)")
        else:
            ok(
                "(2) Orthogonal bulk replication",
                f"GSE159610 matched genes in shortlist: {n_rep} (aneurysm vs control contrast)",
            )

        if "marker_celltype_hint" in sh.columns:
            vc = sh["marker_celltype_hint"].astype(str).value_counts().to_dict()
            ok("(5) Cell-type priors", f"marker_celltype_hint distribution (top keys): {list(vc.items())[:6]}")

    # --- Marker reference ---
    markers = BASE / "data/reference/scrna_aneurysm_celltype_markers.tsv"
    if not markers.is_file():
        fail(f"Missing marker panel: {markers}")
    else:
        mk = pd.read_csv(markers, sep="\t")
        if "gene_symbol" not in mk.columns or "cell_type_priority" not in mk.columns:
            fail("Marker TSV missing gene_symbol or cell_type_priority")
        else:
            ok("(5b) Marker reference file", f"{len(mk)} curated markers")

    # --- eQTL dataset index ---
    eqtl_ids = BASE / "data/external/eqtl_catalogue/blood_ge_dataset_ids.txt"
    if not eqtl_ids.is_file():
        fail(f"Missing {eqtl_ids}")
    else:
        n_ids = len([ln for ln in eqtl_ids.read_text().splitlines() if ln.strip()])
        if n_ids < 1:
            fail("blood_ge_dataset_ids.txt empty")
        else:
            ok("(3b) Blood GE dataset index", f"{n_ids} QTD dataset(s) for Catalogue queries")

    # --- Tier-1 supplement (operationalized shortlist) ---
    tier1 = BASE / "results/phases_g/tier1_supplement.tsv"
    note = BASE / "results/phases_g/tier1_scoring_note.json"
    if tier1.is_file() and note.is_file():
        t1 = pd.read_csv(tier1, sep="\t")
        ok(
            "Tier-1 supplement",
            f"tier1_supplement.tsv ({len(t1)} rows) + tier1_scoring_note.json",
        )
    else:
        warn(
            "Tier-1 files missing; run: python phases_g_tier1_supplement.py --top-n 40"
        )

    # --- Germline PRS path (tooling only; cohort DNA still required for real scores) ---
    prs_w = BASE / "results/phases_g/prs_rsid_weights_from_panel.tsv"
    readme = BASE / "data/genotypes/README_GENOTYPE_INPUT.txt"
    val_py = BASE / "validate_prs_genotype_inputs.py"
    prs_py = BASE / "phases_g_prs_from_genotypes.py"
    ref_py = BASE / "prs_download_phase3_reference_and_score.py"
    need_paths = [prs_w, readme, val_py, prs_py]
    if all(p.is_file() for p in need_paths):
        hdr = pd.read_csv(prs_w, sep="\t", nrows=0)
        if {"rsid", "weight"}.issubset(hdr.columns):
            n_lines = sum(1 for _ in open(prs_w, encoding="utf-8")) - 1
            ok(
                "Germline PRS tooling",
                f"prs_rsid_weights_from_panel.tsv n={n_lines} SNPs; README; validate_prs_genotype_inputs.py; phases_g_prs_from_genotypes.py",
            )
        else:
            warn("prs_rsid_weights_from_panel.tsv missing rsid or weight column")
    else:
        miss = [p.name for p in need_paths if not p.is_file()]
        warn(f"Germline PRS tooling incomplete (missing): {miss}")
    if ref_py.is_file():
        ok(
            "(ref) 1000G Phase3 demo downloader",
            "prs_download_phase3_reference_and_score.py (not PHASES patients)",
        )
    else:
        warn("prs_download_phase3_reference_and_score.py missing")

    val_run = BASE / "phases_g_cohort_validation.py"
    if val_run.is_file():
        ok(
            "Cohort validation orchestrator",
            "phases_g_cohort_validation.py (sklearn ROC/calibration/DCA; optional plink2 clump)",
        )
    else:
        warn("phases_g_cohort_validation.py missing")

    report["summary"] = {
        "all_objectives_met": len(report["errors"]) == 0,
        "n_errors": len(report["errors"]),
        "n_warnings": len(report["warnings"]),
    }

    out_json = BASE / "results/phases_g/OBJECTIVES_VERIFICATION.json"
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(report, indent=2), encoding="utf-8")

    if args.json_only:
        print(json.dumps(report, indent=2))
        return 0 if not report["errors"] else 1

    print("PHASES-G objectives verification")
    print("=" * 60)
    for o in report["objectives"]:
        print(f"  [OK] {o['id']}: {o['detail']}")
    for w in report["warnings"]:
        print(f"  [WARN] {w}")
    for e in report["errors"]:
        print(f"  [FAIL] {e}")
    print("=" * 60)
    print(f"Wrote {out_json}")
    if report["errors"]:
        print("STATUS: INCOMPLETE (fix errors above)")
        return 1
    if report["warnings"]:
        print("STATUS: COMPLETE (with warnings)")
        return 0
    print("STATUS: COMPLETE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
