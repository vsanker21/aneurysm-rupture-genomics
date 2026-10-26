"""Assemble a paper-only code archive for GitHub v1.1.0 / Zenodo software DOI."""
from __future__ import annotations

import shutil
import zipfile
from pathlib import Path

BASE = Path(r"E:\Aneurysms Genomics Project")
OUT = BASE / "paper_code_release"
ZIP_PATH = BASE / "aneurysm-rupture-genomics-v1.1.0.zip"

FILES = [
    ".gitignore",
    "README.md",
    "requirements.txt",
    "project_paths.py",
    "run_pipeline.py",
    "pipeline/__init__.py",
    "pipeline/config.py",
    "pipeline/runner.py",
    "pipeline/stages.py",
    "run_deseq2_gse122897.R",
    "run_deseq2_gse122897.py",
    "run_deseq2_gse122897_21v21_explicit.R",
    "metadata/gse122897_sample_metadata_explicit_21v21.tsv",
    "metadata/PHENOTYPE_DECISION.txt",
    "metadata/DESeq2_count_confirmation.txt",
    "download_phases_g_public_data.py",
    "prepare_cibersortx_mixture_gse122897.py",
    "rewrite_cibersortx_mixture_to_gsm.py",
    "verify_cibersortx_mixture_gse122897.py",
    "immune_deconvolution_lm22_gse122897.py",
    "deg_immune_adjusted_gse122897.py",
    "parse_cibersortx_output_gse122897.py",
    "run_phases_g_full_stack.py",
    "run_phases_g_evidence_pipeline.py",
    "run_phases_g_bridge_refresh.py",
    "export_gse122897_cohort_de_for_bridge.py",
    "link_tissue_de_to_phases_g.py",
    "append_mouse_scrna_stopgap_to_shortlist.py",
    "verify_phases_g_objectives.py",
    "phases_g_tier1_supplement.py",
    "run_full_blood_de_analysis.py",
    "analyze_gse106520_blood_bridge.py",
    "blood_bridge_panel.py",
    "blood_bridge_panel_gene_attribution.py",
    "public_cohort_context_utils.py",
    "analyze_gse54083_external_validation.py",
    "validate_gse13353_key_genes.py",
    "analyze_gse15629_external_context.py",
    "analyze_gse26969_external_context.py",
    "analyze_gse36791_external_context.py",
    "analyze_gse73378_blood_context.py",
    "mpp1_marco_multi_cohort_expression_survey.py",
    "combined_validation_summary.py",
    "rupture_consensus_meta_analysis.py",
    "rerun_stouffer_meta_21v21.py",
    "compare_deseq2_21v21_sensitivity.py",
    "pathway_enrichment_gse122897_de.py",
    "generate_supplementary_materials_complete.py",
    "complete_extended_table_3_annotations.py",
    "generate_main_five_figures_ultra_hd.py",
    "manuscript_revised_vs_metrics.py",
    "manuscript_results_detail.py",
    "CIBERSORTx_Job7_Results.txt",
]


def main() -> None:
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)

    copied = []
    missing = []
    for rel in FILES:
        src = BASE / rel
        if not src.is_file():
            missing.append(rel)
            continue
        dest = OUT / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)
        copied.append(rel)

    (OUT / "RELEASE_NOTES.txt").write_text(
        "Aneurysm rupture genomics analysis code v1.1.0\n"
        "==============================================\n\n"
        "This archive contains the analysis scripts used in the paper.\n"
        "Primary GSE122897 contrast: 21 ruptured vs 21 unruptured, assigned from\n"
        "exact GEO source_name labels. GSM3487892 is excluded. Phenotype assignment\n"
        "does not use grepl(\"Ruptured|ruptured\", ...).\n\n"
        "Raw GEO files are not included. Place them under data/ as described in README.md.\n"
        "Companion processed tables: https://doi.org/10.5281/zenodo.21249929\n",
        encoding="utf-8",
    )
    copied.append("RELEASE_NOTES.txt")

    if ZIP_PATH.exists():
        ZIP_PATH.unlink()
    with zipfile.ZipFile(ZIP_PATH, "w", zipfile.ZIP_DEFLATED) as zf:
        for f in OUT.rglob("*"):
            if f.is_file():
                zf.write(f, Path("aneurysm-rupture-genomics-v1.1.0") / f.relative_to(OUT))

    print("copied", len(copied))
    if missing:
        print("missing", missing)
    print("zip", ZIP_PATH, ZIP_PATH.stat().st_size)


if __name__ == "__main__":
    main()
