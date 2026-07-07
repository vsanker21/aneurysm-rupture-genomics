"""
Pipeline stage definitions for the rupture-genomics analysis workflow.

Each stage runs one or more existing analysis scripts as subprocesses.
Manuscript / Word document generation is intentionally excluded.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable


@dataclass(frozen=True)
class Stage:
    name: str
    description: str
    scripts: tuple[str, ...] = ()
    script_args: dict[str, list[str]] = field(default_factory=dict)
    requires: tuple[str, ...] = ()
    optional: bool = False
    skip_if_missing: tuple[str, ...] = ()


def _all_stages() -> list[Stage]:
    return [
        Stage(
            name="download",
            description="Download public GEO/GWAS/eQTL reference data",
            scripts=("download_phases_g_public_data.py",),
            optional=True,
        ),
        Stage(
            name="deseq2",
            description="DESeq2 ruptured vs unruptured on GSE122897 raw counts (R)",
            scripts=("run_deseq2_gse122897.py",),
            requires=(
                "data/sc_rnaseq/GSE122897/GSE122897_raw_counts_GRCh38.p13_NCBI.tsv.gz",
            ),
        ),
        Stage(
            name="phases_g",
            description="Immune-adjusted DE, eQTL/GWAS evidence, 300-gene shortlist",
            scripts=("run_phases_g_full_stack.py",),
            script_args={
                "run_phases_g_full_stack.py": ["--skip-download"],
            },
            requires=("results/gse122897_deseq2_results.tsv",),
        ),
        Stage(
            name="blood_de",
            description="GSE159610 whole-blood differential expression",
            scripts=("run_full_blood_de_analysis.py",),
            optional=True,
            skip_if_missing=("data/blood_rnaseq/GSE159610",),
        ),
        Stage(
            name="bridge",
            description="Link tissue DE to PHASES-G shortlist (tissue bridge)",
            scripts=("run_phases_g_bridge_refresh.py",),
            requires=("results/phases_g/shortlist_evidence.tsv",),
        ),
        Stage(
            name="blood_bridge",
            description="GSE106520 neutrophil bridge + 8-gene panel derivation",
            scripts=(
                "analyze_gse106520_blood_bridge.py",
                "blood_bridge_panel.py",
            ),
            requires=("results/phases_g/shortlist_evidence.tsv",),
        ),
        Stage(
            name="external",
            description="External cohort validation (runs when cohort data present)",
            scripts=(
                "analyze_gse54083_external_validation.py",
                "validate_gse13353_key_genes.py",
                "analyze_gse15629_external_context.py",
                "analyze_gse26969_external_context.py",
                "analyze_gse36791_external_context.py",
                "analyze_gse73378_blood_context.py",
            ),
            optional=True,
        ),
        Stage(
            name="mpp1_marco",
            description="MPP1/MARCO cross-cohort expression survey",
            scripts=("mpp1_marco_multi_cohort_expression_survey.py",),
            requires=("results/phases_g/shortlist_evidence.tsv",),
        ),
        Stage(
            name="validation_summary",
            description="Combined external-validation summary table and figure",
            scripts=("combined_validation_summary.py",),
            optional=True,
        ),
        Stage(
            name="supplementary",
            description="Supplementary tables S1–S12 and figures S1–S5 (TSV/PNG/PDF only)",
            scripts=("generate_supplementary_materials_complete.py",),
            script_args={
                "generate_supplementary_materials_complete.py": ["--skip-docx"],
            },
            requires=(
                "results/phases_g/shortlist_evidence.tsv",
                "results/blood_bridge_panel/blood_bridge_panel_gene_summary.tsv",
            ),
        ),
        Stage(
            name="shortlist_table",
            description="Complete 300-gene shortlist with HGNC symbols and bridge backfill",
            scripts=("complete_extended_table_3_annotations.py",),
            script_args={
                "complete_extended_table_3_annotations.py": ["--tsv-only"],
            },
            requires=("results/supplementary/tables/Supplementary_Table_S2_300gene_shortlist.tsv",),
        ),
        Stage(
            name="figures",
            description="Main manuscript figures 1–5 (PNG/PDF only)",
            scripts=("generate_main_five_figures_ultra_hd.py",),
            optional=True,
        ),
        Stage(
            name="verify",
            description="Verify PHASES-G objectives and required artifacts",
            scripts=("verify_phases_g_objectives.py",),
            optional=True,
        ),
    ]


STAGES: dict[str, Stage] = {s.name: s for s in _all_stages()}
STAGE_ORDER: tuple[str, ...] = tuple(s.name for s in _all_stages())


def resolve_stage_list(selection: str) -> list[str]:
    if selection == "all":
        return list(STAGE_ORDER)
    names = [x.strip() for x in selection.split(",") if x.strip()]
    unknown = [n for n in names if n not in STAGES]
    if unknown:
        raise ValueError(f"Unknown stage(s): {', '.join(unknown)}")
    return names


def check_requirements(root: Path, stage: Stage) -> list[str]:
    missing: list[str] = []
    for rel in stage.requires:
        if not (root / rel).exists():
            missing.append(rel)
    return missing


def should_skip_stage(root: Path, stage: Stage) -> str | None:
    for rel in stage.skip_if_missing:
        if not (root / rel).exists():
            return f"optional data not found: {rel}"
    return None
