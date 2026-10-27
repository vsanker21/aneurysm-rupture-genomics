# Aneurysm Rupture Genomics — Analysis Pipeline

Reproducible end-to-end analysis for intracranial aneurysm rupture transcriptomics:
tissue discovery (GSE122897), immune-adjusted DE, 300-gene shortlist, blood-bridge
panel derivation, external cohort validation, and supplementary tables/figures.

**This repository is analysis-only.** 

## Requirements

- Python 3.10+
- R 4.x with **DESeq2** (Bioconductor) for GSE122897 count-based DE
- Public GEO access (for download stage)

```bash
pip install -r requirements.txt
```

## Quick start

Clone the repository, place or download data under `data/`, then:

```bash
# List pipeline stages
python run_pipeline.py --list-stages

# Full analysis (skip download if data already present)
python run_pipeline.py --stage all

# Core analysis only (after DESeq2 results exist)
python run_pipeline.py --stage phases_g,bridge,blood_bridge,supplementary,shortlist_table

# Dry run (print commands)
python run_pipeline.py --stage all --dry-run
```

Set a custom root if needed:

```bash
export ANEURYSM_GENOMICS_ROOT=/path/to/repo
python run_pipeline.py --stage all
```

## Pipeline stages

| Stage | Description |
|-------|-------------|
| `download` | Public GEO/GWAS/eQTL reference data |
| `deseq2` | DESeq2 on GSE122897 (R) |
| `phases_g` | Immune-adjusted DE, eQTL/GWAS evidence, shortlist |
| `blood_de` | GSE159610 blood DE (optional) |
| `bridge` | Tissue DE ↔ shortlist bridge |
| `blood_bridge` | GSE106520 bridge + 8-gene panel |
| `external` | External cohort validation (optional) |
| `mpp1_marco` | MPP1/MARCO cross-cohort survey |
| `validation_summary` | Combined validation table/figure |
| `supplementary` | Tables S1–S12, figures S1–S5 (TSV/PNG only) |
| `shortlist_table` | Completed 300-gene shortlist with HGNC symbols |
| `figures` | Main figures 1–5 (PNG/PDF only) |
| `verify` | PHASES-G objectives check |

## Key outputs

```
results/
  gse122897_deseq2_results.tsv          # Tissue DESeq2
  phases_g/shortlist_evidence.tsv       # 300-gene evidence table
  blood_bridge_panel/                   # 8-gene panel performance
  supplementary/tables/                 # Supplementary_Table_S*.tsv
  supplementary/figures/                # SuppFig*.png
  figures/main_manuscript/              # Figure_1..5 (optional stage)
```

## Data layout

Expected inputs (see `download` stage or README files under `data/`):

```
data/
  sc_rnaseq/GSE122897/          # Raw counts + metadata
  blood_rnaseq/GSE159610/      # Whole-blood RNA-seq
  gwas/bakker_2020/downloaded/  # IA GWAS sumstats
  external/eqtl_catalogue/     # Blood eQTL index
```

External validation cohorts under `data/microarray/` and `External Validation/`
are used when present; the `external` stage skips missing cohorts gracefully.

## What is NOT in this pipeline

- Word manuscript assembly (`create_document_with_inferences.py`, etc.)
- Submission docx patching (`Tables_Revised.docx`, `Latest Submission Files/`)
- Legacy multi-omics / PRS tracks (`run_complete_project_pipeline.py`)

Run those separately if needed for publication packaging.

## GSE122897 phenotype (primary 21 vs 21)

Primary discovery contrast is **21 ruptured vs 21 unruptured**. Labels are the
exact GEO `source_name` strings `Ruptured intracranial aneurysm` and
`Unruptured intracranial aneurysm`. GSM3487892 (Unknown rupture status) is
excluded. Cortical artery samples and GSM3487858 are excluded. **18 of 60**
series specimens are excluded; 42 are analysed. Phenotype assignment does
**not** use `grepl("Ruptured|ruptured", ...)`.

The sample table is `metadata/gse122897_sample_metadata_explicit_21v21.tsv`.
A frozen 22 vs 21 as-submitted run (GSM3487892 coded ruptured by substring
match) is retained only as a sensitivity archive.

## Citation / methods summary

Discovery: GSE122897 wall RNA-seq (21 ruptured vs 21 unruptured), DESeq2 +
LM22-adjusted OLS shortlist (padj ≤ 0.05 and immune-adjusted P ≤ 0.05).
Blood bridge: directional concordance in GSE159610 and GSE106520, Stouffer
ranking. Panel: top 8 concordant genes; Phases 2–4 held out in performance
reporting.

Three-cohort Stouffer meta-analysis maps GSE122897 Entrez IDs with
`metadata/hgnc_entrez_to_symbol.tsv`, joins GSE54083 and GSE13353 on HGNC
symbols, and applies BH-FDR across the full overlap. The six-gene wall core
is genes concordant in all three cohorts with meta FDR < 1e-6: TGFBI, ANPEP,
MARCO, MPP1, CLEC5A, CTSL.

Code: GitHub release tag `v1.1.1`; software DOI
https://doi.org/10.5281/zenodo.23700289. Processed tables:
https://doi.org/10.5281/zenodo.21249929.

## License

Add your license here before public GitHub release.

## Author

**Vivek Sanker** — [github.com/vsanker21](https://github.com/vsanker21)

**Live repository:** [github.com/vsanker21/aneurysm-rupture-genomics](https://github.com/vsanker21/aneurysm-rupture-genomics)
