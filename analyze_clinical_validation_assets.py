"""
Derive structured downstream analysis artifacts from downloaded clinical cohort assets.

This script turns public protocol / registry / study-site materials into:

- cohort-level variable coverage matrix,
- evidence snippets supporting each inferred capability,
- a minimum external validation data dictionary,
- and a public-only readiness summary.

It does not fabricate access to participant-level data. Instead, it clarifies
what can already be supported by the downloaded public materials and what still
requires restricted-access cohort transfer.
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path

import pandas as pd
from bs4 import BeautifulSoup, XMLParsedAsHTMLWarning
from pypdf import PdfReader


ROOT = Path(__file__).resolve().parent
SCOUT_DIR = ROOT / "results" / "clinical_cohort_scout"
MANIFEST_PATH = SCOUT_DIR / "clinical_cohort_manifest.tsv"
ASSET_AUDIT_PATH = SCOUT_DIR / "clinical_cohort_asset_audit.tsv"

VARIABLE_MATRIX_PATH = SCOUT_DIR / "clinical_validation_variable_matrix.tsv"
SNIPPETS_PATH = SCOUT_DIR / "clinical_validation_evidence_snippets.tsv"
DATA_DICTIONARY_PATH = SCOUT_DIR / "clinical_validation_minimum_data_dictionary.tsv"
SUMMARY_PATH = SCOUT_DIR / "clinical_validation_analysis_summary.md"
LINKED_RESOURCES_PATH = SCOUT_DIR / "clinical_validation_linked_resources.tsv"


FEATURE_PATTERNS: dict[str, list[str]] = {
    "blood_biospecimen": [
        r"blood sample",
        r"blood samples",
        r"whole blood",
        r"biological specimen bank",
        r"specimen bank",
        r"serum",
        r"plasma",
        r"paxgene",
    ],
    "biomarker_plan": [
        r"biomarker",
        r"biomarkers",
        r"blood sample",
        r"specimen bank",
    ],
    "dna_resource": [
        r"\bdna\b",
        r"whole-genome sequencing",
        r"whole genome sequencing",
        r"genotyping",
        r"cryopreserved",
    ],
    "rna_or_expression_resource": [
        r"\brna\b",
        r"rna-seq",
        r"whole blood rna",
        r"transcript",
    ],
    "prospective_followup": [
        r"prospective",
        r"follow-up",
        r"follow up",
        r"longitudinal",
        r"natural history",
    ],
    "rupture_endpoint": [
        r"primary outcome.*rupture",
        r"rupture resulting in hospital admission or death",
        r"future risk",
        r"rupture risk",
        r"subsequently ruptured",
    ],
    "growth_endpoint": [
        r"aneurysm growth",
        r"uia growth",
        r"risk of ia evolution",
        r"growth will be determined",
    ],
    "phases_named": [
        r"\bphases\b",
        r"population,\s*hypertension,\s*age,\s*size,\s*earlier",
    ],
    "age_named": [r"\bage\b", r"baseline demographics"],
    "hypertension_named": [r"hypertension", r"blood pressure"],
    "size_named": [r"aneurysm size", r"\bsize\b", r"7 mm", r"3 and 7 mm"],
    "site_named": [r"\bsite\b", r"location", r"posterior circulation", r"anterior"],
    "earlier_sah_named": [r"earlier subarachnoid hemorrhage", r"earlier subarachnoid haemorrhage", r"history of subarachnoid"],
    "sex_named": [r"\bsex\b", r"\bgender\b", r"\bmale\b", r"\bfemale\b"],
    "smoking_named": [r"smoking", r"smoking cessation"],
    "treatment_tracking": [r"treatment", r"repair", r"endovascular", r"surgery", r"occlusion treatments"],
    "imaging_tracking": [r"imaging", r"mri", r"hospital admissions", r"aneurysm characteristics", r"hemodynamic", r"morphological"],
}


MINIMUM_DATA_DICTIONARY = [
    ("patient_id", "string", "required", "Stable de-identified participant identifier for joining clinical and assay data."),
    ("aneurysm_id", "string", "required", "Unique aneurysm-level identifier when patients contribute multiple aneurysms."),
    ("age_years", "numeric", "required", "Age at baseline specimen collection or baseline imaging."),
    ("sex", "categorical", "required", "Biological sex or sex assigned at birth used in clinical records."),
    ("hypertension", "binary", "required", "Hypertension status at baseline."),
    ("smoking_status", "categorical", "recommended", "Never / former / current smoking, preferably at baseline."),
    ("prior_sah", "binary", "required", "History of earlier subarachnoid haemorrhage before baseline."),
    ("aneurysm_size_mm", "numeric", "required", "Baseline aneurysm maximum diameter in millimetres."),
    ("aneurysm_site", "categorical", "required", "PHASES-compatible site / location category."),
    ("aneurysm_multiplicity", "binary", "recommended", "Multiple aneurysms present at baseline."),
    ("specimen_type", "categorical", "required", "Whole blood, PAXgene RNA, plasma, serum, DNA, or other specimen."),
    ("specimen_collection_date", "date", "required", "Collection date for baseline specimen."),
    ("specimen_processing_protocol", "string", "recommended", "Extraction, stabilization, and processing protocol."),
    ("storage_temperature", "string", "recommended", "Storage conditions and freeze-thaw constraints."),
    ("rupture_outcome", "binary", "required", "Observed future aneurysm rupture after baseline."),
    ("rupture_date", "date", "required_if_event", "Date of rupture event if observed."),
    ("censor_date", "date", "required", "Administrative censoring or last follow-up date."),
    ("treatment_before_event", "binary", "required", "Whether treatment occurred before rupture or censoring."),
    ("treatment_date", "date", "required_if_treated", "Date of intervention for treated aneurysms."),
    ("followup_imaging_date", "date", "recommended", "Date of follow-up imaging if growth analyses are performed."),
    ("aneurysm_growth", "binary", "recommended", "Growth endpoint for secondary analyses."),
]


def _normalize_text(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def _extract_text(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix == ".pdf":
        reader = PdfReader(str(path))
        pages = [(page.extract_text() or "") for page in reader.pages]
        return _normalize_text(" ".join(pages))
    if suffix in {".html", ".htm"}:
        html = path.read_text(encoding="utf-8", errors="ignore")
        warnings.filterwarnings("ignore", category=XMLParsedAsHTMLWarning)
        soup = BeautifulSoup(html, "lxml")
        return _normalize_text(soup.get_text(" ", strip=True))
    return _normalize_text(path.read_text(encoding="utf-8", errors="ignore"))


def _first_match_snippet(text: str, patterns: list[str], window: int = 160) -> str:
    for pattern in patterns:
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match:
            start = max(0, match.start() - window)
            end = min(len(text), match.end() + window)
            snippet = text[start:end]
            return snippet.strip()
    return ""


def _feature_status(text: str, patterns: list[str]) -> tuple[int, str]:
    snippet = _first_match_snippet(text, patterns)
    return int(bool(snippet)), snippet


def _read_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    manifest = pd.read_csv(MANIFEST_PATH, sep="\t")
    assets = pd.read_csv(ASSET_AUDIT_PATH, sep="\t")
    return manifest, assets


def _build_text_map(assets: pd.DataFrame) -> dict[str, str]:
    text_map: dict[str, str] = {}
    downloadable = assets[assets["status"] == "downloaded"].copy()
    for path_str in downloadable["local_path"].dropna():
        path = ROOT / str(path_str)
        if path.is_file():
            text_map[str(path_str)] = _extract_text(path)
    return text_map


def _build_variable_matrix(manifest: pd.DataFrame, assets: pd.DataFrame, text_map: dict[str, str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    snippet_rows: list[dict[str, object]] = []

    for _, cohort in manifest.iterrows():
        cohort_id = str(cohort["cohort_id"])
        cohort_assets = assets[(assets["cohort_id"] == cohort_id) & (assets["status"] == "downloaded")].copy()
        texts = [text_map.get(str(path), "") for path in cohort_assets["local_path"].dropna()]
        combined_text = _normalize_text(" ".join(texts))

        row: dict[str, object] = {
            "cohort_id": cohort_id,
            "cohort_name": cohort["cohort_name"],
            "public_data_access_level": cohort["public_data_access_level"],
            "baseline_blood_status": cohort["baseline_blood_status"],
            "rupture_followup_status": cohort["rupture_followup_status"],
            "phases_ready_status": cohort["phases_ready_status"],
            "downloadable_public_assets": int(cohort["downloadable_public_assets"]),
            "public_patient_level_data_available": 0,
        }

        for feature, patterns in FEATURE_PATTERNS.items():
            status, snippet = _feature_status(combined_text, patterns)
            row[feature] = status
            if snippet:
                snippet_rows.append(
                    {
                        "cohort_id": cohort_id,
                        "feature": feature,
                        "snippet": snippet,
                    }
                )

        phases_core_count = int(
            row["age_named"]
            + row["hypertension_named"]
            + row["size_named"]
            + row["site_named"]
            + row["earlier_sah_named"]
        )
        row["phases_core_variable_count"] = phases_core_count
        evidence_score = int(
            10 * row["blood_biospecimen"]
            + 8 * row["biomarker_plan"]
            + 8 * row["prospective_followup"]
            + 10 * row["rupture_endpoint"]
            + 4 * row["growth_endpoint"]
            + 6 * row["phases_named"]
            + 3 * phases_core_count
            + 4 * row["treatment_tracking"]
            + 4 * row["imaging_tracking"]
            + 4 * row["smoking_named"]
            + 3 * row["sex_named"]
            + 4 * row["dna_resource"]
            + 4 * row["rna_or_expression_resource"]
        )
        design_priority_score = evidence_score
        blood_status = str(cohort["baseline_blood_status"])
        rupture_status = str(cohort["rupture_followup_status"])
        phases_status = str(cohort["phases_ready_status"])

        if blood_status == "yes":
            design_priority_score += 12
        elif blood_status == "dna_only":
            design_priority_score += 4

        if rupture_status == "yes":
            design_priority_score += 15
        elif rupture_status == "growth_only":
            design_priority_score -= 8

        if phases_status == "yes":
            design_priority_score += 10
        elif phases_status == "partial":
            design_priority_score += 6

        if str(cohort["public_data_access_level"]) == "controlled_access":
            design_priority_score -= 6

        row["public_evidence_score"] = evidence_score
        row["design_priority_score"] = int(design_priority_score)

        if cohort_id in {"CIAP-1", "CIAP-2"}:
            row["recommended_analysis_role"] = "priority_future_external_validation"
        elif cohort_id == "ROAR":
            row["recommended_analysis_role"] = "priority_clinical_outcome_validation_needs_biospecimen_confirmation"
        elif cohort_id == "UCAN":
            row["recommended_analysis_role"] = "secondary_growth_instability_validation"
        elif cohort_id == "phs003072":
            row["recommended_analysis_role"] = "supportive_blood_transcriptome_comparator"
        elif cohort_id == "ROAR-DNA":
            row["recommended_analysis_role"] = "genetic_extension_only"
        elif cohort_id == "ISUIA":
            row["recommended_analysis_role"] = "legacy_natural_history_reference"
        else:
            row["recommended_analysis_role"] = "background_context"

        rows.append(row)

    matrix = (
        pd.DataFrame(rows)
        .sort_values(["design_priority_score", "public_evidence_score", "cohort_id"], ascending=[False, False, True])
        .reset_index(drop=True)
    )
    snippets = pd.DataFrame(snippet_rows).sort_values(["cohort_id", "feature"]).reset_index(drop=True)
    return matrix, snippets


def _build_linked_resources() -> pd.DataFrame:
    rows = [
        {
            "resource_id": "ICAN",
            "resource_name": "French intracranial aneurysm biobank linked via WECAN",
            "evidence_source": "UCAN__wecan_project_page",
            "resource_type": "linked_biobank",
            "public_download_status": "indirect_only",
            "summary": "Publicly downloaded WECAN materials describe the ICAN biobank as including 3400 intracranial aneurysm patients with clinical, imaging, and DNA data.",
            "implication": "Potentially important companion validation / multimodal resource, but not directly downloadable as a standalone cohort from the current scripted pass.",
        }
    ]
    return pd.DataFrame(rows)


def _write_data_dictionary() -> pd.DataFrame:
    df = pd.DataFrame(
        MINIMUM_DATA_DICTIONARY,
        columns=["field_name", "data_type", "priority", "description"],
    )
    df.to_csv(DATA_DICTIONARY_PATH, sep="\t", index=False)
    return df


def _write_summary(matrix: pd.DataFrame, linked: pd.DataFrame) -> None:
    top_cols = [
        "cohort_id",
        "cohort_name",
        "design_priority_score",
        "public_evidence_score",
        "blood_biospecimen",
        "rupture_endpoint",
        "growth_endpoint",
        "phases_core_variable_count",
        "recommended_analysis_role",
    ]
    lines = [
        "# Clinical Validation Analysis Summary",
        "",
        "Derived from the downloaded public cohort materials already present in `results/clinical_cohort_scout/downloads`.",
        "",
        "## Main Conclusions",
        "",
        "- No cohort in the current public-only download set provides participant-level data that can be used immediately for external rupture-risk validation.",
        "- `CIAP-1` and `CIAP-2` remain the strongest public-design matches because the downloaded materials support prospective follow-up, rupture-risk framing, and blood / biomarker collection.",
        "- `ROAR` is the strongest public clinical-outcome framework with explicit PHASES validation and national follow-up linkage, but the downloaded materials do not explicitly confirm banked blood for molecular testing.",
        "- `UCAN` is well supported as a growth / instability cohort rather than a primary hard-rupture validation cohort.",
        "- `phs003072` is the only downloaded cohort that explicitly points to a blood RNA dataset, but it is controlled-access and cross-sectional.",
        "",
        "## Cohort Matrix",
        "",
        matrix[top_cols].to_markdown(index=False),
        "",
        "## Public-Only Downstream Analyses Completed",
        "",
        "1. Structured variable-coverage extraction from HTML and PDF cohort materials.",
        "2. Cohort-level evidence scoring for blood readiness, endpoint type, and PHASES-variable support.",
        "3. A minimum external-validation data dictionary defining the locked fields required for `PHASES + panel` transport testing.",
        "4. Linked-resource capture for companion assets mentioned within downloaded materials.",
        "",
        "## Still Requires Restricted or Collaborative Access",
        "",
        "- Per-patient baseline clinical variables for PHASES calculation.",
        "- Per-patient biospecimen-linked expression or protein measurements.",
        "- Rupture dates, censoring dates, and treatment dates needed for time-to-event or calibration analyses.",
        "",
        "## Linked Resources",
        "",
        linked.to_markdown(index=False),
    ]
    SUMMARY_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    manifest, assets = _read_inputs()
    text_map = _build_text_map(assets)
    matrix, snippets = _build_variable_matrix(manifest, assets, text_map)
    linked = _build_linked_resources()
    _write_data_dictionary()

    matrix.to_csv(VARIABLE_MATRIX_PATH, sep="\t", index=False)
    snippets.to_csv(SNIPPETS_PATH, sep="\t", index=False)
    linked.to_csv(LINKED_RESOURCES_PATH, sep="\t", index=False)
    _write_summary(matrix, linked)

    print(matrix[["cohort_id", "public_evidence_score", "recommended_analysis_role"]].to_string(index=False))
    print(f"\nWrote: {VARIABLE_MATRIX_PATH.relative_to(ROOT)}")
    print(f"Wrote: {SNIPPETS_PATH.relative_to(ROOT)}")
    print(f"Wrote: {DATA_DICTIONARY_PATH.relative_to(ROOT)}")
    print(f"Wrote: {LINKED_RESOURCES_PATH.relative_to(ROOT)}")
    print(f"Wrote: {SUMMARY_PATH.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
