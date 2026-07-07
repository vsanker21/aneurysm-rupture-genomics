"""
Programmatically scout clinical aneurysm cohorts for prospective biomarker validation.

This script does not claim that controlled-access cohort data are publicly downloadable.
Instead, it:

1. Downloads all public-facing protocols / metadata pages that can be reached now.
2. Records when a source is blocked, anti-bot protected, or controlled-access only.
3. Produces a structured manifest ranking cohorts by suitability for a study requiring:
   - baseline blood or other biospecimens,
   - documented aneurysm rupture follow-up,
   - and PHASES-compatible variables.

Outputs are written to `results/clinical_cohort_scout/`:
  - `clinical_cohort_manifest.tsv`
  - `clinical_cohort_manifest.json`
  - `clinical_cohort_scout_summary.md`
  - `downloads/` with fetched public assets
"""

from __future__ import annotations

import json
import re
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "results" / "clinical_cohort_scout"
DOWNLOAD_DIR = OUT_DIR / "downloads"
MANIFEST_TSV = OUT_DIR / "clinical_cohort_manifest.tsv"
MANIFEST_JSON = OUT_DIR / "clinical_cohort_manifest.json"
SUMMARY_MD = OUT_DIR / "clinical_cohort_scout_summary.md"
ASSET_AUDIT_TSV = OUT_DIR / "clinical_cohort_asset_audit.tsv"
DOWNSTREAM_PLAN_MD = OUT_DIR / "clinical_cohort_downstream_steps.md"

REQUEST_HEADERS = {
    "User-Agent": "AneurysmClinicalCohortScout/1.0 (public-access cohort feasibility audit)"
}


@dataclass
class SourceAsset:
    label: str
    url: str
    expected_format: str = "text"  # text | html | pdf | json
    role: str = ""


@dataclass
class AssetAuditRecord:
    cohort_id: str
    label: str
    role: str
    url: str
    expected_format: str
    observed_format: str
    content_type: str
    file_size_bytes: int
    status: str
    local_path: str
    note: str


@dataclass
class CohortRecord:
    cohort_id: str
    cohort_name: str
    cohort_role: str
    region: str
    evidence_urls: list[str] = field(default_factory=list)
    public_assets_downloaded: list[str] = field(default_factory=list)
    public_asset_failures: list[str] = field(default_factory=list)
    blood_or_biospecimen_evidence: str = ""
    rupture_followup_evidence: str = ""
    phases_variable_evidence: str = ""
    baseline_blood_status: str = ""  # yes / no / uncertain / dna_only
    rupture_followup_status: str = ""  # yes / no / growth_only / uncertain
    phases_ready_status: str = ""  # yes / partial / uncertain / no
    public_data_access_level: str = ""  # public_protocol_only / controlled_access / blocked / unknown
    suitability_priority: str = ""  # very_high / high / moderate / low
    suitability_score: int = 0
    rationale: str = ""
    next_action: str = ""
    notes: str = ""
    total_public_assets_checked: int = 0
    downloadable_public_assets: int = 0
    blocked_or_failed_assets: int = 0


def _safe_slug(value: str) -> str:
    return re.sub(r"[^a-zA-Z0-9._-]+", "_", value).strip("_")


def _fetch_url(url: str, timeout: int = 90) -> tuple[int, bytes, dict[str, str]]:
    last_error: Exception | None = None
    for attempt in range(4):
        req = Request(url, headers=REQUEST_HEADERS)
        try:
            with urlopen(req, timeout=timeout) as response:
                headers = {k.lower(): v for k, v in response.info().items()}
                return int(response.status), response.read(), headers
        except HTTPError as exc:
            last_error = exc
            if exc.code not in {429, 500, 502, 503, 504} or attempt == 3:
                raise
            time.sleep(2.0 * (attempt + 1))
        except URLError as exc:
            last_error = exc
            if attempt == 3:
                raise
            time.sleep(2.0 * (attempt + 1))
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            if attempt == 3:
                raise
            time.sleep(2.0 * (attempt + 1))
    raise RuntimeError(f"Failed to fetch {url}: {last_error}")


def _detect_format(payload: bytes, content_type: str) -> str:
    if payload.startswith(b"%PDF"):
        return "pdf"
    prefix = payload[:400].decode("utf-8", errors="ignore").lstrip().lower()
    if "application/json" in (content_type or "").lower() or prefix.startswith("{") or prefix.startswith("["):
        return "json"
    if prefix.startswith("<!doctype html") or prefix.startswith("<html") or "<html" in prefix[:200]:
        return "html"
    return "text"


def _save_asset(cohort_id: str, asset: SourceAsset) -> AssetAuditRecord:
    note = ""
    try:
        status, payload, headers = _fetch_url(asset.url)
        if status != 200:
            return AssetAuditRecord(
                cohort_id=cohort_id,
                label=asset.label,
                role=asset.role,
                url=asset.url,
                expected_format=asset.expected_format,
                observed_format="",
                content_type=headers.get("content-type", ""),
                file_size_bytes=0,
                status=f"http_{status}",
                local_path="",
                note=f"HTTP {status}",
            )
        observed_format = _detect_format(payload, headers.get("content-type", ""))
        suffix = {
            "text": ".txt",
            "html": ".html",
            "pdf": ".pdf",
            "json": ".json",
        }.get(observed_format, ".bin")
        if observed_format != asset.expected_format:
            note = f"Expected {asset.expected_format}, downloaded {observed_format}."
        if observed_format in {"text", "html"}:
            text = payload.decode("utf-8", errors="ignore")
            if "Making sure you're not a bot!" in text or "Proof-of-Work" in text:
                return AssetAuditRecord(
                    cohort_id=cohort_id,
                    label=asset.label,
                    role=asset.role,
                    url=asset.url,
                    expected_format=asset.expected_format,
                    observed_format=observed_format,
                    content_type=headers.get("content-type", ""),
                    file_size_bytes=len(payload),
                    status="blocked",
                    local_path="",
                    note="Blocked by anti-bot challenge.",
                )
            if "Access denied" in text or "Sign in" in text:
                return AssetAuditRecord(
                    cohort_id=cohort_id,
                    label=asset.label,
                    role=asset.role,
                    url=asset.url,
                    expected_format=asset.expected_format,
                    observed_format=observed_format,
                    content_type=headers.get("content-type", ""),
                    file_size_bytes=len(payload),
                    status="auth_wall",
                    local_path="",
                    note="Authentication wall.",
                )
        out_name = f"{_safe_slug(cohort_id)}__{_safe_slug(asset.label)}{suffix}"
        out_path = DOWNLOAD_DIR / out_name
        out_path.write_bytes(payload)
        return AssetAuditRecord(
            cohort_id=cohort_id,
            label=asset.label,
            role=asset.role,
            url=asset.url,
            expected_format=asset.expected_format,
            observed_format=observed_format,
            content_type=headers.get("content-type", ""),
            file_size_bytes=len(payload),
            status="downloaded",
            local_path=out_path.relative_to(ROOT).as_posix(),
            note=note,
        )
    except HTTPError as exc:
        return AssetAuditRecord(
            cohort_id=cohort_id,
            label=asset.label,
            role=asset.role,
            url=asset.url,
            expected_format=asset.expected_format,
            observed_format="",
            content_type="",
            file_size_bytes=0,
            status=f"http_error_{exc.code}",
            local_path="",
            note=f"HTTPError {exc.code}",
        )
    except Exception as exc:  # noqa: BLE001
        return AssetAuditRecord(
            cohort_id=cohort_id,
            label=asset.label,
            role=asset.role,
            url=asset.url,
            expected_format=asset.expected_format,
            observed_format="",
            content_type="",
            file_size_bytes=0,
            status="error",
            local_path="",
            note=f"{type(exc).__name__}: {exc}",
        )


COHORTS: list[tuple[CohortRecord, list[SourceAsset]]] = [
    (
        CohortRecord(
            cohort_id="CIAP-1",
            cohort_name="China Intracranial Aneurysm Project-1",
            cohort_role="Prospective natural-history UIA registry",
            region="China",
            blood_or_biospecimen_evidence="Protocol states that blood samples are collected at first visit and stored in a central Biological Specimen Bank.",
            rupture_followup_evidence="Prospective UIA cohort with aneurysm rupture as the primary outcome and scheduled follow-up.",
            phases_variable_evidence="Age, hypertension, aneurysm size, and site are clearly captured; PHASES-style modeling is discussed, but earlier SAH should still be explicitly verified in the data dictionary.",
            baseline_blood_status="yes",
            rupture_followup_status="yes",
            phases_ready_status="partial",
            public_data_access_level="public_protocol_only",
            suitability_priority="very_high",
            suitability_score=92,
            rationale="Closest publicly documented match to the desired design: baseline blood banking plus prospective rupture follow-up in a large UIA cohort.",
            next_action="Contact CIAP investigators to request data dictionary, biospecimen SOPs, and collaboration / transfer pathway.",
        ),
        [
            SourceAsset(
                label="protocol_pdf",
                url="https://translational-medicine.biomedcentral.com/counter/pdf/10.1186/s12967-019-2092-z.pdf",
                expected_format="pdf",
                role="primary_protocol",
            ),
            SourceAsset(
                label="protocol_article_html",
                url="https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-019-2092-z",
                expected_format="html",
                role="primary_protocol",
            ),
            SourceAsset(
                label="clinicaltrials_page",
                url="https://clinicaltrials.gov/study/NCT03117803",
                expected_format="html",
                role="registry_page",
            ),
        ],
    ),
    (
        CohortRecord(
            cohort_id="CIAP-2",
            cohort_name="China Intracranial Aneurysm Project-2",
            cohort_role="Nested multidimensional rupture-risk model from CIAP-1",
            region="China",
            blood_or_biospecimen_evidence="Protocol explicitly evaluates biomarkers and genetics within a nested 500-patient UIA cohort and references blood sample collection.",
            rupture_followup_evidence="Rupture-risk model is built on the prospective CIAP-1 follow-up cohort.",
            phases_variable_evidence="Clinical factors are integrated alongside imaging, biomarkers, and hemodynamics; PHASES-like variables appear available via CIAP-1 but should be formally confirmed.",
            baseline_blood_status="yes",
            rupture_followup_status="yes",
            phases_ready_status="partial",
            public_data_access_level="public_protocol_only",
            suitability_priority="very_high",
            suitability_score=90,
            rationale="Probably the single best biologically enriched cohort for a locked biomarker validation if collaboration can be secured.",
            next_action="Prioritize investigator outreach because CIAP-2 already frames biomarker-supported rupture prediction.",
        ),
        [
            SourceAsset(
                label="protocol_pdf",
                url="https://translational-medicine.biomedcentral.com/counter/pdf/10.1186/s12967-018-1641-1.pdf",
                expected_format="pdf",
                role="primary_protocol",
            ),
            SourceAsset(
                label="pmc_protocol_html",
                url="https://pmc.ncbi.nlm.nih.gov/articles/PMC5988060/",
                expected_format="html",
                role="primary_protocol",
            ),
            SourceAsset(
                label="clinicaltrials_page",
                url="https://clinicaltrials.gov/study/NCT03133624",
                expected_format="html",
                role="registry_page",
            ),
        ],
    ),
    (
        CohortRecord(
            cohort_id="ROAR",
            cohort_name="Risk of Aneurysm Rupture study",
            cohort_role="Large UK longitudinal UIA cohort for PHASES validation",
            region="United Kingdom",
            blood_or_biospecimen_evidence="Core ROAR protocol does not clearly state a blood biobank in the main study text.",
            rupture_followup_evidence="Longitudinal national-linkage design with rupture resulting in admission or death as the primary endpoint.",
            phases_variable_evidence="ROAR objective 1 is explicit PHASES validation, so PHASES-compatible variables are clearly available.",
            baseline_blood_status="uncertain",
            rupture_followup_status="yes",
            phases_ready_status="yes",
            public_data_access_level="public_protocol_only",
            suitability_priority="high",
            suitability_score=84,
            rationale="Outstanding outcome cohort for testing PHASES and incremental prediction, but biospecimen availability must be confirmed before assay planning.",
            next_action="Ask ROAR investigators whether banked whole blood, PAXgene RNA, plasma, or serum exist for any subset of the cohort.",
        ),
        [
            SourceAsset(
                label="pmc_protocol_html",
                url="https://pmc.ncbi.nlm.nih.gov/articles/PMC10030903/",
                expected_format="html",
                role="primary_protocol",
            ),
            SourceAsset(
                label="bmj_full_html",
                url="https://bmjopen.bmj.com/content/13/3/e070504.full",
                expected_format="html",
                role="journal_page",
            ),
            SourceAsset(
                label="roar_site_page",
                url="https://roarstudy.co.uk/roar-study",
                expected_format="html",
                role="study_site",
            ),
            SourceAsset(
                label="roar_gp_page",
                url="https://roarstudy.co.uk/roar-gp",
                expected_format="html",
                role="study_site",
            ),
        ],
    ),
    (
        CohortRecord(
            cohort_id="ROAR-DNA",
            cohort_name="ROAR-DNA",
            cohort_role="Genetic extension linked to ROAR",
            region="United Kingdom",
            blood_or_biospecimen_evidence="Public site states DNA from 6000 ROAR patients will be collected; this is a DNA resource rather than a transcript or plasma biomarker resource.",
            rupture_followup_evidence="Linked to ROAR long-term follow-up and designed to predict rupture risk genetically.",
            phases_variable_evidence="Inherited from ROAR context, but biomarker use is DNA-focused.",
            baseline_blood_status="dna_only",
            rupture_followup_status="yes",
            phases_ready_status="yes",
            public_data_access_level="public_protocol_only",
            suitability_priority="moderate",
            suitability_score=67,
            rationale="Valuable for future genetics integration, but not a direct substitute for banked expression- or protein-ready blood.",
            next_action="Use only as a complementary access target unless broader biospecimen banking is confirmed.",
        ),
        [
            SourceAsset(
                label="roar_dna_page",
                url="https://roarstudy.co.uk/roar-dna-3",
                expected_format="html",
                role="study_site",
            ),
        ],
    ),
    (
        CohortRecord(
            cohort_id="ISUIA",
            cohort_name="International Study of Unruptured Intracranial Aneurysms",
            cohort_role="Landmark prospective natural-history UIA study",
            region="North America and Europe",
            blood_or_biospecimen_evidence="No public evidence of an accessible blood or biospecimen repository was found in the widely cited ISUIA materials.",
            rupture_followup_evidence="Prospective assessment of unruptured aneurysm natural history with rupture outcomes.",
            phases_variable_evidence="Core predictors later used in PHASES are represented in ISUIA-style analyses.",
            baseline_blood_status="no",
            rupture_followup_status="yes",
            phases_ready_status="yes",
            public_data_access_level="public_protocol_only",
            suitability_priority="moderate",
            suitability_score=63,
            rationale="Excellent outcome dataset historically, but currently unsuitable for rapid biomarker validation unless hidden archived biospecimens can be confirmed by investigators.",
            next_action="Only pursue if an investigator explicitly confirms specimen availability and retrievable identifiers for matching outcomes.",
        ),
        [
            SourceAsset(
                label="pubmed_abstract",
                url="https://pubmed.ncbi.nlm.nih.gov/12867109/",
                expected_format="html",
                role="abstract_page",
            ),
        ],
    ),
    (
        CohortRecord(
            cohort_id="UCAN",
            cohort_name="Unruptured Cerebral Aneurysm: prediction of evolution",
            cohort_role="Prospective multicentre UIA evolution cohort",
            region="France",
            blood_or_biospecimen_evidence="Publicly reachable project pages indicate extensive clinical, biological, and imaging data, but specimen type and biobank details still need direct investigator confirmation.",
            rupture_followup_evidence="Public descriptions emphasize aneurysm growth / evolution; rupture appears secondary or limited.",
            phases_variable_evidence="Not fully confirmed from a downloadable protocol.",
            baseline_blood_status="uncertain",
            rupture_followup_status="growth_only",
            phases_ready_status="uncertain",
            public_data_access_level="public_protocol_only",
            suitability_priority="moderate",
            suitability_score=58,
            rationale="Potentially strong growth / instability resource, but less clearly aligned to hard rupture endpoints than CIAP or ROAR.",
            next_action="Treat as a secondary collaboration target focused on growth, instability, and mechanistic biomarker discovery.",
        ),
        [
            SourceAsset(
                label="hal_project_page",
                url="https://hal.science/hal-02889889v1",
                expected_format="html",
                role="primary_protocol",
            ),
            SourceAsset(
                label="hal_pdf",
                url="https://hal.science/hal-02889889v1/file/NEU-D-20-00179.pdf",
                expected_format="pdf",
                role="primary_protocol",
            ),
            SourceAsset(
                label="umr_team_page",
                url="https://umr1087.univ-nantes.fr/research/research-teams/intracranial-aneurysm",
                expected_format="html",
                role="study_site",
            ),
            SourceAsset(
                label="wecan_project_page",
                url="https://anr.fr/Project-ANR-21-CE17-0006",
                expected_format="html",
                role="related_project",
            ),
        ],
    ),
    (
        CohortRecord(
            cohort_id="phs003072",
            cohort_name="RNA-Seq of Whole Blood from Patients with Intracranial Aneurysms",
            cohort_role="Controlled-access whole-blood RNA cohort",
            region="United States",
            blood_or_biospecimen_evidence="Whole-blood RNA-seq samples are explicitly described.",
            rupture_followup_evidence="Risk groups are cross-sectional, stratified by PHASES, size, and Rupture Resemblance Score rather than prospective rupture follow-up.",
            phases_variable_evidence="PHASES-derived risk strata are explicitly used.",
            baseline_blood_status="yes",
            rupture_followup_status="no",
            phases_ready_status="yes",
            public_data_access_level="controlled_access",
            suitability_priority="moderate",
            suitability_score=55,
            rationale="Useful supporting blood-transcriptome resource, but not the retrospective-prospective rupture-outcome cohort needed for the target paper.",
            next_action="Use only as supportive comparator data or method-development resource, not as the primary clinical validation cohort.",
        ),
        [
            SourceAsset(
                label="ega_page",
                url="https://www.ega-archive.org/studies/phs003072",
                expected_format="html",
                role="controlled_access_metadata",
            ),
        ],
    ),
    (
        CohortRecord(
            cohort_id="FIA",
            cohort_name="Familial Intracranial Aneurysm Linkage Study",
            cohort_role="Familial susceptibility cohort",
            region="North America, Australia, New Zealand",
            blood_or_biospecimen_evidence="Blood samples and cryopreserved white blood cells are explicitly described.",
            rupture_followup_evidence="Designed for susceptibility and linkage mapping, not prospective rupture-outcome biomarker validation.",
            phases_variable_evidence="Not oriented around PHASES-style clinical rupture prediction.",
            baseline_blood_status="yes",
            rupture_followup_status="no",
            phases_ready_status="no",
            public_data_access_level="controlled_access",
            suitability_priority="low",
            suitability_score=34,
            rationale="Important genetics resource, but misaligned with the target clinical blood biomarker validation question.",
            next_action="Keep as a background genetics resource only.",
        ),
        [
            SourceAsset(
                label="ega_page",
                url="https://www.ega-archive.org/studies/phs000293",
                expected_format="html",
                role="controlled_access_metadata",
            ),
        ],
    ),
]


def _build_manifest() -> tuple[pd.DataFrame, pd.DataFrame]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    DOWNLOAD_DIR.mkdir(parents=True, exist_ok=True)

    records: list[CohortRecord] = []
    asset_records: list[AssetAuditRecord] = []

    for cohort, assets in COHORTS:
        cohort.evidence_urls = [asset.url for asset in assets]
        cohort.total_public_assets_checked = len(assets)
        for asset in assets:
            audit = _save_asset(cohort.cohort_id, asset)
            asset_records.append(audit)
            if audit.local_path:
                cohort.public_assets_downloaded.append(audit.local_path)
                cohort.downloadable_public_assets += 1
            else:
                failure_text = f"{asset.label}: {audit.status}"
                if audit.note:
                    failure_text += f" ({audit.note})"
                cohort.public_asset_failures.append(failure_text)
                cohort.blocked_or_failed_assets += 1

        if cohort.public_asset_failures and not cohort.public_assets_downloaded:
            cohort.notes = "No downloadable public assets were retrieved for this cohort during the scripted pass."
        records.append(cohort)

    df = pd.DataFrame(asdict(record) for record in records)
    df = df.sort_values(["suitability_score", "cohort_id"], ascending=[False, True]).reset_index(drop=True)
    asset_df = pd.DataFrame(asdict(record) for record in asset_records)
    asset_df = asset_df.sort_values(["cohort_id", "status", "label"]).reset_index(drop=True)
    return df, asset_df


def _write_summary(df: pd.DataFrame) -> None:
    top = df[
        [
            "cohort_id",
            "cohort_name",
            "baseline_blood_status",
            "rupture_followup_status",
            "phases_ready_status",
            "public_data_access_level",
            "suitability_priority",
            "next_action",
        ]
    ]

    counts = {
        "very_high": int((df["suitability_priority"] == "very_high").sum()),
        "high": int((df["suitability_priority"] == "high").sum()),
        "moderate": int((df["suitability_priority"] == "moderate").sum()),
        "low": int((df["suitability_priority"] == "low").sum()),
    }
    public_only = int((df["public_data_access_level"] == "public_protocol_only").sum())
    controlled = int((df["public_data_access_level"] == "controlled_access").sum())
    blocked = int((df["public_data_access_level"] == "blocked").sum())
    downloadable_assets = int(df["downloadable_public_assets"].sum())
    checked_assets = int(df["total_public_assets_checked"].sum())

    lines = [
        "# Clinical Cohort Scout Summary",
        "",
        "Generated programmatically from public-facing cohort protocols, archive pages, and study websites.",
        "",
        "## Overview",
        "",
        f"- Total target cohorts reviewed: `{len(df)}`",
        f"- `very_high` suitability: `{counts['very_high']}`",
        f"- `high` suitability: `{counts['high']}`",
        f"- `moderate` suitability: `{counts['moderate']}`",
        f"- `low` suitability: `{counts['low']}`",
        f"- Public protocol / metadata only: `{public_only}`",
        f"- Controlled-access study data referenced: `{controlled}`",
        f"- Public source fetches blocked by anti-bot protection: `{blocked}`",
        f"- Public assets downloaded successfully: `{downloadable_assets}` / `{checked_assets}`",
        "",
        "## Ranked Cohorts",
        "",
        top.to_markdown(index=False),
        "",
        "## Key Feasibility Conclusions",
        "",
        "- `CIAP-1` and `CIAP-2` are the closest documented matches to the target design because they combine prospective UIA follow-up with blood banking / biomarker plans.",
        "- `ROAR` is the strongest outcome and PHASES-validation cohort, but biospecimen availability is not explicit in the core protocol and must be confirmed.",
        "- `ISUIA` remains historically important for natural history, but no public biospecimen layer was confirmed.",
        "- `UCAN` is programmatically reachable, but the downloadable materials frame it more as a growth / evolution cohort than a hard rupture-outcome cohort.",
        "- `phs003072` and `FIA` are useful supportive resources, not substitutes for a true retrospective-prospective clinical validation cohort.",
        "",
        "## Immediate Access Strategy",
        "",
        "1. Prioritize collaboration requests to `CIAP-2` and `ROAR`.",
        "2. Ask specifically for data dictionary, biospecimen inventory, specimen type, storage conditions, and permissions for retrospective molecular testing.",
        "3. Treat `ISUIA` as an investigator-confirmation target only if archived biospecimens still exist.",
        "4. Keep `UCAN` as a secondary growth / instability cohort rather than the primary rupture-validation target.",
    ]
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_downstream_plan(df: pd.DataFrame, asset_df: pd.DataFrame) -> None:
    downloadable = df[df["downloadable_public_assets"] > 0][["cohort_id", "cohort_name", "downloadable_public_assets", "public_data_access_level"]]
    lines = [
        "# Clinical Cohort Downstream Steps",
        "",
        "This plan is restricted to what can be supported by the data and public assets fetched programmatically in the current workspace.",
        "",
        "## What Is Downloadable Now",
        "",
        downloadable.to_markdown(index=False),
        "",
        "## What Is Not Downloadable Now",
        "",
        "- Patient-level prospective validation data from `CIAP-1`, `CIAP-2`, `ROAR`, `UCAN`, and `ISUIA` are not openly downloadable from the public endpoints audited here.",
        "- `phs003072` and `FIA` expose study metadata pages, but their underlying participant-level data remain controlled-access.",
        "- No currently fetched source provides an immediately reusable public blood-plus-follow-up validation matrix for the reduced panel.",
        "",
        "## Downstream Steps That Can Proceed Immediately",
        "",
        "1. Lock the cohort priority order for access requests: `CIAP-2`, `CIAP-1`, `ROAR`, `UCAN`, then `ISUIA`.",
        "2. Use the downloaded protocols to build a minimum required data dictionary for external validation.",
        "3. Predefine the analysis package so a future transfer cohort can be processed without changing the model specification after access is granted.",
        "4. Keep `phs003072` as a supportive blood-transcriptome comparator, not as the primary clinical validation cohort.",
        "",
        "## Required External Validation Variables",
        "",
        "- Patient identifiers suitable for de-identified linkage",
        "- Baseline age, sex, hypertension, smoking, prior SAH, aneurysm size, aneurysm site, multiplicity",
        "- Baseline biospecimen type, collection timing, storage temperature, extraction protocol",
        "- Follow-up duration and rupture outcome date or censoring date",
        "- Treatment status and date, because treatment is a competing event for observed rupture",
        "- Imaging dates and growth status where available",
        "",
        "## Planned Analyses Once Access Is Secured",
        "",
        "1. Recompute the reduced blood-panel score in the external cohort using the locked gene set and directionality.",
        "2. Test discrimination for future rupture using AUROC, AUPRC, calibration slope/intercept, Brier score, and bootstrap confidence intervals.",
        "3. Compare `PHASES` alone versus `PHASES + panel` using likelihood-ratio testing, calibration, and decision-curve analysis.",
        "4. Run sensitivity analyses excluding treated-early aneurysms and evaluating time-to-event models where dates are available.",
        "5. Evaluate transportability across sex, aneurysm size bands, and anterior versus posterior circulation subgroups.",
        "",
        "## Asset-Level Audit",
        "",
        asset_df[["cohort_id", "label", "role", "status", "local_path", "note"]].to_markdown(index=False),
    ]
    DOWNSTREAM_PLAN_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    df, asset_df = _build_manifest()
    df.to_csv(MANIFEST_TSV, sep="\t", index=False)
    MANIFEST_JSON.write_text(df.to_json(orient="records", indent=2), encoding="utf-8")
    asset_df.to_csv(ASSET_AUDIT_TSV, sep="\t", index=False)
    _write_summary(df)
    _write_downstream_plan(df, asset_df)
    print(df[["cohort_id", "suitability_priority", "public_data_access_level", "suitability_score"]].to_string(index=False))
    print(f"\nWrote: {MANIFEST_TSV.relative_to(ROOT)}")
    print(f"Wrote: {MANIFEST_JSON.relative_to(ROOT)}")
    print(f"Wrote: {ASSET_AUDIT_TSV.relative_to(ROOT)}")
    print(f"Wrote: {SUMMARY_MD.relative_to(ROOT)}")
    print(f"Wrote: {DOWNSTREAM_PLAN_MD.relative_to(ROOT)}")
    print(f"Downloaded public assets into: {DOWNLOAD_DIR.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
