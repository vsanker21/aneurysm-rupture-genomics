"""
Programmatically discover public aneurysm/SAH transcriptomic cohorts.

This script creates a machine-readable cohort manifest to support a stricter
evidence hierarchy for downstream reprocessing, single-cell-informed cell-state
modelling, transportability analysis, and reduced blood-testable signature work.

Outputs are written to `results/public_cohort_discovery/`:
  - `public_transcriptomic_cohort_manifest.tsv`
  - `public_transcriptomic_cohort_manifest.json`
  - `public_cohort_discovery_summary.md`

It searches:
  - NCBI GEO DataSets (E-utilities `gds` database)
  - EMBL-EBI BioStudies / ArrayExpress

The resulting manifest is intentionally broad. Records are heuristically
classified into tiers so that discovery, orthogonal context, blood context,
single-cell references, and reference-only studies are not conflated.
"""

from __future__ import annotations

import json
import math
import re
import time
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import urlencode
from urllib.error import HTTPError
from urllib.request import Request, urlopen

import pandas as pd


ROOT = Path(__file__).resolve().parent
RESULTS_DIR = ROOT / "results" / "public_cohort_discovery"
INVENTORY_PATH = ROOT / "data" / "complete_dataset_inventory.json"

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BIOSTUDIES_BASE = "https://www.ebi.ac.uk/biostudies/api/v1"
REQUEST_HEADERS = {"User-Agent": "AneurysmCohortDiscovery/1.0 (public transcriptomic cohort audit)"}
REQUEST_DELAY_SEC = 0.15

GEO_SEARCH_PROFILES = [
    {
        "name": "geo_human_expression_general",
        "query": (
            "(intracranial aneurysm OR cerebral aneurysm OR brain aneurysm) "
            "AND Homo sapiens[Organism] "
            'AND ("expression profiling by array"[Filter] '
            'OR "expression profiling by high throughput sequencing"[Filter] '
            'OR "non coding RNA profiling by array"[Filter] '
            'OR "non coding RNA profiling by high throughput sequencing"[Filter])'
        ),
    },
    {
        "name": "geo_human_blood_sah_context",
        "query": (
            "(intracranial aneurysm OR cerebral aneurysm OR subarachnoid hemorrhage OR SAH) "
            "AND (blood OR PBMC OR peripheral blood) "
            "AND Homo sapiens[Organism] "
            'AND ("expression profiling by array"[Filter] '
            'OR "expression profiling by high throughput sequencing"[Filter])'
        ),
    },
    {
        "name": "geo_human_single_cell_reference",
        "query": (
            "(intracranial aneurysm OR cerebral aneurysm) "
            "AND (single cell OR single-cell OR single nucleus OR snRNA OR scRNA OR spatial) "
            "AND Homo sapiens[Organism]"
        ),
    },
    {
        "name": "geo_any_species_single_cell_reference",
        "query": (
            "(intracranial aneurysm OR cerebral aneurysm) "
            "AND (single cell OR single-cell OR single nucleus OR snRNA OR scRNA OR spatial)"
        ),
    },
]

BIOSTUDIES_SEARCH_PROFILES = [
    {"name": "biostudies_intracranial_aneurysm", "query": '"intracranial aneurysm"'},
    {"name": "biostudies_cerebral_aneurysm", "query": '"cerebral aneurysm"'},
    {"name": "biostudies_sah_blood", "query": '"subarachnoid haemorrhage" aneurysm blood transcriptome'},
    {"name": "biostudies_single_cell", "query": '"intracranial aneurysm" "single cell"'},
]

TIER_PRIORITY = {
    "tier1_direct_rupture_tissue": 1,
    "tier2_tissue_context": 2,
    "tier3_blood_rupture_sah_context": 3,
    "tier4_blood_aneurysm_control_context": 4,
    "tier5_human_single_cell_spatial_reference": 5,
    "tier6_nonhuman_or_model_reference": 6,
    "exclude_or_manual_review": 7,
}

LOCAL_PATH_OVERRIDES = {
    "GSE122897": "data/sc_rnaseq/GSE122897",
    "GSE159610": "data/blood_rnaseq/GSE159610",
    "GSE13353": "data/microarray/GSE13353",
    "GSE54083": "data/microarray/GSE54083",
    "GSE75436": "data/microarray/GSE75436",
    "GSE6551": "data/microarray/GSE6551",
    "GSE50867": "data/microarray/GSE50867",
    "GSE26969": "data/microarray/GSE26969",
    "GSE36791": "data/microarray/GSE36791",
    "E-GEOD-15629": "data/microarray/E-GEOD-15629",
    "GSE15629": "data/microarray/E-GEOD-15629",
    "GSE73378": "data/microarray/GSE73378",
    "GSE66238": "data/sc_rnaseq/GSE66238",
    "GSE158558": "data/sc_rnaseq/GSE158558",
    "GSE106520": "data/blood_rnaseq/GSE106520",
    "GSE231922": "data/blood_rnaseq/GSE231922",
}

KNOWN_ACCESSION_OVERRIDES: dict[str, dict[str, str]] = {
    "GSE122897": {
        "sample_source_guess": "aneurysm_wall_tissue",
        "comparison_guess": "ruptured_vs_unruptured",
        "evidence_tier": "tier1_direct_rupture_tissue",
        "recommended_use": "Primary direct rupture tissue discovery cohort",
        "notes": "RNA-seq aneurysm wall dataset already used in the current manuscript pipeline.",
    },
    "GSE54083": {
        "sample_source_guess": "aneurysm_wall_tissue",
        "comparison_guess": "ruptured_vs_unruptured",
        "evidence_tier": "tier1_direct_rupture_tissue",
        "recommended_use": "Independent direct rupture tissue replication cohort",
    },
    "GSE13353": {
        "sample_source_guess": "aneurysm_wall_tissue",
        "comparison_guess": "ruptured_vs_unruptured",
        "evidence_tier": "tier1_direct_rupture_tissue",
        "recommended_use": "Small direct rupture tissue support cohort",
    },
    "GSE26969": {
        "sample_source_guess": "aneurysm_tissue",
        "comparison_guess": "aneurysm_vs_control_vessel",
        "evidence_tier": "tier2_tissue_context",
        "recommended_use": "Orthogonal aneurysm-vs-control vessel tissue context",
    },
    "GSE36791": {
        "sample_source_guess": "blood",
        "comparison_guess": "rupture_sah_vs_control",
        "evidence_tier": "tier3_blood_rupture_sah_context",
        "recommended_use": "Blood post-rupture / SAH context only; not future rupture-risk validation",
    },
    "GSE159610": {
        "sample_source_guess": "blood",
        "comparison_guess": "aneurysm_vs_control",
        "evidence_tier": "tier4_blood_aneurysm_control_context",
        "recommended_use": "Blood aneurysm-vs-control context and bridge to accessible biospecimens",
    },
    "E-GEOD-15629": {
        "sample_source_guess": "aneurysm_tissue",
        "comparison_guess": "mixed_rupture_and_control_context",
        "evidence_tier": "tier2_tissue_context",
        "recommended_use": "Orthogonal mixed tissue context requiring careful phenotype parsing",
    },
    "GSE15629": {
        "sample_source_guess": "aneurysm_tissue",
        "comparison_guess": "mixed_rupture_and_control_context",
        "evidence_tier": "tier2_tissue_context",
        "recommended_use": "Orthogonal mixed tissue context requiring careful phenotype parsing",
    },
    "GSE106520": {
        "sample_source_guess": "blood",
        "comparison_guess": "aneurysm_vs_control",
        "evidence_tier": "tier4_blood_aneurysm_control_context",
        "recommended_use": "Blood neutrophil aneurysm-vs-control context for bridge-signature testing",
    },
}


@dataclass
class CohortRecord:
    accession: str
    title: str
    source_databases: list[str] = field(default_factory=list)
    source_links: list[str] = field(default_factory=list)
    species: str = ""
    assay_type: str = ""
    platform: str = ""
    sample_count: int | None = None
    sample_source_guess: str = ""
    comparison_guess: str = ""
    evidence_tier: str = "exclude_or_manual_review"
    recommended_use: str = ""
    relevance_score: int = 0
    local_status: str = "not_local"
    local_path: str = ""
    raw_or_supplementary_files_hint: str = ""
    publication_date: str = ""
    matched_search_profiles: list[str] = field(default_factory=list)
    matched_terms: list[str] = field(default_factory=list)
    summary_excerpt: str = ""
    notes: str = ""


def _read_json_url(url: str) -> dict[str, Any]:
    last_error: Exception | None = None
    for attempt in range(5):
        request = Request(url, headers=REQUEST_HEADERS)
        try:
            with urlopen(request, timeout=60) as response:
                return json.loads(response.read().decode("utf-8"))
        except HTTPError as exc:
            last_error = exc
            if exc.code not in {429, 500, 502, 503, 504} or attempt == 4:
                raise
            time.sleep(2.0 * (attempt + 1))
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            if attempt == 4:
                raise
            time.sleep(2.0 * (attempt + 1))
    raise RuntimeError(f"Unable to retrieve JSON from {url}: {last_error}")


def _request_pause() -> None:
    time.sleep(REQUEST_DELAY_SEC)


def _clean_text(value: str) -> str:
    return re.sub(r"\s+", " ", value or "").strip()


def _truncate(value: str, max_chars: int = 420) -> str:
    value = _clean_text(value)
    if len(value) <= max_chars:
        return value
    return value[: max_chars - 3].rstrip() + "..."


def _load_local_inventory() -> tuple[dict[str, str], set[str]]:
    local_paths: dict[str, str] = {}
    inventory_accessions: set[str] = set()

    if INVENTORY_PATH.is_file():
        inventory = json.loads(INVENTORY_PATH.read_text(encoding="utf-8"))
        downloaded = inventory.get("downloaded_datasets", {})
        inventory_accessions.update(downloaded.keys())

    for accession, rel_path in LOCAL_PATH_OVERRIDES.items():
        if (ROOT / rel_path).exists():
            local_paths[accession] = rel_path

    for accession in inventory_accessions:
        if accession not in local_paths:
            for root_name in ("data", "results"):
                for match in (ROOT / root_name).rglob(accession):
                    if match.is_dir():
                        local_paths[accession] = match.relative_to(ROOT).as_posix()
                        break
                if accession in local_paths:
                    break

    return local_paths, inventory_accessions


def _geo_esearch(query: str, retmax: int = 200) -> list[str]:
    params = {
        "db": "gds",
        "term": query,
        "retmode": "json",
        "retmax": retmax,
    }
    url = f"{EUTILS_BASE}/esearch.fcgi?{urlencode(params)}"
    data = _read_json_url(url)
    _request_pause()
    return data.get("esearchresult", {}).get("idlist", [])


def _chunked(values: list[str], size: int) -> list[list[str]]:
    return [values[i : i + size] for i in range(0, len(values), size)]


def _geo_esummary(ids: list[str]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for chunk in _chunked(ids, 100):
        params = {
            "db": "gds",
            "id": ",".join(chunk),
            "retmode": "json",
        }
        url = f"{EUTILS_BASE}/esummary.fcgi?{urlencode(params)}"
        data = _read_json_url(url)
        result = data.get("result", {})
        for uid in result.get("uids", []):
            record = result.get(uid)
            if isinstance(record, dict):
                records.append(record)
        _request_pause()
    return records


def _biostudies_search(query: str, page_size: int = 25, max_pages: int = 2, max_results: int = 30) -> list[dict[str, Any]]:
    all_hits: list[dict[str, Any]] = []
    for page in range(1, max_pages + 1):
        params = {"query": query, "pageSize": page_size, "page": page}
        url = f"{BIOSTUDIES_BASE}/search?{urlencode(params)}"
        data = _read_json_url(url)
        hits = data.get("hits", [])
        if not hits:
            break
        all_hits.extend(hits)
        if len(all_hits) >= max_results or page * page_size >= int(data.get("totalHits", 0)):
            break
        _request_pause()
    return all_hits[:max_results]


def _biostudies_fetch(accession: str) -> dict[str, Any]:
    url = f"{BIOSTUDIES_BASE}/studies/{accession}"
    data = _read_json_url(url)
    _request_pause()
    return data


def _extract_biostudies_attr(attributes: list[dict[str, Any]], name: str) -> str:
    for attribute in attributes or []:
        if attribute.get("name") == name:
            return str(attribute.get("value", ""))
    return ""


def _keyword_matches(text: str) -> list[str]:
    lowered = text.lower()
    keyword_patterns = {
        "aneurysm": r"\baneurysm\b",
        "intracranial": r"\bintracranial\b",
        "cerebral": r"\bcerebral\b",
        "rupture": r"\bruptur",
        "unruptured": r"\bunruptured\b",
        "SAH": r"\bsubarachnoid|\bsah\b",
        "blood": r"\bblood\b|\bpbmc\b|\bperipheral blood\b",
        "tissue": r"\btissue\b|\bwall\b|\bartery\b|\bvascular\b",
        "single_cell": r"single cell|single-cell|single nucleus|snrna|scrna",
        "spatial": r"\bspatial\b",
        "mouse": r"\bmus musculus\b|\bmouse\b",
    }
    matches = [name for name, pattern in keyword_patterns.items() if re.search(pattern, lowered)]
    return sorted(matches)


def _has_aneurysm_context(text: str) -> bool:
    lowered = text.lower()
    return bool(
        re.search(
            r"intracranial aneurysm|cerebral aneurysm|brain aneurysm|aneurysm wall|saccular intracranial aneurysm|\bia\b|\bias\b",
            lowered,
        )
    )


def _has_transcriptomic_context(text: str) -> bool:
    lowered = text.lower()
    return bool(
        re.search(
            r"expression profiling|gene expression|transcriptom|rna-seq|rna seq|microarray|scrna|single cell|single-cell|snrna",
            lowered,
        )
    )


def _guess_sample_source(text: str) -> str:
    lowered = text.lower()
    if re.search(r"\bipsc\b|induced pluripotent|cell line|in vitro|organoid|fibroblast-derived", lowered):
        return "cell_model"
    if re.search(r"single cell|single-cell|single nucleus|snrna|scrna", lowered):
        return "single_cell"
    if "spatial" in lowered:
        return "spatial"
    if re.search(r"\bblood\b|\bpbmc\b|\bperipheral blood\b|\bwhole blood\b|\bplasma\b|\bserum\b|\bneutrophil|\bmonocyte|\bleukocyte", lowered):
        return "blood"
    if re.search(r"\baneurysm wall\b|\btissue\b|\bartery\b|\bvascular\b|\bwall\b", lowered):
        return "aneurysm_tissue"
    return "unspecified"


def _guess_comparison(text: str) -> str:
    lowered = text.lower()
    has_aneurysm = _has_aneurysm_context(lowered)
    if re.search(r"\bipsc\b|induced pluripotent|cell line|in vitro|organoid", lowered):
        return "model_reference"
    if has_aneurysm and "ruptured" in lowered and "unruptured" in lowered:
        return "ruptured_vs_unruptured"
    if has_aneurysm and re.search(r"subarachnoid haemorrhage|subarachnoid hemorrhage|\basah\b|\bsah\b", lowered) and re.search(r"control|healthy", lowered):
        return "rupture_sah_vs_control"
    if has_aneurysm and re.search(r"control|normal|healthy", lowered):
        return "aneurysm_vs_control"
    if "growth" in lowered or "progression" in lowered:
        return "aneurysm_progression_or_growth"
    if has_aneurysm and re.search(r"single cell|single-cell|single nucleus|snrna|scrna|spatial", lowered):
        return "cell_state_reference"
    if re.search(r"mouse|murine|rat", lowered):
        return "model_reference"
    return "manual_review"


def _score_relevance(text: str, species: str, sample_count: int | None) -> int:
    score = 0
    lowered = text.lower()
    if _has_aneurysm_context(lowered):
        score += 4
    if "ruptur" in lowered:
        score += 2
    if re.search(r"subarachnoid|sah", lowered):
        score += 1
    if re.search(r"single cell|single-cell|single nucleus|snrna|scrna|spatial", lowered):
        score += 2
    if species.lower() == "homo sapiens":
        score += 3
    elif species:
        score += 1
    if sample_count:
        score += min(4, int(math.log10(max(sample_count, 1))) + 1)
    return score


def _classify_tier(species: str, sample_source_guess: str, comparison_guess: str) -> tuple[str, str]:
    species_l = species.lower()
    human = species_l == "homo sapiens"
    single_or_spatial = sample_source_guess in {"single_cell", "spatial"}
    blood = sample_source_guess == "blood"
    tissue = sample_source_guess == "aneurysm_tissue"
    cell_model = sample_source_guess == "cell_model"

    if human and tissue and comparison_guess == "ruptured_vs_unruptured":
        return "tier1_direct_rupture_tissue", "Direct rupture tissue discovery or replication cohort"
    if human and tissue and comparison_guess in {"aneurysm_vs_control", "manual_review"}:
        return "tier2_tissue_context", "Orthogonal aneurysm tissue context cohort"
    if human and blood and comparison_guess == "rupture_sah_vs_control":
        return "tier3_blood_rupture_sah_context", "Blood post-rupture / SAH context only"
    if human and blood and comparison_guess in {"aneurysm_vs_control", "aneurysm_progression_or_growth"}:
        return "tier4_blood_aneurysm_control_context", "Blood-accessible aneurysm context cohort"
    if human and single_or_spatial and comparison_guess == "cell_state_reference":
        return "tier5_human_single_cell_spatial_reference", "Human single-cell / spatial reference for cell-state modelling"
    if cell_model:
        return "exclude_or_manual_review", "Model-system cohort; keep separate from clinical public transcriptomic tiers"
    if not human and (single_or_spatial or comparison_guess == "model_reference" or "mus musculus" in species_l):
        return "tier6_nonhuman_or_model_reference", "Reference-only nonhuman or model cohort"
    return "exclude_or_manual_review", "Manual review required before downstream use"


def _apply_known_overrides(record: CohortRecord) -> None:
    override = KNOWN_ACCESSION_OVERRIDES.get(record.accession)
    if not override:
        return
    for field_name, value in override.items():
        setattr(record, field_name, value)


def _record_from_geo(summary: dict[str, Any], matched_profile: str) -> CohortRecord | None:
    if summary.get("entrytype") != "GSE":
        return None

    accession = str(summary.get("accession", "")).strip()
    if not accession:
        return None

    species = _clean_text(str(summary.get("taxon", "")))
    title = _clean_text(str(summary.get("title", "")))
    summary_text = _clean_text(str(summary.get("summary", "")))
    text = f"{title} {summary_text}"
    sample_count = len(summary.get("samples", [])) if isinstance(summary.get("samples"), list) else None
    assay_type = _clean_text(str(summary.get("gdstype", "")))
    platform = _clean_text(str(summary.get("gpl", "")))
    matched_terms = _keyword_matches(text)
    sample_source_guess = _guess_sample_source(text)
    comparison_guess = _guess_comparison(text)
    evidence_tier, recommended_use = _classify_tier(species, sample_source_guess, comparison_guess)

    record = CohortRecord(
        accession=accession,
        title=title,
        source_databases=["GEO"],
        source_links=[f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"],
        species=species,
        assay_type=assay_type,
        platform=platform,
        sample_count=sample_count,
        sample_source_guess=sample_source_guess,
        comparison_guess=comparison_guess,
        evidence_tier=evidence_tier,
        recommended_use=recommended_use,
        relevance_score=_score_relevance(text, species, sample_count),
        raw_or_supplementary_files_hint=_clean_text(str(summary.get("suppfile", ""))),
        publication_date=_clean_text(str(summary.get("pdat", ""))),
        matched_search_profiles=[matched_profile],
        matched_terms=matched_terms,
        summary_excerpt=_truncate(summary_text),
    )
    _apply_known_overrides(record)
    return record


def _record_from_biostudies(hit: dict[str, Any], detail: dict[str, Any], matched_profile: str) -> CohortRecord | None:
    accession = str(hit.get("accession", "")).strip()
    if not accession:
        return None

    section = detail.get("section", {}) if isinstance(detail, dict) else {}
    section_attributes = section.get("attributes", []) if isinstance(section, dict) else []
    top_attributes = detail.get("attributes", []) if isinstance(detail, dict) else []

    title = _clean_text(
        str(hit.get("title", ""))
        or _extract_biostudies_attr(section_attributes, "Title")
        or _extract_biostudies_attr(top_attributes, "Title")
    )
    description = _clean_text(_extract_biostudies_attr(section_attributes, "Description") or str(hit.get("content", "")))
    species = _clean_text(_extract_biostudies_attr(section_attributes, "Organism"))
    assay_type = _clean_text(_extract_biostudies_attr(section_attributes, "Study type"))
    release_date = _clean_text(_extract_biostudies_attr(top_attributes, "ReleaseDate") or str(hit.get("release_date", "")))
    text = f"{title} {description}"
    matched_terms = _keyword_matches(text)

    subsection_count = 0
    if isinstance(section, dict):
        subsection_count = len(section.get("subsections", []))

    sample_count = None
    sample_match = re.search(r"\btotal of (\d+)\b|\b(\d+)\s+(?:arrays|samples|specimens)\b", text.lower())
    if sample_match:
        numbers = [group for group in sample_match.groups() if group]
        if numbers:
            sample_count = int(numbers[0])

    sample_source_guess = _guess_sample_source(text)
    comparison_guess = _guess_comparison(text)
    evidence_tier, recommended_use = _classify_tier(species, sample_source_guess, comparison_guess)

    record = CohortRecord(
        accession=accession,
        title=title,
        source_databases=["BioStudies"],
        source_links=[
            f"https://www.ebi.ac.uk/biostudies/studies/{accession}",
            f"https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{accession}" if accession.startswith("E-") else "",
        ],
        species=species,
        assay_type=assay_type,
        platform="",
        sample_count=sample_count,
        sample_source_guess=sample_source_guess,
        comparison_guess=comparison_guess,
        evidence_tier=evidence_tier,
        recommended_use=recommended_use,
        relevance_score=_score_relevance(text, species, sample_count) + min(subsection_count, 2),
        raw_or_supplementary_files_hint=f"BioStudies files={hit.get('files', 0)} links={hit.get('links', 0)}",
        publication_date=release_date,
        matched_search_profiles=[matched_profile],
        matched_terms=matched_terms,
        summary_excerpt=_truncate(description or str(hit.get("content", ""))),
    )
    record.source_links = [link for link in record.source_links if link]
    _apply_known_overrides(record)
    return record


def _looks_relevant_from_preview(hit: dict[str, Any]) -> bool:
    preview = _clean_text(f"{hit.get('title', '')} {hit.get('content', '')}").lower()
    return _has_aneurysm_context(preview) and _has_transcriptomic_context(preview)


def _merge_records(existing: CohortRecord, new: CohortRecord) -> CohortRecord:
    merged = existing
    merged.source_databases = sorted(set(existing.source_databases + new.source_databases))
    merged.source_links = sorted(set(existing.source_links + new.source_links))
    merged.matched_search_profiles = sorted(set(existing.matched_search_profiles + new.matched_search_profiles))
    merged.matched_terms = sorted(set(existing.matched_terms + new.matched_terms))

    for field_name in (
        "title",
        "species",
        "assay_type",
        "platform",
        "sample_source_guess",
        "comparison_guess",
        "evidence_tier",
        "recommended_use",
        "raw_or_supplementary_files_hint",
        "publication_date",
        "summary_excerpt",
        "notes",
    ):
        current = getattr(merged, field_name)
        candidate = getattr(new, field_name)
        if (not current) or (field_name == "summary_excerpt" and len(candidate) > len(current)):
            setattr(merged, field_name, candidate)

    if merged.sample_count is None and new.sample_count is not None:
        merged.sample_count = new.sample_count
    merged.relevance_score = max(existing.relevance_score, new.relevance_score)

    if TIER_PRIORITY.get(new.evidence_tier, 99) < TIER_PRIORITY.get(existing.evidence_tier, 99):
        merged.evidence_tier = new.evidence_tier
        merged.recommended_use = new.recommended_use

    return merged


def _attach_local_status(records: dict[str, CohortRecord], local_paths: dict[str, str], inventory_accessions: set[str]) -> None:
    for accession, record in records.items():
        if accession in local_paths:
            record.local_status = "downloaded"
            record.local_path = local_paths[accession]
        elif accession in inventory_accessions:
            record.local_status = "inventory_only"
        else:
            record.local_status = "not_local"


def _records_to_dataframe(records: list[CohortRecord]) -> pd.DataFrame:
    rows = []
    for record in records:
        row = asdict(record)
        row["source_databases"] = "; ".join(record.source_databases)
        row["source_links"] = "; ".join(record.source_links)
        row["matched_search_profiles"] = "; ".join(record.matched_search_profiles)
        row["matched_terms"] = "; ".join(record.matched_terms)
        row["tier_priority"] = TIER_PRIORITY.get(record.evidence_tier, 99)
        rows.append(row)

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    return df.sort_values(
        by=["tier_priority", "local_status", "relevance_score", "sample_count", "accession"],
        ascending=[True, True, False, False, True],
        na_position="last",
    ).reset_index(drop=True)


def _write_summary(df: pd.DataFrame, out_path: Path) -> None:
    lines = [
        "# Public Cohort Discovery Summary",
        "",
        f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        "",
        "## Overview",
        "",
    ]

    if df.empty:
        lines.append("No public cohorts were discovered.")
        out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return

    tier_counts = Counter(df["evidence_tier"])
    species_counts = Counter(df["species"])
    local_counts = Counter(df["local_status"])

    lines.extend(
        [
            f"- Total discovered accessions: `{len(df)}`",
            f"- Already local: `{local_counts.get('downloaded', 0)}`",
            f"- Human single-cell / spatial references: `{tier_counts.get('tier5_human_single_cell_spatial_reference', 0)}`",
            f"- Direct rupture tissue cohorts: `{tier_counts.get('tier1_direct_rupture_tissue', 0)}`",
            f"- Blood-accessible cohorts: `{tier_counts.get('tier3_blood_rupture_sah_context', 0) + tier_counts.get('tier4_blood_aneurysm_control_context', 0)}`",
            "",
            "## Tier Counts",
            "",
        ]
    )
    for tier, count in sorted(tier_counts.items(), key=lambda item: TIER_PRIORITY.get(item[0], 99)):
        lines.append(f"- `{tier}`: `{count}`")

    lines.extend(["", "## Species Counts", ""])
    for species, count in sorted(species_counts.items(), key=lambda item: (-item[1], item[0])):
        label = species or "Unspecified"
        lines.append(f"- `{label}`: `{count}`")

    lines.extend(
        [
            "",
            "## Recommended Next-Step Cohorts",
            "",
            "| Accession | Tier | Samples | Species | Local | Recommended use |",
            "| --- | --- | ---: | --- | --- | --- |",
        ]
    )

    shortlisted = df[df["evidence_tier"] != "exclude_or_manual_review"].head(20)
    for _, row in shortlisted.iterrows():
        sample_count = "" if pd.isna(row["sample_count"]) else int(row["sample_count"])
        lines.append(
            f"| `{row['accession']}` | `{row['evidence_tier']}` | {sample_count} | `{row['species'] or 'Unspecified'}` | "
            f"`{row['local_status']}` | {row['recommended_use']} |"
        )

    lines.extend(
        [
            "",
            "## Recommended Workstream Mapping",
            "",
            "- `tier1_direct_rupture_tissue`: harmonized discovery / replication differential expression and meta-analysis",
            "- `tier2_tissue_context`: orthogonal vessel-wall context, not primary rupture validation",
            "- `tier3_blood_rupture_sah_context`: post-rupture systemic context only; do not present as future rupture-risk biology",
            "- `tier4_blood_aneurysm_control_context`: bridge to blood-accessible translational signatures",
            "- `tier5_human_single_cell_spatial_reference`: derive cell-state programs and vascular-specific deconvolution priors",
            "- `tier6_nonhuman_or_model_reference`: reference-only mechanistic support",
        ]
    )

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def discover_public_cohorts() -> tuple[pd.DataFrame, dict[str, Any]]:
    local_paths, inventory_accessions = _load_local_inventory()
    records: dict[str, CohortRecord] = {}
    search_log: list[dict[str, Any]] = []

    for profile in GEO_SEARCH_PROFILES:
        ids = _geo_esearch(profile["query"])
        summaries = _geo_esummary(ids)
        search_log.append({"source": "GEO", "profile": profile["name"], "hits": len(ids), "records": len(summaries)})
        for summary in summaries:
            record = _record_from_geo(summary, profile["name"])
            if record is None:
                continue
            existing = records.get(record.accession)
            records[record.accession] = _merge_records(existing, record) if existing else record

    for profile in BIOSTUDIES_SEARCH_PROFILES:
        hits = _biostudies_search(profile["query"])
        search_log.append({"source": "BioStudies", "profile": profile["name"], "hits": len(hits)})
        for hit in hits:
            accession = str(hit.get("accession", "")).strip()
            if not accession:
                continue
            if not _looks_relevant_from_preview(hit):
                continue
            detail = _biostudies_fetch(accession)
            record = _record_from_biostudies(hit, detail, profile["name"])
            if record is None:
                continue
            existing = records.get(record.accession)
            records[record.accession] = _merge_records(existing, record) if existing else record

    _attach_local_status(records, local_paths, inventory_accessions)

    df = _records_to_dataframe(list(records.values()))
    metadata = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "result_count": int(len(df)),
        "search_log": search_log,
        "local_paths_detected": local_paths,
        "inventory_accessions": sorted(inventory_accessions),
    }
    return df, metadata


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    df, metadata = discover_public_cohorts()

    manifest_tsv = RESULTS_DIR / "public_transcriptomic_cohort_manifest.tsv"
    manifest_json = RESULTS_DIR / "public_transcriptomic_cohort_manifest.json"
    summary_md = RESULTS_DIR / "public_cohort_discovery_summary.md"

    df.to_csv(manifest_tsv, sep="\t", index=False)
    manifest_payload = {
        "metadata": metadata,
        "records": df.to_dict(orient="records"),
    }
    manifest_json.write_text(json.dumps(manifest_payload, indent=2), encoding="utf-8")
    _write_summary(df, summary_md)

    print(f"Created: {manifest_tsv}")
    print(f"Created: {manifest_json}")
    print(f"Created: {summary_md}")
    if not df.empty:
        tier_counts = df["evidence_tier"].value_counts().to_dict()
        print("Tier counts:", json.dumps(tier_counts, indent=2))
        print("Top cohorts:")
        cols = ["accession", "evidence_tier", "sample_count", "species", "local_status", "recommended_use"]
        print(df[cols].head(15).to_string(index=False))


if __name__ == "__main__":
    main()
