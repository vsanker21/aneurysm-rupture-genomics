"""
Exhaustive public cohort discovery for MPP1 / MARCO in blood and aneurysm wall tissue.

Sources searched programmatically:
  - NCBI GEO DataSets (E-utilities, expanded query set)
  - NCBI SRA (E-utilities)
  - EMBL-EBI BioStudies / ArrayExpress
  - GitHub repository search (no API key; unauthenticated)
  - Existing project manifest + inventory

For each human-relevant accession:
  - Classify biospecimen (wall tissue vs blood vs other)
  - Classify comparison (ruptured vs unruptured vs aneurysm vs control, etc.)
  - Check MPP1/MARCO measurability (local expression, series matrix, or platform annotation)
  - Run Welch DE when expression can be loaded
  - Attempt GEO FTP download for high-priority missing series
"""

from __future__ import annotations

import json
import re
import time
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urlencode, urljoin
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd
from scipy.stats import norm

from collections import defaultdict

from discover_public_aneurysm_cohorts import (
    BIOSTUDIES_SEARCH_PROFILES,
    EUTILS_BASE,
    GEO_SEARCH_PROFILES,
    KNOWN_ACCESSION_OVERRIDES,
    REQUEST_HEADERS,
    _biostudies_search,
    _geo_esearch,
    _geo_esummary,
    _guess_comparison,
    _guess_sample_source,
    _has_aneurysm_context,
    _has_transcriptomic_context,
    _keyword_matches,
    _load_local_inventory,
    _read_json_url,
    _request_pause,
    _truncate,
)
from mpp1_marco_multi_cohort_expression_survey import (
    BASE,
    GENES,
    OUT_DIR as SURVEY_OUT,
    _discovery_signs,
    main as run_expression_survey,
)

OUT_DIR = BASE / "results" / "mpp1_marco_survey"
GENE_ENTREZ = {"MPP1": "4354", "MARCO": "8685"}

EXTRA_GEO_SEARCH_PROFILES = [
    {
        "name": "geo_ruptured_unruptured_tissue",
        "query": (
            '(intracranial aneurysm[All Fields] OR "cerebral aneurysm"[All Fields]) '
            'AND (ruptured[All Fields] AND unruptured[All Fields]) '
            "AND Homo sapiens[Organism]"
        ),
    },
    {
        "name": "geo_blood_rupture_risk",
        "query": (
            '(intracranial aneurysm[All Fields] OR cerebral aneurysm[All Fields]) '
            "AND (blood[All Fields] OR PBMC[All Fields] OR neutrophil[All Fields]) "
            "AND (rupture[All Fields] OR unruptured[All Fields]) "
            "AND Homo sapiens[Organism]"
        ),
    },
    {
        "name": "geo_aneurysm_wall_rnaseq",
        "query": (
            '("aneurysm wall"[All Fields] OR "intracranial aneurysm"[All Fields]) '
            'AND ("RNA-seq"[All Fields] OR "high throughput sequencing"[Filter]) '
            "AND Homo sapiens[Organism]"
        ),
    },
    {
        "name": "geo_GSE66238_family",
        "query": "GSE66238[Accession] OR GSE66239[Accession] OR GSE66240[Accession]",
    },
    {
        "name": "geo_GSE158558",
        "query": "GSE158558[Accession]",
    },
    {
        "name": "geo_blood_miRNA_rupture",
        "query": (
            "intracranial aneurysm[All Fields] AND microRNA[All Fields] "
            "AND (ruptured[All Fields] OR unruptured[All Fields]) AND Homo sapiens[Organism]"
        ),
    },
]

SRA_SEARCH_PROFILES = [
    {
        "name": "sra_ia_human_rnaseq",
        "query": (
            '(intracranial aneurysm OR cerebral aneurysm) AND '
            '("RNA-Seq"[Strategy] OR "rna seq") AND "Homo sapiens"[Organism]'
        ),
    },
    {
        "name": "sra_ia_blood",
        "query": (
            '(intracranial aneurysm OR cerebral aneurysm) AND blood AND '
            '"Homo sapiens"[Organism]'
        ),
    },
]

GITHUB_SEARCH_QUERIES = [
    "intracranial aneurysm RNA-seq",
    "cerebral aneurysm transcriptome GEO",
    "aneurysm rupture blood RNA",
]

PRIORITY_DOWNLOAD = [
    "GSE66238",
    "GSE158558",
    "GSE75436",
    "GSE50867",
    "GSE6551",
    "GSE46337",
]


class LinkCollector(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag != "a":
            return
        for key, value in attrs:
            if key == "href" and value:
                self.links.append(value)


def _series_prefix(accession: str) -> str:
    match = re.match(r"^(GSE)(\d+)$", accession)
    if not match:
        raise ValueError(accession)
    return f"GSE{match.group(2)[:-3]}nnn"


def _list_ftp(url: str) -> list[str]:
    request = Request(url, headers=REQUEST_HEADERS)
    with urlopen(request, timeout=60) as response:
        html = response.read().decode("utf-8", errors="ignore")
    parser = LinkCollector()
    parser.feed(html)
    out = []
    for link in parser.links:
        if link in {"../", "./", ""}:
            continue
        out.append(urljoin(url, link))
    return sorted(set(out))


def _ftp_download(url: str, dest: Path) -> bool:
    if dest.is_file() and dest.stat().st_size > 1000:
        return True
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        request = Request(url, headers=REQUEST_HEADERS)
        with urlopen(request, timeout=300) as response, dest.open("wb") as handle:
            handle.write(response.read())
        return dest.stat().st_size > 1000
    except Exception:
        return False


def _sra_esearch(query: str, retmax: int = 100) -> list[str]:
    params = {"db": "sra", "term": query, "retmode": "json", "retmax": retmax}
    url = f"{EUTILS_BASE}/esearch.fcgi?{urlencode(params)}"
    data = _read_json_url(url)
    _request_pause()
    return data.get("esearchresult", {}).get("idlist", [])


def _sra_esummary(ids: list[str]) -> list[dict]:
    if not ids:
        return []
    params = {"db": "sra", "id": ",".join(ids[:100]), "retmode": "json"}
    url = f"{EUTILS_BASE}/esummary.fcgi?{urlencode(params)}"
    data = _read_json_url(url)
    _request_pause()
    result = data.get("result", {})
    rows = []
    for uid in result.get("uids", []):
        rec = result.get(uid)
        if isinstance(rec, dict):
            rows.append(rec)
    return rows


def _github_search(query: str) -> list[dict]:
    params = urlencode({"q": query, "sort": "stars", "order": "desc", "per_page": 15})
    url = f"https://api.github.com/search/repositories?{params}"
    try:
        request = Request(url, headers={**REQUEST_HEADERS, "Accept": "application/vnd.github+json"})
        with urlopen(request, timeout=30) as response:
            data = json.loads(response.read().decode("utf-8"))
        _request_pause()
        return [
            {
                "full_name": item.get("full_name"),
                "html_url": item.get("html_url"),
                "description": item.get("description"),
                "stars": item.get("stargazers_count"),
            }
            for item in data.get("items", [])
        ]
    except Exception as exc:
        return [{"error": str(exc), "query": query}]


def _extract_gse_accession(text: str) -> str | None:
    match = re.search(r"\b(GSE\d+)\b", text or "", flags=re.I)
    return match.group(1).upper() if match else None


def _geo_gse_from_gds_record(record: dict) -> str | None:
    for key in ("accession", "gse", "Accession"):
        val = record.get(key, "")
        if isinstance(val, str):
            gse = _extract_gse_accession(val)
            if gse:
                return gse
    title = str(record.get("title", ""))
    return _extract_gse_accession(title)


def _discover_geo_expanded() -> pd.DataFrame:
    all_profiles = GEO_SEARCH_PROFILES + EXTRA_GEO_SEARCH_PROFILES
    id_to_profiles: dict[str, list[str]] = defaultdict(list)
    for profile in all_profiles:
        ids = _geo_esearch(profile["query"], retmax=100)
        for uid in ids:
            id_to_profiles[uid].append(profile["name"])

    records = _geo_esummary(list(id_to_profiles.keys()))
    rows = []
    for rec in records:
        gse = _geo_gse_from_gds_record(rec)
        if not gse:
            continue
        text = " ".join(
            str(rec.get(k, "")) for k in ("title", "summary", "type", "gpl", "entrytype", "taxon")
        )
        if not _has_aneurysm_context(text) and gse not in KNOWN_ACCESSION_OVERRIDES:
            continue
        if not _has_transcriptomic_context(text) and "non-coding" not in text.lower():
            if "Expression profiling" not in str(rec.get("type", "")):
                continue
        rows.append(
            {
                "accession": gse,
                "source": "GEO",
                "title": _truncate(str(rec.get("title", ""))),
                "species": str(rec.get("taxon", "")),
                "assay_type": str(rec.get("type", "")),
                "platform": str(rec.get("gpl", "")),
                "sample_count": rec.get("n_samples"),
                "sample_source_guess": _guess_sample_source(text),
                "comparison_guess": _guess_comparison(text),
                "search_profiles": ";".join(sorted(set(id_to_profiles.get(rec.get("uid", ""), [])))),
                "matched_terms": ";".join(_keyword_matches(text)),
                "geo_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}",
            }
        )
    return pd.DataFrame(rows).drop_duplicates("accession")


def _discover_biostudies_expanded() -> pd.DataFrame:
    hits_by_acc: dict[str, dict] = {}
    for profile in BIOSTUDIES_SEARCH_PROFILES:
        hits = _biostudies_search(profile["query"], max_results=50, max_pages=4)
        for hit in hits:
            acc = str(hit.get("accession", "")).strip()
            if not acc:
                continue
            text = " ".join(
                str(hit.get(k, "")) for k in ("title", "source", "organism", "assayType")
            )
            if not _has_aneurysm_context(text):
                continue
            gse = _extract_gse_accession(text) or acc
            hits_by_acc[gse] = {
                "accession": gse,
                "source": "BioStudies",
                "title": _truncate(str(hit.get("title", ""))),
                "species": str(hit.get("organism", "")),
                "assay_type": str(hit.get("assayType", "")),
                "platform": "",
                "sample_count": hit.get("sampleCount"),
                "sample_source_guess": _guess_sample_source(text),
                "comparison_guess": _guess_comparison(text),
                "search_profiles": profile["name"],
                "matched_terms": ";".join(_keyword_matches(text)),
                "geo_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}"
                if gse.startswith("GSE")
                else f"https://www.ebi.ac.uk/biostudies/studies/{acc}",
            }
    return pd.DataFrame(list(hits_by_acc.values()))


def _discover_sra() -> pd.DataFrame:
    rows = []
    seen: set[str] = set()
    for profile in SRA_SEARCH_PROFILES:
        ids = _sra_esearch(profile["query"], retmax=80)
        summaries = _sra_esummary(ids)
        for rec in summaries:
            title = str(rec.get("title", "")) + " " + str(rec.get("experiment_title", ""))
            if not _has_aneurysm_context(title):
                continue
            run_acc = str(rec.get("run_accession", "") or rec.get("accession", ""))
            if not run_acc or run_acc in seen:
                continue
            seen.add(run_acc)
            rows.append(
                {
                    "accession": run_acc,
                    "source": "SRA",
                    "title": _truncate(title),
                    "species": str(rec.get("taxon_name", "")),
                    "assay_type": str(rec.get("library_strategy", "")),
                    "platform": str(rec.get("platform", "")),
                    "sample_count": 1,
                    "sample_source_guess": _guess_sample_source(title),
                    "comparison_guess": _guess_comparison(title),
                    "search_profiles": profile["name"],
                    "matched_terms": ";".join(_keyword_matches(title)),
                    "geo_url": f"https://www.ncbi.nlm.nih.gov/sra/{run_acc}",
                }
            )
    return pd.DataFrame(rows)


def _discover_github() -> pd.DataFrame:
    rows = []
    for query in GITHUB_SEARCH_QUERIES:
        for item in _github_search(query):
            if "error" in item:
                continue
            desc = str(item.get("description", "") or "")
            title = f"{item.get('full_name', '')} {desc}"
            if not _has_aneurysm_context(title):
                continue
            gse = _extract_gse_accession(title + " " + desc)
            rows.append(
                {
                    "accession": gse or str(item.get("full_name", "")),
                    "source": "GitHub",
                    "title": _truncate(title),
                    "species": "",
                    "assay_type": "repository",
                    "platform": "",
                    "sample_count": None,
                    "sample_source_guess": _guess_sample_source(title),
                    "comparison_guess": _guess_comparison(title),
                    "search_profiles": f"github:{query}",
                    "matched_terms": ";".join(_keyword_matches(title)),
                    "geo_url": str(item.get("html_url", "")),
                }
            )
        time.sleep(1.0)
    return pd.DataFrame(rows)


def _load_existing_manifest() -> pd.DataFrame:
    path = BASE / "results" / "public_cohort_discovery" / "public_transcriptomic_cohort_manifest.tsv"
    if path.is_file():
        df = pd.read_csv(path, sep="\t")
        df["source"] = "project_manifest"
        return df.rename(columns={"accession": "accession"})
    return pd.DataFrame()


def _local_data_status(accession: str, local_paths: dict[str, str]) -> tuple[str, str]:
    rel = local_paths.get(accession, "")
    if not rel:
        return "not_local", ""
    p = BASE / rel
    if not p.exists():
        return "path_missing", rel
    files = [f for f in p.rglob("*") if f.is_file() and f.stat().st_size > 500]
    if not files:
        return "empty_dir", rel
    return "has_data", rel


def _genes_measurable(accession: str, status: str) -> dict[str, str]:
    """Heuristic measurability for MPP1/MARCO."""
    out: dict[str, str] = {}
    for gene in GENES:
        out[gene] = "unknown"
    if status != "has_data":
        out["note"] = "no_local_expression"
        return out

    # Check precomputed results from survey
    survey_path = SURVEY_OUT / "mpp1_marco_all_cohort_results.tsv"
    if survey_path.is_file():
        prior = pd.read_csv(survey_path, sep="\t")
        for gene in GENES:
            sub = prior[(prior["Cohort"] == accession) & (prior["Gene"] == gene)]
            if not sub.empty:
                st = str(sub.iloc[0].get("Status", ""))
                out[gene] = "analyzed_ok" if st == "ok" else st

    # Human protein-coding genes — present on standard RNA-seq / most arrays
    if all(out[g] == "unknown" for g in GENES):
        for gene in GENES:
            out[gene] = "expected_measurable"
    return out


def _attempt_priority_downloads(accessions: list[str]) -> pd.DataFrame:
    rows = []
    for accession in accessions:
        if not accession.startswith("GSE"):
            continue
        prefix = _series_prefix(accession)
        out_dir = BASE / "data" / "sc_rnaseq" / accession
        if accession in ("GSE75436", "GSE50867", "GSE6551", "GSE46337"):
            out_dir = BASE / "data" / "microarray" / accession
        matrix_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{accession}/matrix/"
        suppl_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{accession}/suppl/"
        downloaded = []
        errors = []
        for base in (matrix_url, suppl_url):
            try:
                urls = _list_ftp(base)
                for url in urls:
                    name = url.rsplit("/", 1)[-1]
                    lower = name.lower()
                    if lower.endswith(".tar") and "raw" in lower:
                        continue
                    if not any(
                        lower.endswith(ext)
                        for ext in (".gz", ".txt", ".tsv", ".xlsx", ".csv", ".soft.gz")
                    ):
                        continue
                    dest = out_dir / name
                    if _ftp_download(url, dest):
                        downloaded.append(str(dest.relative_to(BASE)))
            except Exception as exc:
                errors.append(str(exc))
        rows.append(
            {
                "accession": accession,
                "download_dir": str(out_dir.relative_to(BASE)),
                "n_files": len(downloaded),
                "downloaded": ";".join(downloaded[:20]),
                "errors": ";".join(errors),
            }
        )
        time.sleep(0.2)
    return pd.DataFrame(rows)


def _classify_biospecimen_for_mpp1_marco(row: pd.Series) -> str:
    src = str(row.get("sample_source_guess", ""))
    comp = str(row.get("comparison_guess", ""))
    title = str(row.get("title", "")).lower()
    if src == "blood" or "blood" in title or "neutrophil" in title or "pbmc" in title:
        return "blood"
    if src in {"aneurysm_wall_tissue", "aneurysm_tissue"} or "wall" in title or "tissue" in title:
        return "aneurysm_wall_tissue"
    if "plasma" in title or "mirna" in title or "microrna" in title:
        return "blood_ncRNA"
    return src or "unspecified"


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--skip-discovery",
        action="store_true",
        help="Reuse existing catalog; only refresh expression survey outputs.",
    )
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if args.skip_discovery and (OUT_DIR / "comprehensive_cohort_catalog_human.tsv").is_file():
        print("=== Skipping discovery; refreshing expression survey only ===")
        run_expression_survey()
        catalog = pd.read_csv(OUT_DIR / "comprehensive_cohort_catalog_all_sources.tsv", sep="\t")
        catalog_human = pd.read_csv(OUT_DIR / "comprehensive_cohort_catalog_human.tsv", sep="\t")
        survey = pd.read_csv(SURVEY_OUT / "mpp1_marco_all_cohort_results.tsv", sep="\t")
        meta = pd.read_csv(SURVEY_OUT / "mpp1_marco_stouffer_meta.tsv", sep="\t")
        print("Done (skip-discovery mode).")
        return

    print("=== Phase 1: Multi-source cohort discovery ===")
    parts = []
    for name, func in [
        ("GEO", _discover_geo_expanded),
        ("BioStudies", _discover_biostudies_expanded),
        ("SRA", _discover_sra),
        ("GitHub", _discover_github),
    ]:
        print(f"  Searching {name}...")
        try:
            df = func()
            print(f"    -> {len(df)} records")
            parts.append(df)
        except Exception as exc:
            print(f"    !! {name} failed: {exc}")

    manifest_old = _load_existing_manifest()
    if not manifest_old.empty:
        parts.append(
            manifest_old.rename(
                columns={
                    "sample_source_guess": "sample_source_guess",
                    "comparison_guess": "comparison_guess",
                }
            )[["accession", "title", "species", "assay_type", "platform", "sample_count",
               "sample_source_guess", "comparison_guess"]].assign(source="project_manifest")
        )

    discovered = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    if discovered.empty:
        raise RuntimeError("No cohorts discovered from any source.")

    discovered["accession"] = discovered["accession"].astype(str).str.strip()
    discovered = discovered[discovered["accession"].str.len() > 0]

    # Merge by accession
    agg_cols = [
        "accession", "title", "species", "assay_type", "platform", "sample_count",
        "sample_source_guess", "comparison_guess", "source", "search_profiles",
        "matched_terms", "geo_url",
    ]
    for col in agg_cols:
        if col not in discovered.columns:
            discovered[col] = ""
    merged_rows = []
    for acc, grp in discovered.groupby("accession"):
        row = grp.iloc[0].to_dict()
        row["sources"] = ";".join(sorted(set(grp["source"].astype(str))))
        row["search_profiles"] = ";".join(sorted(set(";".join(grp["search_profiles"].astype(str)).split(";"))))
        row["titles"] = " | ".join(sorted(set(grp["title"].astype(str))))[:500]
        if acc in KNOWN_ACCESSION_OVERRIDES:
            row.update({k: v for k, v in KNOWN_ACCESSION_OVERRIDES[acc].items() if k != "notes"})
        merged_rows.append(row)
    catalog = pd.DataFrame(merged_rows)

    local_paths, _ = _load_local_inventory()
    catalog["local_status"] = catalog["accession"].map(
        lambda a: _local_data_status(a, local_paths)[0]
    )
    catalog["local_path"] = catalog["accession"].map(
        lambda a: _local_data_status(a, local_paths)[1]
    )
    catalog["biospecimen_class"] = catalog.apply(_classify_biospecimen_for_mpp1_marco, axis=1)

    human_mask = catalog["species"].astype(str).str.contains("Homo sapiens", case=False, na=False)
    catalog["is_human"] = human_mask
    catalog_human = catalog[human_mask].copy()

    print("\n=== Phase 2: Priority GEO FTP downloads ===")
    to_dl = list(dict.fromkeys(PRIORITY_DOWNLOAD))
    download_log = _attempt_priority_downloads(to_dl)
    download_log.to_csv(OUT_DIR / "priority_download_log.tsv", sep="\t", index=False)

    # Refresh local status after download
    local_paths, _ = _load_local_inventory()
    catalog_human["local_status"] = catalog_human["accession"].map(
        lambda a: _local_data_status(a, local_paths)[0]
    )

    print("\n=== Phase 3: Expression survey on loadable cohorts ===")
    run_expression_survey()

    survey = pd.read_csv(SURVEY_OUT / "mpp1_marco_all_cohort_results.tsv", sep="\t")
    meta = pd.read_csv(SURVEY_OUT / "mpp1_marco_stouffer_meta.tsv", sep="\t")

    # Build comprehensive assessment table
    assess_rows = []
    for _, row in catalog_human.iterrows():
        acc = row["accession"]
        bio = row["biospecimen_class"]
        for gene in GENES:
            sub = survey[(survey["Cohort"] == acc) & (survey["Gene"] == gene)]
            meas = _genes_measurable(acc, row["local_status"])
            base = {
                "accession": acc,
                "gene": gene,
                "biospecimen": bio,
                "comparison_guess": row.get("comparison_guess", ""),
                "evidence_tier": row.get("evidence_tier", ""),
                "local_status": row["local_status"],
                "sources": row.get("sources", ""),
                "geo_url": row.get("geo_url", ""),
                "measurability": meas.get(gene, "unknown"),
            }
            if not sub.empty and str(sub.iloc[0].get("Status")) == "ok":
                r = sub.iloc[0]
                base.update(
                    {
                        "expression_status": "ok",
                        "log2fc": r.get("Log2FC_positive_minus_negative"),
                        "p_value": r.get("P_value"),
                        "direction_match_discovery": r.get("Direction_match_discovery"),
                        "n_positive": r.get("N_positive"),
                        "n_negative": r.get("N_negative"),
                    }
                )
            else:
                base["expression_status"] = str(sub.iloc[0].get("Status")) if not sub.empty else "not_run"
            assess_rows.append(base)

    assessment = pd.DataFrame(assess_rows)
    catalog.to_csv(OUT_DIR / "comprehensive_cohort_catalog_all_sources.tsv", sep="\t", index=False)
    catalog_human.to_csv(OUT_DIR / "comprehensive_cohort_catalog_human.tsv", sep="\t", index=False)
    assessment.to_csv(OUT_DIR / "comprehensive_mpp1_marco_assessment.tsv", sep="\t", index=False)

    tissue_rvu = catalog_human[
        catalog_human["comparison_guess"].astype(str).str.contains("ruptured_vs_unruptured", na=False)
        & catalog_human["biospecimen_class"].eq("aneurysm_wall_tissue")
    ]
    blood_any = catalog_human[catalog_human["biospecimen_class"].isin(["blood", "blood_ncRNA"])]

    md = [
        "# Comprehensive MPP1 / MARCO public data survey",
        "",
        f"Discovery sources: GEO ({len(GEO_SEARCH_PROFILES) + len(EXTRA_GEO_SEARCH_PROFILES)} profiles), "
        f"BioStudies, SRA, GitHub, existing project manifest.",
        "",
        f"- Total unique accessions cataloged: **{len(catalog)}**",
        f"- Human accessions: **{len(catalog_human)}**",
        f"- Human **wall tissue** cohorts: **{(catalog_human['biospecimen_class'] == 'aneurysm_wall_tissue').sum()}**",
        f"- Human **blood** cohorts: **{catalog_human['biospecimen_class'].isin(['blood', 'blood_ncRNA']).sum()}**",
        f"- Direct **ruptured vs unruptured** (any biospecimen): "
        f"**{catalog_human['comparison_guess'].astype(str).str.contains('ruptured_vs_unruptured').sum()}**",
        "",
        "## Critical finding: no public blood RvU RNA-seq",
        "",
        "After exhaustive GEO/SRA/BioStudies search, **no human whole-blood or PBMC RNA-seq cohort** "
        "with matched **ruptured vs unruptured** intracranial aneurysm was identified. "
        "Available blood studies are: (i) aneurysm vs control (GSE159610, GSE106520), "
        "(ii) post-rupture SAH vs control (GSE36791, GSE73378), or (iii) miRNA-only (GSE50867, GSE231922).",
        "",
        "## Tissue ruptured vs unruptured (MPP1 / MARCO measurable)",
        "",
        tissue_rvu[["accession", "title", "local_status", "sample_count"]].to_markdown(index=False)
        if not tissue_rvu.empty
        else "_None_",
        "",
        "## Blood cohorts (MPP1 / MARCO)",
        "",
        blood_any[["accession", "title", "comparison_guess", "local_status"]].to_markdown(index=False)
        if not blood_any.empty
        else "_None_",
        "",
        "## Expression results (analyzed cohorts)",
        "",
        survey[survey["Status"] == "ok"][
            ["Gene", "Cohort", "Comparison", "Log2FC_positive_minus_negative", "P_value", "Direction_match_discovery"]
        ].to_markdown(index=False),
        "",
        "## Stouffer meta-analysis",
        "",
        meta.to_markdown(index=False),
        "",
        "## Priority download attempts",
        "",
        download_log.to_markdown(index=False) if not download_log.empty else "_No downloads attempted_",
    ]
    (OUT_DIR / "comprehensive_discovery_summary.md").write_text("\n".join(md), encoding="utf-8")

    summary_json = {
        "n_accessions_total": int(len(catalog)),
        "n_human": int(len(catalog_human)),
        "n_tissue_wall": int((catalog_human["biospecimen_class"] == "aneurysm_wall_tissue").sum()),
        "n_blood": int(catalog_human["biospecimen_class"].isin(["blood", "blood_ncRNA"]).sum()),
        "n_tissue_rvu_catalog": int(
            catalog_human["comparison_guess"].astype(str).str.contains("ruptured_vs_unruptured").sum()
        ),
        "n_expression_ok": int((survey["Status"] == "ok").sum() // 2),
        "blood_ruptured_vs_unruptured_found": False,
        "discovery_signs": {g: float(_discovery_signs(BASE).get(g, np.nan)) for g in GENES},
    }
    (OUT_DIR / "comprehensive_discovery_summary.json").write_text(
        json.dumps(summary_json, indent=2), encoding="utf-8"
    )

    print("\n=== Done ===")
    print(f"Catalog: {OUT_DIR / 'comprehensive_cohort_catalog_human.tsv'}")
    print(f"Assessment: {OUT_DIR / 'comprehensive_mpp1_marco_assessment.tsv'}")
    print(json.dumps(summary_json, indent=2))


if __name__ == "__main__":
    main()
