"""
Download shortlisted GEO candidate cohorts discovered by the public cohort audit.

Default use:
  python download_geo_candidate_cohorts.py

This reads `results/public_cohort_discovery/public_transcriptomic_cohort_manifest.tsv`
and downloads matrix plus supplementary files for a shortlist of not-yet-local
cohorts. Raw archives are skipped by default to keep the first pass fast.
"""

from __future__ import annotations

import argparse
import json
import re
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urljoin
from urllib.request import Request, urlopen

import pandas as pd


ROOT = Path(__file__).resolve().parent
MANIFEST_PATH = ROOT / "results" / "public_cohort_discovery" / "public_transcriptomic_cohort_manifest.tsv"
SUMMARY_PATH = ROOT / "results" / "public_cohort_discovery" / "downloaded_candidate_cohorts.json"
REQUEST_HEADERS = {"User-Agent": "AneurysmCohortDownloader/1.0 (public cohort acquisition)"}


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


def _read_text_url(url: str) -> str:
    request = Request(url, headers=REQUEST_HEADERS)
    with urlopen(request, timeout=60) as response:
        return response.read().decode("utf-8", errors="ignore")


def _series_prefix(accession: str) -> str:
    match = re.match(r"^(GSE)(\d+)$", accession)
    if not match:
        raise ValueError(f"Unsupported accession format: {accession}")
    digits = match.group(2)
    return f"GSE{digits[:-3]}nnn"


def _list_directory_files(url: str) -> list[str]:
    html = _read_text_url(url)
    parser = LinkCollector()
    parser.feed(html)
    links = []
    for link in parser.links:
        if not link or link.startswith("?") or link.startswith("/"):
            continue
        if link in {"../", "./"}:
            continue
        links.append(urljoin(url, link))
    return sorted(set(links))


def _download(url: str, destination: Path) -> None:
    request = Request(url, headers=REQUEST_HEADERS)
    with urlopen(request, timeout=120) as response, destination.open("wb") as handle:
        handle.write(response.read())


def _choose_output_dir(row: pd.Series) -> Path:
    assay = str(row.get("assay_type", "")).lower()
    sample_source = str(row.get("sample_source_guess", ""))
    accession = str(row["accession"])

    if "array" in assay:
        return ROOT / "data" / "microarray" / accession
    if sample_source == "blood":
        return ROOT / "data" / "blood_rnaseq" / accession
    return ROOT / "data" / "sc_rnaseq" / accession


def _filter_files(urls: list[str], include_raw: bool) -> list[str]:
    selected = []
    for url in urls:
        filename = url.rsplit("/", 1)[-1]
        lower = filename.lower()
        if not include_raw and (lower.endswith(".tar") or "raw" in lower or lower.endswith(".cel.gz")):
            continue
        selected.append(url)
    return selected


def _download_accession(accession: str, row: pd.Series, include_raw: bool) -> dict[str, object]:
    prefix = _series_prefix(accession)
    out_dir = _choose_output_dir(row)
    out_dir.mkdir(parents=True, exist_ok=True)

    matrix_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{accession}/matrix/"
    suppl_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{accession}/suppl/"

    downloaded: list[str] = []
    errors: list[str] = []

    for base_url in (matrix_url, suppl_url):
        try:
            urls = _list_directory_files(base_url)
            urls = _filter_files(urls, include_raw=include_raw)
            for url in urls:
                filename = url.rsplit("/", 1)[-1]
                if filename in {"", "../"}:
                    continue
                destination = out_dir / filename
                if destination.exists():
                    downloaded.append(str(destination.relative_to(ROOT).as_posix()))
                    continue
                _download(url, destination)
                downloaded.append(str(destination.relative_to(ROOT).as_posix()))
        except Exception as exc:
            errors.append(f"{base_url}: {exc}")

    return {
        "accession": accession,
        "output_dir": str(out_dir.relative_to(ROOT).as_posix()),
        "downloaded_files": downloaded,
        "errors": errors,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Download shortlisted GEO candidate cohorts.")
    parser.add_argument(
        "--accessions",
        nargs="*",
        default=["GSE73378", "GSE106520", "GSE231922"],
        help="GEO accessions to download",
    )
    parser.add_argument(
        "--include-raw",
        action="store_true",
        help="Also download large raw archives such as RAW.tar or CEL bundles",
    )
    args = parser.parse_args()

    if not MANIFEST_PATH.is_file():
        raise FileNotFoundError(f"Manifest not found: {MANIFEST_PATH}")

    manifest = pd.read_csv(MANIFEST_PATH, sep="\t")
    manifest = manifest.set_index("accession", drop=False)

    results: list[dict[str, object]] = []
    for accession in args.accessions:
        if accession not in manifest.index:
            results.append(
                {
                    "accession": accession,
                    "output_dir": "",
                    "downloaded_files": [],
                    "errors": ["Accession not present in public cohort manifest"],
                }
            )
            continue
        row = manifest.loc[accession]
        results.append(_download_accession(accession, row, include_raw=args.include_raw))

    SUMMARY_PATH.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"Created: {SUMMARY_PATH}")
    for result in results:
        print(
            f"{result['accession']}: downloaded={len(result['downloaded_files'])} "
            f"errors={len(result['errors'])} output_dir={result['output_dir']}"
        )


if __name__ == "__main__":
    main()
