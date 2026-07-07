"""
Targeted GEO download for legacy completeness cohorts:
  - GSE50867: circulating miRNA (blood pools; ruptured / unruptured / control)
  - GSE6551: tissue microarray (intracranial aneurysm vs intracranial artery control)
"""

from __future__ import annotations

import json
import re
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urljoin
from urllib.request import Request, urlopen

BASE = Path(r"e:\Aneurysms Genomics Project")
OUT_JSON = BASE / "results" / "mpp1_marco_survey" / "gse50867_gse6551_download_log.json"
REQUEST_HEADERS = {"User-Agent": "AneurysmLegacyCohortDownload/1.0"}

COHORTS = {
    "GSE50867": {
        "out_dir": BASE / "data" / "microarray" / "GSE50867",
        "files": [
            "matrix/GSE50867_series_matrix.txt.gz",
            "soft/GSE50867_family.soft.gz",
        ],
        "note": "Non-coding RNA (miRNA) array; MPP1/MARCO mRNA not expected on platform.",
    },
    "GSE6551": {
        "out_dir": BASE / "data" / "microarray" / "GSE6551",
        "files": [
            "matrix/GSE6551-GPL570_series_matrix.txt.gz",
            "matrix/GSE6551-GPL2507_series_matrix.txt.gz",
            "soft/GSE6551_family.soft.gz",
        ],
        "note": "Tissue mRNA microarray (GPL570); aneurysm vs intracranial artery.",
    },
}


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


def _download(url: str, destination: Path, timeout: int = 600) -> dict[str, object]:
    if destination.is_file() and destination.stat().st_size > 1000:
        return {
            "url": url,
            "path": str(destination.relative_to(BASE)),
            "status": "already_present",
            "bytes": destination.stat().st_size,
        }
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = Request(url, headers=REQUEST_HEADERS)
    with urlopen(request, timeout=timeout) as response, destination.open("wb") as handle:
        handle.write(response.read())
    size = destination.stat().st_size
    if size < 1000:
        destination.unlink(missing_ok=True)
        return {"url": url, "path": str(destination), "status": "failed_too_small", "bytes": size}
    return {
        "url": url,
        "path": str(destination.relative_to(BASE)),
        "status": "downloaded",
        "bytes": size,
    }


def main() -> None:
    results: dict[str, object] = {}
    for accession, cfg in COHORTS.items():
        prefix = _series_prefix(accession)
        base_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{accession}/"
        out_dir: Path = cfg["out_dir"]
        file_results = []
        for rel in cfg["files"]:
            url = urljoin(base_url, rel)
            dest = out_dir / rel.rsplit("/", 1)[-1]
            print(f"Downloading {accession}: {dest.name}...")
            try:
                file_results.append(_download(url, dest))
                print(f"  -> {file_results[-1]['status']} ({file_results[-1].get('bytes', 0)} bytes)")
            except Exception as exc:
                file_results.append({"url": url, "path": str(dest), "status": "error", "error": str(exc)})
                print(f"  !! {exc}")
        results[accession] = {"note": cfg["note"], "files": file_results}

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"\nWrote {OUT_JSON}")


if __name__ == "__main__":
    main()
