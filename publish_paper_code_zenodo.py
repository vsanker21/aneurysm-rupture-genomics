"""Create Zenodo v1.1.1 software record (new version of 23676142)."""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import requests

BASE = Path(r"E:\Aneurysms Genomics Project")
ZIP_PATH = BASE / "aneurysm-rupture-genomics-v1.1.1.zip"
API = "https://zenodo.org/api"
PARENT_RECORD = "23676142"
TOKEN_CANDIDATES = [
    BASE / "Latest Submission Files" / "Aneurysm_Genomics_Token.txt",
    BASE / "Nature Files" / "Zenodo.txt",
]


def load_token() -> str:
    token = os.environ.get("ZENODO_ACCESS_TOKEN", "").strip()
    if token:
        return token
    for path in TOKEN_CANDIDATES:
        if path.is_file():
            token = path.read_text(encoding="utf-8").strip().splitlines()[0].strip()
            if token:
                return token
    raise SystemExit("No Zenodo token found")


def software_metadata() -> dict:
    return {
        "title": (
            "Analysis code v1.1.1: concordance-first transcriptomic screening "
            "for intracranial aneurysm rupture (GSE122897 21 vs 21)"
        ),
        "upload_type": "software",
        "description": (
            "Version 1.1.1 of the paper analysis code. GSE122897 Entrez IDs are mapped "
            "to HGNC symbols so the three-cohort Stouffer meta-analysis joins the full "
            "overlap with GSE54083 and GSE13353. BH-FDR is applied across that overlap. "
            "The six-gene wall core is TGFBI, ANPEP, MARCO, MPP1, CLEC5A, and CTSL "
            "(concordant in all three cohorts, meta FDR < 1e-6). Primary GSE122897 "
            "contrast remains 21 ruptured vs 21 unruptured. Companion processed tables: "
            "https://doi.org/10.5281/zenodo.21249929. Source repository: "
            "https://github.com/vsanker21/aneurysm-rupture-genomics (release tag v1.1.1)."
        ),
        "creators": [{"name": "Sanker, Vivek"}],
        "keywords": [
            "intracranial aneurysm",
            "RNA-seq",
            "analysis code",
            "GSE122897",
            "gene panel",
        ],
        "version": "1.1.1",
        "license": "cc-by-4.0",
        "related_identifiers": [
            {
                "identifier": "10.5281/zenodo.21249929",
                "relation": "isSupplementTo",
                "scheme": "doi",
                "resource_type": "dataset",
            },
            {
                "identifier": "https://github.com/vsanker21/aneurysm-rupture-genomics/releases/tag/v1.1.1",
                "relation": "isIdenticalTo",
                "scheme": "url",
                "resource_type": "software",
            },
            {
                "identifier": "10.5281/zenodo.23676142",
                "relation": "isNewVersionOf",
                "scheme": "doi",
                "resource_type": "software",
            },
        ],
    }


def headers(token: str) -> dict:
    return {"Authorization": f"Bearer {token}"}


def new_version(token: str) -> dict:
    r = requests.post(
        f"{API}/deposit/depositions/{PARENT_RECORD}/actions/newversion",
        headers=headers(token),
        timeout=120,
    )
    if r.status_code >= 400:
        print("newversion status", r.status_code, r.text[:500])
        r = requests.post(
            f"{API}/deposit/depositions",
            json={"metadata": software_metadata()},
            headers=headers(token),
            timeout=120,
        )
        r.raise_for_status()
        return r.json()
    dep = r.json()
    draft_url = dep.get("links", {}).get("latest_draft")
    if draft_url:
        d = requests.get(draft_url, headers=headers(token), timeout=120)
        d.raise_for_status()
        return d.json()
    return dep


def replace_files(dep: dict, token: str) -> None:
    for f in list(dep.get("files") or []):
        file_id = f.get("id")
        if not file_id:
            continue
        requests.delete(
            f"{API}/deposit/depositions/{dep['id']}/files/{file_id}",
            headers=headers(token),
            timeout=120,
        )
    bucket = dep["links"]["bucket"]
    with ZIP_PATH.open("rb") as fp:
        up = requests.put(
            f"{bucket}/{ZIP_PATH.name}",
            data=fp,
            headers=headers(token),
            timeout=600,
        )
    up.raise_for_status()


def main() -> None:
    if not ZIP_PATH.is_file():
        raise SystemExit(f"missing {ZIP_PATH}")
    token = load_token()
    dep = new_version(token)
    dep_id = dep["id"]
    print("draft_id", dep_id)
    md = requests.put(
        f"{API}/deposit/depositions/{dep_id}",
        json={"metadata": software_metadata()},
        headers=headers(token),
        timeout=120,
    )
    md.raise_for_status()
    dep = md.json()
    replace_files(dep, token)
    pub = requests.post(
        f"{API}/deposit/depositions/{dep_id}/actions/publish",
        headers=headers(token),
        timeout=180,
    )
    pub.raise_for_status()
    out = pub.json()
    doi = out.get("doi") or (out.get("metadata") or {}).get("doi")
    rec = out.get("links", {}).get("record_html") or out.get("links", {}).get("html")
    print("DOI", doi)
    print("record", rec)
    print("id", out.get("id"))
    Path(r"E:\Aneurysms Genomics Project\metadata\ZENODO_CODE_DOI_v1.1.1.txt").write_text(
        f"DOI: https://doi.org/{doi}\nRecord: {rec}\nDeposition: {out.get('id')}\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    try:
        main()
    except requests.HTTPError as exc:
        print("HTTPError", exc.response.status_code if exc.response is not None else "", file=sys.stderr)
        if exc.response is not None:
            print(exc.response.text[:1000], file=sys.stderr)
        raise
