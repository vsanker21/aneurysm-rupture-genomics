"""
Programmatic downloads for PHASES-G evidence layer (hypothesis generation).

1) ISGC / Bakker et al. 2020 intracranial aneurysm GWAS summary statistics (Figshare).
   Default: Stage_2 meta-analysis (European + East Asian; ruptured + unruptured IA).

2) eQTL Catalogue index: blood + gene-expression (ge) dataset IDs for cis-eQTL queries
   (GTEx blood, Lepik_2017 blood, TwinsUK blood). This proxies constitutive blood
   cis regulation; it is not the full eQTLGen 3.6 GB file — add that separately if
   you require eQTLGen Phase I verbatim.

Outputs:
  data/gwas/bakker_2020/downloaded/IA.GWAS.BakkerMK.2020.sumstats.Stage_2.txt.gz
  data/external/eqtl_catalogue/blood_ge_dataset_ids.txt
  data/sc_rnaseq/GSE122897/GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz  (+ raw counts)
  results/phases_g/download_manifest.json

Usage:
  python download_phases_g_public_data.py
  python download_phases_g_public_data.py --skip-gwas   # only refresh eQTL index
  python download_phases_g_public_data.py --verify-only
"""
from __future__ import annotations

import gzip
import hashlib
import json
import sys
import time
from pathlib import Path

import requests

try:
    from project_paths import get_project_root
    BASE = get_project_root()
except ImportError:
    BASE = Path(__file__).resolve().parent
GWAS_DIR = BASE / "data/gwas/bakker_2020/downloaded"
EQTL_DIR = BASE / "data/external/eqtl_catalogue"
MANIFEST_DIR = BASE / "results/phases_g"
GSE122897_DIR = BASE / "data/sc_rnaseq/GSE122897"
MARKERS_TSV = BASE / "data/reference/scrna_aneurysm_celltype_markers.tsv"

GSE122897_SUPPL_BASE = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE122nnn/GSE122897/suppl"
)
GSE122897_FILES = (
    "GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz",
    "GSE122897_raw_counts_GRCh38.p13_NCBI.tsv.gz",
)
# Relaxed min sizes (bytes) — catch truncated downloads
MIN_BYTES_GWAS_STAGE2 = 20_000_000
MIN_BYTES_TPM = 2_000_000
MIN_BYTES_RAW = 1_000_000

FIGSHARE_ARTICLE = "11303372"
STAGE2_NAME = "IA.GWAS.BakkerMK.2020.sumstats.Stage_2.txt.gz"
FIGSHARE_FILES_API = f"https://api.figshare.com/v2/articles/{FIGSHARE_ARTICLE}/files"

EQTL_V2_DATASETS = "https://www.ebi.ac.uk/eqtl/api/v2/datasets"


def stream_download(
    url: Path | str,
    dest: Path,
    chunk: int = 1 << 20,
    timeout: int = 600,
    retries: int = 4,
) -> str:
    dest.parent.mkdir(parents=True, exist_ok=True)
    last_err: BaseException | None = None
    for attempt in range(retries):
        try:
            h = hashlib.sha256()
            with requests.get(str(url), stream=True, timeout=timeout) as r:
                r.raise_for_status()
                with open(dest, "wb") as f:
                    for blk in r.iter_content(chunk):
                        if blk:
                            f.write(blk)
                            h.update(blk)
            return h.hexdigest()
        except (requests.Timeout, requests.ConnectionError) as e:
            last_err = e
            if dest.exists():
                try:
                    dest.unlink()
                except OSError:
                    pass
            time.sleep(3.0 * (attempt + 1))
        except requests.HTTPError as e:
            code = e.response.status_code if e.response is not None else 0
            if code in (502, 503, 504) and attempt + 1 < retries:
                last_err = e
                time.sleep(3.0 * (attempt + 1))
            else:
                raise
    raise last_err  # type: ignore[misc]


def download_geo_suppl_file(rel_name: str, dest: Path) -> dict:
    url = f"{GSE122897_SUPPL_BASE}/{rel_name}"
    print(f"Downloading GEO {rel_name} ...")
    sha = stream_download(url, dest, timeout=900, retries=5)
    sz = dest.stat().st_size
    print(f"  Wrote {dest} ({sz / 1e6:.2f} MB) sha256={sha[:16]}...")
    return {"file": str(dest.resolve()), "url": url, "bytes": sz, "sha256": sha}


def download_gse122897_expression_files() -> dict:
    GSE122897_DIR.mkdir(parents=True, exist_ok=True)
    out = {}
    mins = (MIN_BYTES_TPM, MIN_BYTES_RAW)
    for name, min_b in zip(GSE122897_FILES, mins):
        dest = GSE122897_DIR / name
        if dest.is_file() and dest.stat().st_size >= min_b:
            print(f"[skip] GSE122897 {name} already present ({dest.stat().st_size} bytes)")
            out[name] = {"file": str(dest.resolve()), "note": "already present"}
            continue
        out[name] = download_geo_suppl_file(name, dest)
    return out


def _gzip_readable(path: Path) -> None:
    with gzip.open(path, "rb") as gz:
        if not gz.read(4096):
            raise ValueError(f"Empty or unreadable gzip: {path}")


def verify_public_phases_g_inputs() -> list[str]:
    """Return list of error messages; empty means all required public inputs OK."""
    err: list[str] = []
    gwas = GWAS_DIR / STAGE2_NAME
    if not gwas.is_file():
        err.append(f"Missing GWAS: {gwas}")
    elif gwas.stat().st_size < MIN_BYTES_GWAS_STAGE2:
        err.append(
            f"GWAS file too small ({gwas.stat().st_size} B); expected >= {MIN_BYTES_GWAS_STAGE2}"
        )
    else:
        try:
            _gzip_readable(gwas)
        except (OSError, ValueError, EOFError) as e:
            err.append(f"GWAS gzip check failed: {gwas} ({e})")

    eqtl_ids = EQTL_DIR / "blood_ge_dataset_ids.txt"
    if not eqtl_ids.is_file():
        err.append(f"Missing eQTL dataset index: {eqtl_ids}")
    else:
        lines = [ln for ln in eqtl_ids.read_text().splitlines() if ln.strip()]
        if len(lines) < 1:
            err.append("blood_ge_dataset_ids.txt has no QTD lines")

    tpm = GSE122897_DIR / GSE122897_FILES[0]
    raw = GSE122897_DIR / GSE122897_FILES[1]
    if not tpm.is_file() or tpm.stat().st_size < MIN_BYTES_TPM:
        err.append(f"Missing or too small TPM matrix: {tpm}")
    else:
        try:
            _gzip_readable(tpm)
        except (OSError, ValueError, EOFError) as e:
            err.append(f"TPM gzip check failed: {tpm} ({e})")

    if not raw.is_file() or raw.stat().st_size < MIN_BYTES_RAW:
        err.append(f"Missing or too small raw counts: {raw}")
    else:
        try:
            _gzip_readable(raw)
        except (OSError, ValueError, EOFError) as e:
            err.append(f"Raw counts gzip check failed: {raw} ({e})")

    if not MARKERS_TSV.is_file():
        err.append(f"Missing cell-type marker TSV (ship with repo): {MARKERS_TSV}")

    return err


def download_bakker_stage2() -> dict:
    GWAS_DIR.mkdir(parents=True, exist_ok=True)
    dest = GWAS_DIR / STAGE2_NAME
    r = requests.get(FIGSHARE_FILES_API, timeout=60)
    r.raise_for_status()
    files = r.json()
    entry = next((f for f in files if f.get("name") == STAGE2_NAME), None)
    if not entry:
        raise RuntimeError(f"Figshare file {STAGE2_NAME} not found in article listing")
    url = entry["download_url"]
    size_exp = entry.get("size")
    print(f"Downloading GWAS ({STAGE2_NAME}, ~{size_exp / 1e6:.1f} MB nominal)...")
    sha = stream_download(url, dest, timeout=900, retries=5)
    meta = {
        "file": str(dest.resolve()),
        "figshare_download_url": url,
        "bytes_on_disk": dest.stat().st_size,
        "sha256": sha,
        "citation": "Bakker et al. Nat Genet 2020; Figshare 11303372",
    }
    print(f"  Wrote {dest} sha256={sha[:16]}...")
    return meta


def build_blood_ge_dataset_index(include_lineage_ge_contexts: bool = False) -> dict:
    """
    Index GE (gene expression) cis-eQTL datasets for whole blood (default).
    Optional lineage contexts (monocyte, neutrophil, LCL, etc.) add recall for
    hypothesis generation when whole-blood GE is sparse in the Catalogue.

    The evidence pipeline queries each dataset with gene_id+dataset_id (v3 API);
    it does not rely on paginated global /associations?gene_id= alone.
    """
    EQTL_DIR.mkdir(parents=True, exist_ok=True)
    out_ids = EQTL_DIR / "blood_ge_dataset_ids.txt"
    out_meta = EQTL_DIR / "blood_ge_datasets.json"
    rows = []
    lineage_keys = (
        "monocyte",
        "neutrophil",
        "t cell",
        "cd4+ t cell",
        "cd4",
        "t-cell",
        "lcl",
        "lymphoblastoid",
    )
    for page in range(0, 40):
        r = requests.get(
            EQTL_V2_DATASETS,
            params={"size": 200, "start": page * 200},
            timeout=60,
        )
        if r.status_code != 200:
            break
        batch = r.json()
        if not batch:
            break
        for x in batch:
            if (x.get("quant_method") or "").lower() != "ge":
                continue
            tl = str(x.get("tissue_label", "")).lower()
            if tl == "blood":
                ok = True
            elif include_lineage_ge_contexts and any(k in tl for k in lineage_keys):
                ok = True
            else:
                ok = False
            if not ok:
                continue
            rows.append(
                {
                    "dataset_id": x["dataset_id"],
                    "study_id": x.get("study_id"),
                    "study_label": x.get("study_label"),
                    "tissue_label": x.get("tissue_label"),
                    "sample_size": x.get("sample_size"),
                }
            )
        time.sleep(0.2)
    ids = sorted({x["dataset_id"] for x in rows})
    out_ids.write_text("\n".join(ids) + "\n", encoding="utf-8")
    out_meta.write_text(json.dumps(rows, indent=2), encoding="utf-8")
    print(f"  Wrote {out_ids} ({len(ids)} GE datasets)")
    return {"blood_ge_dataset_ids_file": str(out_ids.resolve()), "n_datasets": len(ids)}


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-gwas", action="store_true")
    ap.add_argument("--skip-eqtl-index", action="store_true")
    ap.add_argument(
        "--skip-gse122897-expression",
        action="store_true",
        help="Do not fetch GSE122897 TPM/raw GEO supplemental matrices.",
    )
    ap.add_argument(
        "--include-lineage-ge-contexts",
        action="store_true",
        help="Also index monocyte/neutrophil/T-cell/LCL GE datasets (more API calls in evidence step).",
    )
    ap.add_argument(
        "--verify-only",
        action="store_true",
        help="Skip downloads; verify required public files and exit 1 if any check fails.",
    )
    ap.add_argument(
        "--no-verify",
        action="store_true",
        help="Do not run post-download verification (not recommended).",
    )
    args = ap.parse_args()

    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    manifest: dict = {"timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}

    if args.verify_only:
        manifest["mode"] = "verify_only"
    else:
        if not args.skip_gwas:
            dest = GWAS_DIR / STAGE2_NAME
            if dest.exists() and dest.stat().st_size >= MIN_BYTES_GWAS_STAGE2:
                print(f"[skip] GWAS exists: {dest}")
                manifest["gwas"] = {"file": str(dest.resolve()), "note": "already present"}
            else:
                manifest["gwas"] = download_bakker_stage2()
        else:
            manifest["gwas"] = {"skipped": True}

        if not args.skip_eqtl_index:
            manifest["eqtl_catalogue_blood_ge"] = build_blood_ge_dataset_index(
                include_lineage_ge_contexts=args.include_lineage_ge_contexts
            )
        else:
            manifest["eqtl_catalogue_blood_ge"] = {"skipped": True}

        if not args.skip_gse122897_expression:
            manifest["gse122897_geo_suppl"] = download_gse122897_expression_files()
        else:
            manifest["gse122897_geo_suppl"] = {"skipped": True}

    do_verify = args.verify_only or not args.no_verify
    if do_verify:
        problems = verify_public_phases_g_inputs()
        manifest["verification_errors"] = problems
        manifest["verification_ok"] = len(problems) == 0
        if problems:
            print("[verify] FAILED:", file=sys.stderr)
            for p in problems:
                print(f"  - {p}", file=sys.stderr)
            mp = MANIFEST_DIR / "download_manifest.json"
            mp.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
            sys.exit(1)
        print("[verify] All required public inputs present and readable.")

    mp = MANIFEST_DIR / "download_manifest.json"
    mp.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Manifest: {mp}")


if __name__ == "__main__":
    try:
        main()
    except requests.RequestException as e:
        print(f"[error] HTTP: {e}", file=sys.stderr)
        sys.exit(1)
    except (RuntimeError, OSError, ValueError) as e:
        print(f"[error] {e}", file=sys.stderr)
        sys.exit(1)
