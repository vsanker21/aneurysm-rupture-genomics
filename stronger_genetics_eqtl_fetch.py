"""
Fetch blood cis-eQTL association rows from the eQTL Catalogue (v3 per-dataset API)
for each locus in stronger_genetics loci_manifest.json, restricted to the genomic
window already used for regional IA GWAS export.

Outputs: results/phases_g/stronger_genetics/eqtl_region_{locus_id}.tsv

These files are inputs for R coloc / SMR after allele harmonization, and for the
Python z-score correlation screen (stronger_genetics_zscreen.py).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import pandas as pd
import requests

BASE = Path(__file__).resolve().parent
DEFAULT_MANIFEST = BASE / "results/phases_g/stronger_genetics/loci_manifest.json"
DEFAULT_BLOOD_IDS = BASE / "data/external/eqtl_catalogue/blood_ge_dataset_ids.txt"
DEFAULT_TIER1_SCIENCE = BASE / "results/phases_g/tier1_science_supplement.tsv"
OUT_DIR = BASE / "results/phases_g/stronger_genetics"
EQTL_DATASET_ASSOCS = "https://www.ebi.ac.uk/eqtl/api/v3/datasets/{dataset_id}/associations"
# Gene-level paging on a single QTD often returns [] for genes not quantified in that study.
# Global search (gene + chromosome) finds cis signals from other datasets; we prefer blood QTDs when present.
EQTL_GLOBAL_ASSOCS = "https://www.ebi.ac.uk/eqtl/api/v3/associations"


def norm_chr(v) -> str:
    s = str(v).strip().upper().replace("CHR", "")
    if s in ("23", "X"):
        return "X"
    try:
        return str(int(float(s)))
    except ValueError:
        return s


def eqtl_get(url: str, params: dict, timeout: int, retries: int = 5) -> requests.Response:
    last_err: BaseException | None = None
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=timeout)
            r.raise_for_status()
            return r
        except (requests.Timeout, requests.ConnectionError) as e:
            last_err = e
            time.sleep(2.5 * (attempt + 1))
        except requests.HTTPError as e:
            code = e.response.status_code if e.response is not None else 0
            if code in (500, 502, 503, 504):
                last_err = e
                time.sleep(3.0 * (attempt + 1))
            else:
                raise
    raise last_err  # type: ignore[misc]


def load_blood_ids(path: Path) -> list[str]:
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]


def _ensembl_chr(chrom: str) -> str:
    s = norm_chr(chrom)
    return "X" if s == "X" else s


def lift_grch37_interval_to_grch38(chrom: str, lo: int, hi: int, timeout: int) -> tuple[int, int]:
    """Map inclusive GRCh37 coordinates to GRCh38 for eQTL Catalogue (GRCh38)."""
    ch = _ensembl_chr(chrom)
    url = f"https://rest.ensembl.org/map/human/GRCh37/{ch}:{int(lo)}..{int(hi)}:1/GRCh38"
    headers = {"Accept": "application/json", "User-Agent": "stronger_genetics_eqtl_fetch/1.0"}
    r = requests.get(url, headers=headers, timeout=timeout)
    r.raise_for_status()
    data = r.json()
    maps = data.get("mappings") or []
    if not maps:
        return int(lo), int(hi)
    starts: list[int] = []
    ends: list[int] = []
    for m in maps:
        mp = m.get("mapped") or {}
        if str(mp.get("seq_region_name", "")).replace("chr", "") != ch:
            continue
        starts.append(int(mp["start"]))
        ends.append(int(mp["end"]))
    if not starts:
        return int(lo), int(hi)
    return min(starts), max(ends)


def fetch_locus_eqtl(
    ensembl_gene: str,
    chrom: str,
    lo: int,
    hi: int,
    dataset_ids: list[str],
    max_pages: int,
    timeout: int,
) -> pd.DataFrame:
    """Prefer global gene+chromosome paging first (fast, works when blood QTDs omit the gene)."""
    cwant = norm_chr(chrom)
    df_g = fetch_locus_eqtl_global_fallback(
        ensembl_gene,
        chrom,
        lo,
        hi,
        set(dataset_ids),
        max_pages,
        timeout,
    )
    if len(df_g) > 0:
        return df_g

    rows = []
    for ds_id in sorted(dataset_ids):
        start = 0
        dead_pages = 0
        for _ in range(max(1, max_pages)):
            try:
                r = eqtl_get(
                    EQTL_DATASET_ASSOCS.format(dataset_id=ds_id),
                    {
                        "gene_id": ensembl_gene.split(".")[0],
                        "size": 2500,
                        "start": start,
                    },
                    timeout=timeout,
                )
            except requests.RequestException:
                break
            chunk = r.json() or []
            matched_here = 0
            for a in chunk:
                ch = norm_chr(a.get("chromosome", ""))
                pos = a.get("position")
                if pos is None:
                    continue
                pos = int(pos)
                if ch != cwant or not (lo <= pos <= hi):
                    continue
                matched_here += 1
                rows.append(
                    {
                        "rsid": a.get("rsid"),
                        "chromosome": ch,
                        "position": pos,
                        "pvalue": a.get("pvalue"),
                        "beta": a.get("beta"),
                        "se": a.get("se"),
                        "variant": a.get("variant"),
                        "dataset_id": a.get("dataset_id"),
                        "study_id": a.get("study_id"),
                    }
                )
            time.sleep(0.12)
            if len(chunk) < 2500:
                break
            if matched_here == 0:
                dead_pages += 1
                if dead_pages >= 2:
                    break
            else:
                dead_pages = 0
            start += 2500
    if rows:
        df = pd.DataFrame(rows)
        return df.drop_duplicates(subset=["rsid", "position", "dataset_id"], keep="first")
    return pd.DataFrame()


def fetch_locus_eqtl_global_fallback(
    ensembl_gene: str,
    chrom: str,
    lo: int,
    hi: int,
    blood_dataset_ids: set[str],
    max_pages: int,
    timeout: int,
    http_retries: int = 2,
) -> pd.DataFrame:
    cwant = norm_chr(chrom)
    gene = ensembl_gene.split(".")[0]
    rows: list[dict] = []
    start = 0
    # Bound total API calls: chr-wide streams can be long; cis window is usually hit within few pages.
    cap = min(16, max(max_pages * 2, 8))
    for _ in range(cap):
        try:
            r = eqtl_get(
                EQTL_GLOBAL_ASSOCS,
                {
                    "gene_id": gene,
                    "chromosome": cwant,
                    "size": 2500,
                    "start": start,
                },
                timeout=timeout,
                retries=http_retries,
            )
        except requests.RequestException:
            break
        chunk = r.json() or []
        chunk_pos: list[int] = []
        for a in chunk:
            ch = norm_chr(a.get("chromosome", ""))
            pos = a.get("position")
            if pos is None:
                continue
            pos = int(pos)
            if ch == cwant:
                chunk_pos.append(pos)
            if ch != cwant or not (lo <= pos <= hi):
                continue
            rows.append(
                {
                    "rsid": a.get("rsid"),
                    "chromosome": ch,
                    "position": pos,
                    "pvalue": a.get("pvalue"),
                    "beta": a.get("beta"),
                    "se": a.get("se"),
                    "variant": a.get("variant"),
                    "dataset_id": a.get("dataset_id"),
                    "study_id": a.get("study_id"),
                }
            )
        time.sleep(0.12)
        if len(rows) >= 4000:
            break
        if len(chunk) < 2500:
            break
        if chunk_pos and min(chunk_pos) > hi:
            break
        start += 2500
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    df = df.drop_duplicates(subset=["rsid", "position", "dataset_id"], keep="first")
    in_blood = df[df["dataset_id"].isin(blood_dataset_ids)]
    if len(in_blood) > 0:
        return in_blood
    return df


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    ap.add_argument("--blood-ids", type=Path, default=DEFAULT_BLOOD_IDS)
    ap.add_argument("--max-loci", type=int, default=40, help="Cap loci to fetch (API time).")
    ap.add_argument(
        "--tier1-tsv",
        type=Path,
        default=None,
        help="If provided, only fetch eQTL regions for tier-1 loci in this table (in rank order).",
    )
    ap.add_argument(
        "--tier1-top-n",
        type=int,
        default=10,
        help="When --tier1-tsv is used, number of tier-1 rows to fetch.",
    )
    ap.add_argument("--max-pages-per-dataset", type=int, default=6)
    ap.add_argument("--timeout", type=int, default=240)
    ap.add_argument(
        "--only-locus-ids",
        type=str,
        default="",
        help="Comma-separated locus_id values (e.g. FN3K_rs1046896). Ignores --max-loci cap for these.",
    )
    ap.add_argument(
        "--expand-window-kb",
        type=int,
        default=0,
        help="Extend region_start/end by this many kb on each side (refetch / dense loci). "
        "Also run reslice_gwas_region.py with the same kb so GWAS and eQTL windows match for COLOC.",
    )
    ap.add_argument(
        "--no-clear-eqtl-files",
        action="store_true",
        help="Do not delete existing eqtl_region_*.tsv before fetch (for targeted refetch).",
    )
    ap.add_argument(
        "--no-grch38-lift",
        action="store_true",
        help="Do not lift manifest region_start/end from GRCh37 to GRCh38 before Catalogue filtering "
        "(default: lift; Catalogue summary stats are GRCh38).",
    )
    args = ap.parse_args()

    if not args.manifest.is_file():
        print(f"Missing manifest. Run: prepare_stronger_genetics_followup.py", file=sys.stderr)
        return 1
    if not args.blood_ids.is_file():
        print(f"Missing {args.blood_ids}", file=sys.stderr)
        return 1

    loci = json.loads(args.manifest.read_text(encoding="utf-8"))
    ds_ids = load_blood_ids(args.blood_ids)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    only_ids: set[str] | None = None
    if args.only_locus_ids.strip():
        only_ids = {x.strip() for x in args.only_locus_ids.split(",") if x.strip()}

    if not args.no_clear_eqtl_files:
        for stale in OUT_DIR.glob("eqtl_region_*.tsv"):
            try:
                stale.unlink()
            except OSError:
                pass

    tier1_ids: list[str] | None = None
    if args.tier1_tsv:
        if not args.tier1_tsv.is_file():
            print(f"Missing tier1 TSV: {args.tier1_tsv}", file=sys.stderr)
            return 1
        t1 = pd.read_csv(args.tier1_tsv, sep="\t")
        if "symbol" not in t1.columns or "blood_eqtl_top_rsid" not in t1.columns:
            print("tier1 TSV missing required columns: symbol, blood_eqtl_top_rsid", file=sys.stderr)
            return 1
        t1 = t1.head(int(args.tier1_top_n)).copy()
        t1["locus_id"] = (
            t1["symbol"].astype(str).str.strip()
            + "_"
            + t1["blood_eqtl_top_rsid"].astype(str).str.strip()
        )
        tier1_ids = [str(x) for x in t1["locus_id"].tolist() if str(x)]

    done = 0
    summary = []
    extra = int(args.expand_window_kb) * 1000
    for L in loci:
        if "gwas_region_file" not in L or "region_start" not in L:
            continue
        lid = L["locus_id"]
        if only_ids is not None:
            if lid not in only_ids:
                continue
        else:
            if tier1_ids is not None and lid not in tier1_ids:
                continue
            if done >= args.max_loci:
                break
        ens = str(L.get("ensembl_gene", "")).split(".")[0]
        chrom = L["sentinel_chr"]
        lo, hi = int(L["region_start"]), int(L["region_end"])
        if extra > 0:
            lo = max(1, lo - extra)
            hi = hi + extra
        lo_f, hi_f = lo, hi
        if not args.no_grch38_lift:
            try:
                lo_f, hi_f = lift_grch37_interval_to_grch38(str(chrom), lo, hi, min(120, int(args.timeout)))
            except Exception as e:
                print(f"  [lift] failed ({e}); using GRCh37 bounds for filtering", flush=True)
                lo_f, hi_f = lo, hi
            else:
                if (lo_f, hi_f) != (lo, hi):
                    print(
                        f"  [lift] GWAS region GRCh37 chr{chrom}:{lo}-{hi} -> GRCh38 {lo_f}-{hi_f} for Catalogue",
                        flush=True,
                    )
        print(f"[eQTL] {lid} {ens} chr{chrom}:{lo_f}-{hi_f}", flush=True)
        df = fetch_locus_eqtl(
            ens,
            chrom,
            lo_f,
            hi_f,
            ds_ids,
            args.max_pages_per_dataset,
            args.timeout,
        )
        out = OUT_DIR / f"eqtl_region_{lid}.tsv"
        df.to_csv(out, sep="\t", index=False)
        summary.append(
            {
                "locus_id": lid,
                "symbol": L.get("symbol"),
                "n_eqtl_rows": int(len(df)),
                "eqtl_file": str(out.resolve()),
            }
        )
        done += 1
        if only_ids is not None and done >= len(only_ids):
            break
        if tier1_ids is not None and done >= len(tier1_ids):
            break

    sm_path = OUT_DIR / "eqtl_fetch_summary.json"
    if args.no_clear_eqtl_files and sm_path.is_file():
        try:
            prev = json.loads(sm_path.read_text(encoding="utf-8"))
            by_id = {x["locus_id"]: x for x in prev if isinstance(x, dict) and "locus_id" in x}
            for x in summary:
                by_id[x["locus_id"]] = x
            summary = sorted(by_id.values(), key=lambda z: z["locus_id"])
        except Exception:
            pass
    sm_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"Wrote {sm_path} ({len(summary)} loci)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
