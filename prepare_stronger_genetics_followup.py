"""
Prepare inputs and a runbook for stronger genetics on a *small* gene set.

This does NOT run COLOC, SuSiE, or SMR (those are typically R: coloc, coloc.susie, SMR).
It does:
  1) Read tier-1 supplement or shortlist + optional symbol filter
  2) Resolve sentinel rsIDs against Bakker Stage_2 sumstats (CHR, BP, BETA, SE, P, alleles)
  3) Query Ensembl REST for each Ensembl gene ID → TSS (GRCh38)
  4) Define a single cis window per locus (envelope of TSS±W and sentinel±W, same chromosome)
  5) One streaming pass over Bakker gzip → per-locus regional GWAS TSVs (for coloc.abf / susieR)
  6) Write loci_manifest.json + RUNBOOK_stronger_genetics.txt

eQTLGen Phase I full cis summary stats are large (~GB); the runbook points to
https://eqtlgen.org/cis-eqtls.html and tabix-friendly workflows after local download.

Usage:
  python prepare_stronger_genetics_followup.py --print-runbook
  python prepare_stronger_genetics_followup.py --from-tier1 --max-loci 12 --window-kb 500
  python prepare_stronger_genetics_followup.py --symbols TGFBI,MARCO,CTSL --window-kb 1000
  python prepare_stronger_genetics_followup.py --from-shortlist --no-ensembl --max-loci 300
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
import requests

BASE = Path(__file__).resolve().parent
DEFAULT_TIER1 = BASE / "results/phases_g/tier1_supplement.tsv"
DEFAULT_SHORTLIST = BASE / "results/phases_g/shortlist_evidence.tsv"
DEFAULT_GWAS = (
    BASE / "data/gwas/bakker_2020/downloaded/IA.GWAS.BakkerMK.2020.sumstats.Stage_2.txt.gz"
)
OUT_DIR = BASE / "results/phases_g/stronger_genetics"

ENSEMBL_LOOKUP = "https://rest.ensembl.org/lookup/id/{gene_id}"
# Bakker Stage_2 coordinates are GRCh37/hg19 — use GRCh37 Ensembl for TSS / rs positions when proxying.
ENSEMBL_GRCH37_LOOKUP = "https://grch37.rest.ensembl.org/lookup/id/{gene_id}"
ENSEMBL_GRCH37_VAR = "https://grch37.rest.ensembl.org/variation/human/{rsid}"


def _safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s)[:80]


def norm_chr_label(v) -> str:
    """Normalize chromosome labels between Bakker / Ensembl (X vs 23)."""
    s = str(v).strip().upper().replace("CHR", "")
    if s in ("23", "X"):
        return "X"
    try:
        return str(int(float(s)))
    except ValueError:
        return s


def fetch_tss_grch37(ensembl_gene: str, timeout: int = 60) -> tuple[str, int] | None:
    url = ENSEMBL_GRCH37_LOOKUP.format(gene_id=ensembl_gene.split(".")[0])
    r = requests.get(
        url,
        headers={"Content-Type": "application/json"},
        params={"expand": "0"},
        timeout=timeout,
    )
    if r.status_code != 200:
        return None
    j = r.json()
    chrom = str(j.get("seq_region_name", ""))
    if not chrom:
        return None
    strand = j.get("strand")
    if strand == 1:
        tss = int(j["start"])
    elif strand == -1:
        tss = int(j["end"])
    else:
        tss = int(j.get("start", 0))
    return chrom, tss


def fetch_rs_grch37(rsid: str, timeout: int = 60) -> tuple[str, int] | None:
    rid = str(rsid).strip()
    if not rid.startswith("rs"):
        return None
    url = ENSEMBL_GRCH37_VAR.format(rsid=rid)
    r = requests.get(url, headers={"Content-Type": "application/json"}, timeout=timeout)
    if r.status_code != 200:
        return None
    j = r.json()
    maps = j.get("mappings") or []
    for m in maps:
        an = str(m.get("assembly_name", ""))
        if "GRCh37" in an or "hg19" in an.lower():
            ch = str(m.get("seq_region_name", ""))
            pos = int(m.get("start", m.get("end", 0)))
            if ch and pos:
                return ch, pos
    if maps:
        m = maps[0]
        ch = str(m.get("seq_region_name", ""))
        pos = int(m.get("start", m.get("end", 0)))
        if ch and pos:
            return ch, pos
    return None


def bakker_best_rs_in_region(
    gwas_path: Path, chrom: str, lo: int, hi: int
) -> pd.Series | None:
    """Lowest P SNP with rsID in [lo,hi] on chrom (Bakker coordinates = GRCh37)."""
    cwant = norm_chr_label(chrom)
    chr_aliases = {cwant}
    if cwant == "X":
        chr_aliases.add("23")
    elif cwant == "23":
        chr_aliases.add("X")
    best: pd.Series | None = None
    best_p = float("inf")
    for chunk in pd.read_csv(
        gwas_path, sep="\t", compression="gzip", chunksize=400_000, low_memory=False
    ):
        chunk["_chrn"] = chunk["CHR"].map(norm_chr_label)
        sub = chunk[chunk["_chrn"].isin(chr_aliases)]
        if sub.empty:
            continue
        bp = sub["BP"].astype(int)
        rs_ok = sub["SNP"].astype(str).str.startswith("rs", na=False)
        m = (bp >= lo) & (bp <= hi) & rs_ok
        if not m.any():
            continue
        sub2 = sub.loc[m]
        p = pd.to_numeric(sub2["P"], errors="coerce")
        idx = p.idxmin()
        row = sub2.loc[idx]
        pv = float(row["P"])
        if pv < best_p:
            best_p = pv
            best = row
    return best


def fetch_tss_grch38(ensembl_gene: str, timeout: int = 60) -> tuple[str, int] | None:
    url = ENSEMBL_LOOKUP.format(gene_id=ensembl_gene)
    r = requests.get(
        url,
        headers={"Content-Type": "application/json"},
        params={"expand": "0"},
        timeout=timeout,
    )
    if r.status_code != 200:
        return None
    j = r.json()
    chrom = str(j.get("seq_region_name", ""))
    if not chrom:
        return None
    strand = j.get("strand")
    if strand == 1:
        tss = int(j["start"])
    elif strand == -1:
        tss = int(j["end"])
    else:
        tss = int(j.get("start", 0))
    return chrom, tss


def load_bakker_sentinels(path: Path, rsids: set[str]) -> pd.DataFrame:
    need = {r for r in rsids if r and str(r).startswith("rs")}
    hits = []
    for chunk in pd.read_csv(
        path, sep="\t", compression="gzip", chunksize=400_000, low_memory=False
    ):
        sub = chunk[chunk["SNP"].astype(str).isin(need)]
        if len(sub):
            hits.append(sub)
    if not hits:
        return pd.DataFrame()
    return pd.concat(hits, ignore_index=True).drop_duplicates(subset=["SNP"])


def print_runbook(path: Path | None = None) -> str:
    text = """
STRONGER GENETICS RUNBOOK (COLOC / SuSiE / SMR / eQTLGen)
========================================================

A) Regional GWAS files produced by this script
   - Use one TSV per locus: columns compatible with R read.table
   - For coloc (SuSiE): harmonize alleles with eQTL file; match genome build (GRCh38).

B) COLOC (R package coloc) — abf
   - Typical: coloc.abf(dataset1=list(beta=..., varbeta=...^2, snp=..., position=..., type="quant"),
                       dataset2=list(...), MAF=optional)
   - Your GWAS is binary IA (Bakker); eQTL is quantitative — use coloc "cc" vs "quant" with care;
     see coloc docs for case-control + QTL (type="cc" for GWAS with s = SE^2).
   - https://chr1swallis.github.io/coloc/

C) coloc-SuSiE / SuSiE fine-mapping
   - Often run SuSiE separately on GWAS and eQTL, then coloc.susie().
   - Requires LD reference (same ancestry) or in-sample LD; public GWAS + blood eQTL need matched LD.

D) SMR + HEIDI (GCTA SMR)
   - Needs GWAS in .ma format and eQTL as SMR-format exposure; HEIDI tests pleiotropy vs linkage.
   - https://yanglab.westlake.edu.cn/software/smr/

E) eQTLGen cis summary statistics
   - Portal: https://eqtlgen.org/cis-eqtls.html
   - Download "full cis-eQTL" (large; Phase I often via https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/ ).
   - Subset to gene or region with zcat | awk, or index with tabix after bgzip for random access.
   - Harmonize variant IDs with Bakker (chr_pos_ref_alt or rsID join).
   - If your blood eQTL sentinel rsID is missing from Bakker, pick an IA GWAS proxy SNP in LD or expand GWAS source.

F) Deeper Catalogue API scan (optional)
   - python run_phases_g_evidence_pipeline.py --max-genes N --eqtl-pages-per-dataset 4
   - Does not replace colocalization.

G) Suggested order
   1) Subset eQTLGen cis to gene / locus
   2) Harmonize + COLOC abf (quick screen)
   3) If PP4 promising: SMR/HEIDI or SuSiE-coloc with appropriate LD

H) Python pre-stack (no R) — same locus windows as above
   - python run_stronger_genetics_stack.py
     (prepare with --no-ensembl, Catalogue fetch, z-score Spearman screen; not coloc.abf)
   - python verify_stronger_genetics_completion.py
     (checks manifest, eqtl_fetch_summary.json, zscreen_correlation.tsv)
   - R template pointer: r/coloc_abf_template.R

I) R COLOC (exploratory) — Bakker-present tier-1 subset
   - python phases_g_coloc_candidates.py  -> results/phases_g/tier1_coloc_candidates.tsv
   - python stronger_genetics_eqtl_fetch.py --max-loci 10 --tier1-tsv results/phases_g/tier1_coloc_candidates.tsv --tier1-top-n 10 --max-pages-per-dataset 6
   - Rscript r/run_coloc_abf_from_regions.R --tier1-tsv results/phases_g/tier1_coloc_candidates.tsv --top-n 10
   - Output: results/phases_g/stronger_genetics/coloc/coloc_summary.tsv

"""
    text = text.strip() + "\n"
    if path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
    return text


def main() -> None:
    ap = argparse.ArgumentParser(description="Prepare stronger-genetics follow-up files.")
    ap.add_argument("--from-tier1", action="store_true", help=f"Use {DEFAULT_TIER1}")
    ap.add_argument("--from-shortlist", action="store_true", help=f"Use {DEFAULT_SHORTLIST}")
    ap.add_argument("--tier1-path", type=Path, default=DEFAULT_TIER1)
    ap.add_argument("--shortlist-path", type=Path, default=DEFAULT_SHORTLIST)
    ap.add_argument("--symbols", type=str, default="", help="Comma-separated HGNC symbols (subset)")
    ap.add_argument("--max-loci", type=int, default=15, help="Max genes with rsID to process")
    ap.add_argument("--window-kb", type=int, default=500, help="Half-window (kb) around TSS and sentinel")
    ap.add_argument("--gwas", type=Path, default=DEFAULT_GWAS)
    ap.add_argument(
        "--max-snps-per-locus",
        type=int,
        default=150_000,
        help="Safety cap per locus (keeps smallest P if exceeded).",
    )
    ap.add_argument("--print-runbook", action="store_true", help="Print runbook and exit")
    ap.add_argument("--sleep-ensembl", type=float, default=0.25, help="Seconds between Ensembl calls")
    ap.add_argument(
        "--no-ensembl",
        action="store_true",
        help="Do not call Ensembl; cis window = sentinel +/- window-kb only (faster bulk prep).",
    )
    ap.add_argument(
        "--bakker-proxy-if-missing",
        action="store_true",
        default=True,
        help="If blood eQTL sentinel rsID is absent from Bakker, use GRCh37 rs coordinates + "
        "lowest-P Bakker SNP in the cis window (not formal LD). Default: on.",
    )
    ap.add_argument(
        "--no-bakker-proxy",
        action="store_true",
        help="Disable Bakker proxy; skip genes whose sentinel is missing from Bakker.",
    )
    ap.add_argument(
        "--merge-existing-manifest",
        action="store_true",
        help=f"Merge new loci into existing {OUT_DIR}/loci_manifest.json by locus_id (keep other loci).",
    )
    args = ap.parse_args()
    use_proxy = args.bakker_proxy_if_missing and not args.no_bakker_proxy

    if args.print_runbook:
        print(print_runbook(OUT_DIR / "RUNBOOK_stronger_genetics.txt"))
        print(f"Also wrote: {OUT_DIR / 'RUNBOOK_stronger_genetics.txt'}")
        return

    sym_filter = {s.strip().upper() for s in args.symbols.split(",") if s.strip()}

    if args.from_tier1:
        df = pd.read_csv(args.tier1_path, sep="\t")
    elif args.from_shortlist:
        df = pd.read_csv(args.shortlist_path, sep="\t")
    else:
        df = pd.read_csv(args.tier1_path, sep="\t")

    if sym_filter:
        df = df[df["symbol"].astype(str).str.upper().isin(sym_filter)]

    allow_no_eqtl_rsid = bool(sym_filter)
    if not allow_no_eqtl_rsid:
        df = df[df["blood_eqtl_top_rsid"].astype(str).str.startswith("rs", na=False)]
    df = df.head(args.max_loci)
    if len(df) == 0:
        print("No rows after filters; check --symbols / shortlist.", file=sys.stderr)
        sys.exit(1)

    rsids = {
        str(x).strip()
        for x in df["blood_eqtl_top_rsid"].tolist()
        if str(x).strip().startswith("rs")
    }
    if not args.gwas.is_file():
        print(f"Missing GWAS: {args.gwas}", file=sys.stderr)
        sys.exit(1)

    print("Loading sentinel rows from Bakker sumstats (where rsID listed)...")
    sent = load_bakker_sentinels(args.gwas, rsids) if rsids else pd.DataFrame()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print_runbook(OUT_DIR / "RUNBOOK_stronger_genetics.txt")

    loci: list[dict] = []
    w = args.window_kb * 1000
    n_skip_no_sentinel = 0

    for _, row in df.iterrows():
        sym = str(row.get("symbol", "")).strip()
        ens = str(row.get("ensembl_gene", "")).split(".")[0]
        rs_eqtl = str(row.get("blood_eqtl_top_rsid", "") or "").strip()
        g: pd.Series | None = None
        s_chr, s_bp = "", 0
        proxy_note = ""
        effective_rs_for_id = rs_eqtl

        if rs_eqtl.startswith("rs"):
            sub = sent[sent["SNP"].astype(str) == rs_eqtl] if not sent.empty else pd.DataFrame()
            if not sub.empty:
                g = sub.iloc[0]
                try:
                    s_chr = str(int(float(g["CHR"]))) if pd.notna(g["CHR"]) else ""
                except (ValueError, TypeError):
                    s_chr = norm_chr_label(g["CHR"]) if pd.notna(g["CHR"]) else ""
                if s_chr == "23":
                    s_chr = "X"
                s_bp = int(g["BP"])
            elif use_proxy:
                print(f"  [proxy] {sym}: catalogue sentinel {rs_eqtl} not in Bakker; resolving GRCh37 position...")
                pos = fetch_rs_grch37(rs_eqtl)
                time.sleep(args.sleep_ensembl)
                if pos is None:
                    n_skip_no_sentinel += 1
                    print(f"  [skip] {sym}: cannot resolve {rs_eqtl} via Ensembl GRCh37")
                    continue
                s_chr, s_bp = pos
                if norm_chr_label(s_chr) == "23":
                    s_chr = "X"
                lo0 = max(1, s_bp - w)
                hi0 = s_bp + w
                print(f"  [proxy] {sym}: scanning Bakker chr{s_chr}:{lo0}-{hi0} for min-P rsID...")
                g = bakker_best_rs_in_region(args.gwas, s_chr, lo0, hi0)
                if g is None:
                    n_skip_no_sentinel += 1
                    if norm_chr_label(s_chr) == "X":
                        print(
                            f"  [skip] {sym}: no Bakker SNPs on chr X in window "
                            "(Bakker Stage_2 summary stats used here appear to exclude chromosome X)."
                        )
                    else:
                        print(f"  [skip] {sym}: no Bakker SNPs in proxy window")
                    continue
                try:
                    s_chr = str(int(float(g["CHR"]))) if pd.notna(g["CHR"]) else norm_chr_label(g["CHR"])
                except (ValueError, TypeError):
                    s_chr = norm_chr_label(g["CHR"]) if pd.notna(g["CHR"]) else s_chr
                if s_chr == "23":
                    s_chr = "X"
                s_bp = int(g["BP"])
                effective_rs_for_id = str(g["SNP"])
                proxy_note = f"bakker_minP_proxy_for_eqtl_sentinel_{rs_eqtl}"
            else:
                n_skip_no_sentinel += 1
                print(f"  [skip] {sym}: sentinel {rs_eqtl} not in Bakker (omit --no-bakker-proxy to use regional min-P proxy)")
                continue
        else:
            if args.no_ensembl:
                n_skip_no_sentinel += 1
                print(f"  [skip] {sym}: no blood_eqtl_top_rsid in shortlist; rerun without --no-ensembl")
                continue
            print(f"  [TSS] {sym}: no catalogue sentinel; GRCh37 TSS window for Bakker anchor...")
            tss = fetch_tss_grch37(ens)
            time.sleep(args.sleep_ensembl)
            if tss is None:
                n_skip_no_sentinel += 1
                print(f"  [skip] {sym}: Ensembl GRCh37 TSS lookup failed")
                continue
            s_chr, s_bp = tss
            if norm_chr_label(s_chr) == "23":
                s_chr = "X"
            lo0 = max(1, s_bp - w)
            hi0 = s_bp + w
            g = bakker_best_rs_in_region(args.gwas, s_chr, lo0, hi0)
            if g is None:
                n_skip_no_sentinel += 1
                print(f"  [skip] {sym}: no Bakker SNPs near TSS")
                continue
            try:
                s_chr = str(int(float(g["CHR"]))) if pd.notna(g["CHR"]) else norm_chr_label(g["CHR"])
            except (ValueError, TypeError):
                s_chr = norm_chr_label(g["CHR"]) if pd.notna(g["CHR"]) else s_chr
            if s_chr == "23":
                s_chr = "X"
            s_bp = int(g["BP"])
            effective_rs_for_id = str(g["SNP"])
            proxy_note = "bakker_minP_proxy_no_catalogue_sentinel_in_shortlist"

        assert g is not None

        if args.no_ensembl:
            t_chr, t_bp = s_chr, s_bp
            tss_note = "sentinel_centered_no_ensembl"
        else:
            print(f"  Ensembl GRCh37 TSS for {sym} ({ens})...")
            tss = fetch_tss_grch37(ens)
            time.sleep(args.sleep_ensembl)
            if tss is None:
                t_chr, t_bp = s_chr, s_bp
                tss_note = "TSS_fallback_bakker_anchor"
            else:
                t_chr, t_bp = tss
                tss_note = "ensembl_grch37"
        if norm_chr_label(t_chr) != norm_chr_label(s_chr):
            print(
                f"  [warn] {sym}: Ensembl GRCh37 TSS chr {t_chr} != Bakker anchor chr {s_chr}; "
                "using sentinel-centered cis window only (coloc/eQTL still use Ensembl gene id)."
            )
            t_chr, t_bp = s_chr, s_bp
            tss_note = f"{tss_note}|ensembl_chr_mismatch_sentinel_window"

        lo = min(t_bp - w, s_bp - w)
        hi = max(t_bp + w, s_bp + w)
        lo = max(1, lo)
        lid = _safe_name(f"{sym}_{effective_rs_for_id}")
        rec: dict = {
            "locus_id": lid,
            "symbol": sym,
            "ensembl_gene": ens,
            "sentinel_rsid": rs_eqtl if rs_eqtl.startswith("rs") else effective_rs_for_id,
            "sentinel_chr": s_chr,
            "sentinel_bp": s_bp,
            "BETA_IA": float(g["BETA"]),
            "SE_IA": float(g["SE"]),
            "P_IA": float(g["P"]),
            "A_EFF": str(g.get("A_EFF", "")),
            "A_NONEFF": str(g.get("A_NONEFF", "")),
            "gene_tss_chr": t_chr,
            "gene_tss_bp": t_bp,
            "tss_source": tss_note,
            "coloc_region_GRCh37": f"{t_chr}:{lo}-{hi}",
            "region_start": lo,
            "region_end": hi,
            "blood_eqtl_top_pvalue": row.get("blood_eqtl_top_pvalue"),
            "blood_eqtl_dataset_id": row.get("blood_eqtl_dataset_id"),
            "IA_GWAS_P_at_eqtl_snp": row.get("IA_GWAS_P_at_eqtl_snp"),
        }
        if proxy_note:
            rec["bakker_anchor_note"] = proxy_note
            rec["ia_gwas_anchor_snp"] = str(g["SNP"])
        loci.append(rec)

    if not loci:
        print("No loci built after Ensembl / GWAS merge.", file=sys.stderr)
        sys.exit(1)

    manifest_path = OUT_DIR / "loci_manifest.json"
    if args.merge_existing_manifest and manifest_path.is_file():
        try:
            prev = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            prev = []
        by_id = {x["locus_id"]: x for x in prev if isinstance(x, dict) and "locus_id" in x}
        for L in loci:
            by_id[L["locus_id"]] = L
        loci = sorted(by_id.values(), key=lambda z: str(z.get("locus_id", "")))
        print(f"[merge] manifest now has {len(loci)} loci (merged with prior file).")

    # One pass regional extraction
    intervals_by_id = {L["locus_id"]: L for L in loci if "region_start" in L}
    hits: dict[str, list] = defaultdict(list)

    print("Streaming Bakker GWAS for regional subsets (one pass)...")
    for chunk in pd.read_csv(
        args.gwas, sep="\t", compression="gzip", chunksize=400_000, low_memory=False
    ):
        chunk["_chrn"] = chunk["CHR"].map(norm_chr_label)
        for lid, L in intervals_by_id.items():
            cwant = norm_chr_label(L["sentinel_chr"])
            lo, hi = int(L["region_start"]), int(L["region_end"])
            sub = chunk[chunk["_chrn"] == cwant]
            if sub.empty:
                continue
            bp = sub["BP"].astype(int)
            m = (bp >= lo) & (bp <= hi)
            if m.any():
                hits[lid].append(sub.loc[m].drop(columns=["_chrn"], errors="ignore"))

    for lid, parts in hits.items():
        if not parts:
            continue
        reg = pd.concat(parts, ignore_index=True).drop_duplicates(subset=["SNP"])
        if len(reg) > args.max_snps_per_locus:
            reg = reg.nsmallest(args.max_snps_per_locus, "P")
            note = "truncated_to_max_snps_per_locus"
        else:
            note = "full_region"
        fn = OUT_DIR / f"gwas_region_{lid}.tsv"
        reg.to_csv(fn, sep="\t", index=False)
        L = intervals_by_id[lid]
        L["gwas_region_file"] = str(fn.resolve())
        L["gwas_region_n_snps"] = int(len(reg))
        L["gwas_region_note"] = note

    manifest_path.write_text(json.dumps(loci, indent=2), encoding="utf-8")
    pd.DataFrame(loci).to_csv(OUT_DIR / "loci_manifest.csv", index=False)

    print(f"Wrote {manifest_path}")
    print(f"Wrote regional GWAS TSVs under {OUT_DIR}")
    print(f"Wrote {OUT_DIR / 'RUNBOOK_stronger_genetics.txt'}")
    if n_skip_no_sentinel:
        print(
            f"Note: {n_skip_no_sentinel} tier-1/shortlist row(s) skipped because the eQTL sentinel "
            "rsID is absent from Bakker IA sumstats (use LD proxy or a GWAS that includes the variant)."
        )


if __name__ == "__main__":
    main()
