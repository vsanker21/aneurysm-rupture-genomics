"""
Download genotype dosages harmonized to prs_rsid_weights_from_panel.tsv (risk allele),
then run phases_g_prs_from_genotypes.py to produce prs_sample_scores.tsv.

REFERENCE COHORT (not PHASES-labelled): 1000 Genomes Phase 3 via Ensembl REST
individual genotypes, filtered to superpopulation EUR by default using the IGSR
integrated_call_samples panel. Some rsIDs lack Phase-3 genotype rows in Ensembl;
those are QC-reported and omitted from dosage output (fewer SNPs contribute to PRS).

REQUIRES: Network access, stdlib + pandas + numpy (same as phases_g_prs_from_genotypes).

Example:
  python prs_download_phase3_reference_and_score.py --eur-only
  python prs_download_phase3_reference_and_score.py --super-pop AFR \\
      --clinical data/genotypes/example_clinical_placeholder.tsv

Real PHASES-aligned DNA must come from your biobank; this script satisfies
'automated downloads + harmonized doses + runnable PRS pipeline' on public data only.

Outputs use column sample_id. For this reference cohort the values are 1000G IDs
(HG….). Replace them with identifiers aligned to your PHASES bookkeeping when using
your own dosage export—or pass --rename-map (legacy_id → sample_id) before PRS.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import time
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent
WEIGHTS_DEFAULT = BASE / "results/phases_g/prs_rsid_weights_from_panel.tsv"
PANEL_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    "integrated_call_samples_v3.20130502.ALL.panel"
)
ENSEMBL_VAR_BATCH = (
    "https://rest.ensembl.org/variation/human?"
    "content-type=application/json;genotypes=1"
)
DOSAGE_LONG_DEFAULT = BASE / "data/genotypes/phases_panel_1kgp3_phase3_reference_long.tsv"
PRS_SCRIPT = BASE / "phases_g_prs_from_genotypes.py"


def download_text(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"Accept": "*/*"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        dest.write_bytes(resp.read())


def apply_sample_id_rename_map(long_df: pd.DataFrame, map_tsv: Path) -> tuple[pd.DataFrame, dict]:
    """
    Rewrite sample_id column from a two-column TSV (legacy_id → sample_id).
    Unmatched ids are left unchanged.
    """
    mp = pd.read_csv(map_tsv, sep="\t", dtype=str)
    if mp.shape[1] < 2:
        raise ValueError("--rename-map TSV requires at least two columns (legacy → sample_id).")
    legacy_col, new_col = mp.columns[0], mp.columns[1]
    d = dict(
        zip(
            mp[legacy_col].astype(str).str.strip(),
            mp[new_col].astype(str).str.strip(),
        )
    )
    before = long_df["sample_id"].astype(str).str.strip()
    after = before.replace(d)
    out = long_df.copy()
    out["sample_id"] = after.values
    n_changed = int((before != after).sum())
    return out, {"n_rows_renamed": n_changed, "n_keys_in_map": len(d), "rename_map_path": str(map_tsv.resolve())}


def load_phase3_sample_superpop(panel_tsv: Path) -> pd.DataFrame:
    df = pd.read_csv(panel_tsv, sep="\t")
    df.columns = [c.strip().lower() for c in df.columns]
    return df[["sample", "super_pop"]].drop_duplicates()


def _ensembl_batch(rsids: list[str], retries: int = 5) -> dict:
    body = json.dumps({"ids": rsids}).encode("utf-8")
    req = urllib.request.Request(
        ENSEMBL_VAR_BATCH,
        data=body,
        method="POST",
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": "prs_download_phase3_reference/1.0 (phases-g-pipeline)",
        },
    )
    wait = 3.0
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=300) as resp:
                return json.loads(resp.read().decode())
        except urllib.error.HTTPError as e:
            if e.code == 429 and attempt < retries - 1:
                time.sleep(wait)
                wait *= 1.8
                continue
            raise
    raise RuntimeError("unreachable")


def parse_phased_gt(gt: str) -> list[str]:
    gt = (gt or "").strip().upper().replace(",", "/")
    if not gt or "." in gt.lower():
        return []
    sep = "|" if "|" in gt else ("/" if "/" in gt else None)
    if sep is None:
        return []
    return [a.strip() for a in gt.split(sep) if a.strip()]


def harmonized_risk_dosage(
    genotype: str,
    risk_allele: str,
    assessed_allele: str,
    other_allele: str,
) -> tuple[float | None, str | None]:
    """
    Count copies of risk_allele in phased/unphased diploid genotype.
    Genotype nucleotides must match the eQTL **assessed/other** pair (not {risk, other}:
    risk can equal one of those two and both may match the same letter only if the row is bad).
    """
    risk = str(risk_allele).strip().upper()
    ass = str(assessed_allele).strip().upper() if assessed_allele is not None else ""
    oth = str(other_allele).strip().upper() if other_allele is not None else ""

    alleles = parse_phased_gt(genotype)
    if len(alleles) != 2:
        return None, "malformed_or_missing_gt"
    if any(len(a) != 1 for a in alleles):
        return None, "not_biallelic_snp"
    for a in alleles:
        if a not in "ACGT":
            return None, "non_acgt"

    if ass and oth:
        allowed = {ass, oth}
        if len(allowed) < 2:
            return None, "degenerate_assessed_other_in_weights"
        stray = set(alleles) - allowed
        if stray:
            return None, "allele_not_matching_eqltl_assessed_other"

    d = sum(1 for a in alleles if a == risk)
    return float(d), None


def fetch_phase3_gt_long(
    rsids: list[str],
    weights: pd.DataFrame,
    keep_samples: set[str],
    *,
    batch_size: int = 25,
    sleep_s: float = 0.8,
    cache_json: Path | None = None,
) -> tuple[pd.DataFrame, dict]:
    """
    Returns long table sample_id (IGSR codes like HG00096 for this datasource), rsid, dosage harmonized to risk_allele + QC stats.
    """
    risk_map = dict(zip(weights["rsid"].astype(str), weights["risk_allele"].astype(str)))
    other_map = dict(zip(weights["rsid"].astype(str), weights["other_allele"].astype(str)))
    assessed_map = dict(
        zip(weights["rsid"].astype(str), weights["assessed_allele"].astype(str))
    )

    qc: dict = {
        "rsid_phase3_records": defaultdict(int),
        "rsid_harmonize_fail_examples": defaultdict(list),
    }

    batches: dict[str, dict] = {}
    unique_rs = list(dict.fromkeys(rsids))

    for i in range(0, len(unique_rs), batch_size):
        chunk = unique_rs[i : i + batch_size]
        print(
            f"[ensembl] variation POST {len(chunk)} ids ({i + 1}-{i + len(chunk)} / {len(unique_rs)})",
            flush=True,
        )
        batches[str(i)] = _ensembl_batch(chunk)
        time.sleep(sleep_s)

    if cache_json is not None:
        cache_json.parent.mkdir(parents=True, exist_ok=True)
        cache_json.write_text(json.dumps(batches), encoding="utf-8")

    rows: list[dict] = []
    for blob in batches.values():
        if not isinstance(blob, dict):
            continue
        for rsid_raw, vd in blob.items():
            rsid = str(vd.get("name") or rsid_raw).strip()
            risk = risk_map.get(rsid)
            if risk is None:
                continue
            gtl = vd.get("genotypes") or []
            n1 = 0
            for row in gtl:
                samp = row.get("sample") or ""
                if not samp.startswith("1000GENOMES:phase_3:"):
                    continue
                sid = samp.split(":")[-1]
                if keep_samples is not None and sid not in keep_samples:
                    continue
                other = str(other_map.get(rsid, ""))
                ass = str(assessed_map.get(rsid, ""))
                dosage, err = harmonized_risk_dosage(
                    row.get("genotype", ""), risk, ass, other
                )
                if dosage is None and err:
                    eg = qc["rsid_harmonize_fail_examples"][rsid]
                    if len(eg) < 3:
                        eg.append(
                            {
                                "sample": sid,
                                "gt": row.get("genotype"),
                                "risk": risk,
                                "weights_assessed": ass,
                                "weights_other": other,
                                "reason": err,
                            }
                        )
                    continue
                rows.append({"sample_id": sid, "rsid": rsid, "dosage": dosage})
                n1 += 1
            qc["rsid_phase3_records"][rsid] = max(qc["rsid_phase3_records"][rsid], n1)

    qc_summary = {
        "n_rsids_requested": len(unique_rs),
        "n_rsids_with_any_phase3_gt": sum(
            1 for rs in unique_rs if qc["rsid_phase3_records"].get(rs, 0) > 0
        ),
        "n_rsids_zero_phase3_in_ensembl": sum(
            1 for rs in unique_rs if qc["rsid_phase3_records"].get(rs, 0) == 0
        ),
    }
    qc_out = {
        "rsid_phase3_records": {k: int(v) for k, v in qc["rsid_phase3_records"].items()},
        "rsid_harmonize_fail_examples": dict(qc["rsid_harmonize_fail_examples"]),
        "summary": qc_summary,
    }
    return pd.DataFrame(rows), qc_out


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--weights", type=Path, default=WEIGHTS_DEFAULT)
    ap.add_argument(
        "--super-pop",
        default="EUR",
        help="Integrated call super_pop column value (EUR, AFR, EAS, AMR, SAS, ALL). ALL keeps all Phase3.",
    )
    ap.add_argument(
        "--eur-only",
        action="store_true",
        help="Shortcut for --super-pop EUR.",
    )
    ap.add_argument(
        "--dosage-long-out",
        type=Path,
        default=DOSAGE_LONG_DEFAULT,
        help="Harmonized long dosage TSV.",
    )
    ap.add_argument(
        "--qc-json",
        type=Path,
        default=BASE / "results/phases_g/prs_reference_phase3_download_qc.json",
    )
    ap.add_argument(
        "--clinical",
        type=Path,
        default=None,
        help="Forwarded to phases_g_prs_from_genotypes.py.",
    )
    ap.add_argument(
        "--prs-out",
        type=Path,
        default=BASE / "results/phases_g/prs_sample_scores.tsv",
    )
    ap.add_argument("--skip-prs-step", action="store_true")
    ap.add_argument("--batch-size", type=int, default=25)
    ap.add_argument("--sleep", type=float, default=1.0, help="Delay between POST batches (polite).")
    ap.add_argument(
        "--max-snps",
        type=int,
        default=None,
        help="Only first N rsIDs from weights (debug / smoke test).",
    )
    ap.add_argument(
        "--rename-map",
        type=Path,
        default=None,
        help="Two-column TSV (legacy_id, sample_id) to replace reference HG IDs with PHASES-aligned keys.",
    )
    args = ap.parse_args()
    sp = "EUR" if args.eur_only else args.super_pop.strip().upper()

    if not args.weights.is_file():
        print(f"Missing weights {args.weights}; run: python phases_g_prs_from_genotypes.py --build-weights-only")
        sys.exit(1)

    wdf = pd.read_csv(args.weights, sep="\t")
    wdf["risk_allele"] = wdf["risk_allele"].astype(str).str.strip()
    wdf["assessed_allele"] = wdf["assessed_allele"].fillna("").astype(str).str.strip()
    wdf["other_allele"] = wdf["other_allele"].fillna("").astype(str).str.strip()
    rsids = wdf["rsid"].astype(str).unique().tolist()
    if args.max_snps is not None:
        rsids = rsids[: max(0, args.max_snps)]
        wdf = wdf[wdf["rsid"].isin(rsids)]

    panel_cache = BASE / "data/genotypes/cache/integrated_call_samples_v3.20130502.ALL.panel"
    download_text(PANEL_URL, panel_cache)
    meta = load_phase3_sample_superpop(panel_cache)
    if sp == "ALL":
        keep = set(meta["sample"].astype(str))
    else:
        keep = set(meta.loc[meta["super_pop"].astype(str) == sp, "sample"].astype(str))

    print(f"Cohort filter: super_pop={sp} (n_samples={len(keep)})")

    long_df, qc = fetch_phase3_gt_long(
        rsids,
        wdf,
        keep,
        batch_size=args.batch_size,
        sleep_s=args.sleep,
        cache_json=(
            BASE / "data/genotypes/cache/ensembl_variation_phase3_batches.json"
            if args.max_snps is None
            else None
        ),
    )

    rename_meta: dict = {}
    if args.rename_map is not None:
        if not args.rename_map.is_file():
            print(f"Missing --rename-map {args.rename_map}", file=sys.stderr)
            sys.exit(1)
        long_df, rename_meta = apply_sample_id_rename_map(long_df, args.rename_map)
        print(
            f"Applied --rename-map: {rename_meta['n_rows_renamed']} long-table rows rewritten",
            flush=True,
        )

    args.dosage_long_out.parent.mkdir(parents=True, exist_ok=True)
    long_df.sort_values(["sample_id", "rsid"]).to_csv(
        args.dosage_long_out, sep="\t", index=False
    )
    hf = qc.get("rsid_harmonize_fail_examples") or {}
    qc_serial = {
        **qc["summary"],
        "harmonize_fail_rsids_preview": {
            k: hf[k] for k in list(hf.keys())[:15] if hf.get(k)
        },
        "weights": str(args.weights.resolve()),
        "dosage_long": str(args.dosage_long_out.resolve()),
        "cohort_super_pop_filter": sp,
        **({"sample_id_rename": rename_meta} if rename_meta else {}),
    }
    args.qc_json.parent.mkdir(parents=True, exist_ok=True)
    args.qc_json.write_text(json.dumps(qc_serial, indent=2), encoding="utf-8")

    print(
        f"Wrote {args.dosage_long_out}\nQC {args.qc_json}\n"
        f"phase3-covering rsIDs: {qc_serial['n_rsids_with_any_phase3_gt']} / "
        f"{qc_serial['n_rsids_requested']}",
        flush=True,
    )

    if args.skip_prs_step:
        return

    cmd = [sys.executable, str(PRS_SCRIPT)]
    if args.weights.is_file():
        cmd.extend(["--from-weights-tsv", "--weights", str(args.weights)])
    cmd.extend(
        [
            "--dosage-long",
            str(args.dosage_long_out),
            "--out",
            str(args.prs_out),
        ]
    )
    if args.clinical and args.clinical.is_file():
        cmd.extend(["--clinical", str(args.clinical)])
    print("RUN:", " ".join(cmd), flush=True)
    subprocess.check_call(cmd, cwd=str(BASE))
    print(f"Done PRS -> {args.prs_out}", flush=True)


if __name__ == "__main__":
    main()
