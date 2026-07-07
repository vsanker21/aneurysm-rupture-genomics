"""
Pre-flight checks for phases_g_prs_from_genotypes.py inputs.

Exits with code 1 on hard failures (missing columns, duplicate long keys).
Prints coverage and clinical overlap summaries (warnings only unless --strict).

Usage:
  python validate_prs_genotype_inputs.py --dosage-long data/genotypes/your_long.tsv \\
      --weights results/phases_g/prs_rsid_weights_from_panel.tsv
  python validate_prs_genotype_inputs.py --dosage-long ... --clinical ... --strict
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent


def _fatal(msg: str) -> None:
    print(f"[FAIL] {msg}", file=sys.stderr)


def _warn(msg: str) -> None:
    print(f"[WARN] {msg}", file=sys.stderr)


def _ok(msg: str) -> None:
    print(f"[OK] {msg}")


def validate_long(
    path: Path,
    weights: pd.DataFrame,
    strict: bool,
) -> tuple[set[str], set[str]]:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    need = {"sample_id", "rsid", "dosage"}
    miss = need - set(df.columns)
    if miss:
        _fatal(f"{path}: missing columns {sorted(miss)}")
        sys.exit(1)
    dup = df.duplicated(subset=["sample_id", "rsid"], keep=False)
    if dup.any():
        n = int(dup.sum())
        _fatal(f"{path}: {n} rows with duplicate (sample_id, rsid) pairs")
        sys.exit(1)
    try:
        dnum = pd.to_numeric(df["dosage"], errors="coerce")
    except Exception as e:
        _fatal(f"{path}: dosage not numeric ({e})")
        sys.exit(1)
    if dnum.isna().any():
        nna = int(dnum.isna().sum())
        _fatal(f"{path}: {nna} non-numeric dosage values")
        sys.exit(1)
    bad = (dnum < 0) | (dnum > 2)
    if bad.any():
        msg = f"{path}: {int(bad.sum())} dosage values outside [0, 2]"
        if strict:
            _fatal(msg)
            sys.exit(1)
        _warn(msg)

    w_rs = set(weights["rsid"].astype(str))
    d_rs = set(df["rsid"].astype(str))
    samples = set(df["sample_id"].astype(str))
    in_both = w_rs & d_rs
    missing_panel = sorted(w_rs - d_rs)[:20]
    _ok(
        f"long: {len(samples)} samples, {len(d_rs)} distinct rsids in dosage; "
        f"{len(in_both)} / {len(w_rs)} weight rsids have at least one dosage row"
    )
    if missing_panel:
        nmiss = len(w_rs - d_rs)
        _warn(f"{nmiss} weight rsids absent from dosage (first 20): {missing_panel}")
        if strict and nmiss > len(w_rs) * 0.5:
            _fatal("--strict: more than half of weight rsids missing from dosage")
            sys.exit(1)
    extra = d_rs - w_rs
    if extra:
        _warn(f"{len(extra)} dosage rsids not in weights file (ignored by PRS merge)")
    return samples, d_rs


def validate_wide(path: Path, weights: pd.DataFrame, strict: bool) -> tuple[set[str], set[str]]:
    df = pd.read_csv(path, sep="\t", nrows=0)
    if "sample_id" not in df.columns:
        _fatal(f"{path}: missing sample_id column")
        sys.exit(1)
    df = pd.read_csv(path, sep="\t", low_memory=False)
    if df["sample_id"].duplicated().any():
        _fatal(f"{path}: duplicate sample_id rows in wide matrix")
        sys.exit(1)
    rs_cols = [c for c in df.columns if c != "sample_id"]
    w_rs = set(weights["rsid"].astype(str))
    d_rs = set(rs_cols)
    samples = set(df["sample_id"].astype(str))
    in_both = w_rs & d_rs
    _ok(
        f"wide: {len(samples)} samples, {len(d_rs)} rsid columns; "
        f"{len(in_both)} / {len(w_rs)} weight rsids present as columns"
    )
    if len(w_rs - d_rs) > 0:
        _warn(f"{len(w_rs - d_rs)} weight rsids missing as columns in wide matrix")
    return samples, d_rs


def validate_clinical(path: Path, dosage_samples: set[str], strict: bool) -> None:
    c = pd.read_csv(path, sep="\t", low_memory=False)
    if "sample_id" not in c.columns:
        _fatal(f"{path}: clinical file must include sample_id")
        sys.exit(1)
    dup = c["sample_id"].duplicated(keep=False)
    if dup.any():
        _fatal(f"{path}: duplicate sample_id rows ({int(dup.sum())} rows)")
        sys.exit(1)
    clin_ids = set(c["sample_id"].astype(str))
    overlap = clin_ids & dosage_samples
    only_clin = clin_ids - dosage_samples
    only_dose = dosage_samples - clin_ids
    _ok(
        f"clinical: {len(clin_ids)} subjects; overlap with dosage samples: "
        f"{len(overlap)} / {len(dosage_samples)}"
    )
    if only_clin:
        _warn(f"{len(only_clin)} clinical-only sample_ids (merge leaves PRS-only columns NaN for dosage rows)")
    if only_dose:
        _warn(f"{len(only_dose)} dosage-only sample_ids (no clinical row; merged columns NaN)")
    if strict and len(overlap) == 0 and len(clin_ids) > 0 and len(dosage_samples) > 0:
        _fatal("--strict: zero sample_id overlap between clinical and dosage")
        sys.exit(1)


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dosage-long", type=Path, default=None)
    ap.add_argument("--dosage-wide", type=Path, default=None)
    ap.add_argument(
        "--weights",
        type=Path,
        default=BASE / "results/phases_g/prs_rsid_weights_from_panel.tsv",
    )
    ap.add_argument("--clinical", type=Path, default=None)
    ap.add_argument(
        "--strict",
        action="store_true",
        help="Fail if dosages outside [0,2] or (with clinical) zero overlap",
    )
    args = ap.parse_args()

    if args.dosage_long is None and args.dosage_wide is None:
        _fatal("Provide --dosage-long or --dosage-wide")
        sys.exit(1)
    if args.dosage_long is not None and args.dosage_wide is not None:
        _fatal("Provide only one of --dosage-long or --dosage-wide")
        sys.exit(1)

    if not args.weights.is_file():
        _fatal(f"Weights not found: {args.weights}")
        sys.exit(1)

    w = pd.read_csv(args.weights, sep="\t")
    if "rsid" not in w.columns or "weight" not in w.columns:
        _fatal("weights TSV must contain rsid and weight")
        sys.exit(1)

    if args.dosage_long is not None:
        if not args.dosage_long.is_file():
            _fatal(f"Missing {args.dosage_long}")
            sys.exit(1)
        samples, _ = validate_long(args.dosage_long, w, args.strict)
    else:
        if not args.dosage_wide.is_file():
            _fatal(f"Missing {args.dosage_wide}")
            sys.exit(1)
        samples, _ = validate_wide(args.dosage_wide, w, args.strict)

    if args.clinical is not None:
        if not args.clinical.is_file():
            _fatal(f"Missing {args.clinical}")
            sys.exit(1)
        validate_clinical(args.clinical, samples, args.strict)

    _ok("validation finished")
    sys.exit(0)


if __name__ == "__main__":
    main()
