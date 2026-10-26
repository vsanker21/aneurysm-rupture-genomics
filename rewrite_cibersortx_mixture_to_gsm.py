"""
Replace Pat## labels in CIBERSORTx results with GEO GSM IDs using a two-column map.

Use when the portal still returns Pat in the Mixture column even though you intended
to use a GSM-headed mixture (wrong file selected in dropdown, cached upload, etc.).

Map file: results/cibersortx/cibersortx_pat_to_gsm_map.tsv (columns Mixture, Sample)
Built by: python build_cibersortx_pat_gsm_map_gse122897.py --results <your_results.txt>

Example:
  python rewrite_cibersortx_mixture_to_gsm.py \\
    --input results/cibersortx/CIBERSORTx_Job7_Results.txt \\
    --output results/cibersortx/CIBERSORTx_Job7_Results_GSM_mixture.txt
"""
import argparse
from pathlib import Path

import pandas as pd

BASE = Path(r"e:\Aneurysms Genomics Project")
DEFAULT_MAP = BASE / "results/cibersortx/cibersortx_pat_to_gsm_map.tsv"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True, help="CIBERSORTx results .txt/.tsv")
    ap.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output path (default: <input_stem>_GSM_mixture.txt)",
    )
    ap.add_argument("--map", type=Path, default=DEFAULT_MAP, help="Pat -> GSM TSV")
    args = ap.parse_args()
    if not args.input.is_file():
        raise FileNotFoundError(args.input)
    if not args.map.is_file():
        raise FileNotFoundError(args.map)

    out = args.output
    if out is None:
        out = args.input.with_name(args.input.stem + "_GSM_mixture.txt")

    sep = "\t" if args.input.suffix.lower() in (".txt", ".tsv") else ","
    df = pd.read_csv(args.input, sep=sep)
    if df.columns[0] != "Mixture":
        df = df.rename(columns={df.columns[0]: "Mixture"})

    mp = pd.read_csv(args.map, sep="\t")
    c0, c1 = mp.columns[0], mp.columns[1]
    mmap = dict(zip(mp[c0].astype(str).str.strip(), mp[c1].astype(str).str.strip()))

    before = df["Mixture"].astype(str).str.strip()
    df["Mixture"] = before.map(mmap)
    n_bad = int(df["Mixture"].isna().sum())
    if n_bad:
        missing = before[df["Mixture"].isna()].unique().tolist()
        raise SystemExit(
            f"{n_bad} rows could not be mapped (e.g. {missing[:5]}). "
            f"Rebuild --map from this results file: "
            f"python build_cibersortx_pat_gsm_map_gse122897.py --results {args.input}"
        )

    df.to_csv(out, sep="\t", index=False)
    print(f"Wrote {out.resolve()} ({len(df)} rows, Mixture = GSM)")


if __name__ == "__main__":
    main()
