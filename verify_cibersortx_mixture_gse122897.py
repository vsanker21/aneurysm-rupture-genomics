"""
Quick check: GSE122897_CIBERSORTx_mixture_TPM.txt uses GEO GSM column headers
(NAME + GSM3487...) so CIBERSORTx output can be parsed without --sample-map.

Run: python verify_cibersortx_mixture_gse122897.py
"""
import re
from pathlib import Path

BASE = Path(r"e:\Aneurysms Genomics Project")
MIX = BASE / "results/cibersortx/GSE122897_CIBERSORTx_mixture_TPM.txt"
META = BASE / "results/gse122897_sample_metadata_supplementary.tsv"


def main():
    if not MIX.is_file():
        raise SystemExit(f"Missing {MIX}")
    line = MIX.open(encoding="utf-8", errors="replace").readline().rstrip("\n\r")
    cols = line.split("\t")
    assert cols[0] == "NAME", f"First column must be NAME, got {cols[0]!r}"
    bad = [c for c in cols[1:] if not re.fullmatch(r"GSM\d+", c)]
    if bad:
        raise SystemExit(f"Non-GSM sample columns (first few): {bad[:8]}")
    import pandas as pd

    meta = pd.read_csv(META, sep="\t")
    expected = set(meta["Sample"].astype(str))
    got = set(cols[1:])
    if expected != got:
        raise SystemExit(
            f"Column mismatch vs supplementary metadata.\n"
            f"  only in mixture: {sorted(got - expected)[:10]}\n"
            f"  only in meta: {sorted(expected - got)[:10]}"
        )
    print(f"OK: {MIX.name}")
    print(f"  {len(cols) - 1} sample columns, all GSM*; matches supplementary TSV.")
    print("  Upload this file to CIBERSORTx, then:")
    print('  python parse_cibersortx_output_gse122897.py --input "<Results.txt>"')


if __name__ == "__main__":
    main()
