"""
Build a CIBERSORTx-compatible mixture matrix from GSE122897 TPM (ruptured + unruptured only).

CIBERSORTx (https://cibersortx.stanford.edu/) expects:
  - Rows = genes (HGNC symbols), columns = samples
  - RNA-seq: TPM (or FPKM); do NOT log-transform before upload
  - Tab-separated .txt

After running this script, upload the output file in CIBERSORTx "Mixture" mode
(use the appropriate single-cell or bulk signature; enable "RNA-seq" profile if offered).

Outputs:
  results/cibersortx/GSE122897_CIBERSORTx_mixture_TPM.txt
  results/cibersortx/gse122897_mixture_sample_manifest.tsv
"""
from pathlib import Path

import pandas as pd
import requests

BASE = Path(r"e:\Aneurysms Genomics Project")
TPM_PATH = BASE / "data/sc_rnaseq/GSE122897/GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
META_PATH = BASE / "results/gse122897_sample_metadata_supplementary.tsv"
OUT_DIR = BASE / "results/cibersortx"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def entrez_to_symbol_batch(ids, batch_size=400):
    sym = {}
    ids = [str(int(float(x))) for x in ids if str(x).strip() and str(x) != "nan"]
    for i in range(0, len(ids), batch_size):
        batch = ids[i : i + batch_size]
        try:
            r = requests.post(
                "https://mygene.info/v3/gene",
                data={"ids": ",".join(batch), "fields": "symbol,entrezgene"},
                timeout=120,
            )
            if r.status_code != 200:
                continue
            data = r.json()
            if isinstance(data, dict):
                data = [data]
            for row in data:
                q = str(row.get("query", "")).strip()
                s = row.get("symbol")
                if q and s:
                    sym[q] = s
        except Exception as e:
            print(f"[warn] mygene batch {i}: {e}")
    return sym


def main():
    if not TPM_PATH.exists():
        raise FileNotFoundError(TPM_PATH)
    meta = pd.read_csv(META_PATH, sep="\t")
    samples = meta["Sample"].astype(str).tolist()

    print("Loading TPM matrix...")
    tpm = pd.read_csv(TPM_PATH, sep="\t", compression="gzip")
    tpm = tpm.rename(columns={tpm.columns[0]: "GeneID"})
    tpm["GeneID"] = tpm["GeneID"].astype(str).str.replace(".0", "", regex=False)

    missing = [s for s in samples if s not in tpm.columns]
    if missing:
        raise ValueError(f"Samples not in TPM matrix: {missing[:5]} ...")

    sub = tpm[["GeneID"] + samples].copy()
    print(f"Mapping {len(sub)} Entrez IDs to symbols (mygene)...")
    mapping = entrez_to_symbol_batch(sub["GeneID"].unique().tolist())
    sub["GeneSymbol"] = sub["GeneID"].map(mapping)
    before = len(sub)
    sub = sub.dropna(subset=["GeneSymbol"])
    print(f"Dropped {before - len(sub)} rows without symbol mapping")

    # Collapse duplicate symbols (mean TPM across Entrez IDs)
    expr = sub.drop(columns=["GeneID"]).groupby("GeneSymbol", as_index=False).mean()
    # CIBERSORTx: first column gene symbol; header row OK
    out_cols = ["NAME"] + samples
    expr = expr.rename(columns={"GeneSymbol": "NAME"})
    out_path = OUT_DIR / "GSE122897_CIBERSORTx_mixture_TPM.txt"
    expr.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote {out_path} ({len(expr)} genes x {len(samples)} samples)")
    print(
        "  Column headers: NAME + GEO GSM IDs (e.g. GSM3487853). "
        "After CIBERSORTx, parse without --sample-map."
    )

    manifest = meta.copy()
    manifest["in_mixture_file"] = True
    manifest.to_csv(OUT_DIR / "gse122897_mixture_sample_manifest.tsv", sep="\t", index=False)

    print("\n--- Next steps (CIBERSORTx) ---")
    print("See results/cibersortx/CIBERSORTX_UPLOAD_STEPS.txt for sign-in, LM22, and RNA-seq options.")
    print("1. Register at https://cibersortx.stanford.edu/ (academic use).")
    print("2. Upload the mixture file above; choose LM22 (or updated signature) + RNA-seq settings.")
    print("3. Save the job output .txt/.csv and run:")
    print("   python parse_cibersortx_output_gse122897.py --input <CIBERSORTx_output.txt>")
    print("   (optional: --out-dir results\\cibersortx\\<run_name> to keep multiple downloads)")


if __name__ == "__main__":
    main()
