"""
DESeq2 analysis for GSE122897 tissue RNA-seq.

Primary contrast is 21 ruptured vs 21 unruptured using explicit GEO
source_name labels. GSM3487892 (Unknown rupture status) is excluded.
This wrapper runs run_deseq2_gse122897.R and does not rewrite it.
"""
from __future__ import annotations

import glob
import os
import shutil
import subprocess


def find_rscript():
    rscript = shutil.which("Rscript")
    if rscript:
        return rscript
    for base in [r"C:\Program Files\R", os.path.expandvars(r"%LOCALAPPDATA%\Programs\R")]:
        if not base or not os.path.exists(base):
            continue
        for p in glob.glob(os.path.join(base, "R-*", "bin", "*", "Rscript.exe")):
            if os.path.isfile(p):
                return p
        for p in glob.glob(os.path.join(base, "R-*", "bin", "Rscript.exe")):
            if os.path.isfile(p):
                return p
    return None


def check_raw_counts_and_metadata():
    data_dir = "data/sc_rnaseq/GSE122897"
    raw_file = os.path.join(data_dir, "GSE122897_raw_counts_GRCh38.p13_NCBI.tsv.gz")
    meta_files = [
        "metadata/gse122897_sample_metadata_explicit_21v21.tsv",
        os.path.join(data_dir, "sample_metadata.tsv"),
        os.path.join(data_dir, "sample_metadata_with_clinical.tsv"),
    ]
    meta_file = next((m for m in meta_files if os.path.exists(m)), None)
    if not os.path.exists(raw_file):
        print(f"[ERROR] Raw counts not found: {raw_file}")
        return None, None
    if meta_file is None:
        print("[ERROR] No sample metadata found")
        return None, None
    return raw_file, meta_file


def main():
    raw_file, meta_file = check_raw_counts_and_metadata()
    if raw_file is None:
        return

    script_path = "run_deseq2_gse122897.R"
    if not os.path.isfile(script_path):
        print(f"[ERROR] Missing {script_path}")
        return

    print("GSE122897 DESeq2 Analysis (primary 21 vs 21, explicit GEO labels)")
    print("=" * 50)
    print(f"Raw counts: {raw_file}")
    print(f"Metadata: {meta_file}")
    print(f"R script: {script_path}")

    rscript = find_rscript()
    if not rscript:
        print("\n[INFO] Rscript not found. Run: Rscript run_deseq2_gse122897.R")
        return
    result = subprocess.run(
        [rscript, script_path],
        capture_output=True,
        text=True,
        timeout=3600,
        cwd=os.getcwd(),
    )
    if result.returncode == 0:
        print(result.stdout)
        print("\n[OK] DESeq2 analysis completed successfully")
        print("  Results: results/gse122897_deseq2_results.tsv")
    else:
        print("[ERROR] R script failed")
        print(result.stderr or result.stdout)


if __name__ == "__main__":
    main()
