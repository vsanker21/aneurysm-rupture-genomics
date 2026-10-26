"""
Bridge-cohort analysis for GSE106520:
neutrophil RNA-seq in unruptured intracranial aneurysm cases versus controls.

This is not a rupture cohort. We use it as a blood-accessible, cell-type-specific
bridge dataset to ask whether the tissue-derived shortlist contains detectable
signal in circulating innate immune cells.
"""

from __future__ import annotations

import gzip
import io
import json
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd

from public_cohort_context_utils import add_gse122897_reference, compute_gene_level_welch


try:
    from project_paths import get_project_root
    BASE = get_project_root()
except ImportError:
    BASE = Path(__file__).resolve().parent
DATA_DIR = BASE / "data" / "blood_rnaseq" / "GSE106520"
SERIES_MATRIX = DATA_DIR / "GSE106520_series_matrix.txt.gz"
OUT_DIR = BASE / "results" / "external_validation_gse106520"
OUT_DIR.mkdir(parents=True, exist_ok=True)

REQUEST_HEADERS = {"User-Agent": "AneurysmBloodBridge/1.0 (GSE106520 harmonization)"}


def _fetch_text(url: str) -> str:
    request = Request(url, headers=REQUEST_HEADERS)
    with urlopen(request, timeout=120) as response:
        return response.read().decode("utf-8", errors="ignore")


def _download_if_missing(url: str, destination: Path) -> None:
    if destination.is_file() and destination.stat().st_size > 0:
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = Request(url, headers=REQUEST_HEADERS)
    with urlopen(request, timeout=300) as response, destination.open("wb") as handle:
        handle.write(response.read())


def _parse_series_matrix(path: Path) -> pd.DataFrame:
    lines = gzip.open(path, "rt", encoding="utf-8", errors="ignore").read().splitlines()
    begin = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_begin"))
    sample_lines = [line for line in lines[:begin] if line.startswith("!Sample_")]

    sample_ids: list[str] = []
    titles: list[str] = []
    sources: list[str] = []
    supp_files: list[str] = []
    characteristics_rows: list[list[str]] = []

    for line in sample_lines:
        parts = [x.strip().strip('"') for x in line.split("\t")[1:]]
        if line.startswith("!Sample_geo_accession"):
            sample_ids = parts
        elif line.startswith("!Sample_title"):
            titles = parts
        elif line.startswith("!Sample_source_name_ch1"):
            sources = parts
        elif line.startswith("!Sample_supplementary_file_1"):
            supp_files = parts
        elif line.startswith("!Sample_characteristics_ch1"):
            characteristics_rows.append(parts)

    meta = pd.DataFrame(
        {
            "Sample": sample_ids,
            "Title": titles,
            "Source": sources,
            "Supplementary_file": supp_files,
        }
    )

    for row in characteristics_rows:
        if not row:
            continue
        key = row[0].split(":", 1)[0].strip().lower().replace(" ", "_").replace("(", "").replace(")", "")
        values = []
        for value in row:
            if ":" in value:
                values.append(value.split(":", 1)[1].strip())
            else:
                values.append(value.strip())
        meta[key] = values

    meta["Group"] = np.where(
        meta["Title"].str.lower().str.startswith("aneurysm_"),
        "Aneurysm",
        "Control",
    )
    meta["Sample_label"] = meta["Title"].str.replace("_", "-", regex=False)
    meta["Cohort"] = meta["cohort"].fillna("")
    meta["Age"] = pd.to_numeric(meta["age"], errors="coerce")
    meta["Sex"] = meta["sex"].astype(str).str.strip()
    meta["Smoker"] = meta["smoker_y/n"].astype(str).str.strip()
    meta["Hypertension"] = pd.to_numeric(meta["hypertension_presence;_0=no,1=yes"], errors="coerce")
    meta["Hyperlipidemia"] = pd.to_numeric(meta["hyperlipidemia_presence;_0=no,1=yes"], errors="coerce")
    meta["Coronary_artery_disease"] = pd.to_numeric(
        meta["coronary_artery_disease_presence;_0=no,1=yes"], errors="coerce"
    )
    meta["Stroke_history"] = pd.to_numeric(meta["stroke_history_presence;_0=no,1=yes"], errors="coerce")
    meta["Diabetes_mellitus"] = pd.to_numeric(meta["diabetes_mellitus_presence;_0=no,1=yes"], errors="coerce")
    meta["Osteoarthritis"] = pd.to_numeric(meta["osteoarthritis_presence;_0=no,1=yes"], errors="coerce")
    return meta


def _ensure_sample_files(meta: pd.DataFrame) -> None:
    for _, row in meta.iterrows():
        url = str(row["Supplementary_file"]).replace("ftp://", "https://")
        filename = url.rsplit("/", 1)[-1]
        _download_if_missing(url, DATA_DIR / filename)


def _load_sample_expression(meta: pd.DataFrame) -> pd.DataFrame:
    matrices: list[pd.DataFrame] = []
    sample_paths: list[str] = []
    for _, row in meta.iterrows():
        url = str(row["Supplementary_file"]).replace("ftp://", "https://")
        filename = url.rsplit("/", 1)[-1]
        sample_path = DATA_DIR / filename
        if not sample_path.is_file():
            raise FileNotFoundError(sample_path)
        sample_paths.append(sample_path.name)
        with gzip.open(sample_path, "rt", encoding="utf-8", errors="ignore") as handle:
            df = pd.read_csv(handle, sep="\t")
        value_col = df.columns[-1]
        df = df[["gene_short_name", value_col]].rename(columns={"gene_short_name": "Gene", value_col: str(row["Sample"])})
        df["Gene"] = df["Gene"].astype(str).str.strip()
        df = df[(df["Gene"] != "") & (df["Gene"] != "-")].copy()
        df[str(row["Sample"])] = pd.to_numeric(df[str(row["Sample"])], errors="coerce")
        df = df.groupby("Gene", as_index=False)[str(row["Sample"])].mean()
        matrices.append(df)

    merged = matrices[0]
    for df in matrices[1:]:
        merged = merged.merge(df, on="Gene", how="outer")
    gene_expr = merged.set_index("Gene").sort_index()
    gene_expr = np.log2(gene_expr + 1.0)
    meta["Local_sample_file"] = sample_paths
    return gene_expr


def main() -> None:
    if not SERIES_MATRIX.is_file():
        raise FileNotFoundError(SERIES_MATRIX)

    meta = _parse_series_matrix(SERIES_MATRIX)
    _ensure_sample_files(meta)
    gene_expr = _load_sample_expression(meta)

    meta.to_csv(OUT_DIR / "gse106520_sample_metadata.tsv", sep="\t", index=False)
    gene_expr.to_csv(OUT_DIR / "gse106520_gene_level_log2fpkm.tsv", sep="\t", index=True)

    aneurysm_cols = meta.loc[meta["Group"] == "Aneurysm", "Sample"].tolist()
    control_cols = meta.loc[meta["Group"] == "Control", "Sample"].tolist()
    res = compute_gene_level_welch(
        gene_expr,
        group_a_cols=aneurysm_cols,
        group_b_cols=control_cols,
        group_a_label="Aneurysm",
        group_b_label="Control",
        effect_name="Log2FC_AvC",
    )
    merged = add_gse122897_reference(res, base=BASE, effect_col="Log2FC_AvC")

    res.to_csv(OUT_DIR / "gse106520_gene_level_context.tsv", sep="\t", index=False)
    merged.to_csv(OUT_DIR / "gse106520_gene_level_context_with_gse122897.tsv", sep="\t", index=False)
    merged.sort_values("P_value", ascending=True).head(30).to_csv(
        OUT_DIR / "gse106520_top30_genes.tsv", sep="\t", index=False
    )

    shortlist_hits = merged[merged["In_PHAGESG_shortlist"] == True].copy()  # noqa: E712
    if not shortlist_hits.empty:
        shortlist_hits.sort_values("P_value", ascending=True).to_csv(
            OUT_DIR / "gse106520_shortlist_overlap.tsv", sep="\t", index=False
        )

    summary = {
        "dataset": "GSE106520",
        "comparison": "neutrophil RNA-seq, unruptured aneurysm vs control",
        "platform": "HiSeq2500 per-sample FPKM tables merged at gene-symbol level",
        "n_total_samples_loaded": int(len(meta)),
        "n_aneurysm": int(len(aneurysm_cols)),
        "n_control": int(len(control_cols)),
        "n_original_cohort": int((meta["Cohort"] == "Original").sum()),
        "n_replication_cohort": int((meta["Cohort"] == "Replication").sum()),
        "n_genes_tested": int(len(res)),
        "n_genes_p_lt_0_05": int((res["P_value"] < 0.05).sum()),
        "n_genes_fdr_lt_0_05": int((res["FDR_BH"] < 0.05).sum()),
        "n_overlap_with_gse122897_symbols": int(merged["GSE122897_log2FC_compare"].notna().sum()),
        "n_direction_match_with_gse122897": int(
            merged.loc[merged["GSE122897_log2FC_compare"].notna(), "DirectionMatch_GSE122897"].fillna(False).sum()
        ),
        "n_shortlist_genes_present": int(merged["In_PHAGESG_shortlist"].sum()),
        "n_shortlist_direction_match": int(
            merged.loc[merged["In_PHAGESG_shortlist"] == True, "DirectionMatch_GSE122897"].fillna(False).sum()  # noqa: E712
        ),
        "n_shortlist_p_lt_0_05": int(
            ((merged["In_PHAGESG_shortlist"] == True) & (merged["P_value"] < 0.05)).sum()  # noqa: E712
        ),
        "n_shortlist_fdr_lt_0_05": int(
            ((merged["In_PHAGESG_shortlist"] == True) & (merged["FDR_BH"] < 0.05)).sum()  # noqa: E712
        ),
        "bridge_context_note": (
            "This cohort profiles circulating neutrophils from unruptured aneurysm patients and is used as a blood-accessible bridge "
            "rather than as direct rupture validation."
        ),
    }
    (OUT_DIR / "gse106520_context_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
