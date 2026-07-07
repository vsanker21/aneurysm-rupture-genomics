"""
Merge raw and refined metadata tables for GSE159610 and ensure labels are present.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd


RAW_METADATA = Path("data/blood_rnaseq/GSE159610/sample_metadata.tsv")
REFINED_METADATA = Path("data/blood_rnaseq/GSE159610/sample_metadata_refined.tsv")
OUTPUT_METADATA = Path("data/blood_rnaseq/GSE159610/sample_metadata_with_labels.tsv")


def _normalize_sample_ids(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize sample identifiers to simplify joins."""
    df = df.copy()
    df["Sample"] = (
        df["Sample"]
        .astype(str)
        .str.strip()
        .str.upper()
    )
    return df


def _extract_group_mapping(raw_df: pd.DataFrame) -> pd.DataFrame:
    """Extract the best available group column from the raw metadata."""
    def _is_informative(series: pd.Series) -> bool:
        cleaned = series.dropna().astype(str).str.strip()
        if cleaned.empty:
            return False
        unique_vals = cleaned.str.lower().unique()
        if len(unique_vals) == 1:
            val = unique_vals[0]
            if val.startswith("public on") or val in {"unknown"}:
                return False
        return True

    group_col: Optional[str] = None
    for col in raw_df.columns:
        if col.strip().lower() == "group":
            if _is_informative(raw_df[col]):
                group_col = col
                break

    if group_col is None:
        for col in raw_df.columns:
            if col.strip().lower() in {"phenotype"}:
                if _is_informative(raw_df[col]):
                    group_col = col
                    break

    if group_col is None:
        group_df = raw_df[["Sample"]].copy()
        group_df["Group"] = pd.NA
        return group_df

    group_df = (
        raw_df[["Sample", group_col]]
        .rename(columns={group_col: "Group"})
        .copy()
    )
    return group_df


def build_metadata_with_labels(
    raw_path: Path = RAW_METADATA,
    refined_path: Path = REFINED_METADATA,
    output_path: Path = OUTPUT_METADATA,
) -> pd.DataFrame:
    """Create metadata file that always includes explicit labels."""
    if not raw_path.exists():
        raise FileNotFoundError(f"Raw metadata file not found: {raw_path}")

    if not refined_path.exists():
        raise FileNotFoundError(f"Refined metadata file not found: {refined_path}")

    raw_df = pd.read_csv(raw_path, sep="\t", dtype=str)
    refined_df = pd.read_csv(refined_path, sep="\t", dtype=str)

    raw_df = _normalize_sample_ids(raw_df)
    refined_df = _normalize_sample_ids(refined_df)

    group_df = _extract_group_mapping(raw_df).drop_duplicates(subset="Sample")

    merged = refined_df.merge(group_df, on="Sample", how="left")
    if "Group" not in merged.columns:
        merged["Group"] = merged["Status"]
    else:
        missing_group_mask = merged["Group"].isna()
        if missing_group_mask.any():
            merged.loc[missing_group_mask, "Group"] = merged.loc[missing_group_mask, "Status"]

    # Keep the key columns up front
    preferred_cols = [col for col in ["Sample", "Status", "Title", "Group"] if col in merged.columns]
    other_cols = [col for col in merged.columns if col not in preferred_cols]
    ordered_df = merged[preferred_cols + other_cols]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    ordered_df.to_csv(output_path, sep="\t", index=False)

    return ordered_df


def main() -> None:
    merged = build_metadata_with_labels()
    print(f"[OK] Saved merged metadata with labels to {OUTPUT_METADATA}")
    print(f"Total samples: {len(merged)} | Columns: {list(merged.columns)}")


if __name__ == "__main__":
    main()

