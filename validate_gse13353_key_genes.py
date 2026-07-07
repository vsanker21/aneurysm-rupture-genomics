"""
External validation: GSE13353 (ruptured vs unruptured IA walls, GPL570)
against key aneurysm-project genes (SMS, MARCO, TGFBI, COQ2 up; CPXM2, LAMA3, LIMCH1 down in ruptured).
"""
from pathlib import Path
import csv
import gzip

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

BASE = Path(r"e:\Aneurysms Genomics Project")
EXT = BASE / "External Validation"
OUT = BASE / "results" / "external_validation_gse13353"
OUT.mkdir(parents=True, exist_ok=True)

# One representative probe per gene (GPL570.annot.gz; multi-probe genes use primary _at where standard)
GENE_PROBES = {
    "SMS": ["202043_s_at"],
    "MARCO": ["205819_at"],
    "TGFBI": ["201506_at"],
    "COQ2": ["213379_at"],
    "CPXM2": ["226824_at"],
    "LAMA3": ["203726_s_at"],
    "LIMCH1": ["212325_at", "212327_at", "212328_at"],
}

EXPECTED = {
    "SMS": "higher_in_ruptured",
    "MARCO": "higher_in_ruptured",
    "TGFBI": "higher_in_ruptured",
    "COQ2": "higher_in_ruptured",
    "CPXM2": "lower_in_ruptured",
    "LAMA3": "lower_in_ruptured",
    "LIMCH1": "lower_in_ruptured",
}


def parse_series_matrix(path: Path):
    meta = {}
    in_table = False
    header = None
    rows = []
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                in_table = False
                continue
            if in_table:
                if header is None:
                    header = next(csv.reader([line], delimiter="\t", quotechar='"'))
                else:
                    rows.append(next(csv.reader([line], delimiter="\t", quotechar='"')))
                continue
            if line.startswith("!Sample_"):
                parts = next(csv.reader([line], delimiter="\t", quotechar='"'))
                key = parts[0].replace("!Sample_", "")
                meta.setdefault(key, []).append(parts[1:])
    samples = header[1:]
    df = pd.DataFrame(rows, columns=header)
    df = df.rename(columns={df.columns[0]: "ID_REF"})
    for c in samples:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df, samples, meta


def main():
    series_path = EXT / "GSE13353_series_matrix.txt.gz"
    if not series_path.exists():
        raise FileNotFoundError(series_path)

    expr, samples, meta = parse_series_matrix(series_path)
    expr = expr.set_index("ID_REF")

    def col(key, idx=0):
        rows = meta.get(key, [])
        if len(rows) > idx:
            return rows[idx]
        return [""] * len(samples)

    titles = col("title", 0)
    chars = col("characteristics_ch1", 0)
    pheno = []
    for t, ch in zip(titles, chars):
        s = (str(t) + " " + str(ch)).lower()
        # "unruptured" contains substring "ruptured" — classify unruptured first
        if "unruptured" in s:
            pheno.append("unruptured")
        elif "ruptured" in s:
            pheno.append("ruptured")
        else:
            pheno.append("unknown")

    smd = pd.DataFrame({"Sample": samples, "Phenotype": pheno, "Title": titles})
    smd = smd[smd["Phenotype"] != "unknown"].copy()
    r_cols = smd.loc[smd["Phenotype"] == "ruptured", "Sample"].tolist()
    u_cols = smd.loc[smd["Phenotype"] == "unruptured", "Sample"].tolist()

    results = []
    for gene, probes in GENE_PROBES.items():
        present = [p for p in probes if p in expr.index]
        if not present:
            results.append(
                {
                    "Gene": gene,
                    "Probes": ";".join(probes),
                    "Status": "probe_missing",
                    "Mean_ruptured": np.nan,
                    "Mean_unruptured": np.nan,
                    "Log2FC_RvU": np.nan,
                    "P_value": np.nan,
                    "Expected": EXPECTED[gene],
                    "Match": "",
                }
            )
            continue
        sub = expr.loc[present, r_cols + u_cols]
        if len(present) > 1:
            vals_r = sub[r_cols].mean(axis=0).values.astype(float)
            vals_u = sub[u_cols].mean(axis=0).values.astype(float)
        else:
            vals_r = sub[r_cols].values.flatten().astype(float)
            vals_u = sub[u_cols].values.flatten().astype(float)

        mr, mu = np.nanmean(vals_r), np.nanmean(vals_u)
        log2fc = np.log2((mr + 1) / (mu + 1))
        _, pval = ttest_ind(vals_r, vals_u, equal_var=False, nan_policy="omit")

        if EXPECTED[gene] == "higher_in_ruptured":
            match = log2fc > 0
        else:
            match = log2fc < 0

        results.append(
            {
                "Gene": gene,
                "Probes": ";".join(present),
                "Status": "ok",
                "Mean_ruptured": mr,
                "Mean_unruptured": mu,
                "Log2FC_RvU": log2fc,
                "P_value": pval,
                "Expected": EXPECTED[gene],
                "Match": bool(match),
            }
        )

    res_df = pd.DataFrame(results)
    res_df.to_csv(OUT / "gse13353_key_gene_validation.tsv", sep="\t", index=False)
    smd.to_csv(OUT / "gse13353_sample_phenotype.tsv", sep="\t", index=False)

    n_match = int(res_df["Match"].sum())
    summary = {
        "dataset": "GSE13353",
        "platform": "GPL570",
        "n_ruptured": len(r_cols),
        "n_unruptured": len(u_cols),
        "genes_tested": len(GENE_PROBES),
        "direction_matches": n_match,
        "fraction_match": n_match / len(GENE_PROBES),
    }
    pd.Series(summary).to_csv(OUT / "gse13353_validation_summary.tsv", sep="\t", header=False)

    print("GSE13353 external validation (key genes)")
    print(summary)
    print(res_df.to_string(index=False))


if __name__ == "__main__":
    main()
