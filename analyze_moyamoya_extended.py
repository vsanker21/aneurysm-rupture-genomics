import csv
import gzip
import json
import tarfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import seaborn as sns
import statsmodels.api as sm
from scipy.stats import ttest_ind
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.model_selection import GridSearchCV, StratifiedKFold, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from statsmodels.stats.multitest import multipletests


BASE = Path(r"e:\Aneurysms Genomics Project\Moya Moya")
OUT = Path(r"e:\Aneurysms Genomics Project\results\moyamoya")
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)


def parse_series_matrix(series_path: Path):
    meta_rows = {}
    expr_header = None
    expr_data = []
    in_table = False

    with gzip.open(series_path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                in_table = False
                continue

            if in_table:
                if expr_header is None:
                    expr_header = next(csv.reader([line], delimiter="\t", quotechar='"'))
                else:
                    row = next(csv.reader([line], delimiter="\t", quotechar='"'))
                    expr_data.append(row)
                continue

            if line.startswith("!Sample_"):
                parts = next(csv.reader([line], delimiter="\t", quotechar='"'))
                key = parts[0].replace("!Sample_", "")
                vals = parts[1:]
                meta_rows.setdefault(key, []).append(vals)

    samples = expr_header[1:]
    expr_df = pd.DataFrame(expr_data, columns=expr_header)
    expr_df = expr_df.rename(columns={expr_df.columns[0]: "ID_REF"})
    for c in samples:
        expr_df[c] = pd.to_numeric(expr_df[c], errors="coerce")
    expr_df["ID_REF"] = pd.to_numeric(expr_df["ID_REF"], errors="coerce").astype("Int64")
    expr_df = expr_df.dropna(subset=["ID_REF"]).copy()

    sample_md = pd.DataFrame({"Sample": samples})

    def get_row(key, idx=0):
        rows = meta_rows.get(key, [])
        if len(rows) > idx:
            return rows[idx]
        return [""] * len(samples)

    sample_md["Title"] = get_row("title", 0)
    sample_md["Source"] = get_row("source_name_ch1", 0)
    sample_md["Tissue"] = get_row("characteristics_ch1", 0)
    sample_md["Sex_raw"] = get_row("characteristics_ch1", 1)
    sample_md["Age_raw"] = get_row("characteristics_ch1", 2)

    def infer_group(title):
        t = str(title)
        if "[MMD" in t:
            return "MMD"
        if "[IA" in t:
            return "IA"
        if "[EPI" in t:
            return "EPI"
        return "Unknown"

    sample_md["Group"] = sample_md["Title"].map(infer_group)
    sample_md["Sex"] = (
        sample_md["Sex_raw"]
        .str.replace("gender:", "", regex=False)
        .str.strip()
        .replace({"Fenale": "Female"})
    )
    sample_md["Age"] = (
        sample_md["Age_raw"]
        .str.replace("age:", "", regex=False)
        .str.replace("y", "", regex=False)
        .str.strip()
    )
    sample_md["Age"] = pd.to_numeric(sample_md["Age"], errors="coerce")
    sample_md["SexMale"] = (sample_md["Sex"] == "Male").astype(int)
    return expr_df, sample_md


def parse_probe_map(raw_tar: Path):
    with tarfile.open(raw_tar, "r") as t:
        m = t.extractfile("GSM4771802_MMD1.txt.gz")
        txt = gzip.decompress(m.read()).decode("utf-8", "ignore").splitlines()

    feat_header = None
    feat_rows = []
    for ln in txt:
        if ln.startswith("FEATURES\t"):
            feat_header = ln.split("\t")
            continue
        if feat_header and ln.startswith("DATA\t"):
            vals = ln.split("\t")
            if len(vals) == len(feat_header):
                feat_rows.append(vals)

    fdf = pd.DataFrame(feat_rows, columns=feat_header)
    fdf["ProbeUID"] = pd.to_numeric(fdf["ProbeUID"], errors="coerce").astype("Int64")
    probe_map = fdf[
        ["ProbeUID", "ProbeName", "GeneName", "SystematicName", "Description"]
    ].dropna(subset=["ProbeUID"])
    probe_map = probe_map.rename(columns={"ProbeUID": "ID_REF"})
    return probe_map


def curated_symbol(x: str):
    s = str(x).split(";")[0].strip()
    if s in {"", "nan", "None"}:
        return None
    bad_prefix = ("A_", "ENST", "LOC", "lnc-", "XLOC", "TCONS", "AK", "AF")
    if s.startswith(bad_prefix):
        return None
    if not s.replace("-", "").replace(".", "").isalnum():
        return None
    return s


def run_de(expr_matrix, sample_md, probe_map, case_group, ctrl_group, prefix):
    case_samples = sample_md.loc[sample_md["Group"] == case_group, "Sample"].tolist()
    ctrl_samples = sample_md.loc[sample_md["Group"] == ctrl_group, "Sample"].tolist()

    Xa = expr_matrix[case_samples]
    Xb = expr_matrix[ctrl_samples]
    mean_a = Xa.mean(axis=1)
    mean_b = Xb.mean(axis=1)
    log2fc = mean_a - mean_b
    _, pvals = ttest_ind(Xa.values, Xb.values, axis=1, equal_var=False, nan_policy="omit")

    res = pd.DataFrame(
        {
            "ID_REF": expr_matrix.index,
            "Mean_case": mean_a.values,
            "Mean_ctrl": mean_b.values,
            "Log2FC_case_vs_ctrl": log2fc.values,
            "P_value": pvals,
        }
    )
    res = res.replace([np.inf, -np.inf], np.nan).dropna(subset=["P_value"])
    res["FDR_BH"] = multipletests(res["P_value"].values, method="fdr_bh")[1]
    res["P_Bonferroni"] = multipletests(res["P_value"].values, method="bonferroni")[1]
    res = res.sort_values("P_value").merge(probe_map, on="ID_REF", how="left")
    res["GeneSymbol"] = res["GeneName"].map(curated_symbol)
    res.to_csv(OUT / f"{prefix}_de.tsv", sep="\t", index=False)

    top = res[res["GeneSymbol"].notna()].drop_duplicates("GeneSymbol", keep="first")
    top.to_csv(OUT / f"{prefix}_top_curated.tsv", sep="\t", index=False)
    return res, top, case_samples, ctrl_samples


def run_adjusted_lm(expr_matrix, sample_md, de_df, case_group, ctrl_group, prefix):
    sub_md = sample_md[sample_md["Group"].isin([case_group, ctrl_group])].copy()
    sub_md["Y"] = (sub_md["Group"] == case_group).astype(int)
    sub_md = sub_md.dropna(subset=["Age"])

    covars = ["Y", "Age", "SexMale"]
    # Drop constant covariates in subset to avoid singular designs
    covars = [c for c in covars if sub_md[c].nunique() > 1]
    design = sm.add_constant(sub_md[covars], has_constant="add")

    probes = de_df["ID_REF"].head(5000).tolist()
    X = expr_matrix.loc[expr_matrix.index.isin(probes), sub_md["Sample"]]
    pvals = []
    betas = []
    for _, row in X.iterrows():
        y = row.values
        try:
            fit = sm.OLS(y, design.values).fit()
            # coefficient for Y (case indicator)
            y_idx = list(design.columns).index("Y")
            betas.append(fit.params[y_idx])
            pvals.append(fit.pvalues[y_idx])
        except Exception:
            betas.append(np.nan)
            pvals.append(np.nan)

    lm = pd.DataFrame(
        {
            "ID_REF": X.index,
            "Beta_case_indicator": betas,
            "P_value_adj": pvals,
        }
    ).dropna(subset=["P_value_adj"])
    lm["FDR_BH_adj"] = multipletests(lm["P_value_adj"].values, method="fdr_bh")[1]
    lm = lm.merge(de_df[["ID_REF", "GeneSymbol", "GeneName"]], on="ID_REF", how="left")
    lm.to_csv(OUT / f"{prefix}_age_sex_adjusted_lm.tsv", sep="\t", index=False)
    return lm, covars


def run_pathway_strict(top_df, prefix):
    symbols = top_df["GeneSymbol"].dropna().drop_duplicates().tolist()
    if len(symbols) < 20:
        return pd.DataFrame()

    payload = {
        "organism": "hsapiens",
        "query": symbols[:500],
        "sources": ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"],
        "ordered": True,
        "all_results": True,
        "user_threshold": 0.05,
        "no_iea": False,
        "domain_scope": "annotated",
    }
    r = requests.post("https://biit.cs.ut.ee/gprofiler/api/gost/profile", json=payload, timeout=120)
    rows = []
    if r.status_code == 200:
        for it in r.json().get("result", []):
            inter = it.get("intersections", [[]])
            ol = inter[0] if inter else []
            rows.append(
                {
                    "Source": it.get("source", ""),
                    "Term_ID": it.get("native", it.get("id", "")),
                    "Term_Name": it.get("name", ""),
                    "P_value": it.get("p_value", 1.0),
                    "Significant": it.get("significant", False),
                    "Term_Size": it.get("term_size", 0),
                    "Query_Size": it.get("query_size", 0),
                    "Overlap_Count": len(ol),
                    "Overlap_Genes": ",".join(ol[:30]),
                }
            )
    enr = pd.DataFrame(rows)
    if len(enr):
        enr = enr.sort_values("P_value")
        enr["FDR_BH"] = multipletests(enr["P_value"].values, method="fdr_bh")[1]
        enr.to_csv(OUT / f"{prefix}_pathway_strict_ranked.tsv", sep="\t", index=False)
        enr[enr["FDR_BH"] < 0.05].to_csv(
            OUT / f"{prefix}_pathway_strict_fdr05.tsv", sep="\t", index=False
        )
    return enr


def make_embeddings(expr_matrix, sample_md):
    use_samples = sample_md["Sample"].tolist()
    X = expr_matrix[use_samples].T
    Xs = StandardScaler().fit_transform(X)

    pca = PCA(n_components=2, random_state=42)
    p = pca.fit_transform(Xs)
    pca_df = pd.DataFrame(
        {"Sample": use_samples, "PC1": p[:, 0], "PC2": p[:, 1]}
    ).merge(sample_md[["Sample", "Group", "Age", "Sex"]], on="Sample", how="left")
    pca_df.to_csv(OUT / "gse157628_pca_coords.tsv", sep="\t", index=False)

    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="Group", style="Sex", s=90)
    for _, r in pca_df.iterrows():
        plt.text(r["PC1"], r["PC2"], r["Sample"][-3:], fontsize=7, alpha=0.7)
    plt.title("GSE157628 PCA (all samples)")
    plt.tight_layout()
    plt.savefig(FIG / "gse157628_pca.png", dpi=300)
    plt.close()

    tsne = TSNE(n_components=2, random_state=42, perplexity=5)
    t = tsne.fit_transform(Xs)
    tsne_df = pd.DataFrame(
        {"Sample": use_samples, "tSNE1": t[:, 0], "tSNE2": t[:, 1]}
    ).merge(sample_md[["Sample", "Group", "Age", "Sex"]], on="Sample", how="left")
    tsne_df.to_csv(OUT / "gse157628_tsne_coords.tsv", sep="\t", index=False)

    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=tsne_df, x="tSNE1", y="tSNE2", hue="Group", style="Sex", s=90)
    plt.title("GSE157628 t-SNE (all samples)")
    plt.tight_layout()
    plt.savefig(FIG / "gse157628_tsne.png", dpi=300)
    plt.close()

    umap_status = "not_available"
    try:
        import umap  # type: ignore

        reducer = umap.UMAP(n_components=2, random_state=42, n_neighbors=5, min_dist=0.2)
        u = reducer.fit_transform(Xs)
        umap_df = pd.DataFrame(
            {"Sample": use_samples, "UMAP1": u[:, 0], "UMAP2": u[:, 1]}
        ).merge(sample_md[["Sample", "Group", "Age", "Sex"]], on="Sample", how="left")
        umap_df.to_csv(OUT / "gse157628_umap_coords.tsv", sep="\t", index=False)
        plt.figure(figsize=(8, 6))
        sns.scatterplot(data=umap_df, x="UMAP1", y="UMAP2", hue="Group", style="Sex", s=90)
        plt.title("GSE157628 UMAP (all samples)")
        plt.tight_layout()
        plt.savefig(FIG / "gse157628_umap.png", dpi=300)
        plt.close()
        umap_status = "ok"
    except Exception:
        umap_status = "failed_or_missing_package"
    return pca_df, tsne_df, umap_status


def classifier_auc(expr_matrix, sample_md, case_group, ctrl_group, prefix):
    sub = sample_md[sample_md["Group"].isin([case_group, ctrl_group])].copy()
    y = (sub["Group"] == case_group).astype(int).values
    samples = sub["Sample"].tolist()
    X = expr_matrix[samples].T
    k_cap = min(100, X.shape[1])
    pipe = Pipeline(
        [
            ("scaler", StandardScaler()),
            ("kbest", SelectKBest(score_func=f_classif)),
            ("clf", LogisticRegression(max_iter=2000, solver="liblinear")),
        ]
    )

    min_class_n = int(sub["Group"].value_counts().min())
    if min_class_n < 5:
        pd.DataFrame(
            {
                "Metric": [
                    "Nested_CV_AUC_Mean",
                    "Nested_CV_AUC_Std",
                    "NonNested_CV_AUC_Mean",
                    "NonNested_CV_AUC_Std",
                    "N",
                    "N_case",
                    "N_ctrl",
                    "Outer_Folds",
                    "Inner_Folds",
                    "Residualized_On_Age_Sex",
                    "Status",
                ],
                "Value": [
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    len(sub),
                    int((sub["Group"] == case_group).sum()),
                    int((sub["Group"] == ctrl_group).sum()),
                    np.nan,
                    np.nan,
                    0,
                    "insufficient_controls_for_reliable_cv",
                ],
            }
        ).to_csv(OUT / f"{prefix}_classifier_cv.tsv", sep="\t", index=False)
        return np.array([]), np.array([])

    outer_splits = max(2, min(5, min_class_n))
    inner_splits = max(2, min(3, min_class_n - 1)) if min_class_n > 2 else 2
    outer_cv = StratifiedKFold(n_splits=outer_splits, shuffle=True, random_state=42)
    inner_cv = StratifiedKFold(n_splits=inner_splits, shuffle=True, random_state=1)

    k_grid = sorted(set([5, 10, 20, min(50, k_cap), k_cap]))
    param_grid = {
        "kbest__k": [k for k in k_grid if 1 <= k <= X.shape[1]],
        "clf__C": [0.1, 1.0, 10.0],
    }
    search = GridSearchCV(
        estimator=pipe,
        param_grid=param_grid,
        scoring="roc_auc",
        cv=inner_cv,
        n_jobs=-1,
    )
    nested_auc = cross_val_score(search, X, y, cv=outer_cv, scoring="roc_auc")

    # Keep a non-nested estimate for transparency, but prefer nested AUC for claims.
    non_nested_auc = cross_val_score(pipe, X, y, cv=outer_cv, scoring="roc_auc")

    pd.DataFrame(
        {
            "Metric": [
                "Nested_CV_AUC_Mean",
                "Nested_CV_AUC_Std",
                "NonNested_CV_AUC_Mean",
                "NonNested_CV_AUC_Std",
                "N",
                "N_case",
                "N_ctrl",
                "Outer_Folds",
                "Inner_Folds",
                "Residualized_On_Age_Sex",
                "Status",
            ],
            "Value": [
                nested_auc.mean(),
                nested_auc.std(ddof=1) if len(nested_auc) > 1 else 0.0,
                non_nested_auc.mean(),
                non_nested_auc.std(ddof=1) if len(non_nested_auc) > 1 else 0.0,
                len(sub),
                int((sub["Group"] == case_group).sum()),
                int((sub["Group"] == ctrl_group).sum()),
                outer_splits,
                inner_splits,
                0,
                "ok",
            ],
        }
    ).to_csv(OUT / f"{prefix}_classifier_cv.tsv", sep="\t", index=False)
    return nested_auc, non_nested_auc


def main():
    series_path = BASE / "GSE157628_series_matrix.txt.gz"
    raw_tar = BASE / "GSE157628_RAW.tar"

    expr_df, sample_md = parse_series_matrix(series_path)
    probe_map = parse_probe_map(raw_tar)
    sample_md.to_csv(OUT / "gse157628_sample_metadata.tsv", sep="\t", index=False)
    probe_map.to_csv(OUT / "gse157628_probe_annotation.tsv", sep="\t", index=False)

    expr_matrix = expr_df.set_index("ID_REF")

    # Contrasts
    de_mia, top_mia, _, _ = run_de(expr_matrix, sample_md, probe_map, "MMD", "IA", "gse157628_mmd_vs_ia")
    de_mepi, top_mepi, _, _ = run_de(expr_matrix, sample_md, probe_map, "MMD", "EPI", "gse157628_mmd_vs_epi")

    # Adjusted linear models
    lm_mia, cov_mia = run_adjusted_lm(expr_matrix, sample_md, de_mia, "MMD", "IA", "gse157628_mmd_vs_ia")
    lm_mepi, cov_mepi = run_adjusted_lm(expr_matrix, sample_md, de_mepi, "MMD", "EPI", "gse157628_mmd_vs_epi")

    # Stricter pathway analysis on curated symbols
    enr_mia = run_pathway_strict(top_mia, "gse157628_mmd_vs_ia")
    enr_mepi = run_pathway_strict(top_mepi, "gse157628_mmd_vs_epi")

    # Embeddings and confound visualization
    pca_df, tsne_df, umap_status = make_embeddings(expr_matrix, sample_md)

    # Classifier AUC per contrast
    auc_mia_nested, auc_mia_non_nested = classifier_auc(
        expr_matrix, sample_md, "MMD", "IA", "gse157628_mmd_vs_ia"
    )
    auc_mepi_nested, auc_mepi_non_nested = classifier_auc(
        expr_matrix, sample_md, "MMD", "EPI", "gse157628_mmd_vs_epi"
    )

    sex_counts_by_group = (
        sample_md.groupby("Group")["Sex"].value_counts().unstack(fill_value=0).to_dict(orient="index")
    )

    summary = {
        "dataset": "GSE157628",
        "group_counts": sample_md["Group"].value_counts().to_dict(),
        "age_mean_by_group": sample_md.groupby("Group")["Age"].mean().round(2).to_dict(),
        "sex_counts_by_group": sex_counts_by_group,
        "mmd_vs_ia": {
            "de_fdr_005": int((de_mia["FDR_BH"] < 0.05).sum()),
            "de_p_001_absfc2": int(
                ((de_mia["P_value"] < 1e-3) & (np.abs(de_mia["Log2FC_case_vs_ctrl"]) > 1)).sum()
            ),
            "adj_lm_covariates": cov_mia,
            "adj_lm_fdr_005": int((lm_mia["FDR_BH_adj"] < 0.05).sum()),
            "auc_nested_mean": float(auc_mia_nested.mean()) if len(auc_mia_nested) else None,
            "auc_nested_std": float(auc_mia_nested.std(ddof=1)) if len(auc_mia_nested) > 1 else None,
            "auc_non_nested_mean": float(auc_mia_non_nested.mean()) if len(auc_mia_non_nested) else None,
            "auc_non_nested_std": float(auc_mia_non_nested.std(ddof=1))
            if len(auc_mia_non_nested) > 1
            else None,
            "pathway_fdr_005": int((enr_mia["FDR_BH"] < 0.05).sum()) if len(enr_mia) else 0,
        },
        "mmd_vs_epi": {
            "de_fdr_005": int((de_mepi["FDR_BH"] < 0.05).sum()),
            "de_p_001_absfc2": int(
                ((de_mepi["P_value"] < 1e-3) & (np.abs(de_mepi["Log2FC_case_vs_ctrl"]) > 1)).sum()
            ),
            "adj_lm_covariates": cov_mepi,
            "adj_lm_fdr_005": int((lm_mepi["FDR_BH_adj"] < 0.05).sum()),
            "auc_nested_mean": float(auc_mepi_nested.mean()) if len(auc_mepi_nested) else None,
            "auc_nested_std": float(auc_mepi_nested.std(ddof=1)) if len(auc_mepi_nested) > 1 else None,
            "auc_non_nested_mean": float(auc_mepi_non_nested.mean()) if len(auc_mepi_non_nested) else None,
            "auc_non_nested_std": float(auc_mepi_non_nested.std(ddof=1))
            if len(auc_mepi_non_nested) > 1
            else None,
            "pathway_fdr_005": int((enr_mepi["FDR_BH"] < 0.05).sum()) if len(enr_mepi) else 0,
        },
        "embedding_files": {
            "pca": str(OUT / "gse157628_pca_coords.tsv"),
            "tsne": str(OUT / "gse157628_tsne_coords.tsv"),
            "umap_status": umap_status,
            "figure_dir": str(FIG),
        },
    }

    with open(OUT / "gse157628_extended_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
