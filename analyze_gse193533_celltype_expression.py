"""
Exploratory: infer cell-type expression of unresolved genes using mouse scRNA-seq (GSE193533).

Rationale:
  Some Tier-1 genes remain "cell-type unresolved" in human marker-table approaches.
  As a stopgap, we can quantify which broad cell compartment expresses these genes in a
  mouse aneurysm model scRNA dataset (GSE193533). This is NOT equivalent to human IA wall
  cell-of-origin, but can provide supportive evidence for "immune vs vessel-wall" tendency.

Inputs:
  - data/sc_rnaseq/GSE193533/*_filtered_feature_bc_matrix.h5 (3 samples)

Outputs:
  - results/phases_g/gse193533_mouse_celltype_cluster_map.tsv
  - results/phases_g/gse193533_mouse_celltype_gene_expression.tsv

Method:
  - Load the 10x filtered h5 matrices, concatenate
  - Normalize/log1p, PCA, neighbors, Leiden clustering
  - Coarse cluster annotation by mean expression of canonical marker sets:
      EC / SMC / fibroblast / immune / other
  - For a target gene list (symbols), compute mean log-normalized expression per annotated cell type
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


BASE = Path(r"e:\Aneurysms Genomics Project")
IN_DIR = BASE / r"data\sc_rnaseq\GSE193533"
OUT_DIR = BASE / r"results\phases_g"


TARGET_HUMAN = ["MARCO", "PFKFB4", "SHF", "OASL", "COQ2", "RACGAP1", "ADAMTSL2", "SLC38A7", "MAPK13"]

HUMAN_TO_MOUSE_CANDIDATES = {
    # Common case: same symbol except mouse capitalization
    # Special cases / families:
    "OASL": ["Oasl", "Oasl1", "Oasl2"],
}


MARKERS = {
    "EC": ["Pecam1", "Cdh5", "Vwf", "Kdr", "Flt1", "Robo4", "Cldn5"],
    "SMC": ["Acta2", "Tagln", "Myh11", "Cnn1", "Mylk", "Smtn"],
    "FIBROBLAST": ["Dcn", "Lum", "Col1a1", "Col1a2", "Col3a1", "Postn", "Fn1"],
    "IMMUNE": ["Ptprc", "Lst1", "Tyrobp", "Fcer1g", "C1qa", "C1qb", "C1qc", "S100a8", "S100a9"],
}


def _present(var_names: pd.Index, genes: list[str]) -> list[str]:
    s = set(var_names.astype(str))
    return [g for g in genes if g in s]


def annotate_clusters(adata: sc.AnnData, cluster_key: str = "leiden") -> pd.DataFrame:
    """
    Returns a table of cluster -> best cell type + per-type scores.
    """
    # Ensure log-normalized data is in X
    clus = adata.obs[cluster_key].astype(str)
    df = []
    for c in sorted(clus.unique(), key=lambda x: int(x) if x.isdigit() else x):
        idx = clus == c
        scores = {}
        for ct, genes in MARKERS.items():
            present = _present(adata.var_names, genes)
            if not present:
                scores[ct] = np.nan
                continue
            # mean of mean expressions across markers
            sub = adata[idx, present]
            scores[ct] = float(np.asarray(sub.X.mean(axis=0)).ravel().mean())
        best = max(
            [k for k, v in scores.items() if np.isfinite(v)],
            key=lambda k: scores[k],
            default="OTHER",
        )
        df.append({"cluster": c, "best_cell_type": best, **{f"score_{k}": v for k, v in scores.items()}})
    return pd.DataFrame(df)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    h5s = sorted(IN_DIR.glob("*_filtered_feature_bc_matrix.h5"))
    if not h5s:
        raise FileNotFoundError(f"No filtered_feature_bc_matrix.h5 files found in {IN_DIR}")

    adatas = []
    for p in h5s:
        a = sc.read_10x_h5(str(p))
        a.var_names_make_unique()
        a.obs["sample"] = p.stem
        adatas.append(a)

    adata = adatas[0].concatenate(*adatas[1:], batch_key="batch", batch_categories=[a.obs["sample"][0] for a in adatas])

    # Basic QC
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=5)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # Preserve full gene space for downstream expression summaries
    adata.raw = adata

    # Use HVGs for embedding/clustering only (keep full matrix in .raw)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3", subset=False)
    adata_hvg = adata[:, adata.var["highly_variable"].values].copy()

    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, svd_solver="arpack")
    sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30)
    sc.tl.leiden(adata_hvg, resolution=0.6, key_added="leiden")
    adata.obs["leiden"] = adata_hvg.obs["leiden"].astype(str).values

    cmap = annotate_clusters(adata, "leiden")
    # Assign cell type to each cell
    m = cmap.set_index("cluster")["best_cell_type"].to_dict()
    adata.obs["cell_type"] = adata.obs["leiden"].astype(str).map(lambda x: m.get(x, "OTHER"))

    cmap_out = OUT_DIR / "gse193533_mouse_celltype_cluster_map.tsv"
    cmap.to_csv(cmap_out, sep="\t", index=False)

    # Compute per-cell-type mean expression for target genes (try mouse-style symbols)
    # Most mouse gene symbols are TitleCase (e.g. Marco).
    # For gene families, use explicit candidate lists.
    var_set = set(adata.raw.var_names.astype(str))
    target_mouse_primary = {}
    for g in TARGET_HUMAN:
        cands = HUMAN_TO_MOUSE_CANDIDATES.get(g, [g.capitalize()])
        chosen = next((c for c in cands if c in var_set), None)
        target_mouse_primary[g] = chosen
    present = [m for m in target_mouse_primary.values() if m is not None]

    rows = []
    for ct in sorted(adata.obs["cell_type"].unique()):
        idx = adata.obs["cell_type"] == ct
        for g_h in TARGET_HUMAN:
            g_m = target_mouse_primary.get(g_h)
            if not g_m:
                rows.append({"gene_human": g_h, "gene_mouse": g_m, "cell_type": ct, "mean_log1p": np.nan})
                continue
            # Use full gene matrix stored in raw
            val = float(np.asarray(adata.raw[idx, g_m].X.mean()).ravel()[0])
            rows.append({"gene_human": g_h, "gene_mouse": g_m, "cell_type": ct, "mean_log1p": val})

    out = pd.DataFrame(rows)
    out_path = OUT_DIR / "gse193533_mouse_celltype_gene_expression.tsv"
    out.to_csv(out_path, sep="\t", index=False)

    print(f"Wrote {cmap_out}")
    print(f"Wrote {out_path}")
    print("Cell-type counts:", adata.obs["cell_type"].value_counts().to_string())
    print("Target genes present:", present)


if __name__ == "__main__":
    main()

