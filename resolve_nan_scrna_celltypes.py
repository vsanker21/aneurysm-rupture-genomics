"""
Resolve NaN scRNA best-cell-type assignments by adding synonym/alias rows to the
marker-strength table, based on MyGene.info symbol/alias history.

Inputs:
  - results/phases_g/tier1_supplement_with_scrna_specificity_cellmarker.tsv
  - data/reference/scrna_celltype_marker_strengths_cellmarker.tsv

Outputs:
  - data/reference/scrna_celltype_marker_strengths_cellmarker_patched.tsv
  - results/phases_g/nan_scrna_resolution_report.tsv

Strategy:
  1) Identify top40 Tier-1 symbols whose scrna_best_cell_type is NaN.
  2) For each symbol, query MyGene.info for aliases/prev symbols.
  3) If any alias is present in the marker table, copy that alias’ (cell_type, strength)
     rows onto the queried symbol as additional rows (so it can score).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import requests


BASE = Path(r"e:\Aneurysms Genomics Project")
TIER1 = BASE / r"results\phases_g\tier1_supplement_with_scrna_specificity_cellmarker.tsv"
MARKERS = BASE / r"data\reference\scrna_celltype_marker_strengths_cellmarker.tsv"
OUT_MARKERS = BASE / r"data\reference\scrna_celltype_marker_strengths_cellmarker_patched.tsv"
OUT_REPORT = BASE / r"results\phases_g\nan_scrna_resolution_report.tsv"


def mygene_symbol_aliases(symbol: str) -> set[str]:
    url = "https://mygene.info/v3/query"
    params = {
        "q": symbol,
        "scopes": "symbol,alias",
        "fields": "symbol,alias,prev_symbol,ensembl.gene,entrezgene",
        "species": "human",
        "size": 5,
    }
    r = requests.get(url, params=params, timeout=30)
    r.raise_for_status()
    data = r.json()
    hits = data.get("hits", []) if isinstance(data, dict) else []
    out: set[str] = set()
    for h in hits:
        if not isinstance(h, dict):
            continue
        sym = h.get("symbol")
        if sym:
            out.add(str(sym).upper())
        alias = h.get("alias", [])
        if isinstance(alias, str):
            out.add(alias.upper())
        elif isinstance(alias, list):
            out |= {str(a).upper() for a in alias if a}
        prev = h.get("prev_symbol", [])
        if isinstance(prev, str):
            out.add(prev.upper())
        elif isinstance(prev, list):
            out |= {str(a).upper() for a in prev if a}
    out.discard(symbol.upper())
    return out


def main() -> None:
    tier = pd.read_csv(TIER1, sep="\t")
    # Load specificity columns by re-deriving from the appended shortlist is unnecessary here;
    # we just need a stable list of symbols in the top40.
    symbols = tier["symbol"].dropna().astype(str).str.upper().tolist()

    marker = pd.read_csv(MARKERS, sep="\t")
    marker["gene_symbol"] = marker["gene_symbol"].astype(str).str.upper().str.strip()
    marker_genes = set(marker["gene_symbol"])

    # Determine which of these symbols are missing in the marker table
    missing = [s for s in symbols if s not in marker_genes]

    patched_rows = []
    report_rows = []
    # Pre-index marker rows by gene for fast copying
    by_gene = {g: marker[marker["gene_symbol"] == g] for g in marker_genes}

    for sym in missing:
        try:
            aliases = mygene_symbol_aliases(sym)
        except Exception as e:
            report_rows.append({"symbol": sym, "status": "query_failed", "detail": str(e)[:200]})
            continue

        # Find a usable alias present in marker table
        candidates = [a for a in aliases if a in marker_genes]
        if not candidates:
            report_rows.append({"symbol": sym, "status": "no_alias_match", "detail": ";".join(sorted(list(aliases))[:20])})
            continue

        # Choose alias with most rows (most supported cell types)
        best_alias = max(candidates, key=lambda a: len(by_gene[a]))
        src = by_gene[best_alias].copy()
        src["gene_symbol"] = sym
        src["source"] = src["source"].astype(str) + f"|alias_of:{best_alias}"
        patched_rows.append(src)
        report_rows.append({"symbol": sym, "status": "patched", "detail": f"alias={best_alias} rows={len(src)}"})

    patched = pd.concat([marker] + patched_rows, ignore_index=True) if patched_rows else marker
    patched = (
        patched.groupby(["gene_symbol", "cell_type"], as_index=False)["marker_strength"]
        .max()
        .assign(source="CellMarker_aggregated+aliases")
    )
    OUT_MARKERS.parent.mkdir(parents=True, exist_ok=True)
    patched.to_csv(OUT_MARKERS, sep="\t", index=False)

    rep = pd.DataFrame(report_rows)
    rep.to_csv(OUT_REPORT, sep="\t", index=False)

    print(f"Wrote {OUT_MARKERS}")
    print(f"Wrote {OUT_REPORT}")
    print(f"Tier1 symbols missing from marker table: {len(missing)}; patched: {(rep['status']=='patched').sum() if len(rep) else 0}")


if __name__ == "__main__":
    main()

