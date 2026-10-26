"""
Export GSE122897 DESeq2 (ruptured vs unruptured) to a TSV for link_tissue_de_to_phases_g.py.

Columns: symbol, signature_group, tissue_log2FC, tissue_padj
Maps Entrez Gene_ID -> HGNC symbol via local tables first, then mygene.info batches.

Usage:
  python export_gse122897_cohort_de_for_bridge.py
  python export_gse122897_cohort_de_for_bridge.py --max-padj 0.05
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import requests

BASE = Path(__file__).resolve().parent
DESEQ2 = BASE / "results/gse122897_deseq2_results.tsv"
OUT = BASE / "results/phases_g/cohort_de_gse122897_for_phases_g_bridge.tsv"


def load_local_entrez_to_symbol() -> dict[str, str]:
    m: dict[str, str] = {}
    for f in [
        BASE / "results/gene_symbol_to_entrez_mapping.tsv",
        BASE / "results/gene_id_mapping_unified.tsv",
    ]:
        if not f.is_file():
            continue
        df = pd.read_csv(f, sep="\t")
        if "Entrez_ID" in df.columns and "Gene_Symbol" in df.columns:
            for _, r in df.dropna(subset=["Gene_Symbol"]).iterrows():
                try:
                    val = r["Entrez_ID"]
                    if isinstance(val, str) and str(val).startswith("ENSG"):
                        continue
                    eid = str(int(float(val)))
                except (ValueError, TypeError):
                    continue
                sym = str(r["Gene_Symbol"]).strip()
                if sym and sym.lower() != "nan":
                    m[eid] = sym
    sl = BASE / "results/phases_g/shortlist_evidence.tsv"
    if sl.is_file():
        s = pd.read_csv(sl, sep="\t")
        if "Gene_ID" in s.columns and "symbol" in s.columns:
            for _, r in s.iterrows():
                try:
                    eid = str(int(float(r["Gene_ID"])))
                except (ValueError, TypeError):
                    continue
                sym = str(r["symbol"]).strip()
                if sym:
                    m.setdefault(eid, sym)
    return m


def fetch_symbols_mygene(entrez_ids: list[str], batch_size: int = 900) -> dict[str, str]:
    out: dict[str, str] = {}
    url = "https://mygene.info/v3/gene"
    for i in range(0, len(entrez_ids), batch_size):
        batch = entrez_ids[i : i + batch_size]
        try:
            r = requests.post(
                url,
                json={"ids": ",".join(batch), "fields": "symbol,entrezgene"},
                timeout=120,
            )
            r.raise_for_status()
            rows = r.json()
            if not isinstance(rows, list):
                continue
            for row in rows:
                if not isinstance(row, dict):
                    continue
                eid = row.get("entrezgene") or row.get("_id")
                sym = row.get("symbol")
                if eid is not None and sym:
                    out[str(int(eid))] = str(sym)
        except (requests.RequestException, ValueError, TypeError) as e:
            print(f"[warn] mygene batch {i}: {e}", file=sys.stderr)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--deseq2", type=Path, default=DESEQ2)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument(
        "--max-padj",
        type=float,
        default=0.05,
        help="Export genes with DESeq2 padj below this (default 0.05). Use 1 for all genes (large file).",
    )
    ap.add_argument("--skip-mygene", action="store_true", help="Only use local symbol maps.")
    args = ap.parse_args()

    if not args.deseq2.is_file():
        print(f"Missing {args.deseq2}", file=sys.stderr)
        return 1

    df = pd.read_csv(args.deseq2, sep="\t")
    df.columns = [str(c).strip('"') for c in df.columns]
    df["Gene_ID"] = df["Gene_ID"].astype(str).str.replace(".0", "", regex=False)
    df["padj"] = pd.to_numeric(df["padj"], errors="coerce")
    df["log2FoldChange"] = pd.to_numeric(df["log2FoldChange"], errors="coerce")
    df = df.dropna(subset=["padj", "log2FoldChange"])
    if args.max_padj < 1:
        df = df[df["padj"] <= float(args.max_padj)]

    local = load_local_entrez_to_symbol()
    ids = df["Gene_ID"].unique().tolist()
    need = [i for i in ids if i not in local]
    if need and not args.skip_mygene:
        print(f"[export] mygene lookup for {len(need)} Entrez IDs not in local maps...", flush=True)
        fetched = fetch_symbols_mygene(need)
        local.update(fetched)
    elif need and args.skip_mygene:
        print(f"[warn] {len(need)} IDs lack local symbol; omit or use --skip-mygene false", file=sys.stderr)

    rows = []
    for _, r in df.iterrows():
        gid = str(r["Gene_ID"])
        sym = local.get(gid)
        if not sym:
            continue
        lfc = float(r["log2FoldChange"])
        padj = float(r["padj"])
        grp = "upregulated_in_ruptured" if lfc > 0 else "downregulated_in_ruptured"
        rows.append(
            {
                "symbol": sym.upper(),
                "signature_group": grp,
                "tissue_log2FC": lfc,
                "tissue_padj": padj,
            }
        )

    out_df = pd.DataFrame(rows).drop_duplicates(subset=["symbol"], keep="first")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {args.out} ({len(out_df)} genes, padj<={args.max_padj})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
