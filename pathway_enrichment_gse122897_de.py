"""
Pathway enrichment for GSE122897 differentially expressed genes.
Uses g:Profiler API for GO, KEGG, Reactome with FDR control.
Implements: (1) ranked gene list (ordered=TRUE, GSEA-style), (2) threshold-based (padj<0.05).
"""
import pandas as pd
import numpy as np
import requests
import os
import time

def load_entrez_to_symbol():
    """Build Entrez ID -> Gene Symbol mapping from local files."""
    m = {}
    for f in ["results/gene_symbol_to_entrez_mapping.tsv", "results/gene_id_mapping_unified.tsv"]:
        if not os.path.exists(f):
            continue
        df = pd.read_csv(f, sep='\t')
        if 'Entrez_ID' in df.columns and 'Gene_Symbol' in df.columns:
            for _, r in df.dropna(subset=['Entrez_ID', 'Gene_Symbol']).iterrows():
                try:
                    val = r['Entrez_ID']
                    if isinstance(val, str) and val.startswith('ENSG'):
                        continue
                    eid = str(int(float(val))).replace('.0', '')
                except (ValueError, TypeError):
                    continue
                sym = str(r['Gene_Symbol']).strip()
                if sym and sym != 'nan':
                    m[eid] = sym
    return m

def fetch_symbols_mygene_batch(gene_ids, batch_size=400):
    """Fetch gene symbols from mygene.info API in batch (Entrez ID -> symbol)."""
    m = {}
    ids = [str(g).replace('.0', '') for g in gene_ids[:batch_size] if str(g).replace('.0', '').isdigit()]
    if not ids:
        return m
    # mygene query: q=entrezgene:ID1 OR entrezgene:ID2 ...
    q = " OR ".join(f"entrezgene:{i}" for i in ids[:300])
    try:
        r = requests.post(
            'https://mygene.info/v3/query',
            json={
                'q': q,
                'species': 'human',
                'fields': ['symbol', 'entrezgene'],
                'size': len(ids),
            },
            timeout=60
        )
        if r.status_code == 200:
            data = r.json()
            hits = data.get('hits', [])
            for h in hits:
                gid = str(h.get('entrezgene', h.get('_id', '')))
                sym = h.get('symbol', gid)
                if gid and sym:
                    m[gid] = sym
    except Exception as e:
        print(f"  [WARN] mygene batch failed: {e}")
    return m

def load_ranked_de_genes(use_deseq2=True):
    """
    Load full DE results and return ranked gene list (symbols) by -log10(p)*sign(FC).
    Uses DESeq2 for more genes; falls back to Welch.
    """
    e2s = load_entrez_to_symbol()

    # Prefer DESeq2 (raw counts, more power)
    if use_deseq2 and os.path.exists("results/gse122897_deseq2_results.tsv"):
        df = pd.read_csv("results/gse122897_deseq2_results.tsv", sep='\t')
        df.columns = [c.strip('"') for c in df.columns]
        df = df.dropna(subset=['padj', 'log2FoldChange'])
        df['rank_score'] = -np.log10(df['padj'].clip(1e-300)) * np.sign(df['log2FoldChange'])
        df = df.sort_values('rank_score', ascending=False)
        df['Gene_ID'] = df['Gene_ID'].astype(str).str.replace('"', '').str.replace('.0', '')
    else:
        df = pd.read_csv("results/gse122897_differential_expression.tsv", sep='\t')
        df = df.dropna(subset=['P_adj_BH', 'Log2FC'])
        df['rank_score'] = -np.log10(df['P_adj_BH'].clip(1e-300)) * np.sign(df['Log2FC'])
        df = df.sort_values('rank_score', ascending=False)
        df['Gene_ID'] = df['Gene'].astype(str).str.replace('.0', '')

    # Map to symbols
    missing = df[~df['Gene_ID'].isin(e2s)]['Gene_ID'].unique().tolist()[:500]
    if missing:
        api_map = fetch_symbols_mygene_batch(missing)
        e2s.update(api_map)

    df['Gene_Symbol'] = df['Gene_ID'].map(e2s)
    df = df[df['Gene_Symbol'].notna() & (df['Gene_Symbol'] != '')]
    # Dedupe by symbol, keep first (highest rank)
    df = df.drop_duplicates(subset=['Gene_Symbol'], keep='first')
    return df['Gene_Symbol'].tolist()[:500]

def load_threshold_genes(padj_thresh=0.05):
    """Load DE genes passing padj/FDR threshold (for over-representation)."""
    e2s = load_entrez_to_symbol()
    genes = []
    if os.path.exists("results/gse122897_deseq2_results.tsv"):
        df = pd.read_csv("results/gse122897_deseq2_results.tsv", sep='\t')
        df.columns = [c.strip('"') for c in df.columns]
        df = df[df['padj'] < padj_thresh].dropna(subset=['padj'])
        df['Gene_ID'] = df['Gene_ID'].astype(str).str.replace('"', '').str.replace('.0', '')
    else:
        df = pd.read_csv("results/gse122897_top_genes.tsv", sep='\t')
        df = df[df['P_adj_BH'] < padj_thresh].dropna(subset=['P_adj_BH'])
        df['Gene_ID'] = df['Gene'].astype(str).str.replace('.0', '')
    missing = df[~df['Gene_ID'].isin(e2s)]['Gene_ID'].unique().tolist()
    if missing:
        e2s.update(fetch_symbols_mygene_batch(missing))
    df['Gene_Symbol'] = df['Gene_ID'].map(e2s)
    genes = df[df['Gene_Symbol'].notna() & (df['Gene_Symbol'] != '')]['Gene_Symbol'].drop_duplicates().tolist()
    return genes[:500]

def gprofiler_enrichment(gene_list, ordered=True, fdr_threshold=0.05):
    """
    Run pathway enrichment via g:Profiler API.
    ordered=True: GSEA-style ranked analysis. ordered=False: over-representation.
    Sources: GO:BP, GO:MF, GO:CC, KEGG, REAC, WP.
    """
    url = "https://biit.cs.ut.ee/gprofiler/api/gost/profile"
    payload = {
        "organism": "hsapiens",
        "query": gene_list[:500],
        "sources": ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"],
        "user_threshold": fdr_threshold,
        "all_results": True,
        "ordered": ordered,
        "no_iea": False,
        "domain_scope": "annotated",
    }
    try:
        r = requests.post(url, json=payload, timeout=90)
        if r.status_code == 200:
            data = r.json()
            return data.get("result", [])
        return None
    except Exception as e:
        print(f"[ERROR] g:Profiler API failed: {e}")
        return None

def results_to_df(results):
    """Convert g:Profiler results to DataFrame with FDR."""
    rows = []
    for item in results:
        term_id = item.get("native", item.get("id", ""))
        source = item.get("source", "")
        name = item.get("name", "")
        pval = item.get("p_value", 1)
        # g:Profiler uses g:SCS or custom threshold; 'significant' indicates passes
        sig = item.get("significant", False)
        term_size = item.get("term_size", 0)
        query_size = item.get("query_size", 0)
        intersections = item.get("intersections", [[]])
        overlap_genes = intersections[0] if intersections else []
        overlap_str = ",".join(overlap_genes[:20]) if overlap_genes else ""
        if len(overlap_genes) > 20:
            overlap_str += f" (+{len(overlap_genes)-20} more)"
        rows.append({
            "Source": source,
            "Term_ID": term_id,
            "Term_Name": name,
            "P_value": pval,
            "Significant": sig,
            "Term_Size": term_size,
            "Query_Size": query_size,
            "Overlap_Genes": overlap_str,
            "Overlap_Count": len(overlap_genes),
        })
    df = pd.DataFrame(rows)
    # Add Benjamini-Hochberg FDR
    if len(df) > 0:
        from statsmodels.stats.multitest import multipletests
        df = df.sort_values('P_value')
        _, fdr, _, _ = multipletests(df['P_value'].values, method='fdr_bh')
        df['FDR_BH'] = fdr
    return df

def main():
    print("Pathway Enrichment for GSE122897 DE Genes")
    print("=" * 60)

    # 1. Ranked gene list (GSEA-style)
    print("\n[1] Ranked gene list analysis (ordered=TRUE)...")
    ranked = load_ranked_de_genes(use_deseq2=True)
    if not ranked:
        print("[ERROR] No ranked genes. Check DESeq2/Welch results.")
        return
    print(f"  Loaded {len(ranked)} genes (ranked by -log10(padj)*sign(log2FC))")
    symbol_fixes = {'Ccd': 'RUNX2'}
    ranked = [symbol_fixes.get(g, g) for g in ranked if str(g).strip()]

    res_ranked = gprofiler_enrichment(ranked, ordered=True, fdr_threshold=0.05)
    if res_ranked:
        df_ranked = results_to_df(res_ranked)
        out_ranked = "results/gse122897_pathway_enrichment_ranked.tsv"
        df_ranked.to_csv(out_ranked, sep='\t', index=False)
        n_sig = (df_ranked['FDR_BH'] < 0.05).sum() if 'FDR_BH' in df_ranked.columns else (df_ranked['Significant']).sum()
        print(f"  [OK] Saved {len(df_ranked)} pathways to {out_ranked} ({n_sig} at FDR<0.05)")
        for src in ["GO:BP", "KEGG", "REAC"]:
            sub = df_ranked[df_ranked["Source"] == src]
            if len(sub) > 0:
                top = sub.head(3)
                print(f"  Top {src}: {top['Term_Name'].iloc[0][:50]}...")
    else:
        print("  [WARN] No ranked results from g:Profiler")

    # 2. Threshold-based (genes with padj<0.05)
    print("\n[2] Threshold-based analysis (padj<0.05)...")
    thresh_genes = load_threshold_genes(padj_thresh=0.05)
    if len(thresh_genes) > 10:
        thresh_genes = [symbol_fixes.get(g, g) for g in thresh_genes if str(g).strip()]
        res_thresh = gprofiler_enrichment(thresh_genes, ordered=False, fdr_threshold=0.05)
        if res_thresh:
            df_thresh = results_to_df(res_thresh)
            out_thresh = "results/gse122897_pathway_enrichment.tsv"
            df_thresh.to_csv(out_thresh, sep='\t', index=False)
            print(f"  [OK] Saved {len(df_thresh)} pathways to {out_thresh}")
    else:
        print("  [INFO] Few genes at padj<0.05; using ranked results as primary")

    # 3. Create filtered table for manuscript (FDR<0.05 or P<0.05)
    primary = "results/gse122897_pathway_enrichment_ranked.tsv"
    if os.path.exists(primary):
        df = pd.read_csv(primary, sep='\t')
        if 'FDR_BH' in df.columns:
            filt = df[df['FDR_BH'] < 0.05]
        else:
            filt = df[df['P_value'] < 0.05]
        if len(filt) == 0:
            filt = df[df['P_value'] < 0.1].head(20)
        filt.to_csv("results/gse122897_pathway_enrichment_filtered.tsv", sep='\t', index=False)
        print(f"\n[OK] Filtered table: {len(filt)} pathways -> gse122897_pathway_enrichment_filtered.tsv")

if __name__ == "__main__":
    main()
