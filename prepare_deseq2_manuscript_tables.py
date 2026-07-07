"""Prepare DESeq2 and comparison tables for manuscript integration."""
import pandas as pd
import os
import requests
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

def fetch_symbols_mygene(gene_ids):
    """Fetch gene symbols from mygene.info API for missing IDs."""
    m = {}
    for gid in gene_ids:
        gid_str = str(gid).replace('.0', '')
        try:
            r = requests.get(f'https://mygene.info/v3/gene/{gid_str}?fields=symbol', timeout=10)
            if r.status_code == 200:
                d = r.json()
                m[gid_str] = d.get('symbol', gid_str)
            else:
                m[gid_str] = gid_str
        except Exception:
            m[gid_str] = gid_str
        time.sleep(0.2)
    return m

def main():
    e2s = load_entrez_to_symbol()
    os.makedirs("results", exist_ok=True)

    # Top DESeq2 genes annotated
    ds = pd.read_csv("results/gse122897_deseq2_results.tsv", sep='\t')
    ds_top = ds.head(20).copy()
    ds_top['Gene_Symbol'] = ds_top['Gene_ID'].astype(str).str.replace('.0', '').map(e2s)
    # Fill missing symbols via mygene API
    missing = ds_top[ds_top['Gene_Symbol'].isna() | (ds_top['Gene_Symbol'] == '')]['Gene_ID'].astype(str).str.replace('.0', '').tolist()
    if missing:
        api_map = fetch_symbols_mygene(missing)
        for idx, row in ds_top.iterrows():
            gid = str(row['Gene_ID']).replace('.0', '')
            if pd.isna(row['Gene_Symbol']) or row['Gene_Symbol'] == '':
                ds_top.at[idx, 'Gene_Symbol'] = api_map.get(gid, gid)
    ds_top = ds_top[['Gene_ID', 'Gene_Symbol', 'log2FoldChange', 'padj', 'baseMean']]
    ds_top.columns = ['Gene_ID', 'Gene_Symbol', 'Log2FC', 'P_adj', 'baseMean']
    ds_top.to_csv("results/gse122897_deseq2_top_annotated.tsv", sep='\t', index=False)
    print(f"Saved results/gse122897_deseq2_top_annotated.tsv ({len(ds_top)} genes)")

    # Welch vs DESeq2 comparison
    wg = pd.read_csv("results/gse122897_top_genes_annotated.tsv", sep='\t')
    w_ids = set(wg['Gene'].astype(str).str.replace('.0', ''))
    d_sig = ds[ds['padj'] < 0.05].dropna(subset=['padj'])
    d_ids = set(d_sig['Gene_ID'].astype(str).str.replace('.0', ''))
    overlap = w_ids & d_ids
    w_only = w_ids - d_ids
    d_only_count = len(d_ids - w_ids)

    comp = pd.DataFrame([
        {'Metric': 'Welch_t-test_genes_FDR_0.05', 'Value': len(wg[wg['P_adj_BH'] < 0.05]) if 'P_adj_BH' in wg.columns else 'N/A'},
        {'Metric': 'DESeq2_genes_padj_0.05', 'Value': len(d_sig)},
        {'Metric': 'Overlap_top_Welch_in_DESeq2_sig', 'Value': len(overlap)},
        {'Metric': 'Welch_top_10_genes', 'Value': ', '.join(wg.head(10)['Gene_Symbol'].fillna(wg.head(10)['Gene'].astype(str)).tolist())},
        {'Metric': 'DESeq2_top_10_Gene_IDs', 'Value': ', '.join(ds.head(10)['Gene_ID'].astype(str).tolist())},
    ])
    comp.to_csv("results/gse122897_welch_vs_deseq2_comparison.tsv", sep='\t', index=False)
    print(f"Saved results/gse122897_welch_vs_deseq2_comparison.tsv")

    # Pathway enrichment filtered: use ranked results if available (FDR < 0.05)
    if os.path.exists("results/gse122897_pathway_enrichment_ranked.tsv"):
        pw = pd.read_csv("results/gse122897_pathway_enrichment_ranked.tsv", sep='\t')
        if 'FDR_BH' in pw.columns:
            pw_filt = pw[pw['FDR_BH'] < 0.05].head(25)
        else:
            pw_filt = pw[pw['P_value'] < 0.05].head(25)
    else:
        pw = pd.read_csv("results/gse122897_pathway_enrichment.tsv", sep='\t')
        pw_filt = pw[pw['P_value'] < 0.1].head(15)
    pw_filt.to_csv("results/gse122897_pathway_enrichment_filtered.tsv", sep='\t', index=False)
    print(f"Saved results/gse122897_pathway_enrichment_filtered.tsv ({len(pw_filt)} pathways)")

if __name__ == "__main__":
    main()
