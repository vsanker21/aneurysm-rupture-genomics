import os
import pandas as pd
import numpy as np

INFLAMMATION_GENES = ['IL6', 'TNF', 'IL1B', 'CRP', 'CXCL8', 'IL8', 'CXCL1', 'CXCL2']
ECM_GENES = ['MMP9', 'MMP2', 'TIMP1', 'COL3A1', 'COL1A1', 'ELN', 'FBN1', 'LAMA3']
TOP_DE_GENE_IDS = [6611, 124220, 119587, 8685, 3909, 7045]
ENHANCED_EC_MARKERS = ['ANGPT2','ICAM1','VCAM1','THBS1','PERP','CXCL8','CDH5','ENG','KDR','TEK','VWF','PECAM1','LYVE1','SELE','SELP','NOS3','NOS2']

def load_expr():
    expr_file = "data/sc_rnaseq/GSE122897/GSE122897_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
    if not os.path.exists(expr_file):
        raise FileNotFoundError(expr_file)
    return pd.read_csv(expr_file, sep='\t', compression='gzip', index_col=0)

def map_symbols(expr_df, symbols):
    idx_str = expr_df.index.astype(str).tolist()
    hits = [g for g in symbols if g in expr_df.index]
    remaining = [s for s in symbols if s not in hits]
    if remaining:
        try:
            import requests, time
            url = "https://mygene.info/v3/query"
            params = {
                'q': ','.join(remaining),
                'scopes': 'symbol',
                'fields': 'entrezgene',
                'species': 'human'
            }
            r = requests.get(url, params=params, timeout=30)
            if r.status_code == 200:
                data = r.json()
                for item in data:
                    eg = item.get('entrezgene', None)
                    if eg is None:
                        continue
                    eg_str = str(eg)
                    if eg_str in idx_str:
                        hits.append(expr_df.index[idx_str.index(eg_str)])
            time.sleep(0.2)
        except Exception:
            pass
    # deduplicate and ensure present
    seen = set()
    out = []
    for h in hits:
        if h not in seen and str(h) in idx_str:
            out.append(h)
            seen.add(h)
    return out

def main():
    os.makedirs('results', exist_ok=True)
    expr_df = load_expr()
    idx_str_set = set(expr_df.index.astype(str))

    infl_hits = map_symbols(expr_df, INFLAMMATION_GENES)
    ecm_hits = map_symbols(expr_df, ECM_GENES)
    de_hits = [str(gid) for gid in TOP_DE_GENE_IDS if str(gid) in idx_str_set]
    ec_enh_hits = []
    for g in ENHANCED_EC_MARKERS:
        h = map_symbols(expr_df, [g])
        if h:
            ec_enh_hits.append(h[0])

    rows = [
        {'Feature':'Inflammation_Score','Requested':len(INFLAMMATION_GENES),'Mapped':len(infl_hits),'Mapped_Names':','.join(map(str, infl_hits))},
        {'Feature':'ECM_Remodeling_Score','Requested':len(ECM_GENES),'Mapped':len(ecm_hits),'Mapped_Names':','.join(map(str, ecm_hits))},
        {'Feature':'DE_Signature_Score','Requested':len(TOP_DE_GENE_IDS),'Mapped':len(de_hits),'Mapped_Names':','.join(de_hits)},
        {'Feature':'EC_Dysfunction_Score_Enhanced','Requested':len(ENHANCED_EC_MARKERS),'Mapped':len(ec_enh_hits),'Mapped_Names':','.join(map(str, ec_enh_hits))},
    ]

    out = pd.DataFrame(rows)
    out_file = 'results/integration_feature_debug_report.tsv'
    try:
        out.to_csv(out_file, sep='\t', index=False)
        print(f"[OK] Saved mapping diagnostics to {out_file}")
    except Exception as e:
        print(f"[WARNING] Could not save TSV: {e}")

    # Print to stdout
    print("\n=== Feature Mapping Diagnostics (stdout) ===")
    for _, r in out.iterrows():
        print(f"{r['Feature']}: requested={int(r['Requested'])}, mapped={int(r['Mapped'])}")
        print(f"  Mapped IDs: {r['Mapped_Names'][:300]}")

if __name__ == '__main__':
    main()
