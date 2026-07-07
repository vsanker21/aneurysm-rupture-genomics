import os, csv, gzip, tarfile, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
import requests

base = Path(r"e:\Aneurysms Genomics Project\Moya Moya")
out = Path(r"e:\Aneurysms Genomics Project\results\moyamoya")
out.mkdir(parents=True, exist_ok=True)

series_path = base / "GSE157628_series_matrix.txt.gz"
raw_tar = base / "GSE157628_RAW.tar"

meta_rows = {}
expr_header = None
expr_data = []
in_table = False

with gzip.open(series_path, 'rt', encoding='utf-8', errors='ignore') as f:
    for line in f:
        line = line.rstrip('\n')
        if line.startswith('!series_matrix_table_begin'):
            in_table = True
            continue
        if line.startswith('!series_matrix_table_end'):
            in_table = False
            continue
        if in_table:
            if expr_header is None:
                expr_header = next(csv.reader([line], delimiter='\t', quotechar='"'))
            else:
                row = next(csv.reader([line], delimiter='\t', quotechar='"'))
                expr_data.append(row)
            continue
        if line.startswith('!Sample_'):
            parts = next(csv.reader([line], delimiter='\t', quotechar='"'))
            key = parts[0].replace('!Sample_', '')
            vals = parts[1:]
            meta_rows.setdefault(key, []).append(vals)

samples = expr_header[1:]
expr_df = pd.DataFrame(expr_data, columns=expr_header)
expr_df = expr_df.rename(columns={expr_df.columns[0]: 'ID_REF'})
for c in samples:
    expr_df[c] = pd.to_numeric(expr_df[c], errors='coerce')
expr_df['ID_REF'] = pd.to_numeric(expr_df['ID_REF'], errors='coerce').astype('Int64')
expr_df = expr_df.dropna(subset=['ID_REF']).copy()

sample_md = pd.DataFrame({'Sample': samples})

def get_row(key, idx=0):
    rows = meta_rows.get(key, [])
    if len(rows) > idx:
        return rows[idx]
    return [''] * len(samples)

sample_md['Title'] = get_row('title', 0)
sample_md['Source'] = get_row('source_name_ch1', 0)
sample_md['Tissue'] = get_row('characteristics_ch1', 0)
sample_md['Sex_raw'] = get_row('characteristics_ch1', 1)
sample_md['Age_raw'] = get_row('characteristics_ch1', 2)

def infer_group(title):
    t = str(title)
    if '[MMD' in t:
        return 'MMD'
    if '[IA' in t:
        return 'IA'
    if '[EPI' in t:
        return 'EPI'
    return 'Unknown'

sample_md['Group'] = sample_md['Title'].map(infer_group)
sample_md['DiseaseBinary'] = (sample_md['Group'] == 'MMD').astype(int)
sample_md['Sex'] = sample_md['Sex_raw'].str.replace('gender:', '', regex=False).str.strip().replace({'Fenale':'Female'})
sample_md['Age'] = sample_md['Age_raw'].str.replace('age:', '', regex=False).str.replace('y', '', regex=False).str.strip()
sample_md['Age'] = pd.to_numeric(sample_md['Age'], errors='coerce')
sample_md.to_csv(out / 'gse157628_sample_metadata.tsv', sep='\t', index=False)

with tarfile.open(raw_tar, 'r') as t:
    m = t.extractfile('GSM4771802_MMD1.txt.gz')
    txt = gzip.decompress(m.read()).decode('utf-8', 'ignore').splitlines()

feat_header = None
feat_rows = []
for ln in txt:
    if ln.startswith('FEATURES\t'):
        feat_header = ln.split('\t')
        continue
    if feat_header and ln.startswith('DATA\t'):
        vals = ln.split('\t')
        if len(vals) == len(feat_header):
            feat_rows.append(vals)

fdf = pd.DataFrame(feat_rows, columns=feat_header)
fdf['ProbeUID'] = pd.to_numeric(fdf['ProbeUID'], errors='coerce').astype('Int64')
probe_map = fdf[['ProbeUID', 'ProbeName', 'GeneName', 'SystematicName', 'Description']].dropna(subset=['ProbeUID'])
probe_map = probe_map.rename(columns={'ProbeUID':'ID_REF'})
probe_map.to_csv(out / 'gse157628_probe_annotation.tsv', sep='\t', index=False)

mmd_samples = sample_md.loc[sample_md['Group'] == 'MMD', 'Sample'].tolist()
ctrl_samples = sample_md.loc[sample_md['Group'].isin(['IA','EPI']), 'Sample'].tolist()
X = expr_df.set_index('ID_REF')
Xa = X[mmd_samples]
Xb = X[ctrl_samples]
mean_a = Xa.mean(axis=1)
mean_b = Xb.mean(axis=1)
log2fc = mean_a - mean_b
_, pvals = ttest_ind(Xa.values, Xb.values, axis=1, equal_var=False, nan_policy='omit')
res = pd.DataFrame({'ID_REF': X.index, 'Mean_MMD': mean_a.values, 'Mean_Control': mean_b.values, 'Log2FC_MMD_vs_Control': log2fc.values, 'P_value': pvals})
res = res.replace([np.inf, -np.inf], np.nan).dropna(subset=['P_value'])
res['FDR_BH'] = multipletests(res['P_value'].values, method='fdr_bh')[1]
res['P_Bonferroni'] = multipletests(res['P_value'].values, method='bonferroni')[1]
res = res.sort_values('P_value').merge(probe_map, on='ID_REF', how='left')
res.to_csv(out / 'gse157628_mmd_vs_control_de.tsv', sep='\t', index=False)

def valid_symbol(s):
    s = str(s)
    if s in ('nan', 'None', ''):
        return False
    bad_prefix = ('A_', 'ENST', 'LOC', 'lnc-', 'XLOC', 'TCONS')
    return not s.startswith(bad_prefix)

rg = res.copy()
rg['GeneSymbol'] = rg['GeneName'].astype(str).str.split(';').str[0]
rg = rg[rg['GeneSymbol'].map(valid_symbol)]
rg = rg.drop_duplicates(subset=['GeneSymbol'], keep='first')
rg.head(100).to_csv(out / 'gse157628_mmd_vs_control_top_genes.tsv', sep='\t', index=False)

top_probes = res['ID_REF'].head(1000).tolist()
Xclf = X.loc[top_probes, sample_md['Sample']].T
y = sample_md.set_index('Sample').loc[Xclf.index, 'DiseaseBinary'].values
pipe = Pipeline([('scaler', StandardScaler()), ('kbest', SelectKBest(score_func=f_classif, k=min(100, Xclf.shape[1]))), ('clf', LogisticRegression(max_iter=2000, solver='liblinear'))])
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
auc_scores = cross_val_score(pipe, Xclf, y, cv=cv, scoring='roc_auc')
pd.DataFrame({'Metric':['CV_AUC_Mean','CV_AUC_Std','N_Samples','N_MMD','N_Control'],'Value':[auc_scores.mean(),auc_scores.std(ddof=1),len(sample_md),len(mmd_samples),len(ctrl_samples)]}).to_csv(out / 'gse157628_mmd_classifier_cv.tsv', sep='\t', index=False)

enrich_rows = []
pr = rg.copy()
pr['rank_score'] = -np.log10(pr['P_value'].clip(lower=1e-300)) * np.sign(pr['Log2FC_MMD_vs_Control'])
pr = pr.sort_values('rank_score', ascending=False)
ranked_symbols = pr['GeneSymbol'].astype(str).drop_duplicates().head(500).tolist()
if len(ranked_symbols) >= 20:
    payload = {'organism':'hsapiens','query':ranked_symbols,'sources':['GO:BP','GO:MF','GO:CC','KEGG','REAC'],'ordered':True,'all_results':True,'user_threshold':0.05,'no_iea':False,'domain_scope':'annotated'}
    r = requests.post('https://biit.cs.ut.ee/gprofiler/api/gost/profile', json=payload, timeout=120)
    if r.status_code == 200:
        for it in r.json().get('result', []):
            inter = it.get('intersections', [[]])
            ol = inter[0] if inter else []
            enrich_rows.append({'Source':it.get('source',''),'Term_ID':it.get('native',it.get('id','')),'Term_Name':it.get('name',''),'P_value':it.get('p_value',1.0),'Significant':it.get('significant',False),'Term_Size':it.get('term_size',0),'Query_Size':it.get('query_size',0),'Overlap_Count':len(ol),'Overlap_Genes':','.join(ol[:30])})

enr = pd.DataFrame(enrich_rows)
if len(enr):
    enr = enr.sort_values('P_value')
    enr['FDR_BH'] = multipletests(enr['P_value'].values, method='fdr_bh')[1]
    enr.to_csv(out / 'gse157628_mmd_pathway_enrichment_ranked.tsv', sep='\t', index=False)
    enr[enr['FDR_BH'] < 0.05].head(100).to_csv(out / 'gse157628_mmd_pathway_enrichment_fdr05.tsv', sep='\t', index=False)

summary = {
    'dataset':'GSE157628',
    'platform':'GPL16699 (Agilent SurePrint G3 Human GE v2 8x60K)',
    'assay_type':'Expression profiling by array (bulk artery wall tissue; NOT scRNA-seq)',
    'n_total':int(len(sample_md)),
    'group_counts':sample_md['Group'].value_counts().to_dict(),
    'age_mean_by_group':sample_md.groupby('Group')['Age'].mean().round(2).dropna().to_dict(),
    'sex_counts':sample_md['Sex'].value_counts().to_dict(),
    'de_fdr_005':int((res['FDR_BH'] < 0.05).sum()),
    'de_p_001_absfc2':int(((res['P_value'] < 1e-3) & (np.abs(res['Log2FC_MMD_vs_Control']) > 1)).sum()),
    'cv_auc_mean':float(auc_scores.mean()),
    'cv_auc_std':float(auc_scores.std(ddof=1)),
}
if len(enr):
    summary['pathways_fdr_005'] = int((enr['FDR_BH'] < 0.05).sum())
    summary['top_pathways'] = enr.head(10)[['Source','Term_Name','FDR_BH']].to_dict(orient='records')

with open(out / 'gse157628_analysis_summary.json', 'w', encoding='utf-8') as fh:
    json.dump(summary, fh, indent=2)

print('Saved outputs to', out)
print(json.dumps(summary, indent=2))
