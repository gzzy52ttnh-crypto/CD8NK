#!/usr/bin/env python3
"""重新计算早期 T04-T09 中 5 项无法验证的指标

用途：早期探索性代码（T04-T09）已清理，但台账中遗留 5 项无法验证的统计指标。
本脚本从原始数据（GSE243013_T_cells.h5ad + GSE120575）重新计算这些指标并落盘为 CSV。

输出：
  result/T04_T08_recomputed_stats.csv  （SEM/LUSC OR/CD8A/KIR/化疗对照）
  result/T08_GSE120575_recomputed.csv  （GSE120575 NK-like MWU）

运行：/opt/anaconda3/envs/scanpy2/bin/python code/recompute_5_metrics.py
"""

import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats
from scipy.stats import mannwhitneyu
import statsmodels.formula.api as smf
import scipy.sparse as sp

DATA_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ADATA_DIR = os.path.join(DATA_DIR, 'adata')
RESULT_DIR = os.path.join(DATA_DIR, 'result')
OUTPUT_CSV = os.path.join(RESULT_DIR, 'T04_T08_recomputed_stats.csv')
OUTPUT_T08 = os.path.join(RESULT_DIR, 'T08_GSE120575_recomputed.csv')

rows = []

print("=" * 60)
print("重新计算 5 项无法验证的指标")
print("=" * 60)

# ============================================================
# 加载 T_cells h5ad
# ============================================================
print("\n[Loading GSE243013_T_cells.h5ad ...]")
adata = ad.read_h5ad(os.path.join(ADATA_DIR, 'GSE243013_T_cells.h5ad'))

# normalize
if sp.issparse(adata.X):
    counts = np.array(adata.X.sum(axis=1)).flatten()
    scale = 1e4 / np.maximum(counts, 1)
    adata.X = sp.diags(scale) @ adata.X
    adata.X.data = np.log1p(adata.X.data)
else:
    counts = adata.X.sum(axis=1)
    adata.X = adata.X / counts[:, None] * 1e4
    adata.X = np.log1p(adata.X)

# 定义 CD8 / NK-like / Tex
sct = adata.obs['sub_cell_type'].astype(str)
is_cd8 = sct.str.startswith('CD8T')
is_nklike = (sct == 'CD8T_NK-like_FGFBP2')
is_tex = sct.str.contains('Tex')

# 响应分组
resp = adata.obs['pathological_response'].astype(str)

# 大克隆定义 >=5
clone_num = adata.obs['clonotype_number'].astype(int)
is_large = (clone_num >= 5)

# 患者级数据
obs = adata.obs[['sampleID', 'cancer_type', 'chemotherapy', 'anti-PD1_therapy',
                 'pathological_response', 'clonotype', 'clonotype_number',
                 'sub_cell_type']].copy()
obs['is_nklike'] = is_nklike.values
obs['is_large'] = is_large.values

# ============================================================
# 指标 1: T04 SEM（NK-like 大克隆比例的标准误）
# ============================================================
print("\n--- 指标 1: T04 SEM (NK-like 大克隆比例标准误) ---")

patient_stats = []
for sid, grp in obs.groupby('sampleID'):
    large_mask = grp['is_large']
    if large_mask.sum() == 0:
        continue
    nk_large = (large_mask & grp['is_nklike']).sum()
    total_large = large_mask.sum()
    nk_like_large_pct = nk_large / total_large
    pr = grp['pathological_response'].iloc[0]
    patient_stats.append({
        'sampleID': sid,
        'nk_like_large_pct': nk_like_large_pct,
        'total_large_clones': total_large,
        'pathological_response': pr
    })

df_pat = pd.DataFrame(patient_stats)

for grp_name in ['pCR', 'MPR', 'non-MPR']:
    vals = df_pat.loc[df_pat['pathological_response'] == grp_name, 'nk_like_large_pct']
    if len(vals) > 1:
        sem_val = vals.std(ddof=1) / np.sqrt(len(vals))
        mean_val = vals.mean()
        print(f"  {grp_name}: n={len(vals)}, mean={mean_val:.4f}, SEM={sem_val:.6f} ({sem_val*100:.4f}%)")
        rows.append({
            'module': 'T04', 'metric': f'SEM_{grp_name}',
            'value': sem_val, 'value_pct': sem_val * 100,
            'n': len(vals), 'note': f'NK-like大克隆比例标准误 (mean={mean_val:.4f})'
        })

# ============================================================
# 指标 2: T04 cancer_type(LUSC) OR
# ============================================================
print("\n--- 指标 2: T04 cancer_type(LUSC) OR ---")

df_binary = df_pat[df_pat['pathological_response'].isin(['pCR', 'non-MPR'])].copy()
cancer_map = obs.drop_duplicates('sampleID').set_index('sampleID')['cancer_type'].to_dict()
df_binary['cancer_type'] = df_binary['sampleID'].map(cancer_map)
df_binary['resp_bin'] = (df_binary['pathological_response'] == 'pCR').astype(int)
df_binary['is_LUSC'] = (df_binary['cancer_type'] == 'LUSC').astype(int)

print(f"  pCR vs non-MPR: n={len(df_binary)} (pCR={df_binary['resp_bin'].sum()}, non-MPR={(1-df_binary['resp_bin']).sum()})")
print(f"  LUSC={df_binary['is_LUSC'].sum()}, LUAD={(1-df_binary['is_LUSC']).sum()}")

try:
    model = smf.logit('resp_bin ~ is_LUSC', data=df_binary).fit(disp=0)
    or_val = np.exp(model.params['is_LUSC'])
    ci = np.exp(model.conf_int().loc['is_LUSC'])
    p_val = model.pvalues['is_LUSC']
    print(f"  LUSC OR = {or_val:.4f} [95% CI {ci[0]:.4f} - {ci[1]:.4f}], p = {p_val:.4f}")
    rows.append({
        'module': 'T04', 'metric': 'cancer_type_LUSC_OR',
        'value': or_val, 'value_pct': None,
        'n': len(df_binary), 'note': f'LUSC vs LUAD, pCR vs non-MPR, CI=[{ci[0]:.4f},{ci[1]:.4f}], p={p_val:.4f}'
    })
except Exception as e:
    print(f"  Error: {e}")
    from scipy.stats import fisher_exact
    tab = pd.crosstab(df_binary['is_LUSC'], df_binary['resp_bin'])
    print(f"  2x2 table:\n{tab}")
    or_val, p_val = fisher_exact(tab)
    print(f"  Fisher exact: OR = {or_val:.4f}, p = {p_val:.4f}")
    rows.append({
        'module': 'T04', 'metric': 'cancer_type_LUSC_OR',
        'value': or_val, 'value_pct': None,
        'n': len(df_binary), 'note': f'Fisher exact, p={p_val:.4f}'
    })

# 含 MPR 的版本
df_all = df_pat[df_pat['pathological_response'].isin(['pCR', 'MPR', 'non-MPR'])].copy()
df_all['cancer_type'] = df_all['sampleID'].map(cancer_map)
df_all['resp_bin'] = (df_all['pathological_response'] == 'pCR').astype(int)
df_all['is_LUSC'] = (df_all['cancer_type'] == 'LUSC').astype(int)

try:
    model2 = smf.logit('resp_bin ~ is_LUSC', data=df_all).fit(disp=0)
    or_val2 = np.exp(model2.params['is_LUSC'])
    ci2 = np.exp(model2.conf_int().loc['is_LUSC'])
    p_val2 = model2.pvalues['is_LUSC']
    print(f"  [含MPR] LUSC OR = {or_val2:.4f} [95% CI {ci2[0]:.4f} - {ci2[1]:.4f}], p = {p_val2:.4f} (n={len(df_all)})")
    rows.append({
        'module': 'T04', 'metric': 'cancer_type_LUSC_OR_withMPR',
        'value': or_val2, 'value_pct': None,
        'n': len(df_all), 'note': f'含MPR, pCR vs MPR+non-MPR, CI=[{ci2[0]:.4f},{ci2[1]:.4f}], p={p_val2:.4f}'
    })
except Exception as e:
    print(f"  [含MPR] Error: {e}")

# ============================================================
# 指标 3: T04 CD8A/CD8B/KIR 表达值
# ============================================================
print("\n--- 指标 3: T04 CD8A/CD8B/KIR 表达值 ---")

genes_expr = ['CD8A', 'CD8B', 'KIR2DL1', 'KIR2DL3', 'KIR3DL1', 'KIR3DL2', 'KIR2DS4']
available_genes = [g for g in genes_expr if g in adata.var_names]
missing_genes = [g for g in genes_expr if g not in adata.var_names]
print(f"  可用基因: {available_genes}")
if missing_genes:
    print(f"  缺失基因: {missing_genes}")

cd8_mask = is_cd8.values
nklike_mask = is_nklike.values
other_cd8_mask = cd8_mask & ~nklike_mask

print(f"  NK-like cells: {nklike_mask.sum()}")
print(f"  Other CD8T cells: {other_cd8_mask.sum()}")

for gene in available_genes:
    gene_idx = list(adata.var_names).index(gene)
    if sp.issparse(adata.X):
        expr_nk = np.array(adata.X[nklike_mask, gene_idx].todense()).flatten()
        expr_other = np.array(adata.X[other_cd8_mask, gene_idx].todense()).flatten()
    else:
        expr_nk = adata.X[nklike_mask, gene_idx]
        expr_other = adata.X[other_cd8_mask, gene_idx]

    mean_nk = expr_nk.mean()
    mean_other = expr_other.mean()
    try:
        mwu_p = mannwhitneyu(expr_nk, expr_other, alternative='two-sided').pvalue
    except:
        mwu_p = np.nan

    print(f"  {gene}: NK-like mean={mean_nk:.4f}, Other CD8T mean={mean_other:.4f}, MWU p={mwu_p:.4e}")
    rows.append({
        'module': 'T04', 'metric': f'{gene}_expression',
        'value': mean_nk, 'value_pct': mean_other,
        'n': int(nklike_mask.sum()), 'note': f'NK-like mean={mean_nk:.4f}, Other CD8T mean={mean_other:.4f}, MWU p={mwu_p:.4e}'
    })

# ============================================================
# 指标 4: T08 NK-like MWU p (GSE120575 细胞级 + 患者级)
# ============================================================
print("\n--- 指标 4: T08 NK-like MWU p (GSE120575) ---")

nk_genes_path = os.path.join(ADATA_DIR, 'GSE120575', 'GSE120575_nk_genes.csv')
pred_path = os.path.join(ADATA_DIR, 'GSE120575', 'GSE120575_Patient_Predictions.csv')

t08_rows = []

if os.path.exists(nk_genes_path) and os.path.exists(pred_path):
    df_nk = pd.read_csv(nk_genes_path, index_col=0)
    print(f"  GSE120575 NK genes matrix: {df_nk.shape[0]} genes x {df_nk.shape[1]} cells")

    cell_names = df_nk.columns.tolist()
    patient_ids = [c.split('_')[1] for c in cell_names]

    df_pred = pd.read_csv(pred_path)
    print(f"  Predictions file: {len(df_pred)} patients")

    pid_col = 'patient_id' if 'patient_id' in df_pred.columns else df_pred.columns[0]
    resp_col = 'NR_R' if 'NR_R' in df_pred.columns else df_pred.columns[1]

    pred_map = {}
    for _, row in df_pred.iterrows():
        pid = str(row[pid_col])
        resp = str(row[resp_col])
        pred_map[pid] = resp

    cell_response = []
    cell_patient = []
    for c in cell_names:
        pid = c.split('_')[1]
        cell_patient.append(pid)
        cell_response.append(pred_map.get(pid, 'Unknown'))

    nk_score = df_nk.mean(axis=0)

    df_cell = pd.DataFrame({
        'patient': cell_patient,
        'response': cell_response,
        'nk_score': nk_score.values
    })

    df_cell_r = df_cell[df_cell['response'] == 'R']
    df_cell_nr = df_cell[df_cell['response'] == 'NR']

    print(f"  R cells: {len(df_cell_r)}, NR cells: {len(df_cell_nr)}")
    n_r_pat = df_cell_r['patient'].nunique()
    n_nr_pat = df_cell_nr['patient'].nunique()
    print(f"  R patients: {n_r_pat}, NR patients: {n_nr_pat}")

    if len(df_cell_r) > 0 and len(df_cell_nr) > 0:
        mwu_p_cell = mannwhitneyu(df_cell_r['nk_score'], df_cell_nr['nk_score'], alternative='two-sided').pvalue
        mean_r = df_cell_r['nk_score'].mean()
        mean_nr = df_cell_nr['nk_score'].mean()
        print(f"  [细胞级] MWU p = {mwu_p_cell:.6f}, R mean = {mean_r:.4f}, NR mean = {mean_nr:.4f}")
        rows.append({
            'module': 'T08', 'metric': 'NK_like_MWU_cell_level',
            'value': mwu_p_cell, 'value_pct': None,
            'n': len(df_cell_r) + len(df_cell_nr),
            'note': f'细胞级, R mean={mean_r:.4f}, NR mean={mean_nr:.4f}'
        })
        t08_rows.append({
            'module': 'T08', 'metric': 'NK_like_MWU_cell_level',
            'value': mwu_p_cell, 'n': len(df_cell_r) + len(df_cell_nr),
            'note': f'细胞级, R mean={mean_r:.4f}, NR mean={mean_nr:.4f}'
        })

        pb = df_cell.groupby(['patient', 'response'], observed=True)['nk_score'].mean().reset_index()
        pb_r = pb[pb['response'] == 'R']
        pb_nr = pb[pb['response'] == 'NR']
        if len(pb_r) > 0 and len(pb_nr) > 0:
            mwu_p_pat = mannwhitneyu(pb_r['nk_score'], pb_nr['nk_score'], alternative='two-sided').pvalue
            mean_pb_r = pb_r['nk_score'].mean()
            mean_pb_nr = pb_nr['nk_score'].mean()
            print(f"  [患者级pseudobulk] MWU p = {mwu_p_pat:.6f}, R mean = {mean_pb_r:.4f}, NR mean = {mean_pb_nr:.4f}")
            rows.append({
                'module': 'T08', 'metric': 'NK_like_MWU_patient_level',
                'value': mwu_p_pat, 'value_pct': None,
                'n': len(pb_r) + len(pb_nr),
                'note': f'患者级pseudobulk, R mean={mean_pb_r:.4f}, NR mean={mean_pb_nr:.4f}'
            })
            t08_rows.append({
                'module': 'T08', 'metric': 'NK_like_MWU_patient_level',
                'value': mwu_p_pat, 'n': len(pb_r) + len(pb_nr),
                'note': f'患者级pseudobulk, R mean={mean_pb_r:.4f}, NR mean={mean_pb_nr:.4f}'
            })

    df_t08 = pd.DataFrame(t08_rows)
    df_t08.to_csv(OUTPUT_T08, index=False)
    print(f"  T08 结果已保存: {OUTPUT_T08}")
else:
    print(f"  文件不存在: {nk_genes_path} or {pred_path}")

# ============================================================
# 指标 5: T04 化疗对照 n
# ============================================================
print("\n--- 指标 5: T04 化疗对照 n ---")

df_patient = obs.drop_duplicates('sampleID')[['sampleID', 'chemotherapy', 'anti-PD1_therapy', 'pathological_response', 'cancer_type']].copy()

no_chemo = df_patient[df_patient['chemotherapy'] == 'No']
print(f"  chemotherapy='No': n={len(no_chemo)}")
rows.append({
    'module': 'T04', 'metric': 'chemo_No_count',
    'value': len(no_chemo), 'value_pct': None,
    'n': len(df_patient), 'note': f'chemotherapy=No 的患者数'
})

no_antiPD1 = df_patient[df_patient['anti-PD1_therapy'] == 'No']
print(f"  anti-PD1_therapy='No': n={len(no_antiPD1)}")
rows.append({
    'module': 'T04', 'metric': 'antiPD1_No_count',
    'value': len(no_antiPD1), 'value_pct': None,
    'n': len(df_patient), 'note': f'anti-PD1_therapy=No 的患者数'
})

both_no = df_patient[(df_patient['chemotherapy'] == 'No') & (df_patient['anti-PD1_therapy'] == 'No')]
print(f"  chemo=No AND anti-PD1=No: n={len(both_no)}")
rows.append({
    'module': 'T04', 'metric': 'both_No_count',
    'value': len(both_no), 'value_pct': None,
    'n': len(df_patient), 'note': f'chemo=No AND anti-PD1=No 的患者数'
})

# ============================================================
# 保存
# ============================================================
df_out = pd.DataFrame(rows)
df_out.to_csv(OUTPUT_CSV, index=False)
print(f"\n{'=' * 60}")
print(f"结果已保存至: {OUTPUT_CSV}")
print(f"共 {len(df_out)} 项指标")
print(f"{'=' * 60}")
print(df_out[['module', 'metric', 'value', 'note']].to_string(index=False))
