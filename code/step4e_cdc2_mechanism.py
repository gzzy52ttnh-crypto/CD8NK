"""
step4e_cdc2_mechanism.py
cDC2_CD1C 抗原呈递差异深度挖掘

背景：Step 4d 发现 cDC2_CD1C 髓系亚群在 Taxane vs Pemetrexed 间存在差异（p=0.020），
这是唯一支持方案特异性的髓系差异线索。本脚本深入分析：
1. BH-FDR 校正后 cDC2 差异是否仍显著
2. cDC2 抗原呈递功能状态（MHC II、共刺激分子表达）
3. cDC2 比例与 NK-dominant ratio 的相关性
4. cDC2 比例与响应（pCR）的相关性
5. cDC2 与 SPP1+TAM 的关系

输出：cdc2_mechanism_results.csv, FigS8_cdc2_mechanism.png/pdf
"""
import config as cfg
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, spearmanr, chi2_contingency
from scipy.sparse import issparse
from statsmodels.stats.multitest import multipletests
from _common import score_gene_set
import os, warnings
warnings.filterwarnings('ignore')

print("[step4e_cdc2_mechanism] Starting...")
os.makedirs(cfg.RESULT_DIR, exist_ok=True)

results = {}

# ============================================================
# 1. 加载数据与化疗方案分组
# ============================================================
print("\n" + "="*60)
print("1. 加载数据与化疗方案分组")
print("="*60)

myeloid = sc.read_h5ad(cfg.H5AD_MYELOID)
print(f"Myeloid: {myeloid.shape}")

# 化疗方案分组
adata_t = sc.read_h5ad(cfg.H5AD_T)
chemo_pat = adata_t.obs.drop_duplicates('sampleID')[['sampleID', 'chemotherapy']].copy()
chemo_pat.columns = ['patient_id', 'chemotherapy']
del adata_t

from _common import classify_chemo

chemo_pat['chemo_class'] = chemo_pat['chemotherapy'].apply(classify_chemo)
tax_pats = set(chemo_pat[chemo_pat['chemo_class']=='Platinum+Taxane']['patient_id'])
peme_pats = set(chemo_pat[chemo_pat['chemo_class']=='Platinum+Pemetrexed']['patient_id'])
print(f"Taxane: {len(tax_pats)}, Pemetrexed: {len(peme_pats)}")

# 响应标签
pat_resp = chemo_pat.set_index('patient_id')['chemo_class'].copy()

# ============================================================
# 2. BH-FDR 校正验证 cDC2 差异显著性
# ============================================================
print("\n" + "="*60)
print("2. BH-FDR 校正验证")
print("="*60)

mye_obs = myeloid.obs
mye_subtypes = mye_obs[cfg.COL_CELL_TYPE].value_counts().index.tolist()
print(f"Myeloid subtypes: {len(mye_subtypes)}")

# 每患者各髓系亚群比例
pat_mye_comp = []
for p in set(tax_pats) | set(peme_pats):
    p_mask = mye_obs[cfg.COL_SAMPLE] == p
    n_total = p_mask.sum()
    if n_total < 20:
        continue
    grp = 'Taxane' if p in tax_pats else 'Pemetrexed'
    row = {'patient': p, 'group': grp, 'total_mye': n_total}
    for st in mye_subtypes[:10]:
        row[st] = (mye_obs.loc[p_mask, cfg.COL_CELL_TYPE] == st).sum() / n_total
    pat_mye_comp.append(row)

df_mye_comp = pd.DataFrame(pat_mye_comp)
tax_mye = df_mye_comp[df_mye_comp['group']=='Taxane']
peme_mye = df_mye_comp[df_mye_comp['group']=='Pemetrexed']

# 对 top 10 髓系亚群做 MWU + BH-FDR
raw_pvals = []
subtype_names = []
for st in mye_subtypes[:10]:
    if st not in df_mye_comp.columns:
        continue
    tax_vals = tax_mye[st].values
    peme_vals = peme_mye[st].values
    if len(tax_vals) < 3 or len(peme_vals) < 3:
        raw_pvals.append(1.0)
        subtype_names.append(st)
        continue
    stat, pval = mannwhitneyu(tax_vals, peme_vals, alternative='two-sided')
    raw_pvals.append(pval)
    subtype_names.append(st)
    print(f"  {st}: Taxane={np.mean(tax_vals):.4f}, Peme={np.mean(peme_vals):.4f}, raw p={pval:.4f}")

# BH-FDR 校正
_, fdr_pvals, _, _ = multipletests(raw_pvals, alpha=0.05, method='fdr_bh')
print("\nBH-FDR 校正结果:")
for st, raw_p, fdr_p in zip(subtype_names, raw_pvals, fdr_pvals):
    sig = '**' if fdr_p < 0.05 else '*' if raw_p < 0.05 else ''
    print(f"  {st}: raw p={raw_p:.4f}, FDR p={fdr_p:.4f}{sig}")
    results[f'Mye_{st}_raw_p'] = raw_p
    results[f'Mye_{st}_fdr_p'] = fdr_p

# ============================================================
# 3. cDC2_CD1C 抗原呈递功能状态分析
# ============================================================
print("\n" + "="*60)
print("3. cDC2_CD1C 抗原呈递功能状态")
print("="*60)

cdc2_mask = mye_obs[cfg.COL_CELL_TYPE] == 'cDC2_CD1C'
cdc2 = myeloid[cdc2_mask].copy()
print(f"cDC2_CD1C cells: {cdc2.n_obs}")

# 抗原呈递基因集
MHC_II_GENES = ['HLA-DRA', 'HLA-DRB1', 'HLA-DRB5', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DPA1', 'HLA-DPB1']
COSTIM_GENES = ['CD80', 'CD86', 'CD40', 'ICOSLG', 'PDCD1LG2', 'CD274']
ANTIGEN_PROCESSING = ['HLA-A', 'HLA-B', 'HLA-C', 'TAP1', 'TAP2', 'B2M']
DC_MARKERS = ['CD1C', 'CLEC10A', 'FCER1A', 'CD1A']

# 基因集评分（统一使用 _common.score_gene_set）
mhc2_valid = [g for g in MHC_II_GENES if g in cdc2.var_names]
costim_valid = [g for g in COSTIM_GENES if g in cdc2.var_names]
proc_valid = [g for g in ANTIGEN_PROCESSING if g in cdc2.var_names]
dc_valid = [g for g in DC_MARKERS if g in cdc2.var_names]
print(f"  MHC II: {len(mhc2_valid)}/{len(MHC_II_GENES)} genes found")
print(f"  Co-stimulatory: {len(costim_valid)}/{len(COSTIM_GENES)} genes found")
print(f"  Antigen processing: {len(proc_valid)}/{len(ANTIGEN_PROCESSING)} genes found")
print(f"  DC markers: {len(dc_valid)}/{len(DC_MARKERS)} genes found")

mhc2_score = score_gene_set(cdc2, MHC_II_GENES, 'MHC II')
costim_score = score_gene_set(cdc2, COSTIM_GENES, 'Co-stimulatory')
proc_score = score_gene_set(cdc2, ANTIGEN_PROCESSING, 'Antigen processing')
dc_score = score_gene_set(cdc2, DC_MARKERS, 'DC markers')

# 患者水平功能评分
pat_cdc2_func = []
for p in set(tax_pats) | set(peme_pats):
    p_mask = cdc2.obs[cfg.COL_SAMPLE] == p
    n = p_mask.sum()
    if n < 5:
        continue
    grp = 'Taxane' if p in tax_pats else 'Pemetrexed'
    pat_cdc2_func.append({
        'patient': p,
        'group': grp,
        'n_cdc2': n,
        'MHCII_score': np.mean(mhc2_score[p_mask]),
        'Costim_score': np.mean(costim_score[p_mask]),
        'Processing_score': np.mean(proc_score[p_mask]),
        'DC_marker_score': np.mean(dc_score[p_mask]),
    })

df_cdc2_func = pd.DataFrame(pat_cdc2_func)
tax_cdc2 = df_cdc2_func[df_cdc2_func['group']=='Taxane']
peme_cdc2 = df_cdc2_func[df_cdc2_func['group']=='Pemetrexed']

print(f"\nTaxane cDC2 patients: {len(tax_cdc2)}")
print(f"Pemetrexed cDC2 patients: {len(peme_cdc2)}")

func_pvals = []
func_metrics = []
for metric in ['MHCII_score', 'Costim_score', 'Processing_score', 'DC_marker_score']:
    tax_vals = tax_cdc2[metric].values
    peme_vals = peme_cdc2[metric].values
    stat, pval = mannwhitneyu(tax_vals, peme_vals, alternative='two-sided')
    func_pvals.append(pval)
    func_metrics.append(metric)
    results[f'cDC2_{metric}_tax_mean'] = np.mean(tax_vals)
    results[f'cDC2_{metric}_peme_mean'] = np.mean(peme_vals)
    results[f'cDC2_{metric}_pvalue'] = pval

# FDR校正（Benjamini-Hochberg）
_, func_qvals, _, _ = multipletests(func_pvals, method='fdr_bh')
for metric, pval, qval in zip(func_metrics, func_pvals, func_qvals):
    results[f'cDC2_{metric}_qvalue'] = qval
    sig_p = ' *' if pval < 0.05 else ''
    sig_q = ' **' if qval < 0.05 else ''
    tax_mean = results[f'cDC2_{metric}_tax_mean']
    peme_mean = results[f'cDC2_{metric}_peme_mean']
    print(f"  {metric}: Taxane={tax_mean:.4f}, Peme={peme_mean:.4f}, MWU p={pval:.4f}{sig_p}, q={qval:.4f}{sig_q}")

# ============================================================
# 4. cDC2 比例与 NK-dominant ratio 的相关性
# ============================================================
print("\n" + "="*60)
print("4. cDC2 比例与 NK-dominant ratio 的相关性")
print("="*60)

patient_metrics = pd.read_csv(cfg.result_path('per_patient_metrics.csv'))
cdc2_ratios = {}
for p in df_mye_comp['patient']:
    if 'cDC2_CD1C' in df_mye_comp.columns:
        cdc2_ratios[p] = df_mye_comp.loc[df_mye_comp['patient']==p, 'cDC2_CD1C'].iloc[0]
    else:
        cdc2_ratios[p] = 0.0

nk_ratios = {}
for _, row in patient_metrics.iterrows():
    nk_ratios[row['patient_id']] = row['nk_dominant_ratio']

common_pats = set(cdc2_ratios.keys()) & set(nk_ratios.keys())
cdc2_vals = [cdc2_ratios[p] for p in common_pats]
nk_vals = [nk_ratios[p] for p in common_pats]

sr, sp = spearmanr(cdc2_vals, nk_vals)
print(f"Spearman r={sr:.4f}, p={sp:.4f} (n={len(common_pats)})")
results['cDC2_vs_NK_spearman_r'] = sr
results['cDC2_vs_NK_spearman_p'] = sp

# 按化疗方案分层的相关性
tax_common = [p for p in common_pats if p in tax_pats]
peme_common = [p for p in common_pats if p in peme_pats]

tax_cdc2_vals = [cdc2_ratios[p] for p in tax_common]
tax_nk = [nk_ratios[p] for p in tax_common]
sr_tax, sp_tax = spearmanr(tax_cdc2_vals, tax_nk)
print(f"Taxane only: Spearman r={sr_tax:.4f}, p={sp_tax:.4f} (n={len(tax_common)})")
results['cDC2_vs_NK_spearman_r_taxane'] = sr_tax
results['cDC2_vs_NK_spearman_p_taxane'] = sp_tax

peme_cdc2_vals = [cdc2_ratios[p] for p in peme_common]
peme_nk = [nk_ratios[p] for p in peme_common]
sr_peme, sp_peme = spearmanr(peme_cdc2_vals, peme_nk)
print(f"Pemetrexed only: Spearman r={sr_peme:.4f}, p={sp_peme:.4f} (n={len(peme_common)})")
results['cDC2_vs_NK_spearman_r_peme'] = sr_peme
results['cDC2_vs_NK_spearman_p_peme'] = sp_peme

# ============================================================
# 5. cDC2 比例与响应（pCR）的相关性
# ============================================================
print("\n" + "="*60)
print("5. cDC2 比例与响应（pCR）的相关性")
print("="*60)

# 加载响应标签
resp_map = {}
adata_t = sc.read_h5ad(cfg.H5AD_T)
for p in adata_t.obs[cfg.COL_SAMPLE].unique():
    resp = adata_t.obs.loc[adata_t.obs[cfg.COL_SAMPLE]==p, cfg.COL_RESPONSE].iloc[0]
    resp_map[p] = resp
del adata_t

cdc2_pcr = []
cdc2_nmp = []
for p in common_pats:
    if resp_map.get(p, '') == 'pCR':
        cdc2_pcr.append(cdc2_ratios[p])
    elif resp_map.get(p, '') == 'non-MPR':
        cdc2_nmp.append(cdc2_ratios[p])

stat_mwu, p_mwu = mannwhitneyu(cdc2_pcr, cdc2_nmp, alternative='two-sided')
print(f"pCR vs non-MPR: MWU p={p_mwu:.4f}")
print(f"  pCR mean={np.mean(cdc2_pcr):.4f} (n={len(cdc2_pcr)})")
print(f"  non-MPR mean={np.mean(cdc2_nmp):.4f} (n={len(cdc2_nmp)})")
results['cDC2_pCR_mean'] = np.mean(cdc2_pcr)
results['cDC2_nonMPR_mean'] = np.mean(cdc2_nmp)
results['cDC2_pCR_vs_NMP_p'] = p_mwu

# ============================================================
# 6. cDC2 与 SPP1+TAM 的关系
# ============================================================
print("\n" + "="*60)
print("6. cDC2 与 SPP1+TAM 的关系")
print("="*60)

mye_pat = pd.read_csv(cfg.result_path('myeloid_per_patient.csv'))
spp1_map = {}
for _, row in mye_pat.iterrows():
    spp1_map[row['patient_id']] = float(row['SPP1_TAM_ratio']) if pd.notna(row['SPP1_TAM_ratio']) else np.nan

common_spp1 = [p for p in common_pats if p in spp1_map and not np.isnan(spp1_map[p])]
cdc2_spp1_vals = [cdc2_ratios[p] for p in common_spp1]
spp1_vals = [spp1_map[p] for p in common_spp1]

sr_spp1, sp_spp1 = spearmanr(cdc2_spp1_vals, spp1_vals)
print(f"cDC2 ratio vs SPP1+TAM ratio: Spearman r={sr_spp1:.4f}, p={sp_spp1:.4f} (n={len(common_spp1)})")
results['cDC2_vs_SPP1_spearman_r'] = sr_spp1
results['cDC2_vs_SPP1_spearman_p'] = sp_spp1

# ============================================================
# 7. 保存结果 + 绘图
# ============================================================
print("\n" + "="*60)
print("7. 保存结果 + 绘图")
print("="*60)

res_df = pd.DataFrame({'metric': list(results.keys()), 'value': list(results.values())})
res_df.to_csv(cfg.result_path('cdc2_mechanism_results.csv'), index=False)
print(f"cdc2_mechanism_results.csv saved ({len(res_df)} metrics)")

# 绘图：4 个 panel
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel A: cDC2 比例 Taxane vs Pemetrexed
axA = axes[0, 0]
tax_cdc2_ratio = [cdc2_ratios[p] for p in tax_common]
peme_cdc2_ratio = [cdc2_ratios[p] for p in peme_common]
bp_a = axA.boxplot([tax_cdc2_ratio, peme_cdc2_ratio], patch_artist=True, labels=['Taxane', 'Pemetrexed'])
bp_a['boxes'][0].set_facecolor('#3498DB'); bp_a['boxes'][1].set_facecolor('#E74C3C')
stat_a, p_a = mannwhitneyu(tax_cdc2_ratio, peme_cdc2_ratio, alternative='two-sided')
fdr_a = fdr_pvals[subtype_names.index('cDC2_CD1C')] if 'cDC2_CD1C' in subtype_names else 1.0
axA.set_title(f'Panel A: cDC2_CD1C Proportion\nMWU p={p_a:.4f}, FDR p={fdr_a:.4f}')
axA.set_ylabel('cDC2_CD1C proportion')

# Panel B: cDC2 MHC II 评分
axB = axes[0, 1]
bp_b = axB.boxplot([tax_cdc2['MHCII_score'].values, peme_cdc2['MHCII_score'].values],
                   patch_artist=True, labels=['Taxane', 'Pemetrexed'])
bp_b['boxes'][0].set_facecolor('#3498DB'); bp_b['boxes'][1].set_facecolor('#E74C3C')
stat_b, p_b = mannwhitneyu(tax_cdc2['MHCII_score'], peme_cdc2['MHCII_score'])
axB.set_title(f'Panel B: cDC2 MHC II Expression\nMWU p={p_b:.4f}')
axB.set_ylabel('MHC II score')

# Panel C: cDC2 vs NK-dominant ratio scatter
axC = axes[1, 0]
axC.scatter(cdc2_vals, nk_vals, alpha=0.6, c='steelblue', s=20)
axC.set_xlabel('cDC2_CD1C proportion')
axC.set_ylabel('NK-dominant ratio')
axC.set_title(f'Panel C: cDC2 vs NK-dominant\nSpearman r={sr:.4f}, p={sp:.4f}')

# Panel D: cDC2 vs response
axD = axes[1, 1]
bp_d = axD.boxplot([cdc2_pcr, cdc2_nmp], patch_artist=True, labels=['pCR', 'non-MPR'])
bp_d['boxes'][0].set_facecolor('#4CAF50'); bp_d['boxes'][1].set_facecolor('#F44336')
stat_d, p_d = mannwhitneyu(cdc2_pcr, cdc2_nmp)
axD.set_title(f'Panel D: cDC2 Proportion by Response\nMWU p={p_d:.4f}')
axD.set_ylabel('cDC2_CD1C proportion')

plt.tight_layout()
fig.savefig(cfg.result_path('FigS8_cdc2_mechanism.pdf'), dpi=150, bbox_inches='tight')
fig.savefig(cfg.result_path('FigS8_cdc2_mechanism.png'), dpi=150, bbox_inches='tight')
plt.close(fig)
print("FigS8_cdc2_mechanism saved")

print("\n" + "="*60)
print("8. 小结")
print("="*60)
print(f"cDC2_CD1C Taxane vs Pemetrexed: raw p={p_a:.4f}, FDR p={fdr_a:.4f}")
print(f"cDC2 MHC II expression: p={p_b:.4f}")
print(f"cDC2 vs NK-dominant: r={sr:.4f}, p={sp:.4f}")
print(f"cDC2 vs pCR: p={p_d:.4f}")
print("\n[step4e_cdc2_mechanism] Done.")
