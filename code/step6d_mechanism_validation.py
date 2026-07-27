"""
step6d_mechanism_validation.py
Step 6.9-6.12: 机制验证与临床转化补充分析
- Step 6.9: Taxane vs Peme NK-like功能状态比较
- Step 6.10: Taxane vs Peme SPP1+TAM比例比较（髓系闸门）
- Step 6.11: IRS评分分层验证（临床模型方案特异性）
- Step 6.12: Table 1队列基线汇总表
"""
import config as cfg
from _common import NKLIKE_SIGNATURE
import pandas as pd
import numpy as np
import scipy.sparse as sp
from scipy.stats import chi2, mannwhitneyu
from sklearn.metrics import roc_auc_score
from statsmodels.stats.multitest import multipletests
import anndata as ad
import os
import gzip
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

print("[step6d_mechanism_validation] Starting...")
os.makedirs(cfg.RESULT_DIR, exist_ok=True)

p_values_for_fdr = []
p_value_labels = []

# 功能基因集
EFFECTOR_GENES = ['GZMB', 'GZMH', 'GZMA', 'PRF1', 'GNLY', 'IFNG', 'NKG7', 'CTSW']
CHEMOKINE_GENES = ['CX3CR1', 'XCL1', 'XCL2', 'CCL5', 'CXCR3']  # 注意XCL1/2可能不在数据中
ACTIVATION_GENES = ['KLRD1', 'KLRC1', 'KLRC2', 'KLRB1', 'KLRF1', 'FCGR3A', 'CD160', 'CRTAM', 'SH2D1B', 'TYROBP', 'FCER1G']
TRANSCRIPTION_GENES = ['TBX21', 'EOMES', 'ZNF683']

# 髓系相关基因
SPP1_GENE = 'SPP1'
TAM_MARKERS = ['CD68', 'CD163', 'MSR1', 'MRC1']  # TAM标记

# ============================================================
# 读取数据
# ============================================================
print("\n" + "="*70)
print("Step 0: 读取 GSE241934 数据")
print("="*70)

h5ad_path = os.path.join(cfg.RESULT_DIR, 'GSE241934_all.h5ad')
adata = ad.read_h5ad(h5ad_path)
print(f"  数据维度: {adata.shape}")

# 仅保留RWC队列（Taxane vs Pemetrexed对比，IIT只有Taxane）
adata_rwc = adata[adata.obs['cohort'] == 'RWC'].copy()
print(f"  RWC队列: {adata_rwc.shape[0]} cells")

# 过滤掉Unknown和Other方案
adata_rwc = adata_rwc[adata_rwc.obs['chemo_class'].isin(['Taxane', 'Pemetrexed'])].copy()
print(f"  RWC Taxane+Peme: {adata_rwc.shape[0]} cells")
print(f"  Taxane患者: {adata_rwc[adata_rwc.obs['chemo_class']=='Taxane'].obs['sampleID'].nunique()}")
print(f"  Pemetrexed患者: {adata_rwc[adata_rwc.obs['chemo_class']=='Pemetrexed'].obs['sampleID'].nunique()}")

# 患者水平元数据
patient_meta = adata.obs.groupby('sampleID').agg({
    'cohort': 'first',
    'chemo_class': 'first',
    'response_binary': 'first',
    'Pathological Response': 'first',
    'Gender': 'first',
    'Age': 'first',
    'Histology': 'first',
    'EGFR': 'first',
    'Smoking_History': 'first',
    'PD-L1 TPS': 'first',
    'Cycles': 'first',
    'PD1': 'first'
}).reset_index()
print(f"  总患者数: {len(patient_meta)}")

# ============================================================
# Step 6.9: Taxane vs Peme NK-like功能状态比较
# ============================================================
print("\n" + "="*70)
print("Step 6.9: Taxane vs Peme NK-like功能状态比较")
print("="*70)

# 提取NK-like细胞（CD8T_Teff_FGFBP2）
adata_nklike = adata_rwc[adata_rwc.obs['cell.type'] == 'CD8T_Teff_FGFBP2'].copy()
print(f"  NK-like细胞数: {adata_nklike.shape[0]}")

# 计算各功能基因集的平均表达（log1p归一化）
def calc_signature_score(adata, gene_list):
    """计算签名得分：基因平均表达（log1p后）"""
    avail_genes = [g for g in gene_list if g in adata.var_names]
    if len(avail_genes) == 0:
        return np.zeros(adata.shape[0])
    X = adata[:, avail_genes].X
    if sp.issparse(X):
        X = X.toarray()
    return np.mean(X, axis=1)

# 各功能模块得分
print("\n  计算功能模块得分...")
effector_scores = calc_signature_score(adata_nklike, EFFECTOR_GENES)
chemokine_scores = calc_signature_score(adata_nklike, CHEMOKINE_GENES)
activation_scores = calc_signature_score(adata_nklike, ACTIVATION_GENES)
transcription_scores = calc_signature_score(adata_nklike, TRANSCRIPTION_GENES)

# 添加到obs
adata_nklike.obs['effector_score'] = effector_scores
adata_nklike.obs['chemokine_score'] = chemokine_scores
adata_nklike.obs['activation_score'] = activation_scores
adata_nklike.obs['transcription_score'] = transcription_scores

# 患者水平平均得分
patient_func = adata_nklike.obs.groupby('sampleID').agg({
    'chemo_class': 'first',
    'response_binary': 'first',
    'effector_score': 'mean',
    'chemokine_score': 'mean',
    'activation_score': 'mean',
    'transcription_score': 'mean'
}).reset_index()

# Taxane vs Peme 比较（不考虑响应，只看方案差异）
taxane_patients = patient_func[patient_func['chemo_class'] == 'Taxane']
peme_patients = patient_func[patient_func['chemo_class'] == 'Pemetrexed']

print(f"\n  Taxane: {len(taxane_patients)} 例, Pemetrexed: {len(peme_patients)} 例")
print(f"\n  功能状态比较（患者水平均值）:")
print(f"  {'模块':<20} {'Taxane均值':>12} {'Peme均值':>12} {'倍数变化':>10} {'MWU p':>10}")
print("  " + "-"*70)

func_results = []
for mod in ['effector_score', 'chemokine_score', 'activation_score', 'transcription_score']:
    tax_vals = taxane_patients[mod].values
    peme_vals = peme_patients[mod].values
    stat, p = stats.mannwhitneyu(tax_vals, peme_vals, alternative='two-sided')
    fc = np.mean(tax_vals) / np.mean(peme_vals) if np.mean(peme_vals) > 0 else np.nan
    p_values_for_fdr.append(p)
    p_value_labels.append(f'Func_{mod}')
    print(f"  {mod:<20} {np.mean(tax_vals):>12.4f} {np.mean(peme_vals):>12.4f} {fc:>10.2f} {p:>10.4f}")
    func_results.append({
        'module': mod,
        'taxane_mean': np.mean(tax_vals),
        'peme_mean': np.mean(peme_vals),
        'fold_change': fc,
        'mwu_p': p,
        'taxane_n': len(tax_vals),
        'peme_n': len(peme_vals)
    })

# 保存结果
df_func = pd.DataFrame(func_results)
df_func.to_csv(cfg.result_path('GSE241934_functional_comparison.csv'), index=False)
print(f"\n  ✅ 功能比较结果已保存: GSE241934_functional_comparison.csv")

# 可视化：NK-like功能状态雷达图
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Panel A: 箱线图（四模块对比）
ax = axes[0]
modules = ['effector_score', 'chemokine_score', 'activation_score', 'transcription_score']
module_labels = ['Effector', 'Chemokine', 'Activation', 'Transcription']
x_pos = np.arange(len(modules))
width = 0.35

tax_means = [np.mean(taxane_patients[m]) for m in modules]
tax_sems = [stats.sem(taxane_patients[m]) for m in modules]
peme_means = [np.mean(peme_patients[m]) for m in modules]
peme_sems = [stats.sem(peme_patients[m]) for m in modules]

ax.bar(x_pos - width/2, tax_means, width, yerr=tax_sems, 
       label='Taxane', color='#E64B35', alpha=0.8, capsize=3)
ax.bar(x_pos + width/2, peme_means, width, yerr=peme_sems,
       label='Pemetrexed', color='#4DBBD5', alpha=0.8, capsize=3)

ax.set_xticks(x_pos)
ax.set_xticklabels(module_labels, rotation=15, ha='right')
ax.set_ylabel('Mean expression (log1p)')
ax.set_title('NK-like cell functional status\n(Taxane vs Pemetrexed)')
ax.legend()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel B: 各基因表达对比热图样式
ax = axes[1]
all_nk_genes = EFFECTOR_GENES + CHEMOKINE_GENES + ACTIVATION_GENES + TRANSCRIPTION_GENES
all_nk_genes = [g for g in all_nk_genes if g in adata_nklike.var_names]

tax_gene_means = []
peme_gene_means = []
for g in all_nk_genes:
    g_expr = adata_nklike[:, g].X
    if sp.issparse(g_expr):
        g_expr = g_expr.toarray().flatten()
    tax_mask = adata_nklike.obs['chemo_class'] == 'Taxane'
    peme_mask = adata_nklike.obs['chemo_class'] == 'Pemetrexed'
    tax_gene_means.append(np.mean(g_expr[tax_mask]))
    peme_gene_means.append(np.mean(g_expr[peme_mask]))

# 归一化到0-1便于可视化
all_vals = tax_gene_means + peme_gene_means
vmin, vmax = min(all_vals), max(all_vals)
tax_norm = [(v - vmin) / (vmax - vmin) for v in tax_gene_means]
peme_norm = [(v - vmin) / (vmax - vmin) for v in peme_gene_means]

y_pos = np.arange(len(all_nk_genes))
ax.barh(y_pos + 0.2, tax_norm, 0.4, label='Taxane', color='#E64B35', alpha=0.8)
ax.barh(y_pos - 0.2, peme_norm, 0.4, label='Pemetrexed', color='#4DBBD5', alpha=0.8)
ax.set_yticks(y_pos)
ax.set_yticklabels(all_nk_genes, fontsize=8)
ax.set_xlabel('Normalized expression')
ax.set_title('NK-like signature gene expression')
ax.legend()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.invert_yaxis()

plt.tight_layout()
fig_path = cfg.result_path('GSE241934_Fig7_functional_comparison.png')
plt.savefig(fig_path, dpi=150, bbox_inches='tight')
plt.savefig(fig_path.replace('.png', '.pdf'), bbox_inches='tight')
plt.close()
print(f"  ✅ 功能比较图已保存: GSE241934_Fig7_functional_comparison.png")

# ============================================================
# Step 6.10: Taxane vs Peme SPP1+TAM比例比较（髓系闸门）
# ============================================================
print("\n" + "="*70)
print("Step 6.10: Taxane vs Peme SPP1+TAM比例比较（髓系闸门）")
print("="*70)

# 定义SPP1+ TAM：表达SPP1的巨噬细胞
# 方法1：所有巨噬细胞中SPP1+的比例
# 方法2：SPP1+ TAM占所有免疫细胞的比例

# 提取髓系细胞
mye_cell_types = [ct for ct in adata_rwc.obs['cell.type'].unique() 
                  if 'Macro' in str(ct) or 'Mono' in str(ct)]
adata_mye = adata_rwc[adata_rwc.obs['cell.type'].isin(mye_cell_types)].copy()
print(f"  髓系细胞数: {adata_mye.shape[0]}")
print(f"  髓系细胞类型: {sorted(mye_cell_types)}")

# 计算每个细胞的SPP1表达
spp1_expr = adata_mye[:, SPP1_GENE].X
if sp.issparse(spp1_expr):
    spp1_expr = spp1_expr.toarray().flatten()
adata_mye.obs['SPP1_expr'] = spp1_expr
adata_mye.obs['is_SPP1_pos'] = (spp1_expr > 0).astype(int)

# 患者水平统计
patient_mye = []
for sample_id, group in adata_mye.obs.groupby('sampleID'):
    total_mye = len(group)
    spp1_pos_count = group['is_SPP1_pos'].sum()
    spp1_ratio = spp1_pos_count / total_mye if total_mye > 0 else np.nan
    
    # 各巨噬细胞亚群比例
    macro_total = group[group['cell.type'].str.contains('Macro', na=False)].shape[0]
    spp1_macro = group[(group['cell.type'].str.contains('Macro', na=False)) & (group['is_SPP1_pos'] == 1)].shape[0]
    
    patient_mye.append({
        'sampleID': sample_id,
        'chemo_class': group['chemo_class'].iloc[0],
        'response_binary': group['response_binary'].iloc[0],
        'total_myeloid': total_mye,
        'spp1_pos_count': spp1_pos_count,
        'spp1_ratio_in_mye': spp1_ratio,
        'total_macrophages': macro_total,
        'spp1_macro_count': spp1_macro,
        'spp1_macro_ratio_in_macro': spp1_macro / macro_total if macro_total > 0 else np.nan,
        'spp1_macro_ratio_in_mye': spp1_macro / total_mye if total_mye > 0 else np.nan
    })

df_mye = pd.DataFrame(patient_mye)
print(f"\n  患者数: {len(df_mye)}")

# Taxane vs Peme 比较
tax_mye = df_mye[df_mye['chemo_class'] == 'Taxane']
peme_mye = df_mye[df_mye['chemo_class'] == 'Pemetrexed']

print(f"\n  SPP1+ TAM比例比较:")
print(f"  {'指标':<30} {'Taxane均值':>12} {'Peme均值':>12} {'倍数变化':>10} {'MWU p':>10}")
print("  " + "-"*70)

mye_metrics = [
    ('spp1_ratio_in_mye', 'SPP1+ in Myeloid'),
    ('spp1_macro_ratio_in_macro', 'SPP1+ TAM in Macro'),
    ('spp1_macro_ratio_in_mye', 'SPP1+ TAM in Myeloid'),
]

mye_results = []
for col, label in mye_metrics:
    tax_vals = tax_mye[col].dropna().values
    peme_vals = peme_mye[col].dropna().values
    if len(tax_vals) > 0 and len(peme_vals) > 0:
        stat, p = stats.mannwhitneyu(tax_vals, peme_vals, alternative='two-sided')
        fc = np.mean(tax_vals) / np.mean(peme_vals) if np.mean(peme_vals) > 0 else np.nan
        p_values_for_fdr.append(p)
        p_value_labels.append(f'Mye_{label}')
        print(f"  {label:<30} {np.mean(tax_vals):>12.4f} {np.mean(peme_vals):>12.4f} {fc:>10.2f} {p:>10.4f}")
        mye_results.append({
            'metric': label,
            'taxane_mean': np.mean(tax_vals),
            'peme_mean': np.mean(peme_vals),
            'fold_change': fc,
            'mwu_p': p,
            'taxane_n': len(tax_vals),
            'peme_n': len(peme_vals)
        })

# 保存结果
df_mye_results = pd.DataFrame(mye_results)
df_mye_results.to_csv(cfg.result_path('GSE241934_SPP1_TAM_comparison.csv'), index=False)
df_mye.to_csv(cfg.result_path('GSE241934_patient_myeloid_metrics.csv'), index=False)
print(f"\n  ✅ 髓系比较结果已保存: GSE241934_SPP1_TAM_comparison.csv")

# 可视化：SPP1+ TAM比例对比
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Panel A: SPP1+ TAM在巨噬细胞中的比例
ax = axes[0]
box_data = [tax_mye['spp1_macro_ratio_in_macro'].dropna().values, 
            peme_mye['spp1_macro_ratio_in_macro'].dropna().values]
bp = ax.boxplot(box_data, labels=['Taxane', 'Pemetrexed'], patch_artist=True,
                medianprops={'color': 'black'})
colors = ['#E64B35', '#4DBBD5']
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
# 叠加散点
for i, data in enumerate(box_data):
    x = np.random.normal(i + 1, 0.05, size=len(data))
    ax.scatter(x, data, alpha=0.6, color='black', s=30, zorder=3)
ax.set_ylabel('SPP1+ TAM ratio in macrophages')
ax.set_title('SPP1+ TAM proportion\n(in macrophages)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel B: SPP1+ TAM在髓系细胞中的比例
ax = axes[1]
box_data = [tax_mye['spp1_macro_ratio_in_mye'].dropna().values,
            peme_mye['spp1_macro_ratio_in_mye'].dropna().values]
bp = ax.boxplot(box_data, labels=['Taxane', 'Pemetrexed'], patch_artist=True,
                medianprops={'color': 'black'})
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
for i, data in enumerate(box_data):
    x = np.random.normal(i + 1, 0.05, size=len(data))
    ax.scatter(x, data, alpha=0.6, color='black', s=30, zorder=3)
ax.set_ylabel('SPP1+ TAM ratio in myeloid')
ax.set_title('SPP1+ TAM proportion\n(in myeloid cells)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel C: 各巨噬细胞亚群SPP1+比例
ax = axes[2]
macro_types = sorted([ct for ct in mye_cell_types if 'Macro' in str(ct)])
tax_spp1_ratios = []
peme_spp1_ratios = []
for mt in macro_types:
    tax_mask = (adata_mye.obs['cell.type'] == mt) & (adata_mye.obs['chemo_class'] == 'Taxane')
    peme_mask = (adata_mye.obs['cell.type'] == mt) & (adata_mye.obs['chemo_class'] == 'Pemetrexed')
    tax_r = adata_mye.obs.loc[tax_mask, 'is_SPP1_pos'].mean() if tax_mask.sum() > 0 else 0
    peme_r = adata_mye.obs.loc[peme_mask, 'is_SPP1_pos'].mean() if peme_mask.sum() > 0 else 0
    tax_spp1_ratios.append(tax_r)
    peme_spp1_ratios.append(peme_r)

y_pos = np.arange(len(macro_types))
ax.barh(y_pos + 0.2, tax_spp1_ratios, 0.4, label='Taxane', color='#E64B35', alpha=0.8)
ax.barh(y_pos - 0.2, peme_spp1_ratios, 0.4, label='Pemetrexed', color='#4DBBD5', alpha=0.8)
ax.set_yticks(y_pos)
ax.set_yticklabels([mt.replace('Macro_', '') for mt in macro_types])
ax.set_xlabel('SPP1+ proportion')
ax.set_title('SPP1+ proportion by\nmacrophage subset')
ax.legend()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.invert_yaxis()

plt.tight_layout()
fig_path = cfg.result_path('GSE241934_Fig8_SPP1_TAM.png')
plt.savefig(fig_path, dpi=150, bbox_inches='tight')
plt.savefig(fig_path.replace('.png', '.pdf'), bbox_inches='tight')
plt.close()
print(f"  ✅ SPP1+ TAM对比图已保存: GSE241934_Fig8_SPP1_TAM.png")

# ============================================================
# Step 6.11: IRS评分分层验证（临床模型方案特异性）
# ============================================================
print("\n" + "="*70)
print("Step 6.11: IRS评分分层验证（临床模型方案特异性）")
print("="*70)

# 计算IRS：NK-like比例 * (1 - SPP1+TAM比例归一化)
# 先合并NK和髓系数据
patient_nk = adata_rwc.obs.groupby('sampleID').agg({
    'chemo_class': 'first',
    'response_binary': 'first',
    'cell.type': lambda x: (x == 'CD8T_Teff_FGFBP2').sum() / len(x)
}).reset_index()
patient_nk = patient_nk.rename(columns={'cell.type': 'nklike_ratio_all'})

# 合并
df_irs = patient_nk.merge(df_mye[['sampleID', 'spp1_macro_ratio_in_mye', 'total_myeloid']], 
                          on='sampleID', how='inner')
df_irs = df_irs.merge(patient_meta[['sampleID', 'Gender', 'Age', 'Histology', 'EGFR', 'Smoking_History']], 
                      on='sampleID', how='left')

# 归一化SPP1+ TAM比例（Min-Max到0-1）
spp1_min = df_irs['spp1_macro_ratio_in_mye'].min()
spp1_max = df_irs['spp1_macro_ratio_in_mye'].max()
if spp1_max > spp1_min:
    df_irs['spp1_norm'] = (df_irs['spp1_macro_ratio_in_mye'] - spp1_min) / (spp1_max - spp1_min)
else:
    print(f"  ⚠️ SPP1+ TAM比例无变异（min={spp1_min:.4f}, max={spp1_max:.4f}），IRS将退化为纯NK指标")
    df_irs['spp1_norm'] = 0.0

# 计算IRS = NK比例 * (1 - SPP1_norm)
df_irs['IRS'] = df_irs['nklike_ratio_all'] * (1.0 - df_irs['spp1_norm'])

print(f"  有IRS数据的患者数: {len(df_irs)}")
print(f"  IRS范围: [{df_irs['IRS'].min():.4f}, {df_irs['IRS'].max():.4f}]")

# 按方案分组验证IRS的预测效力
print(f"\n  IRS预测能力比较（R vs NR）:")
print(f"  {'方案':<15} {'n':>5} {'R':>3} {'NR':>3} {'R均值':>10} {'NR均值':>10} {'MWU p':>10} {'AUC':>8}")
print("  " + "-"*70)

irs_results = []
for chemo in ['Taxane', 'Pemetrexed']:
    df_sub = df_irs[df_irs['chemo_class'] == chemo].copy()
    r_vals = df_sub[df_sub['response_binary'] == 'R']['IRS'].values
    nr_vals = df_sub[df_sub['response_binary'] == 'NR']['IRS'].values
    
    if len(r_vals) > 0 and len(nr_vals) > 0:
        stat, p = stats.mannwhitneyu(r_vals, nr_vals, alternative='two-sided')
        p_values_for_fdr.append(p)
        p_value_labels.append(f'IRS_{chemo}')
        # AUC
        all_vals = np.concatenate([r_vals, nr_vals])
        all_labels = np.concatenate([np.ones(len(r_vals)), np.zeros(len(nr_vals))])
        from sklearn.metrics import roc_auc_score
        try:
            auc = roc_auc_score(all_labels, all_vals)
        except:
            auc = np.nan
        
        print(f"  {chemo:<15} {len(df_sub):>5} {len(r_vals):>3} {len(nr_vals):>3} "
              f"{np.mean(r_vals):>10.4f} {np.mean(nr_vals):>10.4f} {p:>10.4f} {auc:>8.3f}")
        irs_results.append({
            'chemo_class': chemo,
            'n': len(df_sub),
            'n_R': len(r_vals),
            'n_NR': len(nr_vals),
            'R_mean': np.mean(r_vals),
            'NR_mean': np.mean(nr_vals),
            'mwu_p': p,
            'auc': auc
        })

# 交互效应检验（Logistic回归 + statsmodels 精确p值）
print(f"\n  Logistic回归交互效应（IRS × chemo_class）:")
import statsmodels.api as sm
df_irs_model = df_irs.dropna(subset=['IRS', 'response_binary']).copy()
df_irs_model['y'] = (df_irs_model['response_binary'] == 'R').astype(int)
df_irs_model['chemo_taxane'] = (df_irs_model['chemo_class'] == 'Taxane').astype(int)
df_irs_model['IRS_z'] = (df_irs_model['IRS'] - df_irs_model['IRS'].mean()) / df_irs_model['IRS'].std()
df_irs_model['interaction'] = df_irs_model['IRS_z'] * df_irs_model['chemo_taxane']

y = df_irs_model['y'].values

if len(np.unique(y)) == 2 and len(y) > 5:
    try:
        # 完整模型: y ~ IRS + chemo + IRS:chemo
        X_full = sm.add_constant(df_irs_model[['IRS_z', 'chemo_taxane', 'interaction']])
        model_full = sm.Logit(y, X_full).fit(disp=0)
        or_interact = np.exp(model_full.params['interaction'])
        ci_interact = np.exp(model_full.conf_int().loc['interaction'])
        p_interact = model_full.pvalues['interaction']

        # 似然比检验
        X_noint = sm.add_constant(df_irs_model[['IRS_z', 'chemo_taxane']])
        model_noint = sm.Logit(y, X_noint).fit(disp=0)
        lr_stat = -2 * (model_noint.llf - model_full.llf)
        lr_p = 1 - chi2.cdf(lr_stat, df=1)

        print(f"    交互项 OR = {or_interact:.4f}, 95%CI = [{ci_interact[0]:.4f}, {ci_interact[1]:.4f}], Wald p = {p_interact:.6f}")
        print(f"    似然比检验: χ² = {lr_stat:.4f}, p = {lr_p:.6f}")
        print(f"    IRS_z OR = {np.exp(model_full.params['IRS_z']):.4f}")
        print(f"    chemo_taxane OR = {np.exp(model_full.params['chemo_taxane']):.4f}")

        irs_results.append({
            'chemo_class': 'Interaction',
            'n': len(y),
            'interaction_or': or_interact,
            'interaction_ci_lower': ci_interact[0],
            'interaction_ci_upper': ci_interact[1],
            'interaction_p_wald': p_interact,
            'interaction_lr_stat': lr_stat,
            'interaction_lr_p': lr_p,
            'irs_main_or': np.exp(model_full.params['IRS_z']),
            'chemo_main_or': np.exp(model_full.params['chemo_taxane']),
        })
    except Exception as e:
        print(f"    回归失败: {e}")
        irs_results.append({
            'chemo_class': 'Interaction',
            'n': len(y),
            'error': str(e),
        })

# 保存结果
df_irs_results = pd.DataFrame(irs_results)
df_irs_results.to_csv(cfg.result_path('GSE241934_IRS_validation.csv'), index=False)
df_irs.to_csv(cfg.result_path('GSE241934_patient_IRS.csv'), index=False)
print(f"\n  ✅ IRS验证结果已保存: GSE241934_IRS_validation.csv")

# 可视化：IRS分层验证
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: IRS箱线图（按方案和响应）
ax = axes[0]
tax_r = df_irs[(df_irs['chemo_class'] == 'Taxane') & (df_irs['response_binary'] == 'R')]['IRS'].values
tax_nr = df_irs[(df_irs['chemo_class'] == 'Taxane') & (df_irs['response_binary'] == 'NR')]['IRS'].values
peme_r = df_irs[(df_irs['chemo_class'] == 'Pemetrexed') & (df_irs['response_binary'] == 'R')]['IRS'].values
peme_nr = df_irs[(df_irs['chemo_class'] == 'Pemetrexed') & (df_irs['response_binary'] == 'NR')]['IRS'].values

box_data = [tax_r, tax_nr, peme_r, peme_nr]
labels = ['Taxane\nR', 'Taxane\nNR', 'Peme\nR', 'Peme\nNR']
colors_box = ['#E64B35', '#E64B3580', '#4DBBD5', '#4DBBD580']

bp = ax.boxplot(box_data, labels=labels, patch_artist=True,
                medianprops={'color': 'black'})
for patch, color in zip(bp['boxes'], colors_box):
    patch.set_facecolor(color)
    patch.set_alpha(0.8)
for i, data in enumerate(box_data):
    x = np.random.normal(i + 1, 0.05, size=len(data))
    ax.scatter(x, data, alpha=0.6, color='black', s=30, zorder=3)
ax.set_ylabel('IRS')
ax.set_title('IRS by treatment and response')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel B: ROC曲线对比
ax = axes[1]
for chemo, color in [('Taxane', '#E64B35'), ('Pemetrexed', '#4DBBD5')]:
    df_sub = df_irs[df_irs['chemo_class'] == chemo].dropna(subset=['IRS', 'response_binary'])
    if len(df_sub) > 5 and df_sub['response_binary'].nunique() == 2:
        from sklearn.metrics import roc_curve, roc_auc_score
        y_true = (df_sub['response_binary'] == 'R').astype(int)
        y_score = df_sub['IRS'].values
        fpr, tpr, _ = roc_curve(y_true, y_score)
        auc = roc_auc_score(y_true, y_score)
        ax.plot(fpr, tpr, color=color, label=f'{chemo} (AUC={auc:.3f})', linewidth=2)

ax.plot([0, 1], [0, 1], 'k--', alpha=0.3)
ax.set_xlabel('False Positive Rate')
ax.set_ylabel('True Positive Rate')
ax.set_title('ROC: IRS predicting response')
ax.legend()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
fig_path = cfg.result_path('GSE241934_Fig9_IRS_validation.png')
plt.savefig(fig_path, dpi=150, bbox_inches='tight')
plt.savefig(fig_path.replace('.png', '.pdf'), bbox_inches='tight')
plt.close()
print(f"  ✅ IRS验证图已保存: GSE241934_Fig9_IRS_validation.png")

# ============================================================
# Step 6.12: Table 1队列基线汇总表
# ============================================================
print("\n" + "="*70)
print("Step 6.12: Table 1队列基线汇总表")
print("="*70)

# 构建Table 1
def calc_table1(df, group_col='chemo_class', total_col='sampleID'):
    """生成基线特征表"""
    table_rows = []
    
    # 总样本量
    n_total = df[total_col].nunique()
    table_rows.append({'Characteristic': 'n', 'Total': n_total})
    
    # 分类变量
    cat_vars = [
        ('cohort', 'Cohort'),
        ('response_binary', 'Response (R/NR)'),
        ('Pathological Response', 'Pathological Response'),
        ('Gender', 'Gender'),
        ('Histology', 'Histology'),
        ('EGFR', 'EGFR Mutation'),
        ('Smoking_History', 'Smoking History'),
        ('PD1', 'PD-1 Inhibitor'),
    ]
    
    for col, label in cat_vars:
        if col in df.columns:
            counts = df.groupby(total_col)[col].first().value_counts()
            for cat, count in counts.items():
                pct = count / n_total * 100
                table_rows.append({
                    'Characteristic': f'{label} - {cat}',
                    'Total': f'{count} ({pct:.1f}%)'
                })
    
    # 连续变量
    cont_vars = [
        ('Age', 'Age (years)'),
        ('Cycles', 'Treatment Cycles'),
    ]
    
    for col, label in cont_vars:
        if col in df.columns:
            vals = pd.to_numeric(df.groupby(total_col)[col].first(), errors='coerce').dropna()
            if len(vals) > 0:
                table_rows.append({
                    'Characteristic': f'{label} (mean ± SD)',
                    'Total': f'{vals.mean():.1f} ± {vals.std():.1f}'
                })
                table_rows.append({
                    'Characteristic': f'{label} (median [IQR])',
                    'Total': f'{vals.median():.1f} [{vals.quantile(0.25):.1f}-{vals.quantile(0.75):.1f}]'
                })
    
    return pd.DataFrame(table_rows)

# 分方案生成Table 1（RWC队列）
df_rwc_patients = adata_rwc.obs.groupby('sampleID').first().reset_index()
print(f"  RWC队列患者数: {len(df_rwc_patients)}")

# 总表
table1_total = calc_table1(df_rwc_patients)

# Taxane组
table1_taxane = calc_table1(df_rwc_patients[df_rwc_patients['chemo_class'] == 'Taxane'])
table1_taxane = table1_taxane.rename(columns={'Total': 'Taxane'})

# Pemetrexed组
table1_peme = calc_table1(df_rwc_patients[df_rwc_patients['chemo_class'] == 'Pemetrexed'])
table1_peme = table1_peme.rename(columns={'Total': 'Pemetrexed'})

# 合并
table1_combined = table1_total.merge(table1_taxane, on='Characteristic', how='left')
table1_combined = table1_combined.merge(table1_peme, on='Characteristic', how='left')

# IIT队列表
df_iit_patients = adata[adata.obs['cohort'] == 'IIT'].obs.groupby('sampleID').first().reset_index()
table1_iit = calc_table1(df_iit_patients)
table1_iit = table1_iit.rename(columns={'Total': 'IIT (Taxane)'})

# 全队列表
table1_all = calc_table1(adata.obs.groupby('sampleID').first().reset_index())
table1_all = table1_all.rename(columns={'Total': 'All Patients'})

# 完整Table 1
table1_full = table1_all.merge(table1_iit, on='Characteristic', how='left')
table1_full = table1_full.merge(table1_taxane, on='Characteristic', how='left')
table1_full = table1_full.merge(table1_peme, on='Characteristic', how='left')

# 保存
table1_full.to_csv(cfg.result_path('Table1_baseline_characteristics.csv'), index=False)
print(f"\n  ✅ Table 1已保存: Table1_baseline_characteristics.csv")
print(f"\n  Table 1预览:")
print(table1_full.head(20).to_string())

# ============================================================
# FDR校正汇总
# ============================================================
print("\n" + "="*70)
print("FDR校正结果（Benjamini-Hochberg）")
print("="*70)

if len(p_values_for_fdr) > 0:
    _, q_values, _, _ = multipletests(p_values_for_fdr, method='fdr_bh')
    print(f"\n  总检验数: {len(p_values_for_fdr)}")
    print(f"  {'指标':<40} {'p值':>10} {'q值':>10}")
    print(f"  {'-'*60}")
    for label, p, q in zip(p_value_labels, p_values_for_fdr, q_values):
        sig_p = '*' if p < 0.05 else ''
        sig_q = '*' if q < 0.05 else ''
        print(f"  {label:<40} {p:>10.4f}{sig_p} {q:>10.4f}{sig_q}")
    
    # 保存FDR校正结果
    fdr_df = pd.DataFrame({
        'metric': p_value_labels,
        'p_value': p_values_for_fdr,
        'q_value': q_values
    })
    fdr_df.to_csv(cfg.result_path('GSE241934_fdr_correction_results.csv'), index=False)
    print(f"\n  ✅ FDR校正结果已保存: GSE241934_fdr_correction_results.csv")
else:
    print("  无有效p值用于FDR校正")

# ============================================================
# 汇总
# ============================================================
print("\n" + "="*70)
print("汇总统计")
print("="*70)

summary = {
    'n_patients_total': len(patient_meta),
    'n_patients_iit': patient_meta[patient_meta['cohort'] == 'IIT'].shape[0],
    'n_patients_rwc': patient_meta[patient_meta['cohort'] == 'RWC'].shape[0],
    'n_taxane_total': patient_meta[patient_meta['chemo_class'] == 'Taxane'].shape[0],
    'n_pemetrexed': patient_meta[patient_meta['chemo_class'] == 'Pemetrexed'].shape[0],
    'nklike_cells_total': int((adata.obs['cell.type'] == 'CD8T_Teff_FGFBP2').sum()),
    'myeloid_cells_total': int(adata.obs['cell.type'].str.contains('Macro|Mono', na=False).sum()),
    'functional_modules_tested': 4,
    'spp1_tam_taxane_mean': float(tax_mye['spp1_macro_ratio_in_macro'].mean()),
    'spp1_tam_peme_mean': float(peme_mye['spp1_macro_ratio_in_macro'].mean()),
    'irs_taxane_auc': float(irs_results[0]['auc']) if irs_results and 'auc' in irs_results[0] else np.nan,
    'irs_peme_auc': float(irs_results[1]['auc']) if len(irs_results) > 1 and 'auc' in irs_results[1] else np.nan,
}

print(f"  总患者数: {summary['n_patients_total']}")
print(f"    IIT队列: {summary['n_patients_iit']}")
print(f"    RWC队列: {summary['n_patients_rwc']}")
print(f"    Taxane方案: {summary['n_taxane_total']}")
print(f"    Pemetrexed方案: {summary['n_pemetrexed']}")
print(f"  NK-like细胞数: {summary['nklike_cells_total']}")
print(f"  髓系细胞数: {summary['myeloid_cells_total']}")

df_summary = pd.DataFrame([summary])
df_summary.to_csv(cfg.result_path('GSE241934_mechanism_summary.csv'), index=False)

print("\n" + "="*70)
print("[step6d_mechanism_validation] Completed!")
print("="*70)
