"""
step4d_chemo_mechanism.py
化疗方案特异性机制深度探索

背景：SPP1+ TAM 比例在 Taxane vs Pemetrexed 间无显著差异（p=0.1686），
无法解释方案特异性。需从以下 4 个维度补充分析：
1. SPP1+ TAM 极化状态（M1/M2）—— 功能差异而非比例差异
2. NK-like 克隆大小分布 —— 克隆扩增程度差异
3. 髓系亚群组成比较 —— 寻找其他髓系差异
4. SPP1 表达强度比较 —— 表达量差异

输出：chemo_mechanism_results.csv, FigS7_chemo_mechanism.png/pdf
"""
import config as cfg
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, spearmanr, fisher_exact
from scipy.sparse import issparse
from statsmodels.stats.multitest import multipletests
from _common import score_gene_set
import os, warnings
warnings.filterwarnings('ignore')

print("[step4d_chemo_mechanism] Starting...")
os.makedirs(cfg.RESULT_DIR, exist_ok=True)

# ============================================================
# 1. 加载数据与化疗方案分组
# ============================================================
print("\n" + "="*60)
print("1. 加载数据与化疗方案分组")
print("="*60)

# T 细胞数据（克隆分析用）
adata_t = sc.read_h5ad(cfg.H5AD_T)
print(f"T cells: {adata_t.shape}")

# 髓系数据（SPP1+TAM 极化用）
myeloid = sc.read_h5ad(cfg.H5AD_MYELOID)
print(f"Myeloid: {myeloid.shape}")

# 患者化疗方案
chemo_pat = adata_t.obs.drop_duplicates('sampleID')[['sampleID', 'chemotherapy']].copy()
chemo_pat.columns = ['patient_id', 'chemotherapy']
del adata_t

from _common import classify_chemo

chemo_pat['chemo_class'] = chemo_pat['chemotherapy'].apply(classify_chemo)
tax_pats = set(chemo_pat[chemo_pat['chemo_class']=='Platinum+Taxane']['patient_id'])
peme_pats = set(chemo_pat[chemo_pat['chemo_class']=='Platinum+Pemetrexed']['patient_id'])
print(f"Taxane patients: {len(tax_pats)}, Pemetrexed patients: {len(peme_pats)}")

results = {}
p_values_for_fdr = []
p_value_labels = []

# ============================================================
# 2. SPP1+ TAM 极化状态（M1/M2 基因集评分）
# ============================================================
print("\n" + "="*60)
print("2. SPP1+ TAM 极化状态 (M1/M2)")
print("="*60)

# M1/M2 基因集（经典巨噬细胞极化标记）
M1_GENES = ['CD80', 'CD86', 'IL1B', 'TNF', 'IL6', 'NOS2', 'CXCL9', 'CXCL10', 'HLA-DRA', 'HLA-DRB1']
M2_GENES = ['CD163', 'MSR1', 'MRC1', 'TGFB1', 'IL10', 'ARG1', 'VEGFA', 'CCL2', 'CCL22', 'FN1']
TAM_SUPPRESSIVE = ['SPP1', 'LGALS3', 'TREM2', 'APOE', 'CD9', 'PDCD1LG2', 'CD274']

# 筛选 SPP1+ TAM 细胞（Mac/Mono 中 SPP1>0）
mac_mono_mask = myeloid.obs[cfg.COL_CELL_TYPE].str.contains('Mφ|Mac|Mono', regex=True, na=False).values
spp1_idx = list(myeloid.var_names).index('SPP1') if 'SPP1' in myeloid.var_names else -1
if spp1_idx >= 0:
    spp1_expr = myeloid.X[:, spp1_idx]
    spp1_expr = spp1_expr.toarray().flatten() if issparse(spp1_expr) else np.asarray(spp1_expr).flatten()
    spp1_pos_mask = (spp1_expr > 0) & mac_mono_mask
else:
    spp1_pos_mask = mac_mono_mask

spp1_tam = myeloid[spp1_pos_mask].copy()
print(f"SPP1+ TAM cells: {spp1_tam.n_obs}")

# 计算 M1/M2 评分（基因集平均表达）
m1_valid = [g for g in M1_GENES if g in spp1_tam.var_names]
m2_valid = [g for g in M2_GENES if g in spp1_tam.var_names]
supp_valid = [g for g in TAM_SUPPRESSIVE if g in spp1_tam.var_names]
print(f"  M1: {len(m1_valid)}/{len(M1_GENES)} genes found")
print(f"  M2: {len(m2_valid)}/{len(M2_GENES)} genes found")
print(f"  Suppressive: {len(supp_valid)}/{len(TAM_SUPPRESSIVE)} genes found")

m1_score = score_gene_set(spp1_tam, M1_GENES)
m2_score = score_gene_set(spp1_tam, M2_GENES)
supp_score = score_gene_set(spp1_tam, TAM_SUPPRESSIVE)

# 按患者分组，再按化疗方案比较
spp1_tam.obs['M1_score'] = m1_score
spp1_tam.obs['M2_score'] = m2_score
spp1_tam.obs['Supp_score'] = supp_score
spp1_tam.obs['M1_M2_ratio'] = m1_score / (m2_score + 1e-10)

# 患者水平均值
patient_ids = []
patient_m1 = []
patient_m2 = []
patient_supp = []
patient_m1m2 = []
patient_group = []
patient_ncells = []

for p in spp1_tam.obs[cfg.COL_SAMPLE].unique():
    p_mask = spp1_tam.obs[cfg.COL_SAMPLE] == p
    n = p_mask.sum()
    if n < 5:
        continue
    if p in tax_pats:
        grp = 'Taxane'
    elif p in peme_pats:
        grp = 'Pemetrexed'
    else:
        continue
    patient_ids.append(p)
    patient_group.append(grp)
    patient_ncells.append(n)
    patient_m1.append(np.mean(m1_score[p_mask]))
    patient_m2.append(np.mean(m2_score[p_mask]))
    patient_supp.append(np.mean(supp_score[p_mask]))
    patient_m1m2.append(np.mean(m1_score[p_mask]) / (np.mean(m2_score[p_mask]) + 1e-10))

df_polar = pd.DataFrame({
    'patient': patient_ids,
    'group': patient_group,
    'n_cells': patient_ncells,
    'M1_score': patient_m1,
    'M2_score': patient_m2,
    'M1_M2_ratio': patient_m1m2,
    'Suppressive_score': patient_supp,
})

tax_polar = df_polar[df_polar['group']=='Taxane']
peme_polar = df_polar[df_polar['group']=='Pemetrexed']

print(f"\nTaxane SPP1+TAM patients: {len(tax_polar)}")
print(f"Pemetrexed SPP1+TAM patients: {len(peme_polar)}")

for metric in ['M1_score', 'M2_score', 'M1_M2_ratio', 'Suppressive_score']:
    tax_vals = tax_polar[metric].values
    peme_vals = peme_polar[metric].values
    stat, pval = mannwhitneyu(tax_vals, peme_vals, alternative='two-sided')
    results[f'SPP1_TAM_{metric}_tax_mean'] = np.mean(tax_vals)
    results[f'SPP1_TAM_{metric}_peme_mean'] = np.mean(peme_vals)
    results[f'SPP1_TAM_{metric}_pvalue'] = pval
    p_values_for_fdr.append(pval)
    p_value_labels.append(f'SPP1_TAM_{metric}')
    print(f"  {metric}: Taxane={np.mean(tax_vals):.4f}, Peme={np.mean(peme_vals):.4f}, MWU p={pval:.4f}")

# ============================================================
# 3. NK-like 克隆大小分布（Taxane vs Pemetrexed）
# ============================================================
print("\n" + "="*60)
print("3. NK-like 克隆大小分布")
print("="*60)

# 重新加载 T 细胞数据（克隆分析）
adata_t = sc.read_h5ad(cfg.H5AD_T)
clone_obs = adata_t.obs
nk_types = ['CD8T_NK-like_FGFBP2']
tex_types = ['CD8T_Tex_CXCL13', 'CD8T_terminal_Tex_LAYN']

clone_groups = clone_obs.groupby(cfg.COL_CLONOTYPE)
clone_df = clone_groups.agg(
    clone_size=(cfg.COL_CLONOTYPE_NUM, 'first'),
    nk_count=(cfg.COL_CELL_TYPE, lambda x: x.isin(nk_types).sum()),
    tex_count=(cfg.COL_CELL_TYPE, lambda x: x.isin(tex_types).sum()),
).reset_index()
clone_df['total'] = clone_df['nk_count'] + clone_df['tex_count']
clone_df['nk_ratio'] = clone_df['nk_count'] / clone_df['total'].replace(0, np.nan)
clone_df['nk_ratio'] = clone_df['nk_ratio'].fillna(0)
clone_df['category'] = 'Other'
_nk_thresh = cfg.NK_DOMINANT_RATIO_THRESHOLD
clone_df.loc[clone_df['nk_ratio'] >= _nk_thresh, 'category'] = 'NK-dominant'
clone_df.loc[(clone_df['nk_ratio'] < _nk_thresh) & (clone_df['tex_count'] / clone_df['total'].replace(0, np.nan) >= _nk_thresh), 'category'] = 'Tex-dominant'

# 克隆→患者分配（多数规则）
clone_pat = clone_obs.groupby([cfg.COL_CLONOTYPE, cfg.COL_SAMPLE]).size().reset_index(name='n')
clone_pat = clone_pat.sort_values('n', ascending=False)
clone_pat_map = clone_pat.drop_duplicates(cfg.COL_CLONOTYPE).set_index(cfg.COL_CLONOTYPE)[cfg.COL_SAMPLE]

# 大克隆（≥ BIG_CLONE_THRESHOLD）且 NK-dominant
big_nk_clones = clone_df[(clone_df['clone_size'] >= cfg.BIG_CLONE_THRESHOLD) & (clone_df['category'] == 'NK-dominant')].copy()
big_nk_clones['patient'] = big_nk_clones[cfg.COL_CLONOTYPE].map(clone_pat_map)

# 按患者分组，计算 NK-dominant 大克隆的平均克隆大小
pat_clone_stats = []
for p in set(tax_pats) | set(peme_pats):
    p_clones = big_nk_clones[big_nk_clones['patient'] == p]
    if len(p_clones) == 0:
        continue
    if p in tax_pats:
        grp = 'Taxane'
    else:
        grp = 'Pemetrexed'
    pat_clone_stats.append({
        'patient': p,
        'group': grp,
        'n_big_nk_clones': len(p_clones),
        'mean_clone_size': p_clones['clone_size'].mean(),
        'median_clone_size': p_clones['clone_size'].median(),
        'max_clone_size': p_clones['clone_size'].max(),
        'total_nk_cells': p_clones['clone_size'].sum(),
    })

df_clone = pd.DataFrame(pat_clone_stats)
tax_clone = df_clone[df_clone['group']=='Taxane']
peme_clone = df_clone[df_clone['group']=='Pemetrexed']

print(f"\nTaxane patients with NK-dominant big clones: {len(tax_clone)}")
print(f"Pemetrexed patients with NK-dominant big clones: {len(peme_clone)}")

for metric in ['mean_clone_size', 'median_clone_size', 'max_clone_size', 'n_big_nk_clones', 'total_nk_cells']:
    tax_vals = tax_clone[metric].values
    peme_vals = peme_clone[metric].values
    if len(tax_vals) < 3 or len(peme_vals) < 3:
        print(f"  {metric}: 样本不足 (Taxane={len(tax_vals)}, Peme={len(peme_vals)})")
        continue
    stat, pval = mannwhitneyu(tax_vals, peme_vals, alternative='two-sided')
    results[f'Clone_{metric}_tax_mean'] = np.mean(tax_vals)
    results[f'Clone_{metric}_peme_mean'] = np.mean(peme_vals)
    results[f'Clone_{metric}_pvalue'] = pval
    p_values_for_fdr.append(pval)
    p_value_labels.append(f'Clone_{metric}')
    print(f"  {metric}: Taxane={np.mean(tax_vals):.2f}, Peme={np.mean(peme_vals):.2f}, MWU p={pval:.4f}")

# ============================================================
# 4. 髓系亚群组成比较
# ============================================================
print("\n" + "="*60)
print("4. 髓系亚群组成比较")
print("="*60)

# 每患者髓系亚群比例
mye_obs = myeloid.obs
mye_subtypes = mye_obs[cfg.COL_CELL_TYPE].value_counts().index.tolist()
print(f"  Myeloid subtypes: {len(mye_subtypes)}")

pat_mye_comp = []
for p in set(tax_pats) | set(peme_pats):
    p_mask = mye_obs[cfg.COL_SAMPLE] == p
    n_total = p_mask.sum()
    if n_total < 20:
        continue
    if p in tax_pats:
        grp = 'Taxane'
    else:
        grp = 'Pemetrexed'
    row = {'patient': p, 'group': grp, 'total_mye': n_total}
    for st in mye_subtypes[:8]:  # top 8
        row[st] = (mye_obs.loc[p_mask, cfg.COL_CELL_TYPE] == st).sum() / n_total
    pat_mye_comp.append(row)

df_mye_comp = pd.DataFrame(pat_mye_comp)
tax_mye = df_mye_comp[df_mye_comp['group']=='Taxane']
peme_mye = df_mye_comp[df_mye_comp['group']=='Pemetrexed']

print(f"\nTaxane patients with myeloid data: {len(tax_mye)}")
print(f"Pemetrexed patients with myeloid data: {len(peme_mye)}")

for st in mye_subtypes[:8]:
    tax_vals = tax_mye[st].values
    peme_vals = peme_mye[st].values
    if len(tax_vals) < 3 or len(peme_vals) < 3:
        continue
    stat, pval = mannwhitneyu(tax_vals, peme_vals, alternative='two-sided')
    results[f'Mye_{st}_tax_mean'] = np.mean(tax_vals)
    results[f'Mye_{st}_peme_mean'] = np.mean(peme_vals)
    results[f'Mye_{st}_pvalue'] = pval
    p_values_for_fdr.append(pval)
    p_value_labels.append(f'Mye_{st}')
    sig = ' *' if pval < 0.05 else ''
    print(f"  {st}: Taxane={np.mean(tax_vals):.4f}, Peme={np.mean(peme_vals):.4f}, p={pval:.4f}{sig}")

# ============================================================
# 5. SPP1 表达强度比较
# ============================================================
print("\n" + "="*60)
print("5. SPP1 表达强度比较")
print("="*60)

# SPP1+ TAM 细胞中 SPP1 的平均表达强度（不是比例）
spp1_tam_pat = []
for p in set(tax_pats) | set(peme_pats):
    p_mask = (spp1_tam.obs[cfg.COL_SAMPLE] == p).values
    n = p_mask.sum()
    if n < 5:
        continue
    if p in tax_pats:
        grp = 'Taxane'
    else:
        grp = 'Pemetrexed'
    spp1_vals = spp1_tam.X[p_mask, spp1_idx] if spp1_idx >= 0 else np.zeros(n)
    spp1_vals = spp1_vals.toarray().flatten() if issparse(spp1_vals) else np.asarray(spp1_vals).flatten()
    spp1_tam_pat.append({
        'patient': p,
        'group': grp,
        'n_spp1_cells': n,
        'SPP1_mean_expr': np.mean(spp1_vals),
        'SPP1_median_expr': np.median(spp1_vals),
    })

df_spp1_expr = pd.DataFrame(spp1_tam_pat)
tax_spp1 = df_spp1_expr[df_spp1_expr['group']=='Taxane']
peme_spp1 = df_spp1_expr[df_spp1_expr['group']=='Pemetrexed']

for metric in ['SPP1_mean_expr', 'SPP1_median_expr']:
    tax_vals = tax_spp1[metric].values
    peme_vals = peme_spp1[metric].values
    stat, pval = mannwhitneyu(tax_vals, peme_vals, alternative='two-sided')
    results[f'{metric}_tax_mean'] = np.mean(tax_vals)
    results[f'{metric}_peme_mean'] = np.mean(peme_vals)
    results[f'{metric}_pvalue'] = pval
    p_values_for_fdr.append(pval)
    p_value_labels.append(f'{metric}')
    print(f"  {metric}: Taxane={np.mean(tax_vals):.4f}, Peme={np.mean(peme_vals):.4f}, MWU p={pval:.4f}")

# ============================================================
# 6. FDR校正 + 保存结果 + 绘图
# ============================================================
print("\n" + "="*60)
print("6. FDR校正 + 保存结果 + 绘图")
print("="*60)

# FDR校正（Benjamini-Hochberg）
if len(p_values_for_fdr) > 0:
    _, q_values, _, _ = multipletests(p_values_for_fdr, method='fdr_bh')
    print(f"\n  FDR校正结果 ({len(p_values_for_fdr)} tests):")
    print(f"  {'指标':<40} {'p值':>10} {'q值':>10}")
    print(f"  {'-'*60}")
    for label, p, q in zip(p_value_labels, p_values_for_fdr, q_values):
        sig_p = '*' if p < 0.05 else ''
        sig_q = '*' if q < 0.05 else ''
        print(f"  {label:<40} {p:>10.4f}{sig_p} {q:>10.4f}{sig_q}")
    for label, q in zip(p_value_labels, q_values):
        results[f'{label}_qvalue'] = q
else:
    print("  无有效p值用于FDR校正")

# 保存结果 CSV
res_df = pd.DataFrame({'metric': list(results.keys()), 'value': list(results.values())})
res_df.to_csv(cfg.result_path('chemo_mechanism_results.csv'), index=False)
print(f"chemo_mechanism_results.csv saved ({len(res_df)} metrics)")

# 绘图：4 个 panel
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel A: SPP1+ TAM 极化状态（M1/M2 ratio）
axA = axes[0, 0]
tax_vals_a = tax_polar['M1_M2_ratio'].values
peme_vals_a = peme_polar['M1_M2_ratio'].values
bp_a = axA.boxplot([tax_vals_a, peme_vals_a], patch_artist=True, labels=['Taxane', 'Pemetrexed'])
bp_a['boxes'][0].set_facecolor('#3498DB'); bp_a['boxes'][1].set_facecolor('#E74C3C')
stat_a, p_a = mannwhitneyu(tax_vals_a, peme_vals_a, alternative='two-sided')
axA.set_title(f'Panel A: SPP1+ TAM M1/M2 Ratio\nMWU p={p_a:.4f}')
axA.set_ylabel('M1/M2 ratio (patient-level)')

# Panel B: NK-dominant 大克隆平均大小
axB = axes[0, 1]
tax_vals_b = tax_clone['mean_clone_size'].values
peme_vals_b = peme_clone['mean_clone_size'].values
bp_b = axB.boxplot([tax_vals_b, peme_vals_b], patch_artist=True, labels=['Taxane', 'Pemetrexed'])
bp_b['boxes'][0].set_facecolor('#3498DB'); bp_b['boxes'][1].set_facecolor('#E74C3C')
stat_b, p_b = mannwhitneyu(tax_vals_b, peme_vals_b, alternative='two-sided')
axB.set_title(f'Panel B: NK-dominant Clone Mean Size\nMWU p={p_b:.4f}')
axB.set_ylabel('Mean clone size (big clones ≥5)')

# Panel C: 髓系亚群比例（Top 4）
axC = axes[1, 0]
top4_subtypes = [st for st in mye_subtypes[:4] if st in df_mye_comp.columns]
x = np.arange(len(top4_subtypes))
width = 0.35
tax_means = [tax_mye[st].mean() for st in top4_subtypes]
peme_means = [peme_mye[st].mean() for st in top4_subtypes]
axC.bar(x - width/2, tax_means, width, label='Taxane', color='#3498DB', alpha=0.7)
axC.bar(x + width/2, peme_means, width, label='Pemetrexed', color='#E74C3C', alpha=0.7)
axC.set_xticks(x)
axC.set_xticklabels([s[:20] for s in top4_subtypes], rotation=30, ha='right', fontsize=7)
axC.set_title('Panel C: Myeloid Subtype Proportions')
axC.set_ylabel('Mean proportion')
axC.legend()

# Panel D: SPP1 表达强度
axD = axes[1, 1]
tax_vals_d = tax_spp1['SPP1_mean_expr'].values
peme_vals_d = peme_spp1['SPP1_mean_expr'].values
bp_d = axD.boxplot([tax_vals_d, peme_vals_d], patch_artist=True, labels=['Taxane', 'Pemetrexed'])
bp_d['boxes'][0].set_facecolor('#3498DB'); bp_d['boxes'][1].set_facecolor('#E74C3C')
stat_d, p_d = mannwhitneyu(tax_vals_d, peme_vals_d, alternative='two-sided')
axD.set_title(f'Panel D: SPP1 Expression in SPP1+ TAM\nMWU p={p_d:.4f}')
axD.set_ylabel('Mean SPP1 expression')

plt.tight_layout()
fig.savefig(cfg.result_path('FigS7_chemo_mechanism.pdf'), dpi=150, bbox_inches='tight')
fig.savefig(cfg.result_path('FigS7_chemo_mechanism.png'), dpi=150, bbox_inches='tight')
plt.close(fig)
print("FigS7_chemo_mechanism saved")

# ============================================================
# 7. 小结
# ============================================================
print("\n" + "="*60)
print("7. 小结")
print("="*60)
print(f"SPP1+ TAM 比例: Taxane=0.2098, Peme=0.2454, p=0.1686 (NS)")
print(f"M1/M2 ratio: p={p_a:.4f} {'*' if p_a<0.05 else 'NS'}")
print(f"NK-dominant clone size: p={p_b:.4f} {'*' if p_b<0.05 else 'NS'}")
print(f"SPP1 expression: p={p_d:.4f} {'*' if p_d<0.05 else 'NS'}")
print("\n[step4d_chemo_mechanism] Done.")
