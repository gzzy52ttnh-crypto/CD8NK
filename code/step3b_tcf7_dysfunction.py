"""
Step 3b: 核心机制深钻 — 单细胞+空间双重验证

任务一：单细胞层面 TCF7 功能失调分析
  - TCF7 分层 (q33/q66) × SPP1 环境 (High/Low)
  - 验证：High SPP1 环境下 TCF7_High 但杀伤基因低（分化受阻）

任务二：空间层面"干性-效应解耦"验证
  - Decoupling_Index = z(stem_tcell) - z(nk_cyto_tcell)
  - 与 spp1_rec_interact 的 Spearman 相关
"""
import os, sys, warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, mannwhitneyu
from scipy import stats as ss

sys.path.insert(0, os.path.dirname(__file__))
import config as cfg

RESULT = cfg.RESULT_DIR
os.makedirs(RESULT, exist_ok=True)

# ============================================================
# 任务一：单细胞层面 TCF7 功能失调分析
# ============================================================
print("=" * 60)
print("任务一：单细胞 TCF7 功能失调分析")
print("=" * 60)

nk = sc.read_h5ad(cfg.H5AD_NK)
print(f"Loaded NK-like CD8: {nk.shape}")

# 获取 SPP1+ TAM ratio
mye_pp = pd.read_csv(os.path.join(RESULT, 'myeloid_per_patient.csv'))
spp1_map = dict(zip(mye_pp['patient_id'], mye_pp['SPP1_TAM_ratio']))
print(f"SPP1+ TAM ratio: {len(spp1_map)} patients, median={np.median(list(spp1_map.values())):.4f}")

# SPP1 环境分层：中位数分割
spp1_median = np.median(list(spp1_map.values()))
high_spp1_pats = [p for p, v in spp1_map.items() if v >= spp1_median]
low_spp1_pats = [p for p, v in spp1_map.items() if v < spp1_median]
print(f"  High SPP1 (≥{spp1_median:.4f}): {len(high_spp1_pats)} patients")
print(f"  Low SPP1 (<{spp1_median:.4f}): {len(low_spp1_pats)} patients")

# TCF7 表达
tcf7_idx = list(nk.var_names).index('TCF7')
tcf7_expr = nk.X[:, tcf7_idx]
tcf7_expr = tcf7_expr.toarray().flatten() if hasattr(tcf7_expr, 'toarray') else np.asarray(tcf7_expr).flatten()

# TCF7 分层：q33/q66 — 但 TCF7 在 NK-like 中极度稀疏（>2/3 细胞 = 0）
# 改用非零细胞分位数 + 绝对阈值，确保三组有生物学意义
tcf7_nonzero = tcf7_expr[tcf7_expr > 0]
if len(tcf7_nonzero) > 100:
    # 在非零细胞中取分位数
    q33, q66 = np.quantile(tcf7_nonzero, [0.333, 0.667])
else:
    q33, q66 = 1, 2

# 如果分位数仍为 0，使用绝对阈值
if q33 == 0 and q66 == 0:
    q33, q66 = 1, 2

print(f"TCF7 quantiles (nonzero cells): q33={q33:.4f}, q66={q66:.4f}")
print(f"  TCF7=0 cells: {np.sum(tcf7_expr == 0)} ({np.sum(tcf7_expr == 0)/len(tcf7_expr)*100:.1f}%)")

tcf7_group = np.where(tcf7_expr >= q66, 'TCF7_High',
              np.where(tcf7_expr >= q33, 'TCF7_Mid', 'TCF7_Low'))
nk.obs['tcf7_group'] = tcf7_group
print(f"TCF7 groups: High={np.sum(tcf7_group=='TCF7_High')}, "
      f"Mid={np.sum(tcf7_group=='TCF7_Mid')}, Low={np.sum(tcf7_group=='TCF7_Low')}")

# SPP1 环境标注
patient_arr = nk.obs[cfg.COL_SAMPLE].values
spp1_env = np.array(['Unknown'] * len(patient_arr), dtype=object)
for i, p in enumerate(patient_arr):
    if p in high_spp1_pats:
        spp1_env[i] = 'High_SPP1'
    elif p in low_spp1_pats:
        spp1_env[i] = 'Low_SPP1'
nk.obs['spp1_env'] = spp1_env
print(f"SPP1 env: High={np.sum(spp1_env=='High_SPP1')}, Low={np.sum(spp1_env=='Low_SPP1')}, Unknown={np.sum(spp1_env=='Unknown')}")

# 核心对比：TCF7_High 细胞在 High vs Low SPP1 环境下的杀伤基因表达
cyto_genes = ['GZMB', 'PRF1', 'GNLY']
cyto_data = {}

for gene in cyto_genes:
    gidx = list(nk.var_names).index(gene)
    gexpr = nk.X[:, gidx]
    gexpr = gexpr.toarray().flatten() if hasattr(gexpr, 'toarray') else np.asarray(gexpr).flatten()

    # 只看 TCF7_High 细胞
    tcf7_hi_mask = tcf7_group == 'TCF7_High'
    high_env_mask = (spp1_env == 'High_SPP1') & tcf7_hi_mask
    low_env_mask = (spp1_env == 'Low_SPP1') & tcf7_hi_mask

    cyto_data[gene] = {
        'High_SPP1': gexpr[high_env_mask],
        'Low_SPP1': gexpr[low_env_mask],
    }

    # MWU test
    vals_hi = gexpr[high_env_mask]
    vals_lo = gexpr[low_env_mask]
    if len(vals_hi) > 10 and len(vals_lo) > 10:
        _, mwu_p = mannwhitneyu(vals_hi, vals_lo, alternative='two-sided')
        mean_hi = np.mean(vals_hi)
        mean_lo = np.mean(vals_lo)
        print(f"  {gene} (TCF7_High): High_SPP1 mean={mean_hi:.4f} (n={len(vals_hi)}), "
              f"Low_SPP1 mean={mean_lo:.4f} (n={len(vals_lo)}), MWU p={mwu_p:.4e}")
    else:
        print(f"  {gene}: insufficient cells (High={len(vals_hi)}, Low={len(vals_lo)})")

# 绘图：2x2 布局 — 上排3个小提琴图 + 下排 summary bar chart
fig = plt.figure(figsize=(13, 9))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 0.6])
ax_v1 = fig.add_subplot(gs[0, 0])
ax_v2 = fig.add_subplot(gs[0, 1])
ax_v3 = fig.add_subplot(gs[0, 2])
ax_bar = fig.add_subplot(gs[1, :])
axes = [ax_v1, ax_v2, ax_v3]

for ax, gene in zip(axes, cyto_genes):
    hi_vals = cyto_data[gene]['High_SPP1']
    lo_vals = cyto_data[gene]['Low_SPP1']

    np.random.seed(42)
    sample_n = min(5000, len(hi_vals), len(lo_vals))
    hi_sample = np.random.choice(hi_vals, sample_n, replace=False) if len(hi_vals) > sample_n else hi_vals
    lo_sample = np.random.choice(lo_vals, sample_n, replace=False) if len(lo_vals) > sample_n else lo_vals

    data = [lo_sample, hi_sample]
    labels = ['Low SPP1', 'High SPP1']

    parts = ax.violinplot(data, positions=[1, 2], widths=0.7,
                          showmeans=False, showmedians=False, showextrema=False)
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(['#2E86C1', '#E74C3C'][i])
        pc.set_alpha(0.6)
        pc.set_edgecolor('black')
        pc.set_linewidth(0.8)

    bp = ax.boxplot(data, positions=[1, 2], widths=0.15, patch_artist=True,
                    showfliers=False, zorder=3)
    for patch in bp['boxes']:
        patch.set_facecolor('white')
        patch.set_alpha(0.8)
        patch.set_edgecolor('black')
    for median in bp['medians']:
        median.set_color('black')
        median.set_linewidth(2)

    for i, vals in enumerate([lo_sample, hi_sample]):
        if len(vals) > 200:
            jitter_sample = np.random.choice(vals, 200, replace=False)
        else:
            jitter_sample = vals
        x_jitter = np.random.normal(i+1, 0.05, size=len(jitter_sample))
        ax.scatter(x_jitter, jitter_sample, c='black', s=8, alpha=0.3, zorder=2, edgecolors='none')

    _, mwu_p = mannwhitneyu(hi_vals, lo_vals, alternative='two-sided')
    ax.set_xticks([1, 2])
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel(f'{gene} expression', fontsize=10)
    ax.set_title(f'{gene} (p={mwu_p:.2e})', fontsize=11, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_facecolor('#F8F9FA')
    y_min = max(0, min(lo_sample.min(), hi_sample.min()) - 0.1)
    y_max = max(lo_sample.max(), hi_sample.max()) * 1.1
    ax.set_ylim(y_min, y_max)
    ax.grid(True, alpha=0.3, axis='y')

# 下排：summary bar chart
genes_bar = cyto_genes
high_means = [np.mean(cyto_data[g]['High_SPP1']) for g in genes_bar]
low_means = [np.mean(cyto_data[g]['Low_SPP1']) for g in genes_bar]
high_sems = [np.std(cyto_data[g]['High_SPP1'])/np.sqrt(len(cyto_data[g]['High_SPP1'])) for g in genes_bar]
low_sems = [np.std(cyto_data[g]['Low_SPP1'])/np.sqrt(len(cyto_data[g]['Low_SPP1'])) for g in genes_bar]

x_bar = np.arange(len(genes_bar))
bar_w = 0.35
b1 = ax_bar.bar(x_bar - bar_w/2, low_means, bar_w, yerr=low_sems, color='#2E86C1', alpha=0.8,
                label='Low SPP1 env', edgecolor='black', linewidth=0.5, capsize=4)
b2 = ax_bar.bar(x_bar + bar_w/2, high_means, bar_w, yerr=high_sems, color='#E74C3C', alpha=0.8,
                label='High SPP1 env', edgecolor='black', linewidth=0.5, capsize=4)
ax_bar.set_xticks(x_bar)
ax_bar.set_xticklabels(genes_bar, fontsize=11)
ax_bar.set_ylabel('Mean expression (TCF7$^{High}$ NK-like)', fontsize=10)
ax_bar.set_title('Summary: Cytotoxicity (High vs Low SPP1 env)\n(Lower in High SPP1 = differentiation arrest)', fontsize=11)
ax_bar.legend(loc='upper right', fontsize=10)
ax_bar.spines['top'].set_visible(False); ax_bar.spines['right'].set_visible(False)
ax_bar.set_facecolor('#F8F9FA')
ax_bar.grid(True, alpha=0.3, axis='y')
for bar, v in zip(b1, low_means):
    ax_bar.text(bar.get_x() + bar.get_width()/2, v + max(high_means)*0.02, f'{v:.2f}',
                ha='center', va='bottom', fontsize=9, color='#1A5276')
for bar, v in zip(b2, high_means):
    ax_bar.text(bar.get_x() + bar.get_width()/2, v + max(high_means)*0.02, f'{v:.2f}',
                ha='center', va='bottom', fontsize=9, color='#922B21')

plt.suptitle('TCF7$^{High}$ NK-like cells: Cytotoxicity in High vs Low SPP1 environment\n'
             '(Hypothesis: High SPP1 → TCF7$^+$ but cytotoxic-low = differentiation arrest)',
             fontsize=12, y=0.995, fontweight='bold')
plt.tight_layout()
fig.savefig(os.path.join(RESULT, 'Fig4_TCF7_dysfunction.png'), dpi=150, bbox_inches='tight')
fig.savefig(os.path.join(RESULT, 'Fig4_TCF7_dysfunction.pdf'), dpi=150, bbox_inches='tight')
plt.close(fig)
print("Fig4_TCF7_dysfunction saved")

# 保存统计结果
tcf7_stats = []
for gene in cyto_genes:
    hi_vals = cyto_data[gene]['High_SPP1']
    lo_vals = cyto_data[gene]['Low_SPP1']
    _, mwu_p = mannwhitneyu(hi_vals, lo_vals, alternative='two-sided')
    tcf7_stats.append({
        'gene': gene,
        'tcf7_group': 'TCF7_High',
        'high_spp1_mean': float(np.mean(hi_vals)),
        'low_spp1_mean': float(np.mean(lo_vals)),
        'high_spp1_n': len(hi_vals),
        'low_spp1_n': len(lo_vals),
        'mwu_p': mwu_p,
    })
tcf7_stats_df = pd.DataFrame(tcf7_stats)
tcf7_stats_df.to_csv(os.path.join(RESULT, 'tcf7_dysfunction_stats.csv'), index=False)
print("tcf7_dysfunction_stats.csv saved")

# ============================================================
# 任务二：空间层面"干性-效应解耦"验证
# ============================================================
print("\n" + "=" * 60)
print("任务二：空间层面干性-效应解耦验证")
print("=" * 60)

sp = pd.read_csv(os.path.join(RESULT, 'spatial_roi_scores.csv'))
print(f"Loaded spatial data: {sp.shape}")

# Z-score 标准化
stem_z = (sp['stem_tcell'] - sp['stem_tcell'].mean()) / sp['stem_tcell'].std()
cyto_z = (sp['nk_cyto_tcell'] - sp['nk_cyto_tcell'].mean()) / sp['nk_cyto_tcell'].std()

# 解耦指数：> 0 = 干性高但杀伤低（阻滞状态）
sp['decoupling_index'] = stem_z - cyto_z

# Spearman 相关
r_dec, p_dec = spearmanr(sp['spp1_rec_interact'], sp['decoupling_index'])
print(f"Decoupling Index vs SPP1 interaction:")
print(f"  Spearman r={r_dec:.4f}, p={p_dec:.4f}")
print(f"  Decoupling Index: mean={sp['decoupling_index'].mean():.4f}, "
      f"range=[{sp['decoupling_index'].min():.4f}, {sp['decoupling_index'].max():.4f}]")

# 按 SPP1 分组比较解耦指数
high_spp1_roi = sp[sp['spp1_grp'] == 'High']
low_spp1_roi = sp[sp['spp1_grp'] == 'Low']
if len(high_spp1_roi) > 2 and len(low_spp1_roi) > 2:
    _, mwu_dec = mannwhitneyu(high_spp1_roi['decoupling_index'],
                               low_spp1_roi['decoupling_index'],
                               alternative='greater')
    print(f"  High SPP1 vs Low SPP1 ROI: MWU p={mwu_dec:.4f} (one-sided greater)")
    print(f"    High SPP1 mean DI={high_spp1_roi['decoupling_index'].mean():.4f}")
    print(f"    Low SPP1 mean DI={low_spp1_roi['decoupling_index'].mean():.4f}")

# 绘图：散点图 + 回归线
fig2, ax = plt.subplots(figsize=(7, 6))

x = sp['spp1_rec_interact'].values
y = sp['decoupling_index'].values

# 散点，按 SPP1 分组着色
colors = []
for g in sp['spp1_grp']:
    if g == 'High':
        colors.append('#E74C3C')
    elif g == 'Low':
        colors.append('#2E86C1')
    else:
        colors.append('#95A5A6')

ax.scatter(x, y, c=colors, s=60, alpha=0.7, edgecolors='black', linewidth=0.5)

# 回归线
mask = ~np.isnan(x) & ~np.isnan(y)
if mask.sum() > 3:
    slope, intercept, r_val, p_val, stderr = ss.linregress(x[mask], y[mask])
    x_line = np.linspace(x[mask].min(), x[mask].max(), 100)
    ax.plot(x_line, slope * x_line + intercept, 'k--', lw=1.5, alpha=0.8)
    ax.text(0.05, 0.95, f'Spearman r={r_dec:.3f}, p={p_dec:.4f}\nLinear R²={r_val**2:.3f}',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# 零线
ax.axhline(0, color='gray', ls=':', alpha=0.5)
ax.set_xlabel('SPP1 Interaction Score (spp1_rec_interact)')
ax.set_ylabel('Decoupling Index (z[stem] - z[cytotoxicity])')
ax.set_title('Stemness-Effector Decoupling vs SPP1 Signaling\n'
             '(DI > 0: stemness-high/cytotoxicity-low = arrest state)')

# 图例
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#E74C3C', markersize=10, label='High SPP1 ROI'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#2E86C1', markersize=10, label='Low SPP1 ROI'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#95A5A6', markersize=10, label='Mid SPP1 ROI'),
]
ax.legend(handles=legend_elements, loc='lower right')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
fig2.savefig(os.path.join(RESULT, 'FigS_spatial_decoupling.png'), dpi=150, bbox_inches='tight')
fig2.savefig(os.path.join(RESULT, 'FigS_spatial_decoupling.pdf'), dpi=150, bbox_inches='tight')
plt.close(fig2)
print("FigS_spatial_decoupling saved")

# 保存空间解耦统计
decoupling_stats = pd.DataFrame({
    'metric': ['spearman_r', 'spearman_p', 'high_spp1_mean_DI', 'low_spp1_mean_DI',
               'mwu_high_vs_low_p', 'n_roi', 'n_high_spp1', 'n_low_spp1'],
    'value': [r_dec, p_dec,
              high_spp1_roi['decoupling_index'].mean() if len(high_spp1_roi) > 0 else np.nan,
              low_spp1_roi['decoupling_index'].mean() if len(low_spp1_roi) > 0 else np.nan,
              mwu_dec if len(high_spp1_roi) > 2 and len(low_spp1_roi) > 2 else np.nan,
              len(sp), len(high_spp1_roi), len(low_spp1_roi)]
})
decoupling_stats.to_csv(os.path.join(RESULT, 'spatial_decoupling_stats.csv'), index=False)
print("spatial_decoupling_stats.csv saved")

# 保存含 decoupling_index 的空间数据
sp.to_csv(os.path.join(RESULT, 'spatial_roi_scores.csv'), index=False)
print("spatial_roi_scores.csv updated (with decoupling_index)")

print("\n[step3b_tcf7_dysfunction] Done.")
