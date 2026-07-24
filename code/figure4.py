#!/usr/bin/env python3
"""
Figure 4 — 机制链条验证（分子层面）
  Panel A: 散点图 — SPP1+ TAM 比例 vs NK-Locked Ratio（颜色区分响应组）
  Panel B: 小提琴图 — High-SPP1 vs Low-SPP1 的 KLF2 等干性基因表达
  Panel C: 干性基因表达比较 — NK-like vs Tex 中 KLF2/TCF7/SELL
"""
import os, sys
import numpy as np, pandas as pd, scipy.stats as ss
from statsmodels.stats.multitest import multipletests
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _common import load_h5ad, normalize_log1p, gene_mean_per_cell, embedding

DATA = os.path.dirname(HERE)
ADATA = os.path.join(DATA, 'adata')
RESULT = os.path.join(DATA, 'result')
os.makedirs(RESULT, exist_ok=True)

# ============= 加载数据 =============
print('=== Loading data ===', flush=True)

mp = pd.read_csv(os.path.join(RESULT, 'myeloid_per_patient.csv')).dropna(subset=['spp1_tam_frac'])
ppt = pd.read_csv(os.path.join(RESULT, 'per_patient_metrics.csv'))
print(f' myeloid_per_patient.csv patients: {len(mp)}', flush=True)
print(f' per_patient_metrics.csv patients: {len(ppt)}', flush=True)
merged = mp.merge(ppt[['patient','nk_dominant_ratio', 'nk_pct']], on='patient', how='inner')
merged = merged.dropna(subset=['spp1_tam_frac','nk_dominant_ratio'])
print(f'Merged patients: {len(merged)}', flush=True)
mp_only = set(mp['patient']) - set(ppt['patient'])
ppt_only = set(ppt['patient']) - set(mp['patient'])
if mp_only or ppt_only:
    print(f'NOTE: Patients only in myeloid table: {len(mp_only)} {sorted(mp_only)}', flush=True)
    print(f'NOTE: Patients only in per_patient_metrics: {len(ppt_only)} {sorted(ppt_only)}', flush=True)

# ============= Panel A: 散点图 =============
r_pearson, rp_pearson = ss.pearsonr(merged['spp1_tam_frac'], merged['nk_dominant_ratio'])
rs_spearman, rsp_spearman = ss.spearmanr(merged['spp1_tam_frac'], merged['nk_dominant_ratio'])
print(f'Pearson r={r_pearson:.3f} p={rp_pearson:.3e}  Spearman r={rs_spearman:.3f} p={rsp_spearman:.3e}', flush=True)

# 高低SPP1组
mp_sorted = mp.sort_values('spp1_tam_frac')
n = len(mp_sorted)
high_spp1 = set(mp_sorted.tail(n//3)['patient'])
low_spp1 = set(mp_sorted.head(n//3)['patient'])
print(f'High SPP1: {len(high_spp1)} patients, Low SPP1: {len(low_spp1)} patients', flush=True)

# ============= 加载 T 细胞数据 =============
print('=== Loading and sampling T cells ===', flush=True)
t_adata = load_h5ad(os.path.join(ADATA, 'GSE243013_T_cells.h5ad'))
t_adata.var_names = t_adata.var_names.astype(str)
for c in ['sub_cell_type','sampleID']:
    t_adata.obs[c] = t_adata.obs[c].astype(str)

t_adata.obs['is_cd8'] = t_adata.obs['sub_cell_type'].str.startswith('CD8T')
t_adata.obs['is_nklike'] = t_adata.obs['sub_cell_type'] == 'CD8T_NK-like_FGFBP2'
t_adata.obs['is_tex'] = t_adata.obs['sub_cell_type'].str.contains('exhaust|Tex|Terminal|GZMB|ZNF683', case=False, na=False)

# 采样加速
rng = np.random.default_rng(42)
sample_cells = []
for pid, sub in t_adata.obs[t_adata.obs['is_cd8']].groupby('sampleID'):
    n_sample = min(len(sub), 200)
    if len(sub) > 0:
        sample_cells.extend(list(rng.choice(sub.index, n_sample, replace=False)))
print(f'Sampled {len(sample_cells)} CD8+ T cells', flush=True)

sample_pos = {c: i for i, c in enumerate(t_adata.obs.index)}
sample_idx = [sample_pos[c] for c in sample_cells]
t_sample = t_adata[sample_idx].to_memory()
t_logX_sample = normalize_log1p(t_sample)
t_sample.obs = t_sample.obs.reset_index(drop=True)

stem_genes = ['KLF2','TCF7','LEF1','IL7R','SELL','CCR7']
stem_present = [g for g in stem_genes if g in t_sample.var_names]
print(f'Stem genes present: {stem_present}', flush=True)

# ============= Panel B: High vs Low SPP1 的干性基因表达 =============
print('=== Computing per-patient gene expression (High/Low SPP1) ===', flush=True)

gene_rows = []
obs_names = list(t_sample.obs.index)
name_to_pos = {n: i for i, n in enumerate(obs_names)}
for pid, sub in t_sample.obs.groupby('sampleID'):
    if pid not in high_spp1 and pid not in low_spp1:
        continue
    cd8p = sub[sub['is_cd8']]
    if len(cd8p) == 0:
        continue
    cells_pos = [name_to_pos[n] for n in cd8p.index]
    subX = t_logX_sample[cells_pos].toarray() if hasattr(t_logX_sample[cells_pos], 'toarray') else np.asarray(t_logX_sample[cells_pos])
    r = {'patient': pid, 'group': 'High SPP1+ TAM' if pid in high_spp1 else 'Low SPP1+ TAM'}
    for g in stem_present:
        gi = list(t_sample.var_names).index(g)
        r[g] = float(subX[:, gi].mean())
    gene_rows.append(r)

gene_df = pd.DataFrame(gene_rows)

# 统计检验
gene_stats = []
for g in stem_present:
    a = gene_df[gene_df.group=='High SPP1+ TAM'][g]
    b = gene_df[gene_df.group=='Low SPP1+ TAM'][g]
    pu = ss.mannwhitneyu(a, b, alternative='two-sided').pvalue
    gene_stats.append({'gene': g, 'pvalue': pu, 'mean_high': a.mean(), 'mean_low': b.mean()})
gene_stats_df = pd.DataFrame(gene_stats)
gene_stats_df['FDR'] = multipletests(gene_stats_df['pvalue'], method='fdr_bh')[1]
print(gene_stats_df.to_string(), flush=True)

# ============= Panel C: NK-like vs Tex 的干性基因表达 =============
print('=== Computing NK-like vs Tex gene expression ===', flush=True)

nklike_mask = t_sample.obs['is_nklike'].values
tex_mask = t_sample.obs['is_tex'].values
print(f'NK-like cells: {nklike_mask.sum()}, Tex cells: {tex_mask.sum()}', flush=True)

t_arr = t_logX_sample.toarray() if hasattr(t_logX_sample, 'toarray') else np.asarray(t_logX_sample)

# Panel C 要画的基因
panel_c_genes = [g for g in ['KLF2','TCF7','SELL','LEF1','IL7R','CCR7'] if g in stem_present]
nk_vs_tex_stats = []
for g in panel_c_genes:
    gi = list(t_sample.var_names).index(g)
    nk_expr = t_arr[nklike_mask, gi]
    tex_expr = t_arr[tex_mask, gi]
    pu = ss.mannwhitneyu(nk_expr, tex_expr, alternative='two-sided').pvalue
    nk_vs_tex_stats.append({'gene': g, 'pvalue': pu, 'mean_nklike': nk_expr.mean(), 'mean_tex': tex_expr.mean()})
nk_vs_tex_df = pd.DataFrame(nk_vs_tex_stats)
print(nk_vs_tex_df.to_string(), flush=True)

# ============= 组合绘图 =============
print('=== Composing Figure 4 ===', flush=True)
fig = plt.figure(figsize=(16, 12))
gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3, height_ratios=[0.95, 1])

# ---- Panel A: 散点图 ----
ax1 = fig.add_subplot(gs[0, 0])

colors_map = {'pCR': '#E74C3C', 'non-MPR': '#3498DB'}
scatter_colors = [colors_map.get(resp, 'gray') for resp in merged['response']]

ax1.scatter(merged['spp1_tam_frac'], merged['nk_dominant_ratio'],
            c=scatter_colors, s=80, alpha=0.8, edgecolor='black', linewidth=0.8, zorder=3)

z = np.polyfit(merged['spp1_tam_frac'], merged['nk_dominant_ratio'], 1)
xs = np.linspace(merged['spp1_tam_frac'].min(), merged['spp1_tam_frac'].max(), 100)
ax1.plot(xs, np.polyval(z, xs), 'k--', lw=2, alpha=0.7)

# 置信带
residuals = merged['nk_dominant_ratio'] - np.polyval(z, merged['spp1_tam_frac'])
sse = np.sum(residuals**2)
n_pts = len(merged)
std_err = np.sqrt(sse / (n_pts - 2))
x_mean = merged['spp1_tam_frac'].mean()
ssx = np.sum((merged['spp1_tam_frac'] - x_mean)**2)
conf_int = 1.96 * std_err * np.sqrt(1/n_pts + (xs - x_mean)**2 / ssx)
y_pred = np.polyval(z, xs)
ax1.fill_between(xs, y_pred - conf_int, y_pred + conf_int, color='gray', alpha=0.15)

ax1.set_xlabel('SPP1+ TAM fraction of myeloid cells', fontsize=11)
ax1.set_ylabel('NK-Locked ratio', fontsize=11)
ax1.set_title('Myeloid gate vs Clonal fate locking', fontsize=13, fontweight='bold')

stat_text = f'Pearson r = {r_pearson:.3f}\nSpearman r = {rs_spearman:.3f}\np = {rsp_spearman:.2e}'
ax1.text(0.05, 0.95, stat_text, transform=ax1.transAxes, va='top', ha='left',
         fontsize=10, fontweight='bold', color='#2C3E50',
         bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='lightgray', alpha=0.9))

legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#E74C3C', markersize=10, label='pCR', markeredgecolor='black'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#3498DB', markersize=10, label='non-MPR', markeredgecolor='black'),
]
ax1.legend(handles=legend_elements, loc='upper right', frameon=False, fontsize=9)
ax1.text(-0.15, 1.05, 'A', transform=ax1.transAxes, fontsize=20, fontweight='bold')

# ---- Panel B: 小提琴图 High vs Low SPP1 ----
ax2 = fig.add_subplot(gs[0, 1])

genes_plot_b = [g for g in ['KLF2','TCF7','LEF1','IL7R','SELL'] if g in stem_present]
n_genes_b = len(genes_plot_b)

violin_data_high = []
violin_data_low = []
for g in genes_plot_b:
    violin_data_high.append(gene_df[gene_df.group=='High SPP1+ TAM'][g].values)
    violin_data_low.append(gene_df[gene_df.group=='Low SPP1+ TAM'][g].values)

x_pos = np.arange(n_genes_b)
width = 0.35

for i in range(n_genes_b):
    vp_low = ax2.violinplot([violin_data_low[i]], positions=[x_pos[i] - width/2],
                            widths=width*0.9, showmedians=True, showextrema=False)
    for pc in vp_low['bodies']:
        pc.set_facecolor('#3498DB'); pc.set_alpha(0.7)
    vp_low['cmedians'].set_color('#1A5276'); vp_low['cmedians'].set_linewidth(2)

    vp_high = ax2.violinplot([violin_data_high[i]], positions=[x_pos[i] + width/2],
                             widths=width*0.9, showmedians=True, showextrema=False)
    for pc in vp_high['bodies']:
        pc.set_facecolor('#E74C3C'); pc.set_alpha(0.7)
    vp_high['cmedians'].set_color('#922B21'); vp_high['cmedians'].set_linewidth(2)

ax2.set_xticks(x_pos)
ax2.set_xticklabels(genes_plot_b, fontsize=11, fontweight='bold')
ax2.set_ylabel('Mean log-expression in CD8+ T', fontsize=11)
ax2.set_title('Stemness genes: High vs Low SPP1+ TAM', fontsize=13, fontweight='bold')

legend_elems = [
    Patch(facecolor='#3498DB', alpha=0.7, label='Low SPP1+ TAM'),
    Patch(facecolor='#E74C3C', alpha=0.7, label='High SPP1+ TAM'),
]
ax2.legend(handles=legend_elems, loc='upper right', frameon=False, fontsize=9)

# KLF2 显著性标注（动态获取 p 值）
if 'KLF2' in genes_plot_b:
    klf2_idx = genes_plot_b.index('KLF2')
    klf2_row = gene_stats_df[gene_stats_df.gene == 'KLF2']
    if len(klf2_row) > 0:
        klf2_p = klf2_row['FDR'].values[0]
        sig_text = f'FDR p = {klf2_p:.2e}' if klf2_p >= 0.05 else f'FDR p = {klf2_p:.2e} **'
        y_max = max(max(violin_data_high[klf2_idx]), max(violin_data_low[klf2_idx]))
        ax2.plot([klf2_idx - width/2, klf2_idx + width/2], [y_max * 1.05, y_max * 1.05], 'k-', lw=1.2)
        ax2.text(klf2_idx, y_max * 1.08, sig_text, ha='center', fontsize=8, fontweight='bold', color='#C0392B')

ax2.text(-0.15, 1.05, 'B', transform=ax2.transAxes, fontsize=20, fontweight='bold')

# ---- Panel C: NK-like vs Tex 干性基因表达对比 ----
ax3 = fig.add_subplot(gs[1, :])

n_genes_c = len(panel_c_genes)
x_pos_c = np.arange(n_genes_c)
width_c = 0.35

nk_violin_data = []
tex_violin_data = []
for g in panel_c_genes:
    gi = list(t_sample.var_names).index(g)
    nk_violin_data.append(t_arr[nklike_mask, gi])
    tex_violin_data.append(t_arr[tex_mask, gi])

for i in range(n_genes_c):
    vp_nk = ax3.violinplot([nk_violin_data[i]], positions=[x_pos_c[i] - width_c/2],
                           widths=width_c*0.9, showmedians=True, showextrema=False)
    for pc in vp_nk['bodies']:
        pc.set_facecolor('#E74C3C'); pc.set_alpha(0.7)
    vp_nk['cmedians'].set_color('#922B21'); vp_nk['cmedians'].set_linewidth(2)

    vp_tex = ax3.violinplot([tex_violin_data[i]], positions=[x_pos_c[i] + width_c/2],
                            widths=width_c*0.9, showmedians=True, showextrema=False)
    for pc in vp_tex['bodies']:
        pc.set_facecolor('#F39C12'); pc.set_alpha(0.7)
    vp_tex['cmedians'].set_color('#7D6608'); vp_tex['cmedians'].set_linewidth(2)

# 每个基因标注 p 值（动态）
for i, g in enumerate(panel_c_genes):
    row = nk_vs_tex_df[nk_vs_tex_df.gene == g]
    if len(row) > 0:
        p_val = row['pvalue'].values[0]
        y_max = max(max(nk_violin_data[i]), max(tex_violin_data[i]))
        ax3.plot([x_pos_c[i] - width_c/2, x_pos_c[i] + width_c/2],
                 [y_max * 1.05, y_max * 1.05], 'k-', lw=1.2)
        p_str = f'p = {p_val:.2e}' if p_val >= 0.001 else f'p = {p_val:.1e}'
        ax3.text(x_pos_c[i], y_max * 1.08, p_str, ha='center', fontsize=8,
                 fontweight='bold', color='#2C3E50')

ax3.set_xticks(x_pos_c)
ax3.set_xticklabels(panel_c_genes, fontsize=12, fontweight='bold')
ax3.set_ylabel('Log-expression in single cells', fontsize=11)
ax3.set_title('Stemness gene expression: NK-like vs Exhausted CD8+ T', fontsize=13, fontweight='bold')

legend_elems_c = [
    Patch(facecolor='#E74C3C', alpha=0.7, label=f'NK-like (n={int(nklike_mask.sum())})'),
    Patch(facecolor='#F39C12', alpha=0.7, label=f'Exhausted (n={int(tex_mask.sum())})'),
]
ax3.legend(handles=legend_elems_c, loc='upper right', frameon=False, fontsize=10)

ax3.text(-0.02, 1.05, 'C', transform=ax3.transAxes, fontsize=20, fontweight='bold')

plt.tight_layout()
plt.savefig(os.path.join(RESULT, 'Fig4_mechanism.png'), dpi=300)
plt.savefig(os.path.join(RESULT, 'Fig4_mechanism.pdf'))
plt.close()
print('Fig4 saved to result/Fig4_mechanism.png', flush=True)

# ============= 保存统计 =============
stats = dict(
    pearson_r = r_pearson,
    pearson_p = rp_pearson,
    spearman_r = rs_spearman,
    spearman_p = rsp_spearman,
    n_patients = len(merged),
    n_high_spp1 = len(high_spp1),
    n_low_spp1 = len(low_spp1),
    n_nklike_cells = int(nklike_mask.sum()),
    n_tex_cells = int(tex_mask.sum()),
)
pd.Series(stats).to_csv(os.path.join(RESULT, 'fig4_stats.csv'))
gene_stats_df.to_csv(os.path.join(RESULT, 'fig4_gene_stats.csv'), index=False)
nk_vs_tex_df.to_csv(os.path.join(RESULT, 'fig4_nk_vs_tex_stats.csv'), index=False)

print('=== DONE Fig4 ===', flush=True)
