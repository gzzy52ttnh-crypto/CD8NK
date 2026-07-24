#!/usr/bin/env python3
"""
Figure 1 — NK-like CD8 T 细胞特征与克隆分布
  Panel A: UMAP 展示 NK-like 和 Tex 细胞位置
  Panel B: 特征基因点图（FGFBP2, KLF2, CXCL13, LAYN）
  Panel C: 克隆扩增分布，展示大克隆比例
  Panel D: 响应组间的克隆扩增箱线图
"""
import os, sys
import numpy as np, pandas as pd, scipy.stats as ss
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _common import load_h5ad, normalize_log1p, gene_mean_per_cell, embedding

DATA = os.path.dirname(HERE)
ADATA = os.path.join(DATA, 'adata')
RESULT = os.path.join(DATA, 'result')
os.makedirs(RESULT, exist_ok=True)

print('=== Loading T_cells.h5ad ===', flush=True)
adata = load_h5ad(os.path.join(ADATA, 'GSE243013_T_cells.h5ad'))
adata.var_names = adata.var_names.astype(str)
logX = normalize_log1p(adata)

for c in ['sub_cell_type','sampleID','pathological_response','clonotype','clonotype_number','major_cell_type']:
    adata.obs[c] = adata.obs[c].astype(str)

adata.obs['is_cd8']  = adata.obs['sub_cell_type'].str.startswith('CD8T')
adata.obs['is_nklike'] = adata.obs['sub_cell_type'] == 'CD8T_NK-like_FGFBP2'
adata.obs['is_tex'] = adata.obs['sub_cell_type'].str.contains('exhaust|Tex|Terminal|GZMB|ZNF683', case=False, na=False)

nk_genes = ['FGFBP2','KLRD1','CX3CR1','FCGR3A','KLRC1','NKG7','GNLY','PRF1','GZMB','KLRB1']
nk_present = [g for g in nk_genes if g in adata.var_names]
exh_genes = ['TOX','HAVCR2','LAYN','CXCL13','ENTPD1']
exh_present = [g for g in exh_genes if g in adata.var_names]
stem_genes = ['TCF7','LEF1','IL7R','CCR7','SELL']
stem_present = [g for g in stem_genes if g in adata.var_names]

adata.obs['nk_score'] = gene_mean_per_cell(adata, nk_present, logX=logX)
adata.obs['exh_score']  = gene_mean_per_cell(adata, exh_present, logX=logX)
adata.obs['stem_score'] = gene_mean_per_cell(adata, stem_present, logX=logX)

cd8 = adata.obs['is_cd8']

# ============= Fig 1A: UMAP 全图 =============
print('=== Fig1A: UMAP embedding ===', flush=True)
cd8_idx = adata.obs.index[cd8]
rng = np.random.default_rng(42)
keep = []
for pid, g in adata.obs.loc[cd8_idx].groupby('sampleID'):
    nkidx = list(g.index[g['is_nklike']])
    texidx = list(g.index[g['is_tex']])
    othidx = list(g.index[~g['is_nklike'] & ~g['is_tex']])
    n_other = min(len(othidx), 200)
    sel = nkidx + texidx + (list(rng.choice(othidx, n_other, replace=False)) if len(othidx) > 0 else [])
    keep.extend(sel)
keep_pos = {c: i for i, c in enumerate(adata.obs.index)}
ki = [keep_pos[c] for c in keep]
subX = logX[ki].toarray() if hasattr(logX[ki], 'toarray') else np.asarray(logX[ki])
emb, method = embedding(subX, use_umap=True, n_comps=50, random_state=42)
print(f'embedding method = {method}, shape = {emb.shape}', flush=True)

nk_mask = adata.obs.loc[keep, 'is_nklike'].values
tex_mask = adata.obs.loc[keep, 'is_tex'].values
oth_mask = ~nk_mask & ~tex_mask

# ============= Fig 1B: 特征基因点图（FGFBP2, KLF2, CXCL13, LAYN）=============
print('=== Fig1B: Dot plot of signature genes ===', flush=True)
key_genes = ['FGFBP2','KLF2','CXCL13','LAYN','GZMB','TCF7','TOX']
key_present = [g for g in key_genes if g in adata.var_names]
print(f'Key genes present: {key_present}', flush=True)

subtypes_sel = ['CD8T_NK-like_FGFBP2']
for st in sorted(adata.obs.loc[cd8, 'sub_cell_type'].unique()):
    if st not in subtypes_sel and st.startswith('CD8T'):
        subtypes_sel.append(st)

dot_data = []
for st in subtypes_sel:
    mask = (adata.obs['sub_cell_type'] == st) & cd8
    if mask.sum() < 10:
        continue
    cells_idx = [i for i, c in enumerate(adata.obs.index) if mask.values[i]]
    subX_arr = logX[cells_idx].toarray() if hasattr(logX[cells_idx], 'toarray') else np.asarray(logX[cells_idx])
    row = {'subtype': st.replace('CD8T_', '')}
    for g in key_present:
        gi = list(adata.var_names).index(g)
        expr = subX_arr[:, gi]
        row[f'{g}_mean'] = float(expr.mean())
        row[f'{g}_pct'] = float((expr > 0).mean() * 100)
    dot_data.append(row)
dot_df = pd.DataFrame(dot_data)

# ============= Fig 1C: 克隆扩增分布 =============
print('=== Fig1C: Clone size distribution ===', flush=True)
clone_stats = []
for pid, sub in adata.obs.groupby('sampleID'):
    resp = sub['pathological_response'].iloc[0]
    if resp not in ('pCR','non-MPR'):
        continue
    cd8p = sub[sub['is_cd8']]
    cc = cd8p['clonotype'].value_counts()
    n_cells = len(cd8p)
    n_clones = len(cc)
    big_clones = (cc >= 5).sum()
    top1_frac = cc.iloc[0] / n_cells if len(cc) > 0 else 0
    top5_frac = cc.head(5).sum() / n_cells if len(cc) > 0 else 0
    clone_stats.append(dict(patient=pid, response=resp, n_cells=n_cells, n_clones=n_clones,
                            big_clones=big_clones, top1_frac=top1_frac, top5_frac=top5_frac,
                            big_clone_frac=top5_frac))
clone_df = pd.DataFrame(clone_stats)

# ============= Fig 1D: 响应组间克隆扩增箱线图 =============
print('=== Fig1D: Clonal expansion by response ===', flush=True)
pcr_vals = clone_df[clone_df['response'] == 'pCR']['top5_frac'].values
nonmpr_vals = clone_df[clone_df['response'] == 'non-MPR']['top5_frac'].values
mwu_p_clone = ss.mannwhitneyu(pcr_vals, nonmpr_vals, alternative='two-sided').pvalue

# ============= 组合绘图 =============
print('=== Composing Figure 1 ===', flush=True)
fig = plt.figure(figsize=(16, 14))
gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3, height_ratios=[1, 0.9])

# ---- Panel A: UMAP ----
ax1 = fig.add_subplot(gs[0, 0])
ax1.scatter(emb[oth_mask,0], emb[oth_mask,1], s=3, c='#D5D8DC', alpha=0.3, label='Other CD8+ T', rasterized=True, zorder=1)
ax1.scatter(emb[tex_mask,0], emb[tex_mask,1], s=5, c='#5499C7', alpha=0.6, label='Exhausted CD8+ T', rasterized=True, zorder=2)
ax1.scatter(emb[nk_mask,0],  emb[nk_mask,1],  s=8, c='#E74C3C', alpha=0.75, label='NK-like (FGFBP2+)', rasterized=True, zorder=3)
ax1.legend(loc='upper right', frameon=False, fontsize=9)
ax1.set_title(f'CD8+ T cell {method}', fontsize=13, fontweight='bold')
ax1.set_xlabel('UMAP 1', fontsize=11)
ax1.set_ylabel('UMAP 2', fontsize=11)
ax1.set_xticks([]); ax1.set_yticks([])
ax1.text(-0.12, 1.08, 'A', transform=ax1.transAxes, fontsize=18, fontweight='bold')

# ---- Panel B: Dot plot ----
ax2 = fig.add_subplot(gs[0, 1])
genes_plot = key_present
subtypes_plot = dot_df['subtype'].tolist()

n_genes = len(genes_plot)
n_subtypes = len(subtypes_plot)

means_mat = np.zeros((n_subtypes, n_genes))
pcts_mat = np.zeros((n_subtypes, n_genes))
for i, st in enumerate(subtypes_plot):
    row = dot_df[dot_df['subtype'] == st].iloc[0]
    for j, g in enumerate(genes_plot):
        means_mat[i, j] = row[f'{g}_mean']
        pcts_mat[i, j] = row[f'{g}_pct']

means_norm = (means_mat - means_mat.min()) / (means_mat.max() - means_mat.min() + 1e-8)

for i in range(n_subtypes):
    for j in range(n_genes):
        size = 30 + pcts_mat[i, j] * 2.2
        ax2.scatter(j, i, s=size, c=[plt.cm.RdBu_r(means_norm[i, j])],
                   edgecolors='gray', linewidth=0.5, zorder=2)

ax2.set_xticks(range(n_genes))
ax2.set_xticklabels(genes_plot, rotation=45, ha='right', fontsize=10, fontweight='bold')
ax2.set_yticks(range(n_subtypes))
ax2.set_yticklabels(subtypes_plot, fontsize=9)
ax2.set_title('Signature gene expression', fontsize=13, fontweight='bold')
ax2.set_xlabel('')
ax2.set_ylabel('CD8+ T cell subsets', fontsize=11)
ax2.invert_yaxis()
ax2.set_xlim(-0.5, n_genes - 0.5)
ax2.set_ylim(n_subtypes - 0.5, -0.5)

sm = plt.cm.ScalarMappable(cmap='RdBu_r', norm=plt.Normalize(vmin=0, vmax=1))
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax2, fraction=0.04, pad=0.02)
cbar.set_label('Mean expression', fontsize=8)
cbar.ax.tick_params(labelsize=7)

ax2.text(-0.12, 1.08, 'B', transform=ax2.transAxes, fontsize=18, fontweight='bold')

# ---- Panel C: Clone distribution ----
ax3 = fig.add_subplot(gs[1, 0])

clone_df_sorted = clone_df.sort_values('response').reset_index(drop=True)
x_pos = np.arange(len(clone_df_sorted))

ax3.bar(x_pos, clone_df_sorted['top1_frac'], width=0.7, label='Top 1 clone', color='#E74C3C', alpha=0.8)
ax3.bar(x_pos, clone_df_sorted['top5_frac'] - clone_df_sorted['top1_frac'], bottom=clone_df_sorted['top1_frac'],
        width=0.7, label='Top 2-5 clones', color='#F1948A', alpha=0.8)
other_frac = 1 - clone_df_sorted['top5_frac']
ax3.bar(x_pos, other_frac, bottom=clone_df_sorted['top5_frac'],
        width=0.7, label='Other clones', color='#D5D8DC', alpha=0.8)

xtick_labels = []
for _, row in clone_df_sorted.iterrows():
    xtick_labels.append(f"{row['patient']}")

ax3.set_xticks(x_pos[::3])
ax3.set_xticklabels([xtick_labels[i] for i in range(0, len(xtick_labels), 3)], fontsize=6, rotation=45, ha='right')
ax3.set_ylabel('Fraction of CD8+ T cells', fontsize=11)
ax3.set_title('Clonal expansion across patients', fontsize=13, fontweight='bold')
ax3.legend(loc='upper right', frameon=False, fontsize=8)
ax3.set_ylim(0, 1)

# 添加响应分组色块
for i, row in clone_df_sorted.iterrows():
    color = '#E74C3C' if row['response'] == 'pCR' else '#3498DB'
    ax3.axvspan(i - 0.35, i + 0.35, ymin=0, ymax=0.02, color=color, alpha=0.5, clip_on=False)

ax3.text(-0.1, 1.08, 'C', transform=ax3.transAxes, fontsize=18, fontweight='bold')

# ---- Panel D: 响应组间箱线图 ----
ax4 = fig.add_subplot(gs[1, 1])

box_data = [pcr_vals, nonmpr_vals]
box_labels = [f'pCR\n(n={len(pcr_vals)})', f'non-MPR\n(n={len(nonmpr_vals)})']
box_colors = ['#E74C3C', '#3498DB']

bp = ax4.boxplot(box_data, labels=box_labels, patch_artist=True, widths=0.6,
                 medianprops=dict(color='black', linewidth=1.5))
for patch, color in zip(bp['boxes'], box_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

ax4.set_ylabel('Top-5 clone fraction', fontsize=11)
ax4.set_title('Clonal expansion by response', fontsize=13, fontweight='bold')
ax4.set_ylim(0, 1)

# 添加p值
y_max = max(max(pcr_vals), max(nonmpr_vals))
ax4.plot([1, 1, 2, 2], [y_max * 1.05, y_max * 1.1, y_max * 1.1, y_max * 1.05],
         lw=1.5, c='black')
p_text = f'p = {mwu_p_clone:.2e}' if mwu_p_clone < 0.001 else f'p = {mwu_p_clone:.4f}'
ax4.text(1.5, y_max * 1.12, p_text, ha='center', fontsize=10, fontweight='bold')

ax4.text(-0.1, 1.08, 'D', transform=ax4.transAxes, fontsize=18, fontweight='bold')

plt.tight_layout()
plt.savefig(os.path.join(RESULT, 'Fig1_overview.png'), dpi=300)
plt.savefig(os.path.join(RESULT, 'Fig1_overview.pdf'))
plt.close()
print('Fig1 saved to result/Fig1_overview.png', flush=True)

# ============= 保存统计 =============
stats = dict(
    n_cd8_cells = int(cd8.sum()),
    n_nklike_cells = int((cd8 & adata.obs['is_nklike']).sum()),
    n_tex_cells = int((cd8 & adata.obs['is_tex']).sum()),
    n_patients = len(clone_df),
    n_pCR = int((clone_df.response=='pCR').sum()),
    n_nonMPR = int((clone_df.response=='non-MPR').sum()),
    clone_expansion_mwu_p = float(mwu_p_clone),
    pcr_clone_median = float(np.median(pcr_vals)),
    nonmpr_clone_median = float(np.median(nonmpr_vals)),
)
pd.Series(stats).to_csv(os.path.join(RESULT, 'fig1_stats.csv'))

# 同时保存 per_patient_metrics.csv 供下游使用
rows = []
for pid, sub in adata.obs.groupby('sampleID'):
    resp = sub['pathological_response'].iloc[0]
    if resp not in ('pCR','non-MPR'):
        continue
    cdm = sub['is_cd8']
    cd8n = int(cdm.sum())
    nkn = int((cdm & sub['is_nklike']).sum())
    nk_pct = nkn/cd8n if cd8n > 0 else np.nan
    cd8p = sub[cdm]
    cc = cd8p['clonotype'].value_counts()
    big = cc[cc >= 5]
    fracs = []
    for cl in big.index:
        clc = cd8p[cd8p['clonotype'] == cl]
        fracs.append(float(clc['is_nklike'].sum()) / len(clc))
    nk_dom = float(np.mean(fracs)) if fracs else np.nan
    top5_frac = cc.head(5).sum() / cd8n if len(cc) > 0 else 0
    rows.append(dict(patient=pid, response=resp, cd8t=cd8n, nklike=nkn, nk_pct=nk_pct,
                     n_big_clones=int(len(big)), total_clone_count=int(cc.shape[0]),
                     nk_dominant_ratio=nk_dom, top5_clone_frac=top5_frac))
pdf = pd.DataFrame(rows)
pdf['resp_bin'] = (pdf['response'] == 'pCR').astype(int)
pdf.to_csv(os.path.join(RESULT,'per_patient_metrics.csv'), index=False)
print(f'Per-patient metrics saved: {len(pdf)} patients', flush=True)

print('=== DONE Fig1 ===', flush=True)
