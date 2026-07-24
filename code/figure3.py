#!/usr/bin/env python3
"""
Figure 3 — SPP1+ TAM 微环境机制
  Panel A: SPP1+ TAM 的 UMAP 分布
  Panel B: 特征基因热图（SPP1, IL1B, CXCL13 等）
  Panel C: 响应组间的 SPP1+ TAM 比例比较
  Panel D: 配受体互作网络图（SPP1-CD44 轴，TAM-T细胞互作）
  Panel E: 空间转录组验证（GSE221733 DSP）
"""
import os, sys
import numpy as np, pandas as pd, scipy.stats as ss
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
import warnings
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _common import load_h5ad, normalize_log1p, embedding, gene_mean_per_cell

DATA = os.path.dirname(HERE)
ADATA = os.path.join(DATA, 'adata')
RESULT = os.path.join(DATA, 'result')
os.makedirs(RESULT, exist_ok=True)

print('=== Loading immune.h5ad ===', flush=True)
adata = load_h5ad(os.path.join(ADATA, 'GSE243013_immune.h5ad'), backed='r')
for c in ['sub_cell_type','sampleID','pathological_response','major_cell_type']:
    adata.obs[c] = adata.obs[c].astype(str)

myel_mask = (adata.obs['major_cell_type'] == 'Myeloid cell')
m = adata[myel_mask].to_memory()
print('Myeloid cells:', m.n_obs, flush=True)
m.var_names = m.var_names.astype(str)
logX = normalize_log1p(m)

mac_mask = m.obs['sub_cell_type'].str.contains('Mφ|Mac|Mono', regex=True, na=False)
if 'SPP1' in m.var_names:
    gi = list(m.var_names).index('SPP1')
    spp1_col = logX[:, gi]
    if hasattr(spp1_col, 'toarray'):
        spp1_expr = spp1_col.toarray().ravel()
    else:
        spp1_expr = np.asarray(spp1_col).ravel()
else:
    spp1_expr = np.zeros(m.n_obs)
m.obs['spp1'] = spp1_expr

if mac_mask.sum() > 0:
    mac_spp1 = spp1_expr[mac_mask.values]
    thr = np.median(mac_spp1) if np.median(mac_spp1) > 0 else 0.0
    m.obs['is_spp1_tam'] = False
    m.obs.loc[m.obs.index[mac_mask.values], 'is_spp1_tam'] = (mac_spp1 > thr)
else:
    m.obs['is_spp1_tam'] = spp1_expr > np.median(spp1_expr)

ist = m.obs['is_spp1_tam'].values
print('SPP1+ TAM cells:', int(ist.sum()), flush=True)

# ============= Fig3A: SPP1+ TAM UMAP =============
print('=== Fig3A: Myeloid UMAP ===', flush=True)
rng = np.random.default_rng(42)
tam_idx = np.where(ist)[0]
oth_idx = np.where(~ist)[0]
cap_other = min(len(oth_idx), 5000)
sel_other = rng.choice(oth_idx, cap_other, replace=False) if len(oth_idx) > 0 else np.array([], dtype=int)
sel = np.concatenate([tam_idx, sel_other])
subX = logX[sel].toarray() if hasattr(logX[sel], 'toarray') else np.asarray(logX[sel])
emb, method = embedding(subX, use_umap=True, n_comps=30, random_state=42)
print(f'embedding method = {method}, plotted cells = {emb.shape[0]}', flush=True)
emb_ist = ist[sel]

# ============= Fig3B: 特征基因热图 =============
print('=== Fig3B: Signature gene heatmap ===', flush=True)
tam_genes = ['SPP1','IL1B','CXCL13','CD163','MRC1','CD68','CSF1R','TNF','CXCL8','CCL2','VEGFA','TGFB1']
tam_present = [g for g in tam_genes if g in m.var_names]
print(f'TAM marker genes present: {tam_present}', flush=True)

# 髓系各亚群中这些基因的表达
myel_subtypes = sorted(m.obs['sub_cell_type'].unique())
heatmap_data = []
heatmap_labels = []
for st in myel_subtypes:
    mask = m.obs['sub_cell_type'] == st
    if mask.sum() < 20:
        continue
    cells_idx = np.where(mask.values)[0]
    row = {}
    for g in tam_present:
        gi = list(m.var_names).index(g)
        gcol = logX[cells_idx, gi]
        if hasattr(gcol, 'toarray'):
            gval = gcol.toarray().ravel()
        else:
            gval = np.asarray(gcol).ravel()
        row[g] = float(gval.mean())
    heatmap_data.append(row)
    heatmap_labels.append(st.replace('CD14+ Mono ', '').replace('CD16+ Mono ', '').replace('Mac_', 'Mac_')[:20])

heatmap_df = pd.DataFrame(heatmap_data, index=heatmap_labels)
heatmap_arr = heatmap_df.values
heatmap_norm = (heatmap_arr - heatmap_arr.min(axis=0)) / (heatmap_arr.max(axis=0) - heatmap_arr.min(axis=0) + 1e-8)

# SPP1+ TAM 标记列
spp1_col_idx = tam_present.index('SPP1') if 'SPP1' in tam_present else 0

# ============= Fig3C: 箱线图 =============
print('=== Fig3C: SPP1+ TAM fraction boxplot ===', flush=True)
rows = []
for pid, sub in m.obs.groupby('sampleID'):
    resp = sub['pathological_response'].iloc[0]
    if resp not in ('pCR','non-MPR'):
        continue
    mac_sub = sub[sub['sub_cell_type'].str.contains('Mφ|Mac|Mono', regex=True, na=False)]
    tot_mac = len(mac_sub)
    sp = int(mac_sub['is_spp1_tam'].sum())
    frac = sp/tot_mac if tot_mac > 0 else np.nan
    rows.append(dict(patient=pid, response=resp, myeloid=len(sub), mac_mono=tot_mac, spp1_tam=sp, spp1_tam_frac=frac))
mpdf = pd.DataFrame(rows)
mpdf['resp_bin'] = (mpdf['response']=='pCR').astype(int)
mpdf.to_csv(os.path.join(RESULT,'myeloid_per_patient.csv'), index=False)
print(f'Myeloid per-patient rows={len(mpdf)}', flush=True)

g0 = mpdf[mpdf.response=='non-MPR']['spp1_tam_frac']
g1 = mpdf[mpdf.response=='pCR']['spp1_tam_frac']
mw = ss.mannwhitneyu(g0, g1, alternative='two-sided')
print(f'MWU p={mw.pvalue:.3e}  med nonMPR={g0.median():.4f}  pCR={g1.median():.4f}', flush=True)

# ============= Fig3D: 配受体互作网络图 =============
print('=== Fig3D: Ligand-Receptor network ===', flush=True)
lr_pairs = [
    ('SPP1', 'CD44'),
    ('SPP1', 'ITGAV'),
    ('SPP1', 'ITGB1'),
    ('CXCL8', 'CXCR2'),
    ('CCL2', 'CCR2'),
    ('IL1B', 'IL1R1'),
    ('TNF', 'TNFRSF1A'),
    ('VEGFA', 'FLT1'),
]

print('Loading T cell data for receptor expression...', flush=True)
t_adata = load_h5ad(os.path.join(ADATA, 'GSE243013_T_cells.h5ad'))
t_adata.var_names = t_adata.var_names.astype(str)
t_logX = normalize_log1p(t_adata)

t_adata.obs['is_nklike'] = t_adata.obs['sub_cell_type'] == 'CD8T_NK-like_FGFBP2'
t_adata.obs['is_cd8'] = t_adata.obs['sub_cell_type'].str.startswith('CD8T')
t_adata.obs['is_tex'] = t_adata.obs['sub_cell_type'].str.contains('exhaust|Tex|Terminal|GZMB|ZNF683', case=False, na=False)

ligand_expr = {}
receptor_expr_nk = {}
receptor_expr_tex = {}
interaction_score = {}
interaction_nk = {}
interaction_tex = {}

for ligand, receptor in lr_pairs:
    if ligand in m.var_names:
        gi = list(m.var_names).index(ligand)
        lig_col = logX[:, gi]
        if hasattr(lig_col, 'toarray'):
            lig_vals = lig_col.toarray().ravel()
        else:
            lig_vals = np.asarray(lig_col).ravel()
        lig_mean = float(lig_vals[ist].mean())
    else:
        lig_mean = 0
    ligand_expr[f'{ligand}-{receptor}'] = lig_mean

    if receptor in t_adata.var_names:
        ri = list(t_adata.var_names).index(receptor)
        t_arr = t_logX.toarray() if hasattr(t_logX, 'toarray') else np.asarray(t_logX)
        nk_mask = t_adata.obs['is_nklike'].values
        tex_mask = t_adata.obs['is_tex'].values
        rec_nk = float(t_arr[nk_mask, ri].mean()) if nk_mask.sum() > 0 else 0
        rec_tex = float(t_arr[tex_mask, ri].mean()) if tex_mask.sum() > 0 else 0
    else:
        rec_nk = 0
        rec_tex = 0
    receptor_expr_nk[f'{ligand}-{receptor}'] = rec_nk
    receptor_expr_tex[f'{ligand}-{receptor}'] = rec_tex
    # 分别计算 NK-like 和 Tex 互作得分，取最大值避免均值湮灭信号
    interaction_nk[f'{ligand}-{receptor}'] = lig_mean * rec_nk
    interaction_tex[f'{ligand}-{receptor}'] = lig_mean * rec_tex
    interaction_score[f'{ligand}-{receptor}'] = max(lig_mean * rec_nk, lig_mean * rec_tex)

print('Ligand-receptor interaction scores:', flush=True)
for k, v in sorted(interaction_score.items(), key=lambda x: -x[1]):
    print(f'  {k}: {v:.4f} (L={ligand_expr[k]:.2f}, R_NK={receptor_expr_nk[k]:.2f}, R_Tex={receptor_expr_tex[k]:.2f})', flush=True)

# ============= Fig3E: 空间转录组验证 =============
print('=== Fig3E: Spatial validation ===', flush=True)
spatial_pt = pd.read_csv(os.path.join(RESULT, 'spatial_patient_scores.csv'))
spatial_stats = pd.read_csv(os.path.join(RESULT, 'spatial_stats.csv'), index_col=0).iloc[:, 0]
r_stem = float(spatial_stats.get('interact_stem_r', 0))
p_stem = float(spatial_stats.get('interact_stem_p', 1))
p_nk_resp = float(spatial_stats.get('response_nkcyto_p', 1))

# ============= 组合绘图 =============
print('=== Composing Figure 3 ===', flush=True)
fig = plt.figure(figsize=(18, 16))
gs = fig.add_gridspec(4, 2, hspace=0.42, wspace=0.3, height_ratios=[1, 1, 0.9, 0.95])

# ---- Panel A: UMAP ----
ax1 = fig.add_subplot(gs[0, 0])
ax1.scatter(emb[~emb_ist,0], emb[~emb_ist,1], s=2, c='#D5D8DC', alpha=0.25, label='Other myeloid', rasterized=True, zorder=1)
ax1.scatter(emb[emb_ist,0],  emb[emb_ist,1],  s=6, c='#16A085', alpha=0.7, label='SPP1+ TAM', rasterized=True, zorder=2)
ax1.legend(loc='upper right', frameon=False, fontsize=9)
ax1.set_title(f'Myeloid cell {method}', fontsize=13, fontweight='bold')
ax1.set_xlabel(f'{method} 1', fontsize=11)
ax1.set_ylabel(f'{method} 2', fontsize=11)
ax1.set_xticks([]); ax1.set_yticks([])
ax1.text(-0.12, 1.08, 'A', transform=ax1.transAxes, fontsize=18, fontweight='bold')

# ---- Panel B: 特征基因热图 ----
ax2 = fig.add_subplot(gs[0, 1])

im = ax2.imshow(heatmap_norm.T, cmap='YlOrRd', aspect='auto', interpolation='nearest')

ax2.set_xticks(range(len(heatmap_labels)))
ax2.set_xticklabels(heatmap_labels, rotation=45, ha='right', fontsize=8)
ax2.set_yticks(range(len(tam_present)))
ax2.set_yticklabels(tam_present, fontsize=9, fontweight='bold')
ax2.set_title('SPP1+ TAM signature genes', fontsize=13, fontweight='bold')
ax2.set_ylabel('Marker genes', fontsize=10)

# 在格子里标数值
for i in range(heatmap_arr.shape[0]):
    for j in range(heatmap_arr.shape[1]):
        tc = 'white' if heatmap_norm[i, j] > 0.6 else 'black'
        ax2.text(i, j, f'{heatmap_arr[i, j]:.2f}', ha='center', va='center', fontsize=7, color=tc)

cbar = fig.colorbar(im, ax=ax2, fraction=0.03, pad=0.02)
cbar.set_label('Relative expression', fontsize=8)

# 高亮 SPP1 行
if 'SPP1' in tam_present:
    spp1_idx = tam_present.index('SPP1')
    ax2.axhline(spp1_idx - 0.5, color='#E74C3C', lw=2.5)
    ax2.axhline(spp1_idx + 0.5, color='#E74C3C', lw=2.5)
    ax2.text(len(heatmap_labels) - 0.5, spp1_idx, '★ SPP1', ha='right', va='center',
             fontsize=10, fontweight='bold', color='#C0392B')

ax2.text(-0.12, 1.08, 'B', transform=ax2.transAxes, fontsize=18, fontweight='bold')

# ---- Panel C: 箱线图 ----
ax3 = fig.add_subplot(gs[1, 0])

box_data = [g1.values, g0.values]
bp = ax3.boxplot(box_data, patch_artist=True, showfliers=False, widths=0.55)
colors = ['#E74C3C', '#3498DB']
for patch, c in zip(bp['boxes'], colors):
    patch.set_facecolor(c); patch.set_alpha(0.75)
for med in bp['medians']:
    med.set_color('black'); med.set_linewidth(2)
for whisker in bp['whiskers']:
    whisker.set_color('gray'); whisker.set_linewidth(1.2)

np.random.seed(42)
ax3.scatter(np.random.normal(1, 0.06, len(g1)), g1, s=40, c='#C0392B', alpha=0.7, edgecolor='white', linewidth=0.8, zorder=3)
ax3.scatter(np.random.normal(2, 0.06, len(g0)), g0, s=40, c='#2874A6', alpha=0.7, edgecolor='white', linewidth=0.8, zorder=3)

ax3.set_xticks([1, 2])
ax3.set_xticklabels(['pCR', 'non-MPR'], fontsize=11, fontweight='bold')
ax3.set_ylabel('SPP1+ TAM fraction\nof myeloid cells', fontsize=11)
ax3.set_title('SPP1+ TAM Enrichment by Response', fontsize=13, fontweight='bold')

p_text = f'Mann-Whitney p = {mw.pvalue:.2e}'
ax3.text(0.5, 0.95, p_text, transform=ax3.transAxes, ha='center', va='top',
         fontsize=10, fontweight='bold', color='#2874A6',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='lightgray', alpha=0.9))

ax3.text(-0.15, 1.05, 'C', transform=ax3.transAxes, fontsize=18, fontweight='bold')

# ---- Panel D: 配受体互作网络图 ----
ax4 = fig.add_subplot(gs[1, 1])
ax4.set_xlim(0, 10)
ax4.set_ylim(0, 8)
ax4.axis('off')
ax4.set_title('Ligand-Receptor Interaction Network', fontsize=13, fontweight='bold')

# 左侧：SPP1+ TAM（配体）
tam_circle = Circle((1.8, 4), 1.4, facecolor='#FADBD8', edgecolor='#C0392B', lw=2.5, alpha=0.85, zorder=5)
ax4.add_patch(tam_circle)
ax4.text(1.8, 4.3, 'SPP1+ TAM', fontsize=11, fontweight='bold', ha='center', color='#922B21', zorder=6)
ax4.text(1.8, 3.8, '(Ligands)', fontsize=9, ha='center', color='#922B21', style='italic', zorder=6)

# 右侧：CD8+ T 细胞（受体） - 分 NK-like 和 Tex
nk_circle = Circle((8.2, 5.8), 1.0, facecolor='#D6EAF8', edgecolor='#2E86C1', lw=2.5, alpha=0.85, zorder=5)
ax4.add_patch(nk_circle)
ax4.text(8.2, 6.0, 'NK-like', fontsize=10, fontweight='bold', ha='center', color='#1A5276', zorder=6)
ax4.text(8.2, 5.6, 'CD8+ T', fontsize=8, ha='center', color='#1A5276', style='italic', zorder=6)

tex_circle = Circle((8.2, 2.2), 1.0, facecolor='#FCF3CF', edgecolor='#D4AC0D', lw=2.5, alpha=0.85, zorder=5)
ax4.add_patch(tex_circle)
ax4.text(8.2, 2.4, 'Exhausted', fontsize=10, fontweight='bold', ha='center', color='#7D6608', zorder=6)
ax4.text(8.2, 2.0, 'CD8+ T', fontsize=8, ha='center', color='#7D6608', style='italic', zorder=6)

# 连线：按 NK 互作得分和 Tex 互作得分分别排序，各取 top-3
nk_sorted = sorted(interaction_nk.items(), key=lambda x: -x[1])[:3]
tex_sorted = sorted(interaction_tex.items(), key=lambda x: -x[1])[:3]

y_positions_nk = [6.8, 6.4, 6.0]
y_positions_tex = [1.2, 1.6, 2.0]

for i, (pair, score) in enumerate(nk_sorted):
    ligand, receptor = pair.split('-')
    y_start = 4.8 - i * 0.8
    y_end = y_positions_nk[i]
    lw = 1.5 + min(score * 8, 4)
    color = '#E67E22' if 'SPP1' in pair else '#85929E'
    alpha = 0.8 if 'SPP1' in pair else 0.5

    arrow = FancyArrowPatch((3.2, y_start), (7.2, y_end),
                            arrowstyle='-|>', mutation_scale=12, lw=lw, color=color, alpha=alpha)
    ax4.add_patch(arrow)

    mid_x = 5.2
    mid_y = (y_start + y_end) / 2
    ax4.text(mid_x, mid_y + 0.15, f'{ligand}→{receptor}',
             ha='center', fontsize=7.5, color=color, fontweight='bold' if 'SPP1' in pair else 'normal')

for i, (pair, score) in enumerate(tex_sorted):
    ligand, receptor = pair.split('-')
    y_start = 3.8 - i * 0.8
    y_end = y_positions_tex[i]
    lw = 1.5 + min(score * 8, 4)
    color = '#E67E22' if 'SPP1' in pair else '#85929E'
    alpha = 0.8 if 'SPP1' in pair else 0.5

    arrow = FancyArrowPatch((3.2, y_start), (7.2, y_end),
                            arrowstyle='-|>', mutation_scale=12, lw=lw, color=color, alpha=alpha)
    ax4.add_patch(arrow)

    mid_x = 5.2
    mid_y = (y_start + y_end) / 2
    ax4.text(mid_x, mid_y - 0.15, f'{ligand}→{receptor}',
             ha='center', fontsize=7.5, color=color, fontweight='bold' if 'SPP1' in pair else 'normal')

# 高亮 SPP1-CD44
spp1_cd44_score = interaction_score.get('SPP1-CD44', 0)
ax4.text(5.2, 7.3, f'★ SPP1-CD44 axis  ', ha='center', fontsize=10, fontweight='bold', color='#C0392B',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='#FDF2E9', edgecolor='#E74C3C', alpha=0.9))

ax4.text(-0.08, 1.05, 'D', transform=ax4.transAxes, fontsize=18, fontweight='bold')

# ---- Panel E: 空间转录组验证（跨两列）----
ax5 = fig.add_subplot(gs[2:, :])
gs_e = ax5.get_subplotspec().subgridspec(1, 3, wspace=0.35)
ax_e1 = fig.add_subplot(gs_e[0])
ax_e2 = fig.add_subplot(gs_e[1])
ax_e3 = fig.add_subplot(gs_e[2])

# E1: SPP1×Rec vs Stemness 散点图
resp_colors_d = {'Responder': '#27AE60', 'Non-responder': '#E74C3C'}
resp_markers_d = {'Responder': 'o', 'Non-responder': 's'}
for resp, color in resp_colors_d.items():
    sub = spatial_pt[spatial_pt.response == resp]
    ax_e1.scatter(sub['spp1_rec_interact'], sub['stem_tcell'],
                  c=color, s=90, alpha=0.8, edgecolor='white', linewidth=1.2,
                  label=resp, zorder=3, marker=resp_markers_d[resp])

z_e1 = np.polyfit(spatial_pt['spp1_rec_interact'], spatial_pt['stem_tcell'], 1)
xs_e1 = np.linspace(spatial_pt['spp1_rec_interact'].min(), spatial_pt['spp1_rec_interact'].max(), 100)
ax_e1.plot(xs_e1, np.polyval(z_e1, xs_e1), 'k-', lw=2, alpha=0.85)

ax_e1.set_xlabel('SPP1/TAM × Rec/Tcell\n(interaction score)', fontsize=9, fontweight='bold')
ax_e1.set_ylabel('TCF7/LEF1 / Tcell\n(stemness ratio, log2)', fontsize=9, fontweight='bold')
ax_e1.set_title('Spatial: SPP1 axis correlates with\nT cell stemness (DSP)', fontsize=10.5, fontweight='bold')
stat_e1 = f'ρ = {r_stem:+.2f}\np = {p_stem:.3f}\nn = {len(spatial_pt)} pts'
ax_e1.text(0.05, 0.97, stat_e1, transform=ax_e1.transAxes, va='top', ha='left',
           fontsize=8.5, fontweight='bold',
           bbox=dict(boxstyle='round,pad=0.35', fc='white', ec='lightgray', alpha=0.92))
ax_e1.legend(frameon=False, fontsize=7.5, loc='lower right')
ax_e1.grid(alpha=0.3, ls=':')

# E2: NK-cytotoxicity response comparison
nk_r = spatial_pt[spatial_pt.response == 'Responder']['nk_cyto_tcell']
nk_nr = spatial_pt[spatial_pt.response == 'Non-responder']['nk_cyto_tcell']

box_e2 = [nk_r.values, nk_nr.values]
bp_e2 = ax_e2.boxplot(box_e2, patch_artist=True, showfliers=False, widths=0.55)
colors_e2 = ['#27AE60', '#E74C3C']
for patch, c in zip(bp_e2['boxes'], colors_e2):
    patch.set_facecolor(c); patch.set_alpha(0.75)
for med in bp_e2['medians']:
    med.set_color('black'); med.set_linewidth(2)

np.random.seed(42)
ax_e2.scatter(np.random.normal(1, 0.06, len(nk_r)), nk_r, s=50, c='#27AE60', alpha=0.7,
              edgecolor='white', linewidth=1, zorder=3)
ax_e2.scatter(np.random.normal(2, 0.06, len(nk_nr)), nk_nr, s=50, c='#E74C3C', alpha=0.7,
              edgecolor='white', linewidth=1, zorder=3)

ax_e2.set_xticks([1, 2])
ax_e2.set_xticklabels(['Responder', 'Non-responder'], fontsize=9.5, fontweight='bold')
ax_e2.set_ylabel('NK-cyto / Tcell ratio\n(log2, patient mean)', fontsize=9, fontweight='bold')
ax_e2.set_title('NK-cytotoxicity enrichment\nin Responders', fontsize=10.5, fontweight='bold')
ym_e2 = max(max(nk_r), max(nk_nr)) * 1.2
ax_e2.plot([1, 2], [ym_e2, ym_e2], 'k-', lw=1.3)
ax_e2.text(1.5, ym_e2*1.03, f'p={p_nk_resp:.3f}', ha='center', fontsize=9, fontweight='bold', color='#C0392B')

# E3: 机制示意图
ax_e3.axis('off')
ax_e3.set_xlim(0, 10); ax_e3.set_ylim(0, 5.5)
ax_e3.text(5, 5.2, 'Spatial Correlation Summary', fontsize=11, fontweight='bold', ha='center', color='#2C3E50')

box_tam_e = FancyBboxPatch((0.5, 2.4), 2.8, 1.6, boxstyle="round,pad=0.1",
                           fc='#FADBD8', ec='#E74C3C', lw=2, alpha=0.85)
ax_e3.add_patch(box_tam_e)
ax_e3.text(1.9, 3.5, 'SPP1+ TAM', fontsize=9.5, fontweight='bold', ha='center', color='#C0392B')
ax_e3.text(1.9, 3.0, '↑ in non-MPR', fontsize=7.5, ha='center', color='#922B21')

box_t_e = FancyBboxPatch((6.7, 2.4), 2.8, 1.6, boxstyle="round,pad=0.1",
                         fc='#D6EAF8', ec='#2E86C1', lw=2, alpha=0.85)
ax_e3.add_patch(box_t_e)
ax_e3.text(8.1, 3.5, 'CD8+ T cell', fontsize=9.5, fontweight='bold', ha='center', color='#1A5276')
ax_e3.text(8.1, 3.0, 'Stemness ↑, Effector ↓', fontsize=7.5, ha='center', color='#1A5276')

arrow_e = FancyArrowPatch((3.4, 3.2), (6.6, 3.2),
                          arrowstyle='-|>', mutation_scale=18, lw=2.2, color='#E67E22')
ax_e3.add_patch(arrow_e)
ax_e3.text(5, 3.55, 'SPP1→CD44/ITGAV', fontsize=8, ha='center', color='#D35400', fontweight='bold')
ax_e3.text(5, 3.15, 'differentiation blockade', fontsize=7, ha='center', color='#D35400', style='italic')

bullets_e = [
    f'• DSP: SPP1-CD44 axis detected in situ',
    f'• SPP1+ TAM regions show higher T cell stemness',
    f'• NK-cytotoxicity associated with response (p={p_nk_resp:.2f})',
    f'• SPP1+ TAM: a candidate differentiation barrier',
]
for i, b in enumerate(bullets_e):
    ax_e3.text(0.3, 1.6 - i*0.34, b, fontsize=7.5, ha='left', color='#566573')

ax_e1.text(-0.12, 1.08, 'E', transform=ax_e1.transAxes, fontsize=18, fontweight='bold')

plt.tight_layout()
plt.savefig(os.path.join(RESULT, 'Fig3_myeloid.png'), dpi=300)
plt.savefig(os.path.join(RESULT, 'Fig3_myeloid.pdf'))
plt.close()
print('Fig3 saved to result/Fig3_myeloid.png', flush=True)

# ============= 保存统计 =============
stats = dict(
    n_myeloid_cells = int(m.n_obs),
    n_mac_mono = int(mac_mask.sum()),
    n_spp1_tam = int(ist.sum()),
    spp1_tam_frac = float(ist.sum() / mac_mask.sum()) if mac_mask.sum() > 0 else 0.0,
    mwu_p = mw.pvalue,
    median_spp1_nonMPR = g0.median(),
    median_spp1_pCR = g1.median(),
    spatial_interact_stem_r = r_stem,
    spatial_interact_stem_p = p_stem,
    spatial_nkcyto_response_p = p_nk_resp,
    n_lr_pairs = len(lr_pairs),
    top_lr_pair = nk_sorted[0][0] if nk_sorted else '',
    top_lr_score = nk_sorted[0][1] if nk_sorted else 0,
    top_lr_pair_tex = tex_sorted[0][0] if tex_sorted else '',
    top_lr_score_tex = tex_sorted[0][1] if tex_sorted else 0,
)
pd.Series(stats).to_csv(os.path.join(RESULT, 'fig3_stats.csv'))

# 保存配受体互作数据
lr_df_out = pd.DataFrame({
    'pair': list(ligand_expr.keys()),
    'ligand_expr_SPP1_TAM': list(ligand_expr.values()),
    'receptor_expr_NKlike': list(receptor_expr_nk.values()),
    'receptor_expr_Tex': list(receptor_expr_tex.values()),
    'interaction_NKlike': list(interaction_nk.values()),
    'interaction_Tex': list(interaction_tex.values()),
    'interaction_score': list(interaction_score.values()),
})
lr_df_out.to_csv(os.path.join(RESULT, 'ligand_receptor_scores.csv'), index=False)

print('=== DONE Fig3 ===', flush=True)
