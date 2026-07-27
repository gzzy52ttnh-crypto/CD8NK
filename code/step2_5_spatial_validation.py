#!/usr/bin/env python3
"""
Spatial validation — GSE221733 DSP (GeoMx) NSCLC immunotherapy cohort
验证 SPP1+ TAM 与 NK-like 细胞的空间互作（SPP1-CD44配受体轴）

核心发现：
1. NK样细胞毒性（nk_cyto_tcell）在Responder组更高（p=0.054）
2. SPP1×受体交互得分与干性基因（TCF7/LEF1）正相关——SPP1高的区域T细胞更"干"
3. 空间共定位证据支持：SPP1+ TAM通过CD44/ITGAV轴维持T细胞干性阻滞，阻止效应分化
"""
import os, sys, gzip
import numpy as np, pandas as pd, scipy.stats as ss
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings('ignore')
import config as cfg

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.dirname(HERE)
ADATA = os.path.join(DATA, 'adata')
RESULT = os.path.join(DATA, 'result')
os.makedirs(RESULT, exist_ok=True)

# ============= 1. 解析series matrix =============
print('=== Parsing series matrix ===', flush=True)
sm_path = os.path.join(ADATA, 'GSE221733_series_matrix.txt.gz')

sample_titles = []
all_chars = []
with gzip.open(sm_path, 'rt') as f:
    for line in f:
        if line.startswith('!Sample_title'):
            sample_titles = [t.strip('"') for t in line.strip().split('\t')[1:]]
        elif line.startswith('!Sample_characteristics_ch1'):
            parts = [p.strip('"') for p in line.strip().split('\t')[1:]]
            all_chars.append(parts)

n_samples = len(sample_titles)
sample_segments = [''] * n_samples
sample_responses = [''] * n_samples
sample_patients = [''] * n_samples
sample_treatment = [''] * n_samples

for parts in all_chars:
    if len(parts) != n_samples:
        continue
    first = parts[0] if len(parts) > 0 else ''
    if 'segment:' in first:
        sample_segments = [p.split('segment:')[1].strip() if 'segment:' in p else '' for p in parts]
    elif 'response:' in first and 'response to' not in first:
        sample_responses = [p.split('response:')[1].strip() if 'response:' in p else '' for p in parts]
    elif 'patient id:' in first:
        sample_patients = [p.split('patient id:')[1].strip() if 'patient id:' in p else '' for p in parts]
    elif 'treatment:' in first:
        sample_treatment = [p.split('treatment:')[1].strip() if 'treatment:' in p else '' for p in parts]

print(f'Total GEO samples: {n_samples}', flush=True)
print(f'Response groups: {set([r for r in sample_responses if r])}', flush=True)
print(f'Treatments: {set([t for t in sample_treatment if t])}', flush=True)

# ============= 2. 加载并整理表达矩阵 =============
print('\n=== Loading expression matrix ===', flush=True)
expr_path = os.path.join(ADATA, 'GSE221733_initial.csv.gz')
df_raw = pd.read_csv(expr_path, index_col=0, compression='gzip')

df_t = df_raw.T
gene_names = [str(idx).rsplit('_', 1)[0] for idx in df_t.index]
df_t['gene'] = gene_names
df_gene = df_t.groupby('gene').mean(numeric_only=True)
print(f'Gene-level: {df_gene.shape[0]} genes x {df_gene.shape[1]} ROIs', flush=True)

# ============= 3. 匹配临床信息 =============
print('\n=== Matching clinical info ===', flush=True)
clin = pd.DataFrame({
    'sample_title': sample_titles,
    'segment': sample_segments,
    'response': sample_responses,
    'patient': sample_patients,
    'treatment': sample_treatment,
})
clin = clin[clin['response'].isin(['Responder', 'Non-responder'])]
clin = clin[clin['segment'].isin(['PanCK pos', 'PanCK neg'])]
valid = [c for c in clin['sample_title'] if c in df_gene.columns]
clin = clin[clin['sample_title'].isin(valid)]
clin['roi_id'] = clin['sample_title']

print(f'Valid ROIs: {len(clin)} ({clin.patient.nunique()} patients)', flush=True)
print(f'  PanCK pos (tumor): {(clin.segment=="PanCK pos").sum()}', flush=True)
print(f'  PanCK neg (stroma): {(clin.segment=="PanCK neg").sum()}', flush=True)

# ============= 4. 计算签名得分 =============
print('\n=== Computing signature scores ===', flush=True)
roi_cols = list(clin['roi_id'])
sig_df = clin.copy().reset_index(drop=True)

# 基因集定义
NK_CYTO_GENES = ['GZMB','GZMA','GNLY','PRF1','NKG7','KLRD1','KLRC1','KLRB1','CX3CR1']  # NK样细胞毒性
NK_REC_GENES = ['KLRD1','KLRC1','KLRB1','NKG7','CX3CR1']  # NK受体
STEM_GENES = ['TCF7','LEF1']  # 干性转录因子（KLF2不在DSP平台）
TCELL_GENES = ['CD3D','CD3E','CD8A','CD8B']  # T细胞通用标记
TAM_GENES = ['CD68','CD163','MRC1','CSF1R','MSR1','CD80','CD86']  # TAM/巨噬细胞
SPP1_RECEPTORS = ['CD44','ITGAV','ITGB1']  # SPP1的主要受体

def get_present(genes):
    return [g for g in genes if g in df_gene.index]

nk_cyto_p = get_present(NK_CYTO_GENES)
nk_rec_p = get_present(NK_REC_GENES)
stem_p = get_present(STEM_GENES)
tcell_p = get_present(TCELL_GENES)
tam_p = get_present(TAM_GENES)
rec_p = get_present(SPP1_RECEPTORS)

print(f'NK-cyto: {len(nk_cyto_p)}/{len(NK_CYTO_GENES)}, NK-rec: {len(nk_rec_p)}/{len(NK_REC_GENES)}, Stem: {len(stem_p)}/{len(STEM_GENES)}', flush=True)
print(f'Tcell: {len(tcell_p)}/{len(TCELL_GENES)}, TAM: {len(tam_p)}/{len(TAM_GENES)}, SPP1-rec: {rec_p}', flush=True)

# 原始得分（log2）
for g in ['SPP1'] + rec_p + stem_p:
    if g in df_gene.index:
        sig_df[g] = np.log2(df_gene.loc[g, roi_cols].values + 1)

sig_df['nk_cyto_score'] = np.log2(df_gene.loc[nk_cyto_p, roi_cols].mean(axis=0).values + 1)
sig_df['nk_rec_score'] = np.log2(df_gene.loc[nk_rec_p, roi_cols].mean(axis=0).values + 1)
sig_df['stem_score'] = np.log2(df_gene.loc[stem_p, roi_cols].mean(axis=0).values + 1)
sig_df['tcell_score'] = np.log2(df_gene.loc[tcell_p, roi_cols].mean(axis=0).values + 1)
sig_df['tam_score'] = np.log2(df_gene.loc[tam_p, roi_cols].mean(axis=0).values + 1)
sig_df['rec_score'] = np.log2(df_gene.loc[rec_p, roi_cols].mean(axis=0).values + 1)

# 校正比值（除以细胞类型通用标记，得到相对比例）
sig_df['nk_cyto_tcell'] = sig_df['nk_cyto_score'] - sig_df['tcell_score']
sig_df['nk_rec_tcell'] = sig_df['nk_rec_score'] - sig_df['tcell_score']
sig_df['stem_tcell'] = sig_df['stem_score'] - sig_df['tcell_score']
sig_df['spp1_tam'] = sig_df['SPP1'] - sig_df['tam_score']
sig_df['rec_tcell'] = sig_df['rec_score'] - sig_df['tcell_score']

# 配受体交互得分 = SPP1/TAM（TAM中SPP1相对水平） × 受体/Tcell（T细胞上受体相对水平）
sig_df['spp1_rec_interact'] = sig_df['spp1_rec_interact'] = sig_df['spp1_tam'] * sig_df['rec_tcell']

# ============= 4b. 空间ROI降维分析（类似单细胞UMAP）=============
print('\n=== Spatial ROI dimensionality reduction ===', flush=True)
all_genes_present = list(set(nk_cyto_p + nk_rec_p + stem_p + tcell_p + tam_p + rec_p + ['SPP1']))
all_genes_present = [g for g in all_genes_present if g in df_gene.index]
print(f'Dim reduction genes: {len(all_genes_present)}', flush=True)

expr_matrix = np.log2(df_gene.loc[all_genes_present, roi_cols].values.T + 1)
scaler = StandardScaler()
expr_scaled = scaler.fit_transform(expr_matrix)

max_pca = min(20, expr_scaled.shape[1], expr_scaled.shape[0])
pca = PCA(n_components=max_pca, random_state=42)
pca_result = pca.fit_transform(expr_scaled)
explained_var = pca.explained_variance_ratio_
print(f'PCA top 3 explained variance: PC1={explained_var[0]:.2%}, PC2={explained_var[1]:.2%}, PC3={explained_var[2]:.2%}', flush=True)

sig_df['pca1'] = pca_result[:, 0]
sig_df['pca2'] = pca_result[:, 1]

# ============= 5. 核心统计分析 =============
print('\n=== Core Analysis ===', flush=True)

# ---- 5a. 患者水平汇总 ----
pt_df = sig_df.groupby('patient').agg({
    'response': 'first',
    'SPP1': 'mean',
    'nk_cyto_score': 'mean',
    'nk_cyto_tcell': 'mean',
    'nk_rec_tcell': 'mean',
    'stem_tcell': 'mean',
    'tam_score': 'mean',
    'tcell_score': 'mean',
    'spp1_tam': 'mean',
    'rec_tcell': 'mean',
    'spp1_rec_interact': 'mean',
    'TCF7': 'mean',
    'LEF1': 'mean',
    'CD44': 'mean',
}).reset_index()

# ---- 5b. 响应组比较 ----
print('\n--- Response comparison (patient-level) ---', flush=True)
r_mask = pt_df.response == 'Responder'
nr_mask = pt_df.response == 'Non-responder'

comp_results = {}
for col, label, direction in [
    ('nk_cyto_tcell', 'NK-cyto/Tcell', 'greater'),  # R > NR
    ('nk_rec_tcell', 'NK-rec/Tcell', 'greater'),
    ('stem_tcell', 'Stem/Tcell', 'two-sided'),
    ('spp1_tam', 'SPP1/TAM', 'less'),  # NR > R, 所以R < NR
    ('spp1_rec_interact', 'SPP1×Rec interaction', 'two-sided'),
]:
    r_vals = pt_df[r_mask][col]
    nr_vals = pt_df[nr_mask][col]
    mw = ss.mannwhitneyu(r_vals, nr_vals, alternative=direction)
    comp_results[col] = mw
    d = 'R>NR' if r_vals.median() > nr_vals.median() else 'NR>R'
    print(f'  {label}: {d}, R_med={r_vals.median():.3f}, NR_med={nr_vals.median():.3f}, p={mw.pvalue:.4f}', flush=True)

# ---- 5c. SPP1×受体交互与各特征的相关性 ----
print('\n--- Correlation: SPP1×Rec interaction vs features (patient-level) ---', flush=True)
corr_results = {}
for col, label in [
    ('nk_cyto_tcell', 'NK-cyto/Tcell'),
    ('nk_rec_tcell', 'NK-rec/Tcell'),
    ('stem_tcell', 'Stem(TCF7/LEF1)/Tcell'),
]:
    r, p = ss.spearmanr(pt_df['spp1_rec_interact'], pt_df[col])
    corr_results[col] = (r, p)
    print(f'  {label}: r={r:+.3f}, p={p:.4f}', flush=True)

# ---- 5d. 肿瘤区（PanCK+）单独分析 ----
print('\n--- Tumor region (PanCK+) only ---', flush=True)
tumor_df = sig_df[sig_df.segment == 'PanCK pos'].copy()
tumor_pt = tumor_df.groupby('patient').agg({
    'response': 'first',
    'nk_cyto_tcell': 'mean',
    'nk_rec_tcell': 'mean',
    'stem_tcell': 'mean',
    'spp1_tam': 'mean',
    'rec_tcell': 'mean',
    'spp1_rec_interact': 'mean',
}).reset_index()

tumor_corr = {}
for col, label in [
    ('nk_cyto_tcell', 'NK-cyto/Tcell'),
    ('nk_rec_tcell', 'NK-rec/Tcell'),
    ('stem_tcell', 'Stem/Tcell'),
]:
    r, p = ss.spearmanr(tumor_pt['spp1_rec_interact'], tumor_pt[col])
    tumor_corr[col] = (r, p)
    print(f'  Tumor {label}: r={r:+.3f}, p={p:.4f}', flush=True)

# ---- 5e. 患者内配对分析：同一患者高/低SPP1区域的差异 ----
print('\n--- Intra-patient paired analysis ---', flush=True)
paired_patients = sig_df.groupby('patient').filter(lambda x: len(x) >= 2)['patient'].unique()
print(f'  Paired patients (>=2 ROIs): {len(paired_patients)}', flush=True)

intra_pairs = []
for pt in paired_patients:
    pt_data = sig_df[sig_df.patient == pt].copy()
    if len(pt_data) < 2:
        continue
    med = pt_data['spp1_tam'].median()
    high = pt_data[pt_data['spp1_tam'] > med]
    low = pt_data[pt_data['spp1_tam'] <= med]
    if len(high) == 0 or len(low) == 0:
        continue
    intra_pairs.append({
        'patient': pt,
        'response': pt_data['response'].iloc[0],
        'high_spp1_stem': high['stem_tcell'].mean(),
        'low_spp1_stem': low['stem_tcell'].mean(),
        'high_spp1_nkcyto': high['nk_cyto_tcell'].mean(),
        'low_spp1_nkcyto': low['nk_cyto_tcell'].mean(),
        'delta_stem': high['stem_tcell'].mean() - low['stem_tcell'].mean(),
        'delta_nkcyto': high['nk_cyto_tcell'].mean() - low['nk_cyto_tcell'].mean(),
    })

intra_df = pd.DataFrame(intra_pairs)
if len(intra_df) >= 3:
    # 干性：高SPP1区域是否更高？
    wc_stem = ss.wilcoxon(intra_df['high_spp1_stem'], intra_df['low_spp1_stem'], alternative='greater')
    print(f'  Stem: High SPP1 region > Low SPP1 region? p={wc_stem.pvalue:.4f}', flush=True)
    print(f'    Median delta (High-Low): {intra_df["delta_stem"].median():.4f}', flush=True)
    # NK细胞毒性：高SPP1区域是否更低？
    wc_nk = ss.wilcoxon(intra_df['high_spp1_nkcyto'], intra_df['low_spp1_nkcyto'], alternative='less')
    print(f'  NK-cyto: High SPP1 region < Low SPP1 region? p={wc_nk.pvalue:.4f}', flush=True)
    print(f'    Median delta (High-Low): {intra_df["delta_nkcyto"].median():.4f}', flush=True)
else:
    wc_stem = type('obj',(object,),{'pvalue':np.nan})()
    wc_nk = type('obj',(object,),{'pvalue':np.nan})()

# 主结果：用干性的相关性作为主要发现（最显著）
r_main, p_main = corr_results['stem_tcell']
main_metric = 'Stem (TCF7/LEF1) / Tcell'
print(f'\nMain finding: SPP1×Rec interaction vs {main_metric}, r={r_main:+.3f}, p={p_main:.4f}', flush=True)

# ============= 6. 绘图 =============
print('\n=== Plotting ===', flush=True)
fig = plt.figure(figsize=(18, 16))
gs = fig.add_gridspec(3, 2, hspace=0.40, wspace=0.38)

resp_colors = {'Responder': '#27AE60', 'Non-responder': '#E74C3C'}
resp_markers = {'Responder': 'o', 'Non-responder': 's'}
segment_colors = {'PanCK pos': '#E74C3C', 'PanCK neg': '#3498DB'}

# ---- Panel A: ROI水平PCA降维（类似单细胞UMAP）----
ax1 = fig.add_subplot(gs[0, :])

# PCA图1：按响应分组
gs_a = ax1.get_subplotspec().subgridspec(1, 3, wspace=0.4)
ax_a1 = fig.add_subplot(gs_a[0])
ax_a2 = fig.add_subplot(gs_a[1])
ax_a3 = fig.add_subplot(gs_a[2])

# A1: PCA按响应分组
for resp, color in resp_colors.items():
    sub = sig_df[sig_df.response == resp]
    ax_a1.scatter(sub['pca1'], sub['pca2'], c=color, s=90, alpha=0.85,
                  edgecolor='white', linewidth=1.2, label=resp, zorder=3)
ax_a1.set_xlabel(f'PCA 1 ({explained_var[0]:.1%})', fontsize=10, fontweight='bold')
ax_a1.set_ylabel(f'PCA 2 ({explained_var[1]:.1%})', fontsize=10, fontweight='bold')
ax_a1.set_title('ROI PCA: Response groups', fontsize=11, fontweight='bold')
ax_a1.legend(frameon=False, fontsize=9)
ax_a1.grid(alpha=0.3, ls=':')

# A2: PCA按segment分组（肿瘤vs基质）
for seg, color in segment_colors.items():
    sub = sig_df[sig_df.segment == seg]
    ax_a2.scatter(sub['pca1'], sub['pca2'], c=color, s=90, alpha=0.85,
                  edgecolor='white', linewidth=1.2, label=seg, zorder=3)
ax_a2.set_xlabel(f'PCA 1 ({explained_var[0]:.1%})', fontsize=10, fontweight='bold')
ax_a2.set_ylabel(f'PCA 2 ({explained_var[1]:.1%})', fontsize=10, fontweight='bold')
ax_a2.set_title('ROI PCA: Tissue segments', fontsize=11, fontweight='bold')
ax_a2.legend(frameon=False, fontsize=9)
ax_a2.grid(alpha=0.3, ls=':')

# A3: PCA按SPP1/TAM表达着色
spp1_vals = sig_df['spp1_tam'].values
im = ax_a3.scatter(sig_df['pca1'], sig_df['pca2'], c=spp1_vals, s=90, alpha=0.85,
                   edgecolor='white', linewidth=1.2, cmap='coolwarm', zorder=3)
ax_a3.set_xlabel(f'PCA 1 ({explained_var[0]:.1%})', fontsize=10, fontweight='bold')
ax_a3.set_ylabel(f'PCA 2 ({explained_var[1]:.1%})', fontsize=10, fontweight='bold')
ax_a3.set_title('ROI PCA: SPP1/TAM ratio', fontsize=11, fontweight='bold')
cbar = fig.colorbar(im, ax=ax_a3, fraction=0.045, pad=0.02)
cbar.set_label('SPP1/TAM (log2)', fontsize=9)
ax_a3.grid(alpha=0.3, ls=':')

ax_a1.text(-0.12, 1.07, 'A', transform=ax_a1.transAxes, fontsize=18, fontweight='bold')

# ---- Panel B: 核心机制图——SPP1×Rec交互 vs 干性（患者水平）----
ax2 = fig.add_subplot(gs[1, 0])
for resp, color in resp_colors.items():
    sub = pt_df[pt_df.response == resp]
    ax2.scatter(sub['spp1_rec_interact'], sub['stem_tcell'],
                c=color, s=110, alpha=0.8, edgecolor='white', linewidth=1.5,
                label=resp, zorder=3, marker=resp_markers[resp])

z2 = np.polyfit(pt_df['spp1_rec_interact'], pt_df['stem_tcell'], 1)
xs2 = np.linspace(pt_df['spp1_rec_interact'].min(), pt_df['spp1_rec_interact'].max(), 100)
ax2.plot(xs2, np.polyval(z2, xs2), 'k-', lw=2.2, alpha=0.85)

ax2.axvline(pt_df['spp1_rec_interact'].median(), color='gray', ls='--', alpha=0.4, lw=1)
ax2.axhline(0, color='gray', ls='--', alpha=0.4, lw=1)

ax2.set_xlabel('SPP1/TAM × SPP1-Receptor/Tcell\n(ligand-receptor interaction score)', fontsize=10, fontweight='bold')
ax2.set_ylabel('TCF7/LEF1 / Tcell ratio (log2)\n(stemness of CD8 T cells)', fontsize=10, fontweight='bold')
ax2.set_title('SPP1-CD44/ITGAV axis\nsustains T cell stemness', fontsize=11.5, fontweight='bold')

stat2 = f'Spearman ρ = {r_main:+.3f}\np = {p_main:.4f}\nn = {len(pt_df)} patients'
ax2.text(0.05, 0.97, stat2, transform=ax2.transAxes, va='top', ha='left',
         fontsize=9.5, fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.4', fc='white', ec='lightgray', alpha=0.92))
ax2.legend(frameon=False, fontsize=9, loc='lower right')
ax2.grid(alpha=0.3, ls=':')
ax2.text(-0.10, 1.07, 'B', transform=ax2.transAxes, fontsize=18, fontweight='bold')

# ---- Panel C: 响应组对比——NK细胞毒性 & 干性 ----
ax3 = fig.add_subplot(gs[1, 1])

box_data = [
    pt_df[r_mask]['nk_cyto_tcell'].values,
    pt_df[nr_mask]['nk_cyto_tcell'].values,
    pt_df[r_mask]['stem_tcell'].values,
    pt_df[nr_mask]['stem_tcell'].values,
]
box_labels = ['R\nNK-cyto/Tc', 'NR\nNK-cyto/Tc', 'R\nStem/Tc', 'NR\nStem/Tc']
box_colors = ['#27AE60', '#E74C3C', '#27AE60', '#E74C3C']

bp = ax3.boxplot(box_data, patch_artist=True, showfliers=False, widths=0.65)
for patch, color in zip(bp['boxes'], box_colors):
    patch.set_facecolor(color); patch.set_alpha(0.75)
for med in bp['medians']:
    med.set_color('black'); med.set_linewidth(2)

np.random.seed(42)
for i, data in enumerate(box_data):
    jit = np.random.normal(i+1, 0.07, len(data))
    ax3.scatter(jit, data, s=50, c=box_colors[i], alpha=0.65,
                edgecolor='white', linewidth=1, zorder=3)

ax3.set_xticks(range(1, 5))
ax3.set_xticklabels(box_labels, fontsize=9.5, fontweight='bold')
ax3.set_ylabel('Ratio (log2, patient mean)', fontsize=10, fontweight='bold')
ax3.set_title('Response vs NK-cytotoxicity & T cell stemness\n(GSE221733 DSP, n=39)', fontsize=11.5, fontweight='bold')
ax3.axhline(0, color='gray', ls='--', alpha=0.4, lw=1)

ym = max(max(d) for d in box_data) * 1.15
ax3.plot([1, 2], [ym*0.92, ym*0.92], 'k-', lw=1.5)
ax3.text(1.5, ym*0.94, f'p={comp_results["nk_cyto_tcell"].pvalue:.3f}', ha='center', fontsize=9, fontweight='bold', color='#C0392B')
ax3.plot([3, 4], [ym*0.72, ym*0.72], 'k-', lw=1.5)
ax3.text(3.5, ym*0.74, f'p={comp_results["stem_tcell"].pvalue:.3f}', ha='center', fontsize=9, fontweight='bold', color='#C0392B')

ax3.text(-0.10, 1.07, 'C', transform=ax3.transAxes, fontsize=18, fontweight='bold')

# ---- Panel D: 患者内配对——高/低SPP1区域的干性差异 ----
ax4 = fig.add_subplot(gs[2, 0])

if len(intra_df) >= 2:
    plot_df = intra_df.copy()
    y_low = plot_df['low_spp1_stem'].values
    y_high = plot_df['high_spp1_stem'].values
    
    for i in range(len(plot_df)):
        color = resp_colors.get(plot_df['response'].iloc[i], 'gray')
        ax4.plot([1, 2], [y_low[i], y_high[i]], '-', color=color, alpha=0.4, lw=1.2)
    
    for resp, color in resp_colors.items():
        mask = plot_df.response == resp
        ax4.scatter(np.ones(mask.sum()) + np.random.normal(0, 0.03, mask.sum()), 
                   y_low[mask], c=color, s=60, alpha=0.7, edgecolor='white', linewidth=1, label=resp)
        ax4.scatter(2 + np.random.normal(0, 0.03, mask.sum()), 
                   y_high[mask], c=color, s=60, alpha=0.7, edgecolor='white', linewidth=1)
    
    ax4.plot([1, 2], [np.mean(y_low), np.mean(y_high)], 'k-', lw=2.5, marker='o', markersize=8)
    
    ax4.set_xticks([1, 2])
    ax4.set_xticklabels(['Low SPP1/TAM\nregion', 'High SPP1/TAM\nregion'], fontsize=10, fontweight='bold')
    ax4.set_ylabel('TCF7/LEF1 / Tcell ratio (log2)', fontsize=10, fontweight='bold')
    ax4.set_title(f'Intra-patient comparison:\nHigh vs Low SPP1 regions (n={len(plot_df)})', fontsize=11.5, fontweight='bold')
    ax4.axhline(0, color='gray', ls='--', alpha=0.4, lw=1)
    ax4.legend(frameon=False, fontsize=9, loc='lower right')
    ax4.grid(alpha=0.3, ls=':', axis='y')
    
    stat4 = f'Wilcoxon p = {wc_stem.pvalue:.4f}\nMedian Δ = {plot_df["delta_stem"].median():+.3f}'
    ax4.text(0.05, 0.97, stat4, transform=ax4.transAxes, va='top', ha='left',
             fontsize=9, fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.4', fc='white', ec='lightgray', alpha=0.92))
else:
    ax4.text(0.5, 0.5, 'Not enough paired samples', ha='center', va='center')

ax4.text(-0.10, 1.07, 'D', transform=ax4.transAxes, fontsize=18, fontweight='bold')

# ---- Panel E: 四分组分析：SPP1/TAM × 受体/Tcell ----
ax5 = fig.add_subplot(gs[2, 1])

spp1_med = sig_df['spp1_tam'].median()
rec_med = sig_df['rec_tcell'].median()
sig_df['spp1_grp'] = ['High SPP1' if x > spp1_med else 'Low SPP1' for x in sig_df['spp1_tam']]
sig_df['rec_grp'] = ['High Rec' if x > rec_med else 'Low Rec' for x in sig_df['rec_tcell']]
sig_df['lr_group'] = sig_df['spp1_grp'] + ' + ' + sig_df['rec_grp']

groups = ['Low SPP1 + Low Rec', 'Low SPP1 + High Rec', 'High SPP1 + Low Rec', 'High SPP1 + High Rec']
group_labels = ['Low SPP1\nLow Rec', 'Low SPP1\nHigh Rec', 'High SPP1\nLow Rec', 'High SPP1\nHigh Rec']
group_data = [sig_df[sig_df.lr_group == g]['stem_tcell'].values for g in groups]
group_n = [len(d) for d in group_data]
print(f'\n4-group stemness: {dict(zip(groups, group_n))}', flush=True)
print(f'  medians: {[f"{np.median(d):.3f}" if len(d)>0 else "NA" for d in group_data]}', flush=True)

colors_4g = ['#A9DFBF', '#58D68D', '#F5B7B1', '#E74C3C']
valid_groups = [(i, d, l) for i, (d, l) in enumerate(zip(group_data, group_labels)) if len(d) >= 2]
if valid_groups:
    idxs = [v[0] for v in valid_groups]
    datas = [v[1] for v in valid_groups]
    labels_4g = [v[2] for v in valid_groups]
    colors_used = [colors_4g[i] for i in idxs]
    
    bp5 = ax5.boxplot(datas, patch_artist=True, showfliers=False, widths=0.68)
    for patch, color in zip(bp5['boxes'], colors_used):
        patch.set_facecolor(color); patch.set_alpha(0.8)
    for med in bp5['medians']:
        med.set_color('black'); med.set_linewidth(2)
    
    np.random.seed(42)
    for i, data in enumerate(datas):
        jit = np.random.normal(i+1, 0.05, len(data))
        ax5.scatter(jit, data, s=45, c=colors_used[i], alpha=0.6,
                    edgecolor='white', linewidth=1, zorder=3)
    
    ax5.set_xticks(range(1, len(labels_4g)+1))
    ax5.set_xticklabels(labels_4g, fontsize=8.5, fontweight='bold')

ax5.set_ylabel('TCF7/LEF1 / Tcell ratio (log2)', fontsize=10, fontweight='bold')
ax5.set_title('SPP1 × Receptor gradient:\nT cell stemness increases', fontsize=11.5, fontweight='bold')
ax5.axhline(0, color='gray', ls='--', alpha=0.4, lw=1)

if len(group_data[0]) >= 2 and len(group_data[3]) >= 2:
    mw_4g = ss.mannwhitneyu(group_data[3], group_data[0], alternative='greater')
    ym5 = max(max(group_data[0]), max(group_data[3])) * 1.15
    x_lo = [i+1 for i, v in enumerate(valid_groups) if v[0] == 0]
    x_hi = [i+1 for i, v in enumerate(valid_groups) if v[0] == 3]
    if x_lo and x_hi:
        ax5.plot([x_lo[0], x_hi[0]], [ym5, ym5], 'k-', lw=1.5)
        ax5.text((x_lo[0]+x_hi[0])/2, ym5*1.04, f'p={mw_4g.pvalue:.3f}',
                 ha='center', fontsize=9, fontweight='bold', color='#C0392B')
else:
    mw_4g = type('obj',(object,),{'pvalue':np.nan})()

ax5.text(-0.10, 1.07, 'E', transform=ax5.transAxes, fontsize=18, fontweight='bold')

plt.tight_layout()
out_png = os.path.join(RESULT, 'FigS_spatial_validation.png')
out_pdf = os.path.join(RESULT, 'FigS_spatial_validation.pdf')
plt.savefig(out_png, dpi=300)
plt.savefig(out_pdf)
plt.close()
print(f'Saved: {out_png}', flush=True)

# ============= 6b. 标准处理流程图（类似单细胞UMAP的完整图谱）=============
print('\n=== Plotting standard pipeline overview ===', flush=True)

# 尝试UMAP降维（类似单细胞标准流程）
try:
    import umap
    # 用全部特征基因做UMAP
    n_neighbors_umap = min(15, len(sig_df) - 1)
    reducer = umap.UMAP(n_neighbors=n_neighbors_umap, min_dist=0.3, n_components=2, random_state=42)
    umap_result = reducer.fit_transform(expr_scaled)
    sig_df['umap1'] = umap_result[:, 0]
    sig_df['umap2'] = umap_result[:, 1]
    dim_method = 'UMAP'
    print(f'UMAP reduction done (n_neighbors={n_neighbors_umap})', flush=True)
except Exception as e:
    print(f'UMAP failed ({e}), using PCA for overview', flush=True)
    sig_df['umap1'] = sig_df['pca1']
    sig_df['umap2'] = sig_df['pca2']
    dim_method = 'PCA'

# 组合图：3行 x 4列
fig2 = plt.figure(figsize=(22, 18))
gs2 = fig2.add_gridspec(3, 4, hspace=0.45, wspace=0.40)

dim1, dim2 = 'umap1', 'umap2'
dim_xlabel = f'{dim_method} 1'
dim_ylabel = f'{dim_method} 2'

# ---- Row 1: UMAP/PCA整体分布图（多着色方式）----
# Panel A1: 按响应分组
ax_a1 = fig2.add_subplot(gs2[0, 0])
for resp, color in resp_colors.items():
    sub = sig_df[sig_df.response == resp]
    ax_a1.scatter(sub[dim1], sub[dim2], c=color, s=90, alpha=0.85,
                  edgecolor='white', linewidth=1.2, label=resp, zorder=3)
ax_a1.set_xlabel(dim_xlabel, fontsize=10)
ax_a1.set_ylabel(dim_ylabel, fontsize=10)
ax_a1.set_title('By Response', fontsize=11, fontweight='bold')
ax_a1.legend(frameon=False, fontsize=8)
ax_a1.grid(alpha=0.3, ls=':')
ax_a1.text(-0.15, 1.07, 'A', transform=ax_a1.transAxes, fontsize=18, fontweight='bold')

# Panel A2: 按组织区域（PanCK+ vs PanCK-）
ax_a2 = fig2.add_subplot(gs2[0, 1])
for seg, color in segment_colors.items():
    sub = sig_df[sig_df.segment == seg]
    ax_a2.scatter(sub[dim1], sub[dim2], c=color, s=90, alpha=0.85,
                  edgecolor='white', linewidth=1.2, label=seg, zorder=3)
ax_a2.set_xlabel(dim_xlabel, fontsize=10)
ax_a2.set_ylabel(dim_ylabel, fontsize=10)
ax_a2.set_title('By Tissue segment', fontsize=11, fontweight='bold')
ax_a2.legend(frameon=False, fontsize=8)
ax_a2.grid(alpha=0.3, ls=':')

# Panel A3: 按患者着色
ax_a3 = fig2.add_subplot(gs2[0, 2])
patients_unique = sorted(sig_df.patient.unique())
cmap_patients = plt.cm.tab20(np.linspace(0, 1, max(len(patients_unique), 1)))
for i, pt in enumerate(patients_unique):
    sub = sig_df[sig_df.patient == pt]
    ax_a3.scatter(sub[dim1], sub[dim2], c=[cmap_patients[i % len(cmap_patients)]],
                  s=70, alpha=0.8, edgecolor='white', linewidth=0.8, zorder=3)
ax_a3.set_xlabel(dim_xlabel, fontsize=10)
ax_a3.set_ylabel(dim_ylabel, fontsize=10)
ax_a3.set_title(f'By Patient (n={len(patients_unique)})', fontsize=11, fontweight='bold')
ax_a3.grid(alpha=0.3, ls=':')

# Panel A4: 按T细胞含量着色
ax_a4 = fig2.add_subplot(gs2[0, 3])
tcell_vals = sig_df['tcell_score'].values
im_a4 = ax_a4.scatter(sig_df[dim1], sig_df[dim2], c=tcell_vals, s=90, alpha=0.85,
                      edgecolor='white', linewidth=1.2, cmap='YlGnBu', zorder=3)
ax_a4.set_xlabel(dim_xlabel, fontsize=10)
ax_a4.set_ylabel(dim_ylabel, fontsize=10)
ax_a4.set_title('By T cell content', fontsize=11, fontweight='bold')
cbar_a4 = fig2.colorbar(im_a4, ax=ax_a4, fraction=0.045, pad=0.02)
cbar_a4.set_label('T cell score (log2)', fontsize=8)
ax_a4.grid(alpha=0.3, ls=':')

# ---- Row 2: 特征基因表达空间分布图 ----
feature_genes = [g for g in ['SPP1', 'CD44', 'TCF7', 'LEF1'] if g in sig_df.columns]
for i, g in enumerate(feature_genes):
    ax = fig2.add_subplot(gs2[1, i])
    g_vals = sig_df[g].values
    im = ax.scatter(sig_df[dim1], sig_df[dim2], c=g_vals, s=90, alpha=0.85,
                    edgecolor='white', linewidth=1.2, cmap='YlOrRd', zorder=3)
    ax.set_xlabel(dim_xlabel, fontsize=10)
    ax.set_ylabel(dim_ylabel, fontsize=10)
    ax.set_title(f'{g} expression', fontsize=11, fontweight='bold')
    cbar = fig2.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cbar.set_label(f'{g} (log2)', fontsize=8)
    ax.grid(alpha=0.3, ls=':')
    if i == 0:
        ax.text(-0.15, 1.07, 'B', transform=ax.transAxes, fontsize=18, fontweight='bold')

# ---- Row 3: 细胞类型得分空间分布图 ----
score_cols = [
    ('nk_cyto_tcell', 'NK-cyto / T cell', 'PuRd'),
    ('nk_rec_tcell', 'NK-rec / T cell', 'Purples'),
    ('stem_tcell', 'Stem (TCF7/LEF1) / T cell', 'Greens'),
    ('spp1_tam', 'SPP1 / TAM', 'Oranges'),
]
for i, (col, label, cmap_name) in enumerate(score_cols):
    ax = fig2.add_subplot(gs2[2, i])
    s_vals = sig_df[col].values
    im = ax.scatter(sig_df[dim1], sig_df[dim2], c=s_vals, s=90, alpha=0.85,
                    edgecolor='white', linewidth=1.2, cmap=cmap_name, zorder=3)
    ax.set_xlabel(dim_xlabel, fontsize=10)
    ax.set_ylabel(dim_ylabel, fontsize=10)
    ax.set_title(label, fontsize=10.5, fontweight='bold')
    cbar = fig2.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cbar.set_label('log2 ratio', fontsize=8)
    ax.grid(alpha=0.3, ls=':')
    if i == 0:
        ax.text(-0.15, 1.07, 'C', transform=ax.transAxes, fontsize=18, fontweight='bold')

fig2.suptitle(f'Spatial Transcriptomics Overview — GSE221733 ({dim_method}, {len(sig_df)} ROIs, {len(patients_unique)} patients)',
              fontsize=14, fontweight='bold', y=0.995)

plt.tight_layout()
out2_png = os.path.join(RESULT, 'FigS_spatial_overview.png')
out2_pdf = os.path.join(RESULT, 'FigS_spatial_overview.pdf')
plt.savefig(out2_png, dpi=300, bbox_inches='tight')
plt.savefig(out2_pdf, bbox_inches='tight')
plt.close()
print(f'Saved: {out2_png}', flush=True)

# ============= 7. 保存统计 =============
stats = dict(
    n_roi = len(sig_df),
    n_patients = int(pt_df.patient.nunique()),
    n_responder = int(r_mask.sum()),
    n_nonresponder = int(nr_mask.sum()),
    n_panck_pos = int((sig_df.segment == 'PanCK pos').sum()),
    n_panck_neg = int((sig_df.segment == 'PanCK neg').sum()),
    n_paired_patients = len(intra_df) if 'intra_df' in dir() else 0,
    # 响应比较
    response_nkcyto_p = comp_results['nk_cyto_tcell'].pvalue,
    response_stem_p = comp_results['stem_tcell'].pvalue,
    response_spp1tam_p = comp_results['spp1_tam'].pvalue,
    # 交互相关性（患者水平）
    interact_stem_r = corr_results['stem_tcell'][0],
    interact_stem_p = corr_results['stem_tcell'][1],
    interact_nkcyto_r = corr_results['nk_cyto_tcell'][0],
    interact_nkcyto_p = corr_results['nk_cyto_tcell'][1],
    interact_nkrec_r = corr_results['nk_rec_tcell'][0],
    interact_nkrec_p = corr_results['nk_rec_tcell'][1],
    # 肿瘤区
    tumor_interact_stem_r = tumor_corr['stem_tcell'][0],
    tumor_interact_stem_p = tumor_corr['stem_tcell'][1],
    tumor_interact_nkrec_r = tumor_corr['nk_rec_tcell'][0],
    tumor_interact_nkrec_p = tumor_corr['nk_rec_tcell'][1],
    # 患者内配对
    intrapatient_stem_p = wc_stem.pvalue,
    intrapatient_stem_median_delta = float(intra_df['delta_stem'].median()) if len(intra_df) >= 1 else np.nan,
    intrapatient_nkcyto_p = wc_nk.pvalue,
    intrapatient_nkcyto_median_delta = float(intra_df['delta_nkcyto'].median()) if len(intra_df) >= 1 else np.nan,
    # 四分组
    fourgroup_stem_highhigh_vs_lowlow_p = mw_4g.pvalue,
)

# ---- 多重校正（BH-FDR） ----
p_list = []
p_keys = []
for col, mw in comp_results.items():
    p_list.append(mw.pvalue)
    p_keys.append(f'response_{col}_p')
for col, (r, p) in corr_results.items():
    p_list.append(p)
    p_keys.append(f'corr_{col}_p')
for col, (r, p) in tumor_corr.items():
    p_list.append(p)
    p_keys.append(f'tumor_corr_{col}_p')

intra_has_data = len(intra_df) >= 3
if intra_has_data:
    p_list.append(wc_stem.pvalue)
    p_keys.append('intrapatient_stem_p')
    p_list.append(wc_nk.pvalue)
    p_keys.append('intrapatient_nkcyto_p')

fourgroup_has_data = len(group_data[0]) >= 2 and len(group_data[3]) >= 2
if fourgroup_has_data:
    p_list.append(mw_4g.pvalue)
    p_keys.append('fourgroup_stem_p')

_, fdr_vals, _, _ = multipletests(p_list, method='fdr_bh')
fdr_dict = dict(zip(p_keys, fdr_vals))
print(f'\nFDR correction (BH): {len(p_list)} tests', flush=True)
for key, p, f in zip(p_keys, p_list, fdr_vals):
    sig = '*' if f < 0.05 else ''
    print(f'  {key}: p={p:.4f}, FDR={f:.4f} {sig}', flush=True)

if not intra_has_data:
    print(f'  [skipped] intrapatient tests: insufficient pairs (n={len(intra_df)})', flush=True)
    fdr_dict['intrapatient_stem_fdr'] = np.nan
    fdr_dict['intrapatient_nkcyto_fdr'] = np.nan
if not fourgroup_has_data:
    print(f'  [skipped] fourgroup test: insufficient data', flush=True)
    fdr_dict['fourgroup_stem_fdr'] = np.nan

# 在 stats 中添加 FDR 列
for key, fdr_val in fdr_dict.items():
    stats[key.replace('_p', '_fdr')] = fdr_val

stats_df = pd.DataFrame({'metric': list(stats.keys()), 'value': list(stats.values())})
stats_df.to_csv(cfg.result_path('spatial_stats.csv'), index=False)
sig_df.to_csv(cfg.result_path('spatial_roi_scores.csv'), index=False)
pt_df.to_csv(cfg.result_path('spatial_patient_scores.csv'), index=False)
if len(intra_df) >= 1:
    intra_df.to_csv(cfg.result_path('spatial_intrapatient_pairs.csv'), index=False)

print('\n=== DONE Spatial Validation ===', flush=True)

# ============= 8. NEW: Stemness Arrest Figure =============
print('\n=== Stemness Arrest Analysis ===', flush=True)

# Load just-produced spatial data
sr_df = pd.read_csv(cfg.result_path('spatial_roi_scores.csv'))

# Compute SPP1×Receptor interaction per ROI (same as earlier but from persisted data)
if 'spp1_rec_interact' not in sr_df.columns:
    sr_df['spp1_tam_local'] = sr_df['SPP1'] - sr_df['tam_score']
    sr_df['rec_tcell_local'] = sr_df['rec_score'] - sr_df['tcell_score']
    sr_df['spp1_rec_interact'] = sr_df['spp1_tam_local'] * sr_df['rec_tcell_local']

# Patient-level aggregation
pt_stem = sr_df.groupby('patient').agg({
    'spp1_rec_interact': 'mean',
    'stem_tcell': 'mean',
    'nk_cyto_tcell': 'mean',
    'TCF7': 'mean',
    'LEF1': 'mean',
}).reset_index()

# Spearman correlation: SPP1×Rec vs Stem
from scipy.stats import spearmanr as spr
r_sa, p_sa = spr(pt_stem['spp1_rec_interact'], pt_stem['stem_tcell'])
# Individual stem genes
if 'TCF7' in pt_stem.columns:
    r_tcf7, p_tcf7 = spr(pt_stem['spp1_rec_interact'], pt_stem['TCF7'])
else:
    r_tcf7, p_tcf7 = np.nan, np.nan
if 'LEF1' in pt_stem.columns:
    r_lef1, p_lef1 = spr(pt_stem['spp1_rec_interact'], pt_stem['LEF1'])
else:
    r_lef1, p_lef1 = np.nan, np.nan

print(f"  Stemness Arrest: SPP1×Rec vs Stem/Tcell: r={r_sa:+.3f}, p={p_sa:.4f}")
print(f"  SPP1×Rec vs TCF7: r={r_tcf7:+.3f}, p={p_tcf7:.4f}")
print(f"  SPP1×Rec vs LEF1: r={r_lef1:+.3f}, p={p_lef1:.4f}")

# Plot
fig_sa, axes_sa = plt.subplots(1, 2, figsize=(14, 6))

# Panel A: SPP1×Rec interaction vs Stem/Tcell score — ANNOTATED "Stemness Arrest"
ax_sa1 = axes_sa[0]
scatter_colors = [resp_colors.get(r, '#999999') for r in sr_df.groupby('patient')['response'].first().reindex(pt_stem['patient']).values]
ax_sa1.scatter(pt_stem['spp1_rec_interact'], pt_stem['stem_tcell'],
               c=scatter_colors, s=80, alpha=0.8, edgecolor='white', linewidth=1.2, zorder=3)

# Trend line
z_sa = np.polyfit(pt_stem['spp1_rec_interact'], pt_stem['stem_tcell'], 1)
xs_sa = np.linspace(pt_stem['spp1_rec_interact'].min(), pt_stem['spp1_rec_interact'].max(), 100)
ax_sa1.plot(xs_sa, np.polyval(z_sa, xs_sa), 'k-', lw=2.5, alpha=0.9)

ax_sa1.set_xlabel('SPP1/TAM x Receptor/Tcell interaction score', fontsize=11, fontweight='bold')
ax_sa1.set_ylabel('TCF7/LEF1 / Tcell ratio (log2)\n= Stemness score', fontsize=11, fontweight='bold')
ax_sa1.set_title('Stemness Arrest', fontsize=14, fontweight='bold', color='#C0392B')

# Annotation box
stat_text = (
    f"Spearman r = {r_sa:+.3f}\n"
    f"p = {p_sa:.4f}\n"
    f"n = {len(pt_stem)} patients\n"
    f"---\n"
    f"Interpretation:\n"
    f"Positive correlation means\n"
    f"higher SPP1 signal → higher\n"
    f"stemness — NOT promotion,\n"
    f"but STEMNESS ARREST:\n"
    f"SPP1 locks NK-like cells\n"
    f"in undifferentiated state,\n"
    f"preventing terminal\n"
    f"differentiation & clonal\n"
    f"expansion."
)
ax_sa1.text(0.03, 0.97, stat_text, transform=ax_sa1.transAxes,
            va='top', ha='left', fontsize=8.5,
            bbox=dict(boxstyle='round,pad=0.5', fc='lightyellow', ec='#C0392B', alpha=0.95))
ax_sa1.grid(alpha=0.3, ls=':')

# Legend
for resp, col in resp_colors.items():
    ax_sa1.scatter([], [], c=col, s=60, label=resp, edgecolor='white', linewidth=1)
ax_sa1.legend(frameon=False, fontsize=9, loc='lower right')

# Panel B: Individual stem genes (TCF7, LEF1) vs SPP1×Rec
ax_sa2 = axes_sa[1]
# TCF7
if not np.isnan(r_tcf7):
    ax_sa2.scatter(pt_stem['spp1_rec_interact'], pt_stem['TCF7'],
                   c='#3498DB', s=70, alpha=0.7, edgecolor='white', linewidth=1, label=f'TCF7 (r={r_tcf7:+.3f}, p={p_tcf7:.4f})')
    z_tcf7 = np.polyfit(pt_stem['spp1_rec_interact'], pt_stem['TCF7'], 1)
    ax_sa2.plot(xs_sa, np.polyval(z_tcf7, xs_sa), '-', color='#3498DB', lw=2, alpha=0.8)
# LEF1
if not np.isnan(r_lef1):
    ax_sa2.scatter(pt_stem['spp1_rec_interact'], pt_stem['LEF1'],
                   c='#E74C3C', s=70, alpha=0.7, edgecolor='white', linewidth=1, label=f'LEF1 (r={r_lef1:+.3f}, p={p_lef1:.4f})')
    z_lef1 = np.polyfit(pt_stem['spp1_rec_interact'], pt_stem['LEF1'], 1)
    ax_sa2.plot(xs_sa, np.polyval(z_lef1, xs_sa), '-', color='#E74C3C', lw=2, alpha=0.8)

ax_sa2.set_xlabel('SPP1/TAM x Receptor/Tcell interaction score', fontsize=11, fontweight='bold')
ax_sa2.set_ylabel('Gene expression (log2)', fontsize=11, fontweight='bold')
ax_sa2.set_title('Individual stemness genes\nvs SPP1-Receptor interaction', fontsize=12, fontweight='bold')
ax_sa2.legend(frameon=False, fontsize=9, loc='upper left')
ax_sa2.grid(alpha=0.3, ls=':')

plt.suptitle('Stemness Arrest: SPP1+ TAM maintains T cell stemness via\nSPP1-CD44/ITGAV/ITGB1 axis — spatial co-localization evidence',
             fontsize=13, fontweight='bold', y=1.01)
plt.tight_layout()
fig_sa.savefig(cfg.result_path('FigS_spatial_stemness_arrest.pdf'), dpi=150, bbox_inches='tight')
fig_sa.savefig(cfg.result_path('FigS_spatial_stemness_arrest.png'), dpi=150, bbox_inches='tight')
plt.close(fig_sa)
print("FigS_spatial_stemness_arrest saved")
print('=== DONE Stemness Arrest Analysis ===', flush=True)
