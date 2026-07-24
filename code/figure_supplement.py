#!/usr/bin/env python3
"""
Supplementary Figures
  Fig S1: 数据集概览与质控
    Panel A: 患者临床特征表（响应分布、样本量）
    Panel B: 细胞质量过滤标准（nGene, nUMI 分布）
    Panel C: TCR 测序覆盖度分析
  Fig S2: 敏感性分析汇总
    Panel A: NK-like ≥ 10 过滤后 OR 变化
    Panel B: 不同样本量亚组的 NK-Locked Ratio 分布
  Fig S3: 额外验证队列分析
    Panel A: GSE135222 详细结果
    Panel B: 各队列基因匹配率对比
"""
import os, sys
import numpy as np, pandas as pd, scipy.stats as ss
import statsmodels.formula.api as smf
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# 从公共配置加载响应标签（与 figure5.py 统一，避免重复定义）
from config import GSE135222_RESPONSE
import warnings
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _common import load_h5ad, normalize_log1p

DATA = os.path.dirname(HERE)
ADATA = os.path.join(DATA, 'adata')
RESULT = os.path.join(DATA, 'result')
os.makedirs(RESULT, exist_ok=True)

# ============= Fig S1: 数据集概览与质控 =============
print('=== Fig S1: Dataset overview & QC ===', flush=True)
adata = load_h5ad(os.path.join(ADATA, 'GSE243013_T_cells.h5ad'), backed='r')
for c in ['sub_cell_type','sampleID','pathological_response','clonotype']:
    adata.obs[c] = adata.obs[c].astype(str)

# Panel A: 患者临床特征
patient_summary = adata.obs.groupby('sampleID').agg(
    n_cells=('sub_cell_type', 'count'),
    response=('pathological_response', 'first'),
    n_clonotypes=('clonotype', 'nunique'),
).reset_index()
patient_summary = patient_summary[patient_summary['response'].isin(['pCR','non-MPR','unknown','MPR'])]

# Panel B: nGene/nUMI 分布（从 obs 中提取，如果有的话）
obs_cols = list(adata.obs.columns)
qc_cols = [c for c in obs_cols if any(k in c.lower() for k in ['ngene','n_gene','numi','n_umi','total_counts','n_counts'])]
print(f'QC columns found: {qc_cols}', flush=True)

# 如果没有 QC 列，用 n_genes_per_cell 和 total_counts
if not qc_cols:
    # 尝试从 var 中计算
    has_qc = False
    if hasattr(adata, 'X'):
        try:
            sample_idx = np.random.default_rng(42).choice(adata.n_obs, min(5000, adata.n_obs), replace=False)
            sub = adata[sample_idx].to_memory()
            subX = sub.X.toarray() if hasattr(sub.X, 'toarray') else np.asarray(sub.X)
            n_genes_per_cell = (subX > 0).sum(axis=1)
            total_counts = subX.sum(axis=1)
            has_qc = True
        except Exception:
            has_qc = False
else:
    has_qc = True
    for c in qc_cols:
        if 'gene' in c.lower():
            n_genes_per_cell = adata.obs[c].values
        elif 'count' in c.lower() or 'umi' in c.lower():
            total_counts = adata.obs[c].values

# Panel C: TCR 覆盖度
tcr_stats = adata.obs.groupby('sampleID').agg(
    n_cells=('clonotype', 'count'),
    n_clonotypes=('clonotype', 'nunique'),
    response=('pathological_response', 'first'),
).reset_index()
tcr_stats = tcr_stats[tcr_stats['response'].isin(['pCR','non-MPR'])]
tcr_stats['tcr_diversity'] = tcr_stats['n_clonotypes'] / tcr_stats['n_cells']

# 绘图
fig = plt.figure(figsize=(18, 6))
gs = fig.add_gridspec(1, 3, wspace=0.35)

# Panel A: 临床特征表
ax1 = fig.add_subplot(gs[0, 0])
ax1.axis('off')
ax1.set_title('Patient Clinical Characteristics', fontsize=13, fontweight='bold')

resp_counts = patient_summary['response'].value_counts()
table_data = [
    ['Response', 'N patients', 'Mean cells/patient', 'Mean clones/patient'],
]
for resp in ['pCR', 'non-MPR', 'MPR', 'unknown']:
    if resp in resp_counts.index:
        sub = patient_summary[patient_summary['response'] == resp]
        table_data.append([
            resp, str(len(sub)),
            f'{sub["n_cells"].mean():.0f}',
            f'{sub["n_clonotypes"].mean():.0f}',
        ])

table = ax1.table(cellText=table_data[1:], colLabels=table_data[0],
                  loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 1.8)
for i in range(len(table_data[0])):
    table[0, i].set_facecolor('#2C3E50')
    table[0, i].set_text_props(color='white', fontweight='bold')
for i in range(1, len(table_data)):
    if table_data[i][0] == 'pCR':
        for j in range(len(table_data[0])):
            table[i, j].set_facecolor('#FADBD8')
    elif table_data[i][0] == 'non-MPR':
        for j in range(len(table_data[0])):
            table[i, j].set_facecolor('#D6EAF8')

ax1.text(-0.05, 1.05, 'A', transform=ax1.transAxes, fontsize=18, fontweight='bold')

# Panel B: QC 分布
ax2 = fig.add_subplot(gs[0, 1])
if has_qc:
    ax2.hist(n_genes_per_cell, bins=50, color='#5DADE2', alpha=0.7, edgecolor='white')
    ax2.axvline(np.median(n_genes_per_cell), color='#E74C3C', ls='--', lw=2, label=f'Median={np.median(n_genes_per_cell):.0f}')
    ax2.set_xlabel('Genes per cell', fontsize=11)
    ax2.set_ylabel('Number of cells', fontsize=11)
    ax2.set_title('Cell quality: Gene detection', fontsize=13, fontweight='bold')
    ax2.legend(frameon=False, fontsize=9)
else:
    ax2.text(0.5, 0.5, 'QC data not available\nin h5ad file', ha='center', va='center', fontsize=12)
    ax2.set_title('Cell quality', fontsize=13, fontweight='bold')
ax2.text(-0.15, 1.05, 'B', transform=ax2.transAxes, fontsize=18, fontweight='bold')

# Panel C: TCR 覆盖度
ax3 = fig.add_subplot(gs[0, 2])
colors_tcr = {'pCR': '#E74C3C', 'non-MPR': '#3498DB'}
for resp, color in colors_tcr.items():
    sub = tcr_stats[tcr_stats['response'] == resp]
    ax3.scatter(sub['n_cells'], sub['n_clonotypes'], c=color, s=50, alpha=0.7,
                edgecolor='white', linewidth=0.8, label=f'{resp} (n={len(sub)})')

ax3.set_xlabel('Number of CD8+ T cells', fontsize=11)
ax3.set_ylabel('Number of clonotypes', fontsize=11)
ax3.set_title('TCR coverage per patient', fontsize=13, fontweight='bold')
ax3.legend(frameon=False, fontsize=9)
ax3.grid(alpha=0.3, ls=':')

# 标注多样性
tcr_pcr = tcr_stats[tcr_stats['response'] == 'pCR']['tcr_diversity']
tcr_nonmpr = tcr_stats[tcr_stats['response'] == 'non-MPR']['tcr_diversity']
mwu_tcr = ss.mannwhitneyu(tcr_pcr, tcr_nonmpr, alternative='two-sided').pvalue
ax3.text(0.05, 0.95, f'TCR diversity\npCR median={tcr_pcr.median():.3f}\nnon-MPR median={tcr_nonmpr.median():.3f}\nMWU p={mwu_tcr:.3f}',
         transform=ax3.transAxes, va='top', fontsize=8,
         bbox=dict(boxstyle='round', facecolor='white', edgecolor='lightgray', alpha=0.9))

ax3.text(-0.15, 1.05, 'C', transform=ax3.transAxes, fontsize=18, fontweight='bold')

plt.tight_layout()
plt.savefig(os.path.join(RESULT, 'FigS1_dataset_overview.png'), dpi=300)
plt.savefig(os.path.join(RESULT, 'FigS1_dataset_overview.pdf'))
plt.close()
print('FigS1 saved', flush=True)

# ============= Fig S2: 敏感性分析 =============
print('=== Fig S2: Sensitivity analysis ===', flush=True)
ppt = pd.read_csv(os.path.join(RESULT, 'per_patient_metrics.csv'))
ppt = ppt.dropna(subset=['nk_dominant_ratio'])

# Panel A: 不同 NK-like 阈值下的 OR
thresholds = [0, 5, 10, 20, 50]
or_results = []
for thr in thresholds:
    sub = ppt[ppt['nklike'] >= thr]
    if len(sub) < 10 or sub['resp_bin'].nunique() < 2:
        or_results.append({'threshold': thr, 'or': np.nan, 'p': np.nan, 'n': len(sub)})
        continue
    try:
        m = smf.logit('resp_bin ~ nk_dominant_ratio', data=sub).fit(disp=0)
        or_val = float(np.exp(m.params['nk_dominant_ratio']))
        p_val = float(m.pvalues['nk_dominant_ratio'])
        or_results.append({'threshold': thr, 'or': or_val, 'p': p_val, 'n': len(sub)})
    except Exception:
        or_results.append({'threshold': thr, 'or': np.nan, 'p': np.nan, 'n': len(sub)})

or_df = pd.DataFrame(or_results)
print(or_df.to_string(), flush=True)

# Panel B: 样本量 vs NK-Locked Ratio
fig = plt.figure(figsize=(14, 6))
gs2 = fig.add_gridspec(1, 2, wspace=0.3)

ax_s2a = fig.add_subplot(gs2[0, 0])
valid = or_df.dropna(subset=['or'])
ax_s2a.bar(range(len(valid)), valid['or'], color=['#E74C3C' if p < 0.05 else '#85929E' for p in valid['p']],
           alpha=0.8, edgecolor='white', linewidth=1.5)
ax_s2a.set_xticks(range(len(valid)))
ax_s2a.set_xticklabels([f'≥{int(t)} NK-like' for t in valid['threshold']], fontsize=9)
ax_s2a.set_ylabel('Odds Ratio for pCR', fontsize=11)
ax_s2a.set_title('Sensitivity: OR by NK-like cell threshold', fontsize=13, fontweight='bold')
ax_s2a.axhline(1, color='gray', ls='--', lw=1.5)
for i, row in valid.iterrows():
    ax_s2a.text(i, row['or'] + 0.5, f'OR={row["or"]:.1f}\np={row["p"]:.3f}\nn={int(row["n"])}',
                ha='center', fontsize=7, fontweight='bold')
ax_s2a.text(-0.15, 1.05, 'A', transform=ax_s2a.transAxes, fontsize=18, fontweight='bold')

# Panel B: NK-Locked Ratio 分布（按响应组，按样本量分层）
ax_s2b = fig.add_subplot(gs2[0, 1])
for resp, color in [('pCR', '#E74C3C'), ('non-MPR', '#3498DB')]:
    sub = ppt[ppt['response'] == resp]
    ax_s2b.scatter(sub['cd8t'], sub['nk_dominant_ratio'], c=color, s=50, alpha=0.7,
                   edgecolor='white', linewidth=0.8, label=f'{resp} (n={len(sub)})')

ax_s2b.set_xlabel('CD8+ T cells per patient', fontsize=11)
ax_s2b.set_ylabel('NK-Locked ratio', fontsize=11)
ax_s2b.set_title('NK-Locked ratio vs sample size', fontsize=13, fontweight='bold')
ax_s2b.legend(frameon=False, fontsize=9)
ax_s2b.grid(alpha=0.3, ls=':')

# 相关性
r_sample, p_sample = ss.spearmanr(ppt['cd8t'], ppt['nk_dominant_ratio'])
ax_s2b.text(0.05, 0.95, f'Spearman r={r_sample:.3f}\np={p_sample:.3f}',
            transform=ax_s2b.transAxes, va='top', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='white', edgecolor='lightgray', alpha=0.9))
ax_s2b.text(-0.15, 1.05, 'B', transform=ax_s2b.transAxes, fontsize=18, fontweight='bold')

plt.tight_layout()
plt.savefig(os.path.join(RESULT, 'FigS2_sensitivity.png'), dpi=300)
plt.savefig(os.path.join(RESULT, 'FigS2_sensitivity.pdf'))
plt.close()
print('FigS2 saved', flush=True)
or_df.to_csv(os.path.join(RESULT, 'figS2_sensitivity_stats.csv'), index=False)

# ============= Fig S3: 额外验证队列 =============
print('=== Fig S3: Additional validation ===', flush=True)
fig = plt.figure(figsize=(14, 6))
gs3 = fig.add_gridspec(1, 2, wspace=0.3)

# Panel A: GSE135222 详细结果
ax_s3a = fig.add_subplot(gs3[0, 0])
s2_data = pd.read_csv(os.path.join(RESULT, 'sig_GSE135222.csv'), index_col=0)

# GSE135222_RESPONSE 已从 config.py 统一加载（与 figure5.py 共享）

r_vals = []
nr_vals = []
# sig_GSE135222.csv 是 Series 保存格式：index=样本名，第一列为签名值
sig_col = s2_data.columns[0]
for sample in s2_data.index:
    resp = GSE135222_RESPONSE.get(sample)
    if resp == 'R':
        r_vals.append(float(s2_data.loc[sample, sig_col]))
    elif resp == 'NR':
        nr_vals.append(float(s2_data.loc[sample, sig_col]))

p_s3a = ss.mannwhitneyu(nr_vals, r_vals, alternative='two-sided').pvalue

bp = ax_s3a.boxplot([nr_vals, r_vals], patch_artist=True, showfliers=False, widths=0.55)
for patch, c in zip(bp['boxes'], ['#E74C3C', '#27AE60']):
    patch.set_facecolor(c); patch.set_alpha(0.75)
for med in bp['medians']:
    med.set_color('black'); med.set_linewidth(2)

np.random.seed(42)
ax_s3a.scatter(np.random.normal(1, 0.06, len(nr_vals)), nr_vals, s=50, c='#C0392B', alpha=0.8, edgecolor='white', linewidth=1, zorder=3)
ax_s3a.scatter(np.random.normal(2, 0.06, len(r_vals)), r_vals, s=50, c='#1E8449', alpha=0.8, edgecolor='white', linewidth=1, zorder=3)

ax_s3a.set_xticks([1, 2])
ax_s3a.set_xticklabels([f'NR\n(n={len(nr_vals)})', f'R\n(n={len(r_vals)})'], fontsize=11, fontweight='bold')
ax_s3a.set_ylabel('NK-like signature score', fontsize=11)
ax_s3a.set_title('GSE135222 (NSCLC, PFS-based)\nDetailed Results', fontsize=13, fontweight='bold')

y_max = max(max(nr_vals), max(r_vals))
y_range = max(nr_vals + r_vals) - min(nr_vals + r_vals)
ax_s3a.plot([1, 2], [y_max + y_range*0.1, y_max + y_range*0.1], 'k-', lw=1.5)
ax_s3a.text(1.5, y_max + y_range*0.12, f'p = {p_s3a:.4f}', ha='center', fontsize=10, fontweight='bold')
ax_s3a.text(0.5, 0.02, 'Note: R defined as PFS=0 (no progression)\nR sample size limited (n=6)',
            transform=ax_s3a.transAxes, fontsize=8, style='italic', color='gray', va='bottom')
ax_s3a.text(-0.15, 1.05, 'A', transform=ax_s3a.transAxes, fontsize=18, fontweight='bold')

# Panel B: 各队列基因匹配率（从 fig5_gene_match.csv 动态读取，避免硬编码）
ax_s3b = fig.add_subplot(gs3[0, 1])
gene_match = pd.read_csv(os.path.join(RESULT, 'fig5_gene_match.csv'))
datasets = gene_match['dataset'].tolist()
matched = gene_match['genes_matched'].tolist()
total = gene_match['genes_total'].tolist()
dtypes = gene_match['dtype'].tolist()
match_rates = gene_match['match_rate'].tolist()
colors_bar = ['#E74C3C', '#F39C12', '#3498DB', '#2E86C1'][:len(datasets)]
print(f'Gene match rates loaded from fig5_gene_match.csv:\n{gene_match.to_string()}', flush=True)

bars = ax_s3b.bar(range(len(datasets)), matched, color=colors_bar, alpha=0.8, edgecolor='white', linewidth=1.5)
ax_s3b.axhline(24, color='gray', ls='--', lw=1, alpha=0.5)
ax_s3b.set_xticks(range(len(datasets)))
ax_s3b.set_xticklabels(datasets, fontsize=10, rotation=15)
ax_s3b.set_ylabel('NK-like genes matched (out of 24)', fontsize=11)
ax_s3b.set_title('Gene matching rate across cohorts', fontsize=13, fontweight='bold')
ax_s3b.set_ylim(0, 28)

for i, (m, t) in enumerate(zip(matched, total)):
    pct = m / t * 100
    ax_s3b.text(i, m + 0.5, f'{m}/{t}\n({pct:.0f}%)', ha='center', fontsize=9, fontweight='bold')

# 标注数据类型（dtypes 已从 fig5_gene_match.csv 动态读取）
for i, dt in enumerate(dtypes):
    ax_s3b.text(i, -2, dt, ha='center', fontsize=8, color='gray', style='italic')

ax_s3b.text(-0.15, 1.05, 'B', transform=ax_s3b.transAxes, fontsize=18, fontweight='bold')

plt.tight_layout()
plt.savefig(os.path.join(RESULT, 'FigS3_additional_validation.png'), dpi=300)
plt.savefig(os.path.join(RESULT, 'FigS3_additional_validation.pdf'))
plt.close()
print('FigS3 saved', flush=True)

# 保存统计（直接使用从 fig5_gene_match.csv 读取的数据）
gene_match.to_csv(os.path.join(RESULT, 'figS3_gene_match_stats.csv'), index=False)

print('=== DONE Supplementary Figures ===', flush=True)
