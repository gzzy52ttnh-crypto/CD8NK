"""
step6e_tcr_diversity.py
Step 6.13-6.15: TCR克隆深度分析
- Step 6.13: 克隆谱系共享分析
- Step 6.14: 克隆多样性指数比较
- Step 6.15: GSE179994 配对检验
"""
import config as cfg
import _common
import pandas as pd
import numpy as np
import scipy.sparse as sp
import anndata as ad
import os
import gzip
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import entropy as shannon_entropy
from statsmodels.stats.multitest import multipletests

print("[step6e_tcr_diversity] Starting...")
os.makedirs(cfg.RESULT_DIR, exist_ok=True)

# ============================================================
# 读取 TCR 数据（从已解压的目录）
# ============================================================
print("\n" + "="*70)
print("Step 0: 读取 TCR 数据")
print("="*70)

tcr_dir = os.path.join(cfg.ADATA_DIR, 'GSE241934_TCR')
import glob
tcr_files = glob.glob(os.path.join(tcr_dir, '*filtered_contig*.csv.gz'))

# 读取所有患者的 contig 数据
all_contigs = []
for fpath in sorted(tcr_files):
    fname = os.path.basename(fpath)
    parts = fname.split('_')
    if len(parts) >= 3:
        patient_id = parts[1]
        df = pd.read_csv(fpath, compression='gzip')
        df['patient_id'] = patient_id
        all_contigs.append(df)

df_contigs = pd.concat(all_contigs, ignore_index=True)
print(f"  Total contigs: {len(df_contigs)}")
print(f"  Patients: {df_contigs['patient_id'].nunique()}")

# 构建细胞-克隆型映射（每个细胞取主要克隆型，即raw_clonotype_id）
cell_clono = df_contigs[['barcode', 'raw_clonotype_id', 'patient_id']].drop_duplicates()
cell_clono = cell_clono.rename(columns={'raw_clonotype_id': 'clonotype_id'})
cell_clono['cellID'] = cell_clono['patient_id'] + '_' + cell_clono['barcode']

# 读取 h5ad 获取临床信息和细胞类型
h5ad_path = os.path.join(cfg.RESULT_DIR, 'GSE241934_all.h5ad')
adata = ad.read_h5ad(h5ad_path)

# 合并 TCR 与 scRNA-seq
tcr_cells = set(cell_clono['cellID'])
sc_cells = set(adata.obs_names)
common_cells = tcr_cells & sc_cells
print(f"  TCR cells: {len(tcr_cells)}, scRNA cells: {len(sc_cells)}, matched: {len(common_cells)}")

# 仅保留有 TCR 的细胞
adata_tcr = adata[list(common_cells)].copy()
tcr_info = cell_clono.set_index('cellID').loc[adata_tcr.obs_names]
adata_tcr.obs['clonotype_id'] = tcr_info['clonotype_id'].values
adata_tcr.obs['tcr_patient'] = tcr_info['patient_id'].values

print(f"  分析细胞数: {adata_tcr.shape[0]}")

# 患者水平临床信息
patient_meta = adata.obs.groupby('sampleID').agg({
    'cohort': 'first',
    'chemo_class': 'first',
    'response_binary': 'first',
}).reset_index()
patient_meta = patient_meta.rename(columns={'sampleID': 'patient_id'})

# ============================================================
# Step 6.13: 克隆谱系共享分析
# ============================================================
print("\n" + "="*70)
print("Step 6.13: 克隆谱系共享分析")
print("="*70)

# 计算每个克隆在多少患者中出现
clone_patient_counts = adata_tcr.obs.groupby('clonotype_id')['tcr_patient'].nunique()
n_total_clones = len(clone_patient_counts)
n_shared_clones = (clone_patient_counts > 1).sum()
print(f"  总克隆数: {n_total_clones}")
print(f"  跨患者共享克隆: {n_shared_clones} ({n_shared_clones/n_total_clones*100:.2f}%)")

# 克隆大小分布
clone_sizes = adata_tcr.obs.groupby('clonotype_id').size()
print(f"  克隆大小中位数: {clone_sizes.median():.0f}")
print(f"  最大克隆大小: {clone_sizes.max():.0f}")

# 按响应组分析克隆共享模式
# 方法：R组患者之间、NR组患者之间、跨组的共享克隆数
r_patients = patient_meta[(patient_meta['response_binary'] == 'R') & 
                          (patient_meta['chemo_class'] == 'Taxane')]['patient_id'].tolist()
nr_patients = patient_meta[(patient_meta['response_binary'] == 'NR') & 
                           (patient_meta['chemo_class'] == 'Taxane')]['patient_id'].tolist()

print(f"\n  Taxane方案: R组 {len(r_patients)} 例, NR组 {len(nr_patients)} 例")

# 提取Taxane队列
adata_taxane = adata_tcr[adata_tcr.obs['chemo_class'] == 'Taxane'].copy()
print(f"  Taxane队列细胞数: {adata_taxane.shape[0]}")

# 计算每对患者之间的共享克隆数
def calc_pairwise_sharing(adata_sub, patients):
    """计算患者对之间的Jaccard相似系数和共享克隆数"""
    # 构建患者-克隆矩阵
    patient_clones = {}
    for p in patients:
        mask = adata_sub.obs['tcr_patient'] == p
        clones = set(adata_sub.obs.loc[mask, 'clonotype_id'].unique())
        patient_clones[p] = clones
    
    results = []
    for i, p1 in enumerate(patients):
        for j, p2 in enumerate(patients):
            if j <= i:
                continue
            c1 = patient_clones[p1]
            c2 = patient_clones[p2]
            intersection = len(c1 & c2)
            union = len(c1 | c2)
            jaccard = intersection / union if union > 0 else 0
            results.append({
                'patient1': p1,
                'patient2': p2,
                'shared_clones': intersection,
                'union_clones': union,
                'jaccard': jaccard
            })
    return pd.DataFrame(results)

# 三组配对：R-R, NR-NR, R-NR
print("\n  计算克隆共享模式...")
if len(r_patients) >= 2 and len(nr_patients) >= 2:
    rr_sharing = calc_pairwise_sharing(adata_taxane, r_patients)
    nrnr_sharing = calc_pairwise_sharing(adata_taxane, nr_patients)
    rnr_sharing = calc_pairwise_sharing(adata_taxane, r_patients + nr_patients)
    
    # 过滤掉 R-R 和 NR-NR 配对，只保留跨组的
    rnr_sharing = rnr_sharing[
        rnr_sharing['patient1'].isin(r_patients) & rnr_sharing['patient2'].isin(nr_patients) |
        rnr_sharing['patient2'].isin(r_patients) & rnr_sharing['patient1'].isin(nr_patients)
    ].copy()
    
    print(f"\n  克隆共享比较（Taxane方案）:")
    print(f"  {'组对':<15} {'配对数':>8} {'共享克隆中位数':>14} {'Jaccard中位数':>14} {'MWU p':>10}")
    print("  " + "-"*70)
    
    # 统计检验
    rr_jaccard = rr_sharing['jaccard'].values
    nrnr_jaccard = nrnr_sharing['jaccard'].values
    rnr_jaccard = rnr_sharing['jaccard'].values
    
    stat_rr_nr, p_rr_nr = stats.mannwhitneyu(rr_jaccard, nrnr_jaccard, alternative='two-sided')
    stat_rr_rnr, p_rr_rnr = stats.mannwhitneyu(rr_jaccard, rnr_jaccard, alternative='two-sided')
    
    print(f"  {'R-R':<15} {len(rr_sharing):>8} {np.median(rr_sharing['shared_clones']):>14.0f} {np.median(rr_jaccard):>14.4f} {'-':>10}")
    print(f"  {'NR-NR':<15} {len(nrnr_sharing):>8} {np.median(nrnr_sharing['shared_clones']):>14.0f} {np.median(nrnr_jaccard):>14.4f} {p_rr_nr:>10.4f}")
    print(f"  {'R-NR':<15} {len(rnr_sharing):>8} {np.median(rnr_sharing['shared_clones']):>14.0f} {np.median(rnr_jaccard):>14.4f} {p_rr_rnr:>10.4f}")
    
    sharing_results = [
        {'group_pair': 'R-R', 'n_pairs': len(rr_sharing), 
         'median_shared': np.median(rr_sharing['shared_clones']),
         'median_jaccard': np.median(rr_jaccard)},
        {'group_pair': 'NR-NR', 'n_pairs': len(nrnr_sharing),
         'median_shared': np.median(nrnr_sharing['shared_clones']),
         'median_jaccard': np.median(nrnr_jaccard)},
        {'group_pair': 'R-NR', 'n_pairs': len(rnr_sharing),
         'median_shared': np.median(rnr_sharing['shared_clones']),
         'median_jaccard': np.median(rnr_jaccard)},
    ]
    df_sharing = pd.DataFrame(sharing_results)
    df_sharing.to_csv(cfg.result_path('GSE241934_clone_sharing.csv'), index=False)
    print(f"\n  ✅ 克隆共享分析结果已保存")
else:
    print(f"  患者数不足，跳过配对共享分析")
    df_sharing = pd.DataFrame()

# NK-dominant 克隆的共享模式
print(f"\n  NK-dominant 克隆的共享模式:")
# 获取克隆分类信息
# NK-like 24 基因签名统一引用 _common.NKLIKE_SIGNATURE
nk_genes = _common.NKLIKE_SIGNATURE

# 计算每个克隆的 NK-like 细胞比例
def is_nklike_cell(cell_type):
    return 'CD8T_Teff_FGFBP2' in str(cell_type)

def is_tex_cell(cell_type):
    return 'CD8T_Tex_CXCL13' in str(cell_type)

clone_stats = []
for clono, group in adata_taxane.obs.groupby('clonotype_id'):
    total = len(group)
    nk_count = group['cell.type'].apply(is_nklike_cell).sum()
    tex_count = group['cell.type'].apply(is_tex_cell).sum()
    nk_ratio = nk_count / total
    tex_ratio = tex_count / total
    
    _nk_thresh = cfg.NK_DOMINANT_RATIO_THRESHOLD
    if nk_ratio >= _nk_thresh:
        category = 'NK-dominant'
    elif tex_ratio >= _nk_thresh:
        category = 'Tex-dominant'
    elif nk_ratio + tex_ratio >= _nk_thresh:
        category = 'Mixed'
    else:
        category = 'Other'
    
    n_patients = group['tcr_patient'].nunique()
    clone_stats.append({
        'clonotype_id': clono,
        'clone_size': total,
        'nk_count': nk_count,
        'tex_count': tex_count,
        'nk_ratio': nk_ratio,
        'tex_ratio': tex_ratio,
        'category': category,
        'n_patients': n_patients,
        'is_shared': n_patients > 1
    })

df_clone_stats = pd.DataFrame(clone_stats)
_big_thresh = cfg.BIG_CLONE_THRESHOLD
print(f"  大克隆(>={_big_thresh}): {(df_clone_stats['clone_size']>=_big_thresh).sum()}")
print(f"  NK-dominant克隆: {(df_clone_stats['category']=='NK-dominant').sum()}")
print(f"  Tex-dominant克隆: {(df_clone_stats['category']=='Tex-dominant').sum()}")

# 各类型克隆的共享率
big_clones = df_clone_stats[df_clone_stats['clone_size'] >= _big_thresh]
sharing_by_cat = big_clones.groupby('category')['is_shared'].agg(['mean', 'count', 'sum'])
print(f"\n  各类型大克隆的共享率:")
print(sharing_by_cat.to_string())

# 保存克隆统计
df_clone_stats.to_csv(cfg.result_path('GSE241934_clone_stats_with_sharing.csv'), index=False)

# 可视化：克隆共享分析
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Panel A: R-R / NR-NR / R-NR Jaccard 箱线图
ax = axes[0]
if len(df_sharing) > 0 and len(rr_jaccard) > 0 and len(nrnr_jaccard) > 0 and len(rnr_jaccard) > 0:
    box_data = [rr_jaccard, nrnr_jaccard, rnr_jaccard]
    labels = ['R-R', 'NR-NR', 'R-NR']
    bp = ax.boxplot(box_data, labels=labels, patch_artist=True,
                    medianprops={'color': 'black'})
    colors = ['#E64B35', '#4DBBD5', '#999999']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    for i, data in enumerate(box_data):
        x = np.random.normal(i + 1, 0.05, size=len(data))
        ax.scatter(x, data, alpha=0.6, color='black', s=20, zorder=3)
    ax.set_ylabel('Jaccard similarity')
    ax.set_title('Clone sharing\n(by response group)')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
else:
    ax.text(0.5, 0.5, 'Insufficient samples', ha='center', va='center')
    ax.set_title('Clone sharing')

# Panel B: 各克隆类型的共享率柱状图
ax = axes[1]
cats = ['NK-dominant', 'Tex-dominant', 'Mixed', 'Other']
share_rates = []
counts = []
for cat in cats:
    cat_clones = big_clones[big_clones['category'] == cat]
    if len(cat_clones) > 0:
        share_rates.append(cat_clones['is_shared'].mean())
        counts.append(len(cat_clones))
    else:
        share_rates.append(0)
        counts.append(0)

bars = ax.bar(cats, share_rates, color=['#E64B35', '#4DBBD5', '#91D1C2', '#CCCCCC'], alpha=0.8)
for bar, count in zip(bars, counts):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
            f'n={count}', ha='center', va='bottom', fontsize=9)
ax.set_ylabel('Shared clone proportion')
ax.set_title('Clone sharing rate\nby clone category')
ax.set_ylim(0, max(share_rates) * 1.3 if max(share_rates) > 0 else 0.1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.setp(ax.get_xticklabels(), rotation=15, ha='right')

# Panel C: 克隆大小分布（共享vs私有）
ax = axes[2]
shared_sizes = big_clones[big_clones['is_shared']]['clone_size'].values
private_sizes = big_clones[~big_clones['is_shared']]['clone_size'].values
box_data = [shared_sizes, private_sizes]
labels = [f'Shared\n(n={len(shared_sizes)})', f'Private\n(n={len(private_sizes)})']
bp = ax.boxplot(box_data, labels=labels, patch_artist=True,
                medianprops={'color': 'black'})
bp['boxes'][0].set_facecolor('#E64B35')
bp['boxes'][1].set_facecolor('#4DBBD5')
for patch in bp['boxes']:
    patch.set_alpha(0.7)
ax.set_ylabel('Clone size (cells)')
ax.set_title('Clone size: shared vs private')
ax.set_yscale('log')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
fig_path = cfg.result_path('GSE241934_FigS10_clone_sharing.png')
plt.savefig(fig_path, dpi=150, bbox_inches='tight')
plt.savefig(fig_path.replace('.png', '.pdf'), bbox_inches='tight')
plt.close()
print(f"\n  ✅ 克隆共享分析图已保存: GSE241934_FigS10_clone_sharing.png")

# ============================================================
# Step 6.14: 克隆多样性指数比较
# ============================================================
print("\n" + "="*70)
print("Step 6.14: 克隆多样性指数比较")
print("="*70)

# 计算每个患者的克隆多样性指标
diversity_results = []
for patient_id in patient_meta['patient_id']:
    mask = adata_tcr.obs['tcr_patient'] == patient_id
    if mask.sum() == 0:
        continue
    
    clone_counts = adata_tcr.obs.loc[mask, 'clonotype_id'].value_counts().values
    total_cells = clone_counts.sum()
    n_clones = len(clone_counts)
    
    if n_clones < 2:
        shannon = 0.0
        simpson = 0.0
        clonality = 0.0
    else:
        p_i = clone_counts / total_cells
        # Shannon entropy
        shannon = float(shannon_entropy(p_i, base=np.e))
        # Simpson index (1 - sum(p_i^2))
        simpson = 1 - np.sum(p_i ** 2)
        # Clonality = 1 - (Shannon / log(N))
        clonality = 1 - (shannon / np.log(n_clones))
    
    # 大克隆比例（前10%克隆占总细胞的比例）
    top_n = max(1, int(n_clones * 0.1))
    top_clone_cells = clone_counts[:top_n].sum()
    top_clone_frac = top_clone_cells / total_cells
    
    # 临床信息
    meta_row = patient_meta[patient_meta['patient_id'] == patient_id].iloc[0]
    
    diversity_results.append({
        'patient_id': patient_id,
        'cohort': meta_row['cohort'],
        'chemo_class': meta_row['chemo_class'],
        'response': meta_row['response_binary'],
        'total_tcr_cells': total_cells,
        'n_clones': n_clones,
        'shannon_entropy': shannon,
        'simpson_index': simpson,
        'clonality': clonality,
        'top10_clone_fraction': top_clone_frac,
        'clonality_1p': sum(clone_counts >= 10) / n_clones  # 大克隆比例（>=10细胞）
    })

df_diversity = pd.DataFrame(diversity_results)
print(f"  有多样性数据的患者: {len(df_diversity)}")

# Taxane方案中 R vs NR 比较
taxane_div = df_diversity[df_diversity['chemo_class'] == 'Taxane'].copy()
print(f"\n  Taxane方案 R vs NR 克隆多样性比较:")
print(f"  {'指标':<25} {'R均值':>10} {'NR均值':>10} {'变化方向':>10} {'MWU p':>10}")
print("  " + "-"*65)

metrics = [
    ('shannon_entropy', 'Shannon Entropy'),
    ('simpson_index', 'Simpson Index'),
    ('clonality', 'Clonality'),
    ('n_clones', 'Number of Clones'),
    ('top10_clone_fraction', 'Top10% Clone Fraction'),
]

diversity_pvals = []
p_values_for_fdr = []
for col, label in metrics:
    r_vals = taxane_div[taxane_div['response'] == 'R'][col].values
    nr_vals = taxane_div[taxane_div['response'] == 'NR'][col].values
    if len(r_vals) >= 2 and len(nr_vals) >= 2:
        stat, p = stats.mannwhitneyu(r_vals, nr_vals, alternative='two-sided')
        p_values_for_fdr.append(p)
        direction = 'R > NR' if np.mean(r_vals) > np.mean(nr_vals) else 'R < NR'
        print(f"  {label:<25} {np.mean(r_vals):>10.3f} {np.mean(nr_vals):>10.3f} {direction:>10} {p:>10.4f}")
        diversity_pvals.append({'metric': label, 'R_mean': np.mean(r_vals), 'NR_mean': np.mean(nr_vals),
                               'direction': direction, 'mwu_p': p, 'n_R': len(r_vals), 'n_NR': len(nr_vals)})

# FDR校正（Benjamini-Hochberg）
if len(p_values_for_fdr) > 0:
    _, q_values, _, _ = multipletests(p_values_for_fdr, method='fdr_bh')
    print(f"\n  FDR校正结果 ({len(p_values_for_fdr)} tests):")
    print(f"  {'指标':<25} {'p值':>10} {'q值':>10}")
    print(f"  {'-'*45}")
    for i, (dp, q) in enumerate(zip(diversity_pvals, q_values)):
        dp['mwu_q'] = q
        sig_p = '*' if dp['mwu_p'] < 0.05 else ''
        sig_q = '*' if q < 0.05 else ''
        print(f"  {dp['metric']:<25} {dp['mwu_p']:>10.4f}{sig_p} {q:>10.4f}{sig_q}")

# 保存
df_diversity.to_csv(cfg.result_path('GSE241934_clone_diversity.csv'), index=False)
pd.DataFrame(diversity_pvals).to_csv(cfg.result_path('GSE241934_diversity_comparison.csv'), index=False)
print(f"\n  ✅ 克隆多样性结果已保存")

# 可视化：多样性比较
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Panel A: Shannon entropy
ax = axes[0, 0]
r_vals = taxane_div[taxane_div['response'] == 'R']['shannon_entropy'].values
nr_vals = taxane_div[taxane_div['response'] == 'NR']['shannon_entropy'].values
bp = ax.boxplot([r_vals, nr_vals], labels=['R', 'NR'], patch_artist=True,
                medianprops={'color': 'black'})
bp['boxes'][0].set_facecolor('#E64B35')
bp['boxes'][1].set_facecolor('#4DBBD5')
for patch in bp['boxes']:
    patch.set_alpha(0.7)
for i, data in enumerate([r_vals, nr_vals]):
    x = np.random.normal(i + 1, 0.05, size=len(data))
    ax.scatter(x, data, alpha=0.6, color='black', s=30, zorder=3)
ax.set_ylabel('Shannon entropy')
ax.set_title('Shannon entropy\n(Taxane, R vs NR)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel B: Clonality
ax = axes[0, 1]
r_vals = taxane_div[taxane_div['response'] == 'R']['clonality'].values
nr_vals = taxane_div[taxane_div['response'] == 'NR']['clonality'].values
bp = ax.boxplot([r_vals, nr_vals], labels=['R', 'NR'], patch_artist=True,
                medianprops={'color': 'black'})
bp['boxes'][0].set_facecolor('#E64B35')
bp['boxes'][1].set_facecolor('#4DBBD5')
for patch in bp['boxes']:
    patch.set_alpha(0.7)
for i, data in enumerate([r_vals, nr_vals]):
    x = np.random.normal(i + 1, 0.05, size=len(data))
    ax.scatter(x, data, alpha=0.6, color='black', s=30, zorder=3)
ax.set_ylabel('Clonality')
ax.set_title('Clonality\n(Taxane, R vs NR)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel C: Number of clones
ax = axes[1, 0]
r_vals = taxane_div[taxane_div['response'] == 'R']['n_clones'].values
nr_vals = taxane_div[taxane_div['response'] == 'NR']['n_clones'].values
bp = ax.boxplot([r_vals, nr_vals], labels=['R', 'NR'], patch_artist=True,
                medianprops={'color': 'black'})
bp['boxes'][0].set_facecolor('#E64B35')
bp['boxes'][1].set_facecolor('#4DBBD5')
for patch in bp['boxes']:
    patch.set_alpha(0.7)
for i, data in enumerate([r_vals, nr_vals]):
    x = np.random.normal(i + 1, 0.05, size=len(data))
    ax.scatter(x, data, alpha=0.6, color='black', s=30, zorder=3)
ax.set_ylabel('Number of clones')
ax.set_title('Number of clones\n(Taxane, R vs NR)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel D: Simpson index
ax = axes[1, 1]
r_vals = taxane_div[taxane_div['response'] == 'R']['simpson_index'].values
nr_vals = taxane_div[taxane_div['response'] == 'NR']['simpson_index'].values
bp = ax.boxplot([r_vals, nr_vals], labels=['R', 'NR'], patch_artist=True,
                medianprops={'color': 'black'})
bp['boxes'][0].set_facecolor('#E64B35')
bp['boxes'][1].set_facecolor('#4DBBD5')
for patch in bp['boxes']:
    patch.set_alpha(0.7)
for i, data in enumerate([r_vals, nr_vals]):
    x = np.random.normal(i + 1, 0.05, size=len(data))
    ax.scatter(x, data, alpha=0.6, color='black', s=30, zorder=3)
ax.set_ylabel('Simpson index')
ax.set_title('Simpson index\n(Taxane, R vs NR)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
fig_path = cfg.result_path('GSE241934_FigS11_clone_diversity.png')
plt.savefig(fig_path, dpi=150, bbox_inches='tight')
plt.savefig(fig_path.replace('.png', '.pdf'), bbox_inches='tight')
plt.close()
print(f"  ✅ 克隆多样性图已保存: GSE241934_FigS11_clone_diversity.png")

# ============================================================
# Step 6.15: GSE179994 配对检验（治疗前后动态变化）
# ============================================================
print("\n" + "="*70)
print("Step 6.15: GSE179994 配对检验（治疗前后动态变化）")
print("="*70)

# 检查GSE179994数据是否存在
gse17994_dir = os.path.join(cfg.ADATA_DIR, 'GSE179994')
response_scores_path = os.path.join(gse17994_dir, 'GSE179994_patient_timepoint_scores.csv')

if os.path.exists(response_scores_path):
    df_179994 = pd.read_csv(response_scores_path)
    print(f"  GSE179994 数据已加载: {df_179994.shape}")
    print(f"  列名: {list(df_179994.columns)[:10]}")
    
    # 检查是否有配对数据
    patient_col = 'patient' if 'patient' in df_179994.columns else ('patient_id' if 'patient_id' in df_179994.columns else None)
    print(f"\n  患者数: {df_179994[patient_col].nunique() if patient_col else 'unknown'}")
    print(f"  时间点: {df_179994['timepoint'].unique() if 'timepoint' in df_179994.columns else 'unknown'}")
    
    # 查找有两个时间点的患者（配对）
    if 'timepoint' in df_179994.columns and patient_col is not None:
        patient_counts = df_179994.groupby(patient_col).size()
        paired_patients = patient_counts[patient_counts == 2].index.tolist()
        single_patients = patient_counts[patient_counts == 1].index.tolist()
        
        print(f"  有2个时间点的患者（配对）: {len(paired_patients)}")
        print(f"  有1个时间点的患者: {len(single_patients)}")
        
        # 确定哪个是pre哪个是on-treatment
        # 查看有两个时间点的患者
        print(f"\n  配对患者时间点详情:")
        for p in sorted(paired_patients)[:5]:
            rows = df_179994[df_179994[patient_col] == p]
            print(f"    {p}: {list(rows['timepoint'].values)} = {list(rows['mean_score'].round(4).values)}")
        
        # 方法：对于有pre和unknown的患者，假设pre是治疗前，unknown是治疗后（on-treatment）
        # 因为GSE179994是新辅助治疗，通常有pre和on-treatment时间点
        
        if len(paired_patients) >= 3:
            # 提取配对数据
            nk_col = 'mean_score'  # 使用均值得分
            print(f"\n  使用 NK-like 指标: {nk_col}")
            
            pre_vals = []
            post_vals = []
            paired_ids = []
            
            for p in sorted(paired_patients):
                rows = df_179994[df_179994[patient_col] == p]
                timepoints = list(rows['timepoint'].values)
                scores = list(rows[nk_col].values)
                
                # 尝试识别pre和post
                pre_val = None
                post_val = None
                
                if 'pre' in timepoints and 'unknown' in timepoints:
                    pre_val = rows[rows['timepoint'] == 'pre'][nk_col].values[0]
                    post_val = rows[rows['timepoint'] == 'unknown'][nk_col].values[0]
                elif len(timepoints) == 2:
                    # 假设第一个是pre，第二个是post
                    pre_val = scores[0]
                    post_val = scores[1]
                
                if pre_val is not None and post_val is not None:
                    pre_vals.append(pre_val)
                    post_vals.append(post_val)
                    paired_ids.append(p)
            
            if len(pre_vals) >= 3:
                pre_vals = np.array(pre_vals)
                post_vals = np.array(post_vals)
                
                # Wilcoxon 配对检验
                stat, p = stats.wilcoxon(pre_vals, post_vals)
                print(f"\n  配对 Wilcoxon 检验:")
                print(f"    Pre均值: {np.mean(pre_vals):.4f}")
                print(f"    On-treatment均值: {np.mean(post_vals):.4f}")
                print(f"    变化倍数: {np.mean(post_vals)/np.mean(pre_vals):.2f}")
                print(f"    上升比例: {np.mean(post_vals > pre_vals)*100:.1f}%")
                print(f"    p值: {p:.4f}")
                
                # 保存
                df_paired = pd.DataFrame({
                    'patient_id': paired_ids,
                    'pre_score': pre_vals,
                    'post_score': post_vals,
                    'change': post_vals - pre_vals,
                    'fold_change': post_vals / pre_vals
                })
                df_paired.to_csv(cfg.result_path('GSE179994_paired_analysis.csv'), index=False)
                print(f"  ✅ GSE179994 配对分析结果已保存")
                
                # 可视化
                fig, ax = plt.subplots(1, 1, figsize=(5, 5))
                for i in range(len(pre_vals)):
                    color = '#E64B35' if post_vals[i] > pre_vals[i] else '#4DBBD5'
                    ax.plot([0, 1], [pre_vals[i], post_vals[i]], 'o-', 
                            color=color, alpha=0.7, linewidth=1, markersize=5)
                ax.plot([0, 1], [np.mean(pre_vals), np.mean(post_vals)], 'o-', 
                        color='black', linewidth=3, markersize=10, label='Mean')
                ax.set_xticks([0, 1])
                ax.set_xticklabels(['Pre-treatment', 'On-treatment'])
                ax.set_ylabel('NK-like signature score (mean)')
                ax.set_title(f'GSE179994: NK-like score dynamics\n(paired, n={len(pre_vals)}, p={p:.3f})')
                ax.legend()
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                
                plt.tight_layout()
                fig_path = cfg.result_path('GSE179994_FigS12_paired.png')
                plt.savefig(fig_path, dpi=150, bbox_inches='tight')
                plt.savefig(fig_path.replace('.png', '.pdf'), bbox_inches='tight')
                plt.close()
                print(f"  ✅ GSE179994 配对分析图已保存")
            else:
                print(f"  有效配对数不足 ({len(pre_vals)})，跳过检验")
        else:
            print(f"  配对患者不足，跳过配对检验")
    else:
        print(f"  数据格式不匹配，跳过配对检验")
else:
    print(f"  GSE179994 数据文件不存在，跳过配对检验")

# ============================================================
# 汇总
# ============================================================
print("\n" + "="*70)
print("汇总统计")
print("="*70)

summary = {
    'n_patients_with_tcr': adata_tcr.obs['tcr_patient'].nunique(),
    'n_total_clones': n_total_clones,
    'n_shared_clones': n_shared_clones,
    'shared_clone_percent': n_shared_clones / n_total_clones * 100,
    'n_taxane_patients': len(taxane_div),
    'n_taxane_R': (taxane_div['response'] == 'R').sum(),
    'n_taxane_NR': (taxane_div['response'] == 'NR').sum(),
}

print(f"  TCR患者数: {summary['n_patients_with_tcr']}")
print(f"  总克隆数: {summary['n_total_clones']}")
print(f"  共享克隆数: {summary['n_shared_clones']} ({summary['shared_clone_percent']:.2f}%)")

pd.DataFrame([summary]).to_csv(cfg.result_path('GSE241934_tcr_deep_summary.csv'), index=False)

print("\n" + "="*70)
print("[step6e_tcr_diversity] Completed!")
print("="*70)
