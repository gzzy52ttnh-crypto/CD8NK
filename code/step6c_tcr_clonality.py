"""
step6c_tcr_clonality.py
Step 6.7: GSE241934 TCR 克隆型分析
- 从 RAW.tar 中提取 TCR 克隆型数据
- 与 scRNA-seq 细胞类型注释结合
- 计算 NK-dominant 大克隆比例（与主队列方法对齐）
- 用克隆水平指标重新做 Taxane 方案验证
"""
import config as cfg
import pandas as pd
import numpy as np
import os
import tarfile
import gzip
import io
import warnings
warnings.filterwarnings('ignore')

print("[step6c_tcr_clonality] Starting...")
os.makedirs(cfg.RESULT_DIR, exist_ok=True)

# ============================================================
# Step 1: 从 RAW.tar 读取 TCR 数据（不解压到磁盘，内存中处理）
# ============================================================
print("\n" + "="*70)
print("Step 1: 读取 TCR 克隆型数据")
print("="*70)

tcr_dir = os.path.join(cfg.ADATA_DIR, 'GSE241934_TCR')
print(f"  TCR 数据目录: {tcr_dir}")

all_contigs = []
all_clonotypes = []

# 查找所有患者文件
import glob
tcr_files = glob.glob(os.path.join(tcr_dir, '*.csv.gz'))
print(f"  总文件数: {len(tcr_files)}")

# 按患者分组
patient_files = {}
for fpath in tcr_files:
    fname = os.path.basename(fpath)
    parts = fname.split('_')
    if len(parts) >= 3:
        patient_id = parts[1]  # P343
        if patient_id not in patient_files:
            patient_files[patient_id] = {}
        if 'filtered_contig' in fname:
            patient_files[patient_id]['filtered_contig'] = fpath
        elif 'clonotypes' in fname:
            patient_files[patient_id]['clonotypes'] = fpath

print(f"  患者数: {len(patient_files)}")

# 逐患者读取
for pid, files in sorted(patient_files.items()):
    print(f"  处理患者 {pid}...")
    
    # 读取 filtered contig annotations
    if 'filtered_contig' in files:
        try:
            df_contig = pd.read_csv(files['filtered_contig'], compression='gzip')
            df_contig['patient_id'] = pid
            all_contigs.append(df_contig)
            print(f"    contigs: {len(df_contig)}")
        except Exception as e:
            print(f"    contig 读取失败: {e}")
    
    # 读取 clonotypes
    if 'clonotypes' in files:
        try:
            df_clono = pd.read_csv(files['clonotypes'], compression='gzip')
            df_clono['patient_id'] = pid
            all_clonotypes.append(df_clono)
            print(f"    clonotypes: {len(df_clono)}")
        except Exception as e:
            print(f"    clonotype 读取失败: {e}")

df_all_contigs = pd.concat(all_contigs, ignore_index=True)
df_all_clonotypes = pd.concat(all_clonotypes, ignore_index=True)

print(f"\n  Total contigs: {len(df_all_contigs)}")
print(f"  Total clonotypes: {len(df_all_clonotypes)}")
print(f"  患者数: {df_all_contigs['patient_id'].nunique()}")

# ============================================================
# Step 2: TCR barcode 前缀化，与 scRNA-seq cellID 对齐
# ============================================================
print("\n" + "="*70)
print("Step 2: TCR barcode 与 scRNA-seq cellID 对齐")
print("="*70)

# 读取患者得分文件获取患者信息
patient_scores_path = os.path.join(cfg.RESULT_DIR, 'GSE241934_patient_nklike_scores.csv')
df_patient = pd.read_csv(patient_scores_path, index_col=0)
print(f"  有临床信息的患者数: {len(df_patient)}")

# 为 TCR barcode 添加患者前缀: {patient_id}_{barcode}
df_all_contigs['cellID'] = df_all_contigs['patient_id'] + '_' + df_all_contigs['barcode']
print(f"  TCR 细胞数 (unique barcode): {df_all_contigs['cellID'].nunique()}")

# 检查匹配率（读取 h5ad 获取 scRNA-seq cellID）
h5ad_path = os.path.join(cfg.RESULT_DIR, 'GSE241934_all.h5ad')
import anndata as ad
if os.path.exists(h5ad_path):
    adata = ad.read_h5ad(h5ad_path)
    scrna_cells = set(adata.obs_names)
    tcr_cells = set(df_all_contigs['cellID'].unique())
    matched = scrna_cells & tcr_cells
    print(f"  scRNA-seq 细胞数: {len(scrna_cells)}")
    print(f"  TCR 细胞数: {len(tcr_cells)}")
    print(f"  匹配细胞数: {len(matched)}")
    print(f"  匹配率 (TCR中): {len(matched)/len(tcr_cells):.1%}")
    print(f"  匹配率 (scRNA中): {len(matched)/len(scrna_cells):.1%}")
else:
    print(f"  h5ad 不存在，跳过匹配率检查")

# ============================================================
# Step 3: 克隆型定义与细胞类型映射
# ============================================================
print("\n" + "="*70)
print("Step 3: 克隆型-细胞类型映射")
print("="*70)

# 从 contig 中提取 clonotype_id 与细胞的映射
# 注意：10x 的 clonotype_id 是 per-cell 的，同一克隆型的细胞共享相同的 clonotype_id
tcr_cell_clono = df_all_contigs[['cellID', 'raw_clonotype_id', 'patient_id']].drop_duplicates()
tcr_cell_clono = tcr_cell_clono.rename(columns={'raw_clonotype_id': 'clonotype_id'})
print(f"  细胞-克隆型对: {len(tcr_cell_clono)}")

# 获取 scRNA-seq 的细胞类型注释（仅 CD8+ T 细胞）
adata_cd8 = adata[adata.obs['cell.type'].str.contains('CD8', na=False)].copy()
print(f"  CD8+ T 细胞数: {len(adata_cd8)}")

# NK-like 细胞定义（与 Step 6.2 一致）
nklike_celltypes = [ct for ct in adata_cd8.obs['cell.type'].unique() 
                    if 'FGFBP2' in str(ct) and 'CD8' in str(ct)]
print(f"  NK-like 细胞类型: {nklike_celltypes}")

# Tex 细胞定义（与主队列对齐）
tex_celltypes = [ct for ct in adata_cd8.obs['cell.type'].unique() 
                 if 'Tex' in str(ct) or 'exhaust' in str(ct).lower()]
print(f"  Tex 细胞类型: {tex_celltypes}")

# 构建细胞类型表
df_celltype = pd.DataFrame({
    'cellID': adata_cd8.obs_names,
    'cell.type': adata_cd8.obs['cell.type'].values,
    'sampleID': adata_cd8.obs['sampleID'].values,
    'is_nklike': adata_cd8.obs['cell.type'].isin(nklike_celltypes).astype(int).values,
    'is_tex': adata_cd8.obs['cell.type'].isin(tex_celltypes).astype(int).values,
})

# 合并 TCR 克隆型与细胞类型
df_merged = tcr_cell_clono.merge(df_celltype, on='cellID', how='inner')
print(f"\n  合并后细胞数（有TCR + CD8+注释）: {len(df_merged)}")
print(f"  患者数: {df_merged['patient_id'].nunique()}")

# ============================================================
# Step 4: 克隆水平统计（NK-dominant 定义）
# ============================================================
print("\n" + "="*70)
print("Step 4: 克隆水平统计与分类")
print("="*70)

# 按患者+克隆型分组
clone_stats = []

for (pid, clono_id), group in df_merged.groupby(['patient_id', 'clonotype_id']):
    n_total = len(group)
    n_nk = group['is_nklike'].sum()
    n_tex = group['is_tex'].sum()
    nk_ratio = n_nk / n_total if n_total > 0 else 0
    tex_ratio = n_tex / n_total if n_total > 0 else 0
    
    # 克隆分类（与主队列 Step 2 一致，阈值统一引用 config）
    _nk_thresh = cfg.NK_DOMINANT_RATIO_THRESHOLD
    if nk_ratio >= _nk_thresh:
        category = 'NK-dominant'
    elif tex_ratio >= _nk_thresh:
        category = 'Tex-dominant'
    elif nk_ratio + tex_ratio >= _nk_thresh:
        category = 'Mixed'
    else:
        category = 'Other'
    
    clone_stats.append({
        'patient_id': pid,
        'clonotype_id': clono_id,
        'clone_size': n_total,
        'n_nklike': n_nk,
        'n_tex': n_tex,
        'nk_ratio': nk_ratio,
        'tex_ratio': tex_ratio,
        'category': category,
    })

df_clone = pd.DataFrame(clone_stats)
print(f"  总克隆数: {len(df_clone)}")
print(f"  分类分布:")
print(df_clone['category'].value_counts())

# 保存克隆水平统计
clone_out = os.path.join(cfg.RESULT_DIR, 'GSE241934_tcr_clone_stats.csv')
df_clone.to_csv(clone_out, index=False)
print(f"\n  ✅ 克隆统计已保存: {clone_out}")

# ============================================================
# Step 5: 患者水平 NK-dominant 比例（大克隆阈值=5，与主队列一致）
# ============================================================
print("\n" + "="*70)
print("Step 5: 患者水平 NK-dominant 大克隆比例")
print("="*70)

BIG_CLONE_THRESHOLD = cfg.BIG_CLONE_THRESHOLD  # 与主队列 Step 2 一致（统一引用 config）

patient_clone_metrics = []

for pid in sorted(df_clone['patient_id'].unique()):
    df_p = df_clone[df_clone['patient_id'] == pid]
    
    # 大克隆（clone_size >= threshold）
    df_big = df_p[df_p['clone_size'] >= BIG_CLONE_THRESHOLD]
    n_big_clones = len(df_big)
    
    if n_big_clones == 0:
        patient_clone_metrics.append({
            'patient_id': pid,
            'n_total_clones': len(df_p),
            'n_big_clones': 0,
            'nk_dominant_ratio': np.nan,
            'nk_dominant_count': 0,
            'tex_dominant_count': 0,
            'mixed_count': 0,
            'other_count': 0,
        })
        continue
    
    nk_dominant = (df_big['category'] == 'NK-dominant').sum()
    tex_dominant = (df_big['category'] == 'Tex-dominant').sum()
    mixed = (df_big['category'] == 'Mixed').sum()
    other = (df_big['category'] == 'Other').sum()
    nk_ratio = nk_dominant / n_big_clones
    
    patient_clone_metrics.append({
        'patient_id': pid,
        'n_total_clones': len(df_p),
        'n_big_clones': n_big_clones,
        'nk_dominant_ratio': nk_ratio,
        'nk_dominant_count': nk_dominant,
        'tex_dominant_count': tex_dominant,
        'mixed_count': mixed,
        'other_count': other,
    })

df_patient_clone = pd.DataFrame(patient_clone_metrics)
print(f"  有大克隆的患者数: {df_patient_clone['n_big_clones'].gt(0).sum()} / {len(df_patient_clone)}")

# 合并临床信息
df_patient_clone = df_patient_clone.merge(
    df_patient[['cohort', 'chemo_class', 'response_binary', 'analysis_subset']],
    left_on='patient_id', right_index=True, how='left'
)

# 保存患者水平克隆指标
patient_clone_out = os.path.join(cfg.RESULT_DIR, 'GSE241934_patient_clone_metrics.csv')
df_patient_clone.to_csv(patient_clone_out, index=False)
print(f"  ✅ 患者克隆指标已保存: {patient_clone_out}")

# ============================================================
# Step 6: 核心验证 - 用 NK-dominant 克隆比例做 Taxane 方案验证
# ============================================================
print("\n" + "="*70)
print("Step 6: NK-dominant 克隆比例验证（Taxane 方案）")
print("="*70)

from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm

validation_results = []
mwu_p_values = []

# 仅使用有大克隆的患者
df_valid = df_patient_clone[df_patient_clone['n_big_clones'] > 0].copy()
print(f"  有效患者数（有大克隆）: {len(df_valid)}")

for subset_key, subset_label in [
    ('A_IIT_Taxane', 'IIT Taxane (EGFR-MT)'),
    ('B_RWC_Taxane', 'RWC Taxane (EGFR-WT)'),
    ('Taxane_combined', 'Taxane Combined'),
    ('C_RWC_Pemetrexed', 'RWC Pemetrexed'),
]:
    if subset_key == 'Taxane_combined':
        sub = df_valid[df_valid['chemo_class'] == 'Taxane']
    else:
        sub = df_valid[df_valid['analysis_subset'] == subset_key]
    
    if len(sub) < 5:
        print(f"  {subset_label}: n={len(sub)}，样本过少跳过")
        continue
    
    r_vals = sub[sub['response_binary'] == 'R']['nk_dominant_ratio'].values
    nr_vals = sub[sub['response_binary'] == 'NR']['nk_dominant_ratio'].values
    
    if len(r_vals) == 0 or len(nr_vals) == 0:
        print(f"  {subset_label}: R={len(r_vals)}, NR={len(nr_vals)}，跳过")
        continue
    
    # MWU 检验
    stat, mwu_p = mannwhitneyu(r_vals, nr_vals, alternative='two-sided')
    mwu_p_values.append(mwu_p)
    
    # AUC
    try:
        auc = roc_auc_score(
            [1]*len(r_vals) + [0]*len(nr_vals),
            list(r_vals) + list(nr_vals)
        )
    except:
        auc = np.nan
    
    print(f"  {subset_label}: n={len(sub)} (R={len(r_vals)}, NR={len(nr_vals)})")
    print(f"    R 均值: {np.nanmean(r_vals):.4f}, NR 均值: {np.nanmean(nr_vals):.4f}")
    print(f"    MWU p = {mwu_p:.4f}, AUC = {auc:.4f}")
    
    validation_results.append({
        'subset': subset_key,
        'n_total': len(sub),
        'n_R': len(r_vals),
        'n_NR': len(nr_vals),
        'mean_R': np.nanmean(r_vals),
        'mean_NR': np.nanmean(nr_vals),
        'mwu_p': mwu_p,
        'auc': auc,
        'metric': 'nk_dominant_ratio (clone-based)',
    })

# FDR校正（Benjamini-Hochberg）
if len(mwu_p_values) > 0:
    _, q_values, _, _ = multipletests(mwu_p_values, method='fdr_bh')
    print(f"\n  FDR校正结果 ({len(mwu_p_values)} tests):")
    for i, (res, q) in enumerate(zip(validation_results, q_values)):
        res['mwu_q'] = q
        print(f"    {res['subset']}: p={res['mwu_p']:.4f}, q={q:.4f}")

df_val_clone = pd.DataFrame(validation_results)
val_clone_out = os.path.join(cfg.RESULT_DIR, 'GSE241934_tcr_validation_results.csv')
df_val_clone.to_csv(val_clone_out, index=False)
print(f"\n  ✅ TCR 克隆验证结果已保存: {val_clone_out}")

# ============================================================
# Step 7: 与细胞比例方法对比（方法学一致性验证）
# ============================================================
print("\n" + "="*70)
print("Step 7: 克隆比例 vs 细胞比例方法对比")
print("="*70)

# 合并两种方法的患者水平指标
df_compare = df_patient_clone.merge(
    df_patient[['frac_nklike_of_cd8', 'mean_nklike_score_cd8']],
    left_on='patient_id', right_index=True, how='inner'
)

# 计算相关性
from scipy.stats import spearmanr

valid_compare = df_compare[df_compare['n_big_clones'] > 0].copy()

corr_r, corr_p = spearmanr(
    valid_compare['nk_dominant_ratio'],
    valid_compare['frac_nklike_of_cd8'],
    nan_policy='omit'
)
print(f"  NK-dominant克隆比例 vs NK-like细胞比例 (Spearman):")
print(f"    r = {corr_r:.4f}, p = {corr_p:.4f}")
print(f"    n = {len(valid_compare)}")

# 保存对比结果
compare_out = os.path.join(cfg.RESULT_DIR, 'GSE241934_clone_vs_cell_fraction_comparison.csv')
df_compare.to_csv(compare_out, index=False)
print(f"  ✅ 方法对比结果已保存: {compare_out}")

# ============================================================
# Step 8: 可视化
# ============================================================
print("\n" + "="*70)
print("Step 8: 可视化")
print("="*70)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['pdf.fonttype'] = 42

color_R = '#E64B35'
color_NR = '#4DBBD5'
color_taxane = '#E64B35'
color_pemetrexed = '#4DBBD5'

fig, axes = plt.subplots(2, 2, figsize=(14, 11))
fig.suptitle('GSE241934: TCR Clonotype Analysis (NK-dominant clones)', 
             fontsize=14, fontweight='bold', y=1.02)

# Panel A: 各子集 NK-dominant 比例箱线图
ax = axes[0, 0]
subset_info = [
    ('A_IIT_Taxane', 'IIT\nTaxane', '#91D1C2'),
    ('B_RWC_Taxane', 'RWC\nTaxane', '#F39B7F'),
    ('Taxane_combined', 'Taxane\nCombined', color_taxane),
    ('C_RWC_Pemetrexed', 'RWC\nPemetrexed', color_pemetrexed),
]

box_data_r = []
box_data_nr = []
labels = []
for sk, label, _ in subset_info:
    if sk == 'Taxane_combined':
        sub = df_valid[df_valid['chemo_class'] == 'Taxane']
    else:
        sub = df_valid[df_valid['analysis_subset'] == sk]
    if len(sub) < 3:
        continue
    r_vals = sub[sub['response_binary'] == 'R']['nk_dominant_ratio'].values
    nr_vals = sub[sub['response_binary'] == 'NR']['nk_dominant_ratio'].values
    if len(r_vals) > 0 and len(nr_vals) > 0:
        box_data_nr.append(nr_vals)
        box_data_r.append(r_vals)
        labels.append(label)

# 绘制分组箱线图（每组两个箱子：NR和R）
x_positions = []
all_data = []
all_colors = []
for i in range(len(labels)):
    x_positions.extend([i*3, i*3 + 1])
    all_data.append(box_data_nr[i])
    all_data.append(box_data_r[i])
    all_colors.extend([color_NR, color_R])

bp = ax.boxplot(all_data, positions=x_positions, patch_artist=True,
                widths=0.7, medianprops={'color': 'black', 'linewidth': 1.5})
for patch, color in zip(bp['boxes'], all_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

# 配对 MWU p 值标注
for i in range(len(labels)):
    r_data = box_data_r[i]
    nr_data = box_data_nr[i]
    _, p_val = mannwhitneyu(r_data, nr_data, alternative='two-sided')
    y_max = max(np.nanmax(r_data), np.nanmax(nr_data))
    p_text = f'p={p_val:.3f}'
    ax.text(i*3 + 0.5, y_max * 1.1, p_text, ha='center', fontsize=9, fontweight='bold')

ax.set_xticks([i*3 + 0.5 for i in range(len(labels))])
ax.set_xticklabels(labels, fontsize=10)
ax.set_ylabel('NK-dominant clone ratio')
ax.set_title('Panel A: NK-dominant clone ratio by response', fontsize=11, fontweight='bold')
ax.set_ylim(bottom=0)

# Panel B: 克隆分类饼图（Taxane 合并）
ax = axes[0, 1]
df_taxane_big = df_clone[
    df_clone['clone_size'] >= BIG_CLONE_THRESHOLD
].merge(
    df_patient[['chemo_class']], left_on='patient_id', right_index=True, how='left'
)
df_taxane_big = df_taxane_big[df_taxane_big['chemo_class'] == 'Taxane']
cat_counts = df_taxane_big['category'].value_counts()
colors_pie = ['#E64B35', '#4DBBD5', '#91D1C2', '#999999']
ax.pie(cat_counts.values, labels=cat_counts.index, autopct='%1.1f%%',
       colors=colors_pie[:len(cat_counts)], startangle=90)
ax.set_title(f'Panel B: Clone categories (Taxane, n={len(df_taxane_big)})', 
             fontsize=11, fontweight='bold')

# Panel C: 克隆大小分布
ax = axes[1, 0]
ax.hist(df_clone['clone_size'], bins=50, color='#3498DB', alpha=0.7, edgecolor='white')
ax.axvline(x=BIG_CLONE_THRESHOLD, color='red', linestyle='--', linewidth=2, 
           label=f'Threshold = {BIG_CLONE_THRESHOLD}')
ax.set_xlabel('Clone size (cell count)')
ax.set_ylabel('Number of clones')
ax.set_title('Panel C: Clone size distribution', fontsize=11, fontweight='bold')
ax.set_yscale('log')
ax.legend()

# Panel D: 克隆比例 vs 细胞比例散点图
ax = axes[1, 1]
ax.scatter(valid_compare['frac_nklike_of_cd8'], valid_compare['nk_dominant_ratio'],
           c=valid_compare['response_binary'].map({'R': color_R, 'NR': color_NR}),
           alpha=0.7, s=50, edgecolors='white', linewidth=0.5)
ax.set_xlabel('NK-like cell fraction (of CD8+)')
ax.set_ylabel('NK-dominant clone ratio')
ax.set_title(f'Panel D: Clone-based vs cell-based (Spearman r={corr_r:.3f})', 
             fontsize=11, fontweight='bold')

# 添加图例
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=color_R, markersize=10, label='Responder'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=color_NR, markersize=10, label='Non-responder'),
]
ax.legend(handles=legend_elements, loc='upper left', fontsize=9)

plt.tight_layout()
fig_path = os.path.join(cfg.RESULT_DIR, 'GSE241934_Fig6_TCR_clonality.png')
fig.savefig(fig_path, dpi=300, bbox_inches='tight')
fig_pdf = os.path.join(cfg.RESULT_DIR, 'GSE241934_Fig6_TCR_clonality.pdf')
fig.savefig(fig_pdf, bbox_inches='tight')
print(f"  ✅ TCR 克隆分析图已保存: {fig_path}")
plt.close(fig)

# ============================================================
# Step 9: 汇总统计
# ============================================================
print("\n" + "="*70)
print("Step 9: 汇总统计")
print("="*70)

summary = {
    'n_patients_with_tcr': df_all_contigs['patient_id'].nunique(),
    'n_cells_with_tcr_cd8': len(df_merged),
    'n_total_clones': len(df_clone),
    'n_big_clones': (df_clone['clone_size'] >= BIG_CLONE_THRESHOLD).sum(),
    'big_clone_threshold': BIG_CLONE_THRESHOLD,
    'nk_dominant_clones': (df_clone['category'] == 'NK-dominant').sum(),
    'tex_dominant_clones': (df_clone['category'] == 'Tex-dominant').sum(),
    'mixed_clones': (df_clone['category'] == 'Mixed').sum(),
    'other_clones': (df_clone['category'] == 'Other').sum(),
    'patients_with_big_clones': df_valid['patient_id'].nunique(),
    'spearman_r_clone_vs_fraction': corr_r,
    'spearman_p_clone_vs_fraction': corr_p,
}

print("\n  TCR 克隆分析汇总:")
for k, v in summary.items():
    if isinstance(v, float):
        print(f"    {k}: {v:.4f}")
    else:
        print(f"    {k}: {v}")

df_summary = pd.DataFrame([summary])
summary_out = os.path.join(cfg.RESULT_DIR, 'GSE241934_tcr_summary.csv')
df_summary.to_csv(summary_out, index=False)
print(f"\n  ✅ 汇总统计已保存: {summary_out}")

print("\n" + "="*70)
print("[step6c_tcr_clonality] Completed!")
print("="*70)
