"""
step4b_chemo_dynamics.py
化疗对 NK-like CD8+ T 细胞预测效力的影响分析 (GSE179994, Liu 2022 Nature Cancer)

数据集真实队列结构 (核对自论文 Supplementary Table 1):
- GSE179994: NSCLC, 150,849 T 细胞, 36 患者, T 细胞-only (无髓系)
- 子队列 1 (Treatment-naïve): 25 例 (P2-9, P11, P12, P14-18, P20-28, P34)
  * 未接受免疫治疗, 仅单时间点活检, 无响应标签
  * 不纳入化疗影响分析
- 子队列 2 (Treatment A): 8 例 (P1, P10, P13, P19, P29, P30, P33, P35)
  * Pembro + Carbo + Pemetrexed, pre+post 配对, 全部 Responder
- 子队列 3 (Treatment B): 3 例 (P36, P37, P38)
  * Pembro + Carbo + 白蛋白结合紫杉醇, 仅 post, 全部 Non-responder

核心数据限制 (硬性限制, 非提取问题):
- NR 组 (P36/37/38) 仅有 post 样本, 无 pre 配对
- 因此无法进行 R vs NR 的化疗动力学比较
- 8 例配对患者全部为 R, 仅能分析"化疗在响应者中的动力学"
- Treatment A vs B 化疗方案不同 (Pemetrexed vs 白蛋白结合紫杉醇), 亦为混杂因素

分析框架 (受数据限制, 仅能验证部分假设):
H1: 化疗直接杀伤 NK-like 细胞 → NK-like 比例 pre→post 下降 (仅 R 组)
H2: 化疗改变 T 细胞分化路径 → Tex 比例与 cluster 构成变化 (仅 R 组)
H3: 化疗导致免疫重塑 → NK-like 签名得分 pre→post 动力学 (仅 R 组)
H4: 髓系闸门 (SPP1+ TAM) → 跳过 (T-only 数据无髓系)
注: 无法回答"化疗如何差异化影响 R vs NR" (NR 无配对)

统计检验:
1. 配对 Wilcoxon signed-rank (NK-like 得分 pre vs post, 仅 R 组 n=8)
2. 配对 Wilcoxon signed-rank (NK-like high 比例 pre vs post, 仅 R 组 n=8)
3. 配对 Wilcoxon signed-rank (Tex 比例 pre vs post, 仅 R 组 n=8)
4. 卡方检验 (cluster 构成 pre vs post, pooled cells, 仅 R 组)
5. 配对 t-test (敏感性分析, 与 Wilcoxon 对照)

输出:
- Fig_chemo_dynamics.png/pdf (4 panels)
- chemo_dynamics_results.csv
- chemo_paired_metrics.csv
"""
import config as cfg
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon, mannwhitneyu, chi2_contingency, ttest_rel
from statsmodels.stats.multitest import multipletests
import os
import warnings
warnings.filterwarnings('ignore')

print("[step4b_chemo_dynamics] Starting chemo-immune dynamics analysis...")

os.makedirs(cfg.RESULT_DIR, exist_ok=True)

# ============================================================
# 数据加载
# ============================================================
META_PATH = os.path.join(cfg.ROOT_DIR, 'GSE179994_Tcell.metadata.tsv')
SCORES_PATH = os.path.join(cfg.ADATA_DIR, 'GSE179994', 'GSE179994_all_signature_scores.csv')

if not os.path.exists(META_PATH) or not os.path.exists(SCORES_PATH):
    raise FileNotFoundError(f"GSE179994 data not found. Paths:\n  {META_PATH}\n  {SCORES_PATH}")

meta = pd.read_csv(META_PATH, sep='\t')
scores = pd.read_csv(SCORES_PATH)
df = meta.merge(scores, left_on='cellid', right_on='cell', how='inner')
print(f"  Loaded {len(df):,} cells, {df['patient'].nunique()} patients")

# 时间点: sample 列使用 'P1.pre' / 'P1.post.1' 格式
df['timepoint'] = df['sample'].apply(
    lambda x: 'pre' if '.pre' in str(x)
    else ('post' if '.post' in str(x) else 'unknown'))

# CD8 T 细胞
df_cd8 = df[df['celltype'] == 'CD8'].copy()
print(f"  CD8+ T cells: {len(df_cd8):,}")
print(f"  CD8 by timepoint: pre={len(df_cd8[df_cd8['timepoint']=='pre']):,}, "
      f"post={len(df_cd8[df_cd8['timepoint']=='post']):,}")

# 配对患者 (pre+post 都有)
pre_pats = set(df_cd8[df_cd8['timepoint'] == 'pre']['patient'].unique())
post_pats = set(df_cd8[df_cd8['timepoint'] == 'post']['patient'].unique())
paired_pats = sorted(pre_pats & post_pats)
print(f"  Paired patients (pre+post): {len(paired_pats)} -> {paired_pats}")

# 响应标签 + 治疗方案 (核对自论文 Supplementary Table 1)
# 子队列结构:
#   - 25 例 Treatment-naïve (P2-9, P11, P12, P14-18, P20-28, P34): 无响应标签, 不纳入
#   - 8 例 Treatment A (P1, P10, P13, P19, P29, P30, P33, P35): Pembro+Carbo+Pemetrexed, 全 R
#   - 3 例 Treatment B (P36, P37, P38): Pembro+Carbo+白蛋白结合紫杉醇, 全 NR, 仅 post
G17_PATIENT_RESPONSE = {
    'P1': 'Responder',     # P1.post.1=R, P1.post.2=R, P1.post.3=NR (混合, dominant R)
    'P10': 'Responder',
    'P13': 'Responder',    # P13.post.1=R, P13.post.3=NR (混合, dominant R)
    'P19': 'Responder',
    'P29': 'Responder',
    'P30': 'Responder',
    'P33': 'Responder',
    'P35': 'Responder',
}
G17_PATIENT_REGIMEN = {
    'P1': 'Pembro+Carbo+Pemetrexed', 'P10': 'Pembro+Carbo+Pemetrexed',
    'P13': 'Pembro+Carbo+Pemetrexed', 'P19': 'Pembro+Carbo+Pemetrexed',
    'P29': 'Pembro+Carbo+Pemetrexed', 'P30': 'Pembro+Carbo+Pemetrexed',
    'P33': 'Pembro+Carbo+Pemetrexed', 'P35': 'Pembro+Carbo+Pemetrexed',
}
# 核心限制: NR 组 (P36/37/38, Treatment B) 仅有 post 样本, 无法配对
# 因此本分析仅反映"化疗在响应者(Treatment A)中的动力学", 不能外推到 NR 组

# NK-like 签名得分阈值 (q75 of all CD8 cells, 数据驱动)
NKLlike_THRESHOLD = df_cd8['score'].quantile(0.75)
print(f"  NK-like threshold (q75 of CD8 score): {NKLlike_THRESHOLD:.4f}")

# ============================================================
# 患者水平指标计算
# ============================================================
paired_df = df_cd8[df_cd8['patient'].isin(paired_pats)].copy()
paired_df['is_nklike_high'] = (paired_df['score'] > NKLlike_THRESHOLD).astype(int)

# 患者×时间点 汇总
patient_tp_metrics = paired_df.groupby(['patient', 'timepoint']).agg(
    mean_score=('score', 'mean'),
    median_score=('score', 'median'),
    n_cells=('score', 'count'),
    nklike_high_prop=('is_nklike_high', 'mean'),
).reset_index()

# 加入 cluster 比例
cluster_pivot = paired_df.groupby(['patient', 'timepoint', 'cluster']).size().unstack(fill_value=0)
for cl in ['Non-exhausted', 'Tex', 'Prolif.']:
    if cl not in cluster_pivot.columns:
        cluster_pivot[cl] = 0
cluster_pivot['total'] = cluster_pivot.sum(axis=1)
for cl in ['Non-exhausted', 'Tex', 'Prolif.']:
    cluster_pivot[f'{cl}_prop'] = cluster_pivot[cl] / cluster_pivot['total']
cluster_pivot = cluster_pivot.reset_index()

patient_tp_metrics = patient_tp_metrics.merge(cluster_pivot[['patient', 'timepoint', 'Non-exhausted_prop', 'Tex_prop', 'Prolif._prop']],
                                              on=['patient', 'timepoint'], how='left')
patient_tp_metrics['response'] = patient_tp_metrics['patient'].map(G17_PATIENT_RESPONSE)

print("\n=== Per-patient paired metrics ===")
print(patient_tp_metrics.to_string(index=False))

# 保存
patient_tp_metrics.to_csv(cfg.result_path('chemo_paired_metrics.csv'), index=False)
print(f"\n  Saved: chemo_paired_metrics.csv")

# Pivot 为 pre/post 配对格式
pre_vals = patient_tp_metrics[patient_tp_metrics['timepoint'] == 'pre'].set_index('patient')
post_vals = patient_tp_metrics[patient_tp_metrics['timepoint'] == 'post'].set_index('patient')

# ============================================================
# 统计检验
# ============================================================
results = []

# --- Test 1: NK-like signature score pre vs post (Wilcoxon signed-rank) ---
pre_score = pre_vals.loc[paired_pats, 'mean_score'].values
post_score = post_vals.loc[paired_pats, 'mean_score'].values
n_up = int(np.sum(post_score > pre_score))
n_down = int(np.sum(post_score < pre_score))
try:
    stat_w1, p_w1 = wilcoxon(post_score, pre_score, alternative='two-sided')
    stat_w1_alt_less, p_w1_less = wilcoxon(post_score - pre_score, alternative='less')
    stat_w1_alt_greater, p_w1_greater = wilcoxon(post_score - pre_score, alternative='greater')
except ValueError as e:
    stat_w1, p_w1 = np.nan, np.nan
    p_w1_less, p_w1_greater = np.nan, np.nan
    print(f"  [Warning] Wilcoxon failed: {e}")

delta_score = post_score - pre_score
results.append({
    'test': 'NK-like signature score (paired Wilcoxon signed-rank)',
    'n_paired': len(paired_pats),
    'pre_mean': round(float(np.mean(pre_score)), 4),
    'post_mean': round(float(np.mean(post_score)), 4),
    'mean_delta': round(float(np.mean(delta_score)), 4),
    'n_up': n_up,
    'n_down': n_down,
    'statistic': round(float(stat_w1), 4) if not np.isnan(stat_w1) else 'NaN',
    'p_two_sided': round(float(p_w1), 4) if not np.isnan(p_w1) else 'NaN',
    'p_less (post<pre)': round(float(p_w1_less), 4) if not np.isnan(p_w1_less) else 'NaN',
    'p_greater (post>pre)': round(float(p_w1_greater), 4) if not np.isnan(p_w1_greater) else 'NaN',
    'interpretation': 'Chemotherapy reduces NK-like signature' if np.mean(delta_score) < 0 else 'Chemotherapy increases NK-like signature'
})
print(f"\n  Test 1 (NK-like score pre vs post): mean Δ={np.mean(delta_score):.4f}, "
      f"n_up={n_up}, n_down={n_down}, Wilcoxon p={p_w1:.4f}")

# --- Test 2: NK-like high proportion pre vs post (Wilcoxon signed-rank) ---
pre_prop = pre_vals.loc[paired_pats, 'nklike_high_prop'].values
post_prop = post_vals.loc[paired_pats, 'nklike_high_prop'].values
try:
    stat_w2, p_w2 = wilcoxon(post_prop, pre_prop, alternative='two-sided')
    stat_w2_less, p_w2_less = wilcoxon(post_prop, pre_prop, alternative='less')
    stat_w2_greater, p_w2_greater = wilcoxon(post_prop, pre_prop, alternative='greater')
except ValueError:
    stat_w2, p_w2 = np.nan, np.nan
    p_w2_less, p_w2_greater = np.nan, np.nan
delta_prop = post_prop - pre_prop
results.append({
    'test': 'NK-like high proportion (paired Wilcoxon signed-rank)',
    'n_paired': len(paired_pats),
    'pre_mean': round(float(np.mean(pre_prop)), 4),
    'post_mean': round(float(np.mean(post_prop)), 4),
    'mean_delta': round(float(np.mean(delta_prop)), 4),
    'n_up': int(np.sum(delta_prop > 0)),
    'n_down': int(np.sum(delta_prop < 0)),
    'statistic': round(float(stat_w2), 4) if not np.isnan(stat_w2) else np.nan,
    'p_two_sided': round(float(p_w2), 4) if not np.isnan(p_w2) else np.nan,
    'p_less (post<pre)': round(float(p_w2_less), 4) if not np.isnan(p_w2_less) else np.nan,
    'p_greater (post>pre)': round(float(p_w2_greater), 4) if not np.isnan(p_w2_greater) else np.nan,
    'interpretation': 'Chemotherapy reduces NK-like cell proportion' if np.mean(delta_prop) < 0 else 'Chemotherapy increases NK-like cell proportion'
})
print(f"  Test 2 (NK-like high prop pre vs post): mean Δ={np.mean(delta_prop):.4f}, "
      f"Wilcoxon p={p_w2:.4f}")

# --- Test 3: Tex proportion pre vs post (Wilcoxon signed-rank) ---
pre_tex = pre_vals.loc[paired_pats, 'Tex_prop'].values
post_tex = post_vals.loc[paired_pats, 'Tex_prop'].values
try:
    stat_w3, p_w3 = wilcoxon(post_tex, pre_tex, alternative='two-sided')
    stat_w3_less, p_w3_less = wilcoxon(post_tex, pre_tex, alternative='less')
    stat_w3_greater, p_w3_greater = wilcoxon(post_tex, pre_tex, alternative='greater')
except ValueError:
    stat_w3, p_w3 = np.nan, np.nan
    p_w3_less, p_w3_greater = np.nan, np.nan
delta_tex = post_tex - pre_tex
results.append({
    'test': 'Tex proportion (paired Wilcoxon signed-rank)',
    'n_paired': len(paired_pats),
    'pre_mean': round(float(np.mean(pre_tex)), 4),
    'post_mean': round(float(np.mean(post_tex)), 4),
    'mean_delta': round(float(np.mean(delta_tex)), 4),
    'n_up': int(np.sum(delta_tex > 0)),
    'n_down': int(np.sum(delta_tex < 0)),
    'statistic': round(float(stat_w3), 4) if not np.isnan(stat_w3) else np.nan,
    'p_two_sided': round(float(p_w3), 4) if not np.isnan(p_w3) else np.nan,
    'p_less (post<pre)': round(float(p_w3_less), 4) if not np.isnan(p_w3_less) else np.nan,
    'p_greater (post>pre)': round(float(p_w3_greater), 4) if not np.isnan(p_w3_greater) else np.nan,
    'interpretation': 'Chemotherapy reduces Tex proportion' if np.mean(delta_tex) < 0 else 'Chemotherapy increases Tex proportion'
})
print(f"  Test 3 (Tex prop pre vs post): mean Δ={np.mean(delta_tex):.4f}, "
      f"Wilcoxon p={p_w3:.4f}")

# --- Test 4: CD8 cluster composition pre vs post (pooled chi-squared) ---
# Pool cells across all paired patients by timepoint
cluster_counts = paired_df.groupby(['timepoint', 'cluster']).size().unstack(fill_value=0)
# Ensure all 3 clusters present
for cl in ['Non-exhausted', 'Tex', 'Prolif.']:
    if cl not in cluster_counts.columns:
        cluster_counts[cl] = 0
cluster_counts = cluster_counts[['Non-exhausted', 'Tex', 'Prolif.']]
print(f"\n  Pooled cluster counts (pre vs post):\n{cluster_counts}")

chi2, p_chi, dof, expected = chi2_contingency(cluster_counts.values)
results.append({
    'test': 'CD8 cluster composition (pooled chi-squared)',
    'n_paired': len(paired_pats),
    'pre_mean': f"Non-exh={int(cluster_counts.loc['pre','Non-exhausted'])}, "
                f"Tex={int(cluster_counts.loc['pre','Tex'])}, "
                f"Prolif={int(cluster_counts.loc['pre','Prolif.'])}",
    'post_mean': f"Non-exh={int(cluster_counts.loc['post','Non-exhausted'])}, "
                 f"Tex={int(cluster_counts.loc['post','Tex'])}, "
                 f"Prolif={int(cluster_counts.loc['post','Prolif.'])}",
    'mean_delta': np.nan,
    'n_up': np.nan,
    'n_down': np.nan,
    'statistic': round(float(chi2), 4),
    'p_two_sided': round(float(p_chi), 6),
    'p_less (post<pre)': np.nan,
    'p_greater (post>pre)': np.nan,
    'interpretation': 'Cluster composition significantly changed' if p_chi < 0.05 else 'No significant cluster composition change'
})
print(f"  Test 4 (cluster composition chi-squared): chi2={chi2:.2f}, df={dof}, p={p_chi:.6f}")

# --- Test 5: Sensitivity - paired t-test on NK-like score ---
try:
    stat_t, p_t = ttest_rel(post_score, pre_score)
except Exception:
    stat_t, p_t = np.nan, np.nan
results.append({
    'test': 'NK-like signature score (paired t-test, sensitivity)',
    'n_paired': len(paired_pats),
    'pre_mean': round(float(np.mean(pre_score)), 4),
    'post_mean': round(float(np.mean(post_score)), 4),
    'mean_delta': round(float(np.mean(delta_score)), 4),
    'n_up': n_up,
    'n_down': n_down,
    'statistic': round(float(stat_t), 4) if not np.isnan(stat_t) else 'NaN',
    'p_two_sided': round(float(p_t), 4) if not np.isnan(p_t) else 'NaN',
    'p_less (post<pre)': 'NaN',
    'p_greater (post>pre)': 'NaN',
    'interpretation': 'Consistent with Wilcoxon' if (p_t < 0.05) == (p_w1 < 0.05) else 'Divergent from Wilcoxon'
})

# FDR校正（Benjamini-Hochberg）
p_values_wilcoxon = []
for i, r in enumerate(results[:3]):
    p = r['p_two_sided']
    if isinstance(p, float) and not np.isnan(p):
        p_values_wilcoxon.append(p)
    else:
        p_values_wilcoxon.append(np.nan)

if len(p_values_wilcoxon) > 0 and not all(np.isnan(p) for p in p_values_wilcoxon):
    valid_mask = ~np.isnan(p_values_wilcoxon)
    valid_p = np.array(p_values_wilcoxon)[valid_mask]
    if len(valid_p) > 0:
        _, q_values_valid, _, _ = multipletests(valid_p, method='fdr_bh')
        q_values = np.full(len(p_values_wilcoxon), np.nan)
        q_values[valid_mask] = q_values_valid
        for i in range(3):
            results[i]['q_two_sided'] = round(float(q_values[i]), 4) if not np.isnan(q_values[i]) else 'NaN'
else:
    for i in range(3):
        results[i]['q_two_sided'] = 'NaN'

# 保存结果
results_df = pd.DataFrame(results)
results_df.to_csv(cfg.result_path('chemo_dynamics_results.csv'), index=False)
print(f"\n  Saved: chemo_dynamics_results.csv")
print("\n=== Statistical results summary (with FDR correction) ===")
print(results_df[['test', 'n_paired', 'mean_delta', 'p_two_sided', 'q_two_sided', 'interpretation']].to_string(index=False))

# ============================================================
# 可视化 (4-panel figure)
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 11))

COLOR_PRE = '#2196F3'   # Blue for pre-treatment
COLOR_POST = '#FF9800'  # Orange for post-treatment
COLOR_UP = '#4CAF50'
COLOR_DOWN = '#F44336'

# --- Panel A: NK-like signature score pre→post paired dynamics ---
axA = axes[0, 0]
x_pre = np.zeros(len(paired_pats))
x_post = np.ones(len(paired_pats))
for i, pat in enumerate(paired_pats):
    pre_v = pre_vals.loc[pat, 'mean_score']
    post_v = post_vals.loc[pat, 'mean_score']
    color_line = COLOR_UP if post_v > pre_v else COLOR_DOWN
    axA.plot([0, 1], [pre_v, post_v], '-', color=color_line, alpha=0.5, lw=1.2)
    axA.scatter([0, 1], [pre_v, post_v], c=[COLOR_PRE, COLOR_POST], s=80, zorder=5,
                edgecolors='black', linewidth=0.5)
    # 标注患者 ID
    axA.annotate(pat, (1.02, post_v), fontsize=7, va='center')

# 群体均值
axA.plot([0, 1], [np.mean(pre_score), np.mean(post_score)], 'k--', lw=2.5, alpha=0.7,
         label=f'Group mean (Δ={np.mean(delta_score):.4f})')
axA.scatter([0, 1], [np.mean(pre_score), np.mean(post_score)], c='black', s=150, marker='D',
            zorder=10, edgecolors='white', linewidth=1.5)

axA.set_xticks([0, 1])
axA.set_xticklabels(['Pre-treatment', 'Post-treatment'], fontsize=10)
axA.set_ylabel('NK-like signature score (mean per patient)', fontsize=10)
axA.set_xlim(-0.3, 1.5)
axA.set_title(f'Panel A: NK-like Score Dynamics (n={len(paired_pats)} paired)\n'
              f'Wilcoxon p={p_w1:.4f} | ↑{n_up} / ↓{n_down} | Δ={np.mean(delta_score):+.4f}',
              fontsize=11)
axA.legend(loc='upper right', fontsize=8)
axA.grid(alpha=0.3, ls=':')
axA.spines['top'].set_visible(False); axA.spines['right'].set_visible(False)

# --- Panel B: NK-like high proportion pre→post paired dynamics ---
axB = axes[0, 1]
for i, pat in enumerate(paired_pats):
    pre_v = pre_vals.loc[pat, 'nklike_high_prop']
    post_v = post_vals.loc[pat, 'nklike_high_prop']
    color_line = COLOR_UP if post_v > pre_v else COLOR_DOWN
    axB.plot([0, 1], [pre_v, post_v], '-', color=color_line, alpha=0.5, lw=1.2)
    axB.scatter([0, 1], [pre_v, post_v], c=[COLOR_PRE, COLOR_POST], s=80, zorder=5,
                edgecolors='black', linewidth=0.5)
    axB.annotate(pat, (1.02, post_v), fontsize=7, va='center')

axB.plot([0, 1], [np.mean(pre_prop), np.mean(post_prop)], 'k--', lw=2.5, alpha=0.7,
         label=f'Group mean (Δ={np.mean(delta_prop):+.4f})')
axB.scatter([0, 1], [np.mean(pre_prop), np.mean(post_prop)], c='black', s=150, marker='D',
            zorder=10, edgecolors='white', linewidth=1.5)

axB.set_xticks([0, 1])
axB.set_xticklabels(['Pre-treatment', 'Post-treatment'], fontsize=10)
axB.set_ylabel(f'NK-like high proportion (score > q75={NKLlike_THRESHOLD:.3f})', fontsize=10)
axB.set_xlim(-0.3, 1.5)
axB.set_title(f'Panel B: NK-like High Cell Proportion Dynamics\n'
              f'Wilcoxon p={p_w2:.4f} | ↑{int(np.sum(delta_prop>0))} / ↓{int(np.sum(delta_prop<0))}',
              fontsize=11)
axB.legend(loc='upper right', fontsize=8)
axB.grid(alpha=0.3, ls=':')
axB.spines['top'].set_visible(False); axB.spines['right'].set_visible(False)

# --- Panel C: CD8 cluster composition pre vs post (pooled) ---
axC = axes[1, 0]
cluster_pct = cluster_counts.div(cluster_counts.sum(axis=1), axis=0) * 100
clusters_order = ['Non-exhausted', 'Tex', 'Prolif.']
colors_cluster = ['#64B5F6', '#E53935', '#FFB300']

bar_width = 0.35
bottom_pre = 0
bottom_post = 0
for i, cl in enumerate(clusters_order):
    pre_pct = cluster_pct.loc['pre', cl]
    post_pct = cluster_pct.loc['post', cl]
    axC.bar(0, pre_pct, bar_width, bottom=bottom_pre, color=colors_cluster[i], label=cl, edgecolor='white')
    axC.bar(1, post_pct, bar_width, bottom=bottom_post, color=colors_cluster[i], edgecolor='white')
    if pre_pct > 3:
        axC.text(0, bottom_pre + pre_pct/2, f'{pre_pct:.1f}%', ha='center', va='center', fontsize=8, color='white', fontweight='bold')
    if post_pct > 3:
        axC.text(1, bottom_post + post_pct/2, f'{post_pct:.1f}%', ha='center', va='center', fontsize=8, color='white', fontweight='bold')
    bottom_pre += pre_pct
    bottom_post += post_pct

axC.set_xticks([0, 1])
axC.set_xticklabels([f'Pre-treatment\n(n={int(cluster_counts.loc["pre"].sum()):,} cells)',
                     f'Post-treatment\n(n={int(cluster_counts.loc["post"].sum()):,} cells)'], fontsize=10)
axC.set_ylabel('Cluster proportion (%)', fontsize=10)
axC.set_ylim(0, 100)
axC.set_title(f'Panel C: CD8+ Cluster Composition (pooled)\n'
              f'Chi-squared: χ²={chi2:.2f}, df={dof}, p={p_chi:.6f}', fontsize=11)
axC.legend(loc='upper right', fontsize=8, framealpha=0.9)
axC.spines['top'].set_visible(False); axC.spines['right'].set_visible(False)

# --- Panel D: Per-patient Tex proportion pre→post dynamics ---
axD = axes[1, 1]
for i, pat in enumerate(paired_pats):
    pre_v = pre_vals.loc[pat, 'Tex_prop']
    post_v = post_vals.loc[pat, 'Tex_prop']
    color_line = COLOR_UP if post_v > pre_v else COLOR_DOWN
    axD.plot([0, 1], [pre_v, post_v], '-', color=color_line, alpha=0.5, lw=1.2)
    axD.scatter([0, 1], [pre_v, post_v], c=[COLOR_PRE, COLOR_POST], s=80, zorder=5,
                edgecolors='black', linewidth=0.5)
    axD.annotate(pat, (1.02, post_v), fontsize=7, va='center')

axD.plot([0, 1], [np.mean(pre_tex), np.mean(post_tex)], 'k--', lw=2.5, alpha=0.7,
         label=f'Group mean (Δ={np.mean(delta_tex):+.4f})')
axD.scatter([0, 1], [np.mean(pre_tex), np.mean(post_tex)], c='black', s=150, marker='D',
            zorder=10, edgecolors='white', linewidth=1.5)

axD.set_xticks([0, 1])
axD.set_xticklabels(['Pre-treatment', 'Post-treatment'], fontsize=10)
axD.set_ylabel('Tex proportion (within CD8+ T cells)', fontsize=10)
axD.set_xlim(-0.3, 1.5)
axD.set_title(f'Panel D: Tex Proportion Dynamics\n'
              f'Wilcoxon p={p_w3:.4f} | ↑{int(np.sum(delta_tex>0))} / ↓{int(np.sum(delta_tex<0))}',
              fontsize=11)
axD.legend(loc='upper right', fontsize=8)
axD.grid(alpha=0.3, ls=':')
axD.spines['top'].set_visible(False); axD.spines['right'].set_visible(False)

# 添加整体标题与说明
fig.suptitle('Chemotherapy Impact on NK-like CD8+ T Cell Dynamics (GSE179994, NSCLC anti-PD-1+chemo)\n'
             f'8 paired patients (all Responders) | T-cell only dataset | Limitation: no myeloid cells for SPP1+ TAM analysis',
             fontsize=12, y=1.00)

plt.tight_layout()
fig.savefig(cfg.result_path('Fig_chemo_dynamics.pdf'), dpi=150, bbox_inches='tight')
fig.savefig(cfg.result_path('Fig_chemo_dynamics.png'), dpi=150, bbox_inches='tight')
plt.close(fig)
print(f"\n  Saved: Fig_chemo_dynamics.png/pdf")

# ============================================================
# 输出结论
# ============================================================
print("\n" + "="*70)
print("CHEMO DYNAMICS ANALYSIS SUMMARY")
print("="*70)
print(f"""
Dataset: GSE179994 (Liu et al. 2022, Nature Cancer)
- NSCLC, 150,849 T cells, 36 patients (T-cell only, no myeloid)

True cohort structure (verified from Supp Table 1):
- 25 Treatment-naïve (P2-9, P11, P12, P14-18, P20-28, P34): no treatment, no response label
- 8 Treatment A (P1,P10,P13,P19,P29,P30,P33,P35): Pembro+Carbo+Pemetrexed, paired pre+post, ALL Responder
- 3 Treatment B (P36,P37,P38): Pembro+Carbo+白蛋白结合紫杉醇, post-only, ALL Non-responder

This analysis uses ONLY the 8 Treatment A paired patients (all Responder):
{paired_pats}

HARD DATA LIMITATIONS (not extraction issues):
- NR group (P36/37/38, Treatment B) has NO pre-treatment samples → cannot do R vs NR paired dynamics
- 8 paired patients are all Responder → cannot extrapolate to NR
- Treatment A (R) vs Treatment B (NR) use different chemo regimens → confounded
- T-cell only → no SPP1+ TAM myeloid gate analysis (H4 skipped)

Key findings (ONLY in Responders, Treatment A):
1. NK-like signature score: Δmean = {np.mean(delta_score):+.4f}, Wilcoxon p = {p_w1:.4f}
   Direction: {n_up} up, {n_down} down ({'chemotherapy REDUCES' if np.mean(delta_score) < 0 else 'chemotherapy INCREASES'} NK-like score)
2. NK-like high cell proportion: Δmean = {np.mean(delta_prop):+.4f}, Wilcoxon p = {p_w2:.4f}
3. Tex proportion: Δmean = {np.mean(delta_tex):+.4f}, Wilcoxon p = {p_w3:.4f}
4. Cluster composition (pooled chi-squared): chi2={chi2:.2f}, p={p_chi:.6f}
   {'Significant' if p_chi < 0.05 else 'Non-significant'} reshaping of CD8 compartment

Interpretation (limited to Responders only):
- {'Chemotherapy significantly reduces NK-like signature, supporting H1 (direct killing)' if np.mean(delta_score) < 0 and p_w1 < 0.05 else 'No significant chemotherapy effect on NK-like signature (trend: Δ='+f'{np.mean(delta_score):+.4f})'}
- {'Cluster composition significantly reshaped, supporting H2 (differentiation pathway shift)' if p_chi < 0.05 else 'No significant cluster composition change'}
- CANNOT directly explain GSE179994 AUC=0.4722 (NR has no paired pre)
- Hypothesis (needs validation): chemo reshapes CD8 compartment, may reduce NK-like signature predictive power
""")

print("[step4b_chemo_dynamics] Done.")
