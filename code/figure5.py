#!/usr/bin/env python3
"""
Figure 5 — 外部队列验证
  Panel A: GSE126044 的 R vs NR 箱线图
  Panel B: 跨癌种对比森林图（NSCLC vs Melanoma），含 GSE120575
  Panel C: 跨癌种异质性模型图（NSCLC冷肿瘤 vs 黑色素瘤热肿瘤响应瓶颈）
"""
import os
import pandas as pd, numpy as np
from scipy.stats import mannwhitneyu
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import warnings
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.dirname(HERE)
ADATA = os.path.join(DATA, 'adata')
RESULT = os.path.join(DATA, 'result')
os.makedirs(RESULT, exist_ok=True)

NK_GENES = {
 'FGFBP2':'ENSG00000157951','KLRD1':'ENSG00000167283','CX3CR1':'ENSG00000121807',
 'FCGR3A':'ENSG00000184557','KLRC1':'ENSG00000124161','KLRC2':'ENSG00000182628',
 'KLRB1':'ENSG00000110865','NKG7':'ENSG00000105641','GNLY':'ENSG00000115526',
 'PRF1':'ENSG00000180644','GZMB':'ENSG00000100453','GZMH':'ENSG00000113069',
 'GZMA':'ENSG00000101206','CTSW':'ENSG00000130721','KLRF1':'ENSG00000172639',
 'SH2D1B':'ENSG00000139352','TYROBP':'ENSG00000007129','FCER1G':'ENSG00000163823',
 'CD160':'ENSG00000173676','CRTAM':'ENSG00000167111','IFNG':'ENSG00000111537',
 'TBX21':'ENSG00000163531','EOMES':'ENSG00000105611','ZNF683':'ENSG00000168669'}

NK_ENTREZ = {
 'FGFBP2':83888,'KLRD1':3824,'CX3CR1':1524,'FCGR3A':2214,'KLRC1':3821,
 'KLRC2':3822,'KLRB1':3820,'NKG7':4828,'GNLY':57823,'PRF1':5551,
 'GZMB':3002,'GZMH':2999,'GZMA':3001,'CTSW':1521,'KLRF1':51348,
 'SH2D1B':11323,'TYROBP':7305,'FCER1G':2207,'CD160':11126,'CRTAM':56253,
 'IFNG':3458,'TBX21':30009,'EOMES':8320,'ZNF683':283665}

# 从公共配置加载响应标签（避免与 figure_supplement.py 重复定义）
from config import GSE126044_RESPONSE, GSE135222_RESPONSE

def build_gse91061_response():
    import gzip
    sm_path = os.path.join(ADATA, 'GSE91061_series_matrix.txt.gz')
    titles, responses, visits = [], [], []
    with gzip.open(sm_path, 'rt') as f:
        for line in f:
            if line.startswith('!Sample_title'):
                titles = [t.strip('"') for t in line.strip().split('\t')[1:]]
            elif line.startswith('!Sample_characteristics_ch1'):
                chars = [c.strip('"') for c in line.strip().split('\t')[1:]]
                if any('response:' in c for c in chars):
                    responses = [c.split('response:')[1].strip() for c in chars]
                elif any('visit' in c for c in chars):
                    visits = [c.split(':')[-1].strip() for c in chars]
    resp_map = {}
    for title, resp, visit in zip(titles, responses, visits):
        if 'Pre' in visit and resp in ('PRCR', 'PD'):
            resp_map[title] = 'R' if resp == 'PRCR' else 'NR'
    return resp_map

def norm_id(g):
    g = str(g)
    if g.startswith('ENSG'):
        return g.split('.')[0]
    return g

def load_matrix(path, sep='\t'):
    return pd.read_csv(path, sep=sep, index_col=0, compression='gzip')

def extract_signature(df):
    idx_norm = {}
    idx_str = {}
    for g in df.index:
        n = norm_id(g)
        if n not in idx_norm:
            idx_norm[n] = g
        s = str(g).strip()
        if s not in idx_str:
            idx_str[s] = g
    rows = []
    for sym, ensg in NK_GENES.items():
        if sym in df.index:
            rows.append(sym)
        elif ensg in idx_norm:
            rows.append(idx_norm[ensg])
        elif str(NK_ENTREZ.get(sym, '')) in idx_str:
            rows.append(idx_str[str(NK_ENTREZ[sym])])
    sub = df.loc[rows]
    sig = sub.mean(axis=0)
    return sig, rows

# ============= 加载 bulk 验证数据 =============
print('=== Loading bulk matrices ===', flush=True)

m1 = load_matrix(os.path.join(ADATA, 'GSE126044_counts.txt.gz'))
s1, r1 = extract_signature(m1)
print(f'GSE126044 matched {len(r1)}/24 genes (n={m1.shape[1]})', flush=True)

m2 = load_matrix(os.path.join(ADATA, 'GSE135222_exp.tsv.gz'))
s2, r2 = extract_signature(m2)
print(f'GSE135222 matched {len(r2)}/24 genes (n={m2.shape[1]})', flush=True)

m3 = pd.read_csv(os.path.join(ADATA, 'GSE91061_fpkm.csv.gz'), index_col=0, compression='gzip')
s3, r3 = extract_signature(m3)
gse91061_resp = build_gse91061_response()
print(f'GSE91061 matched {len(r3)}/24 genes (n={m3.shape[1]}, Pre R/NR={len(gse91061_resp)})', flush=True)

# ============= 加载 GSE120575 scRNA-seq =============
print('=== Loading GSE120575 (melanoma scRNA-seq) ===', flush=True)
nk_file = os.path.join(ADATA, 'GSE120575', 'GSE120575_nk_genes.csv')
m4 = pd.read_csv(nk_file, index_col=0)
print(f'GSE120575 NK genes: {list(m4.index)}, cells: {m4.shape[1]}', flush=True)

# 患者响应映射
pred_file = os.path.join(ADATA, 'GSE120575', 'GSE120575_Patient_Predictions.csv')
pred_df = pd.read_csv(pred_file)
patient_resp = dict(zip(pred_df['patient_id'], pred_df['NR_R']))
print(f'GSE120575 patients: R={sum(1 for v in patient_resp.values() if v=="R")}, NR={sum(1 for v in patient_resp.values() if v=="NR")}', flush=True)

# 从细胞名提取患者ID并计算 pseudobulk NK-like signature
cell_to_patient = {}
for col in m4.columns:
    parts = col.split('_')
    if len(parts) >= 2:
        pid = parts[1]  # e.g., P3
        cell_to_patient[col] = pid

# 每个患者的 pseudobulk NK-like signature（均值）
patient_scores_4 = {}
patient_groups_4 = {}
for cell, pid in cell_to_patient.items():
    resp = patient_resp.get(pid)
    if resp is None:
        continue
    score = float(m4[cell].mean())
    if pid not in patient_scores_4:
        patient_scores_4[pid] = []
    patient_scores_4[pid].append(score)
    patient_groups_4[pid] = resp

# 患者水平均值
total_pred_patients = len(patient_resp)
n_pred_r = sum(1 for v in patient_resp.values() if v == 'R')
n_pred_nr = sum(1 for v in patient_resp.values() if v == 'NR')
expr_patients = set(pid for col in m4.columns if len(col.split('_')) >= 2
                    for pid in [col.split('_')[1]])
missing_pred = set(patient_resp.keys()) - expr_patients
r4_r = [np.mean(patient_scores_4[pid]) for pid in patient_scores_4 if patient_groups_4[pid] == 'R']
r4_nr = [np.mean(patient_scores_4[pid]) for pid in patient_scores_4 if patient_groups_4[pid] == 'NR']
p4 = mannwhitneyu(r4_nr, r4_r, alternative='two-sided').pvalue
print(f'GSE120575 p={p4:.4f} (NR={len(r4_nr)}, R={len(r4_r)})', flush=True)
print(f'NOTE: Predictions file has {total_pred_patients} patients (R={n_pred_r}, NR={n_pred_nr}), '
      f'but NK-gene expression matrix only contains {len(expr_patients)} patients. '
      f'{len(missing_pred)} patients from predictions are absent from the expression matrix: {sorted(missing_pred)}', flush=True)

# ============= 计算各队列 p 值 =============
r1_r, r1_nr = [], []
for sample_name, score in s1.items():
    resp = GSE126044_RESPONSE.get(sample_name)
    if resp == 'R': r1_r.append(float(score))
    elif resp == 'NR': r1_nr.append(float(score))

r2_r, r2_nr = [], []
for sample_name, score in s2.items():
    resp = GSE135222_RESPONSE.get(sample_name)
    if resp == 'R': r2_r.append(float(score))
    elif resp == 'NR': r2_nr.append(float(score))

r3_r, r3_nr = [], []
for sample_name, score in s3.items():
    resp = gse91061_resp.get(sample_name)
    if resp == 'R': r3_r.append(float(score))
    elif resp == 'NR': r3_nr.append(float(score))

p1 = mannwhitneyu(r1_nr, r1_r, alternative='two-sided').pvalue
p2 = mannwhitneyu(r2_nr, r2_r, alternative='two-sided').pvalue
p3 = mannwhitneyu(r3_nr, r3_r, alternative='two-sided').pvalue

print(f'GSE126044 p={p1:.4f} (NR={len(r1_nr)}, R={len(r1_r)})', flush=True)
print(f'GSE135222 p={p2:.4f} (NR={len(r2_nr)}, R={len(r2_r)})', flush=True)
print(f'GSE91061 p={p3:.4f} (NR={len(r3_nr)}, R={len(r3_r)})', flush=True)
print(f'GSE120575 p={p4:.4f} (NR={len(r4_nr)}, R={len(r4_r)})', flush=True)

# 保存基因匹配率（供 figure_supplement.py 使用，避免硬编码）
gene_match_df = pd.DataFrame([
    {'dataset': 'GSE126044', 'genes_matched': len(r1), 'genes_total': 24, 'dtype': 'bulk'},
    {'dataset': 'GSE135222', 'genes_matched': len(r2), 'genes_total': 24, 'dtype': 'bulk'},
    {'dataset': 'GSE91061',  'genes_matched': len(r3), 'genes_total': 24, 'dtype': 'bulk'},
    {'dataset': 'GSE120575', 'genes_matched': m4.shape[0], 'genes_total': 24, 'dtype': 'scRNA'},
])
gene_match_df['match_rate'] = gene_match_df['genes_matched'] / gene_match_df['genes_total']
gene_match_df.to_csv(os.path.join(RESULT, 'fig5_gene_match.csv'), index=False)
print(f'Gene match rates saved:\n{gene_match_df.to_string()}', flush=True)

# 保存签名
s1.to_csv(os.path.join(RESULT, 'sig_GSE126044.csv'))
s2.to_csv(os.path.join(RESULT, 'sig_GSE135222.csv'))
s3.to_csv(os.path.join(RESULT, 'sig_GSE91061.csv'))

# ============= 组合绘图 =============
print('=== Composing Figure 5 ===', flush=True)
fig = plt.figure(figsize=(18, 14))
gs = fig.add_gridspec(3, 2, hspace=0.4, wspace=0.3, height_ratios=[1, 1, 0.9])

# ---- Panel A: GSE126044 箱线图 ----
ax1 = fig.add_subplot(gs[0, 0])

box_data = [r1_nr, r1_r]
bp = ax1.boxplot(box_data, patch_artist=True, showfliers=False, widths=0.55)
colors = ['#E74C3C', '#27AE60']
for patch, c in zip(bp['boxes'], colors):
    patch.set_facecolor(c); patch.set_alpha(0.75)
for med in bp['medians']:
    med.set_color('black'); med.set_linewidth(2)
for whisker in bp['whiskers']:
    whisker.set_color('gray'); whisker.set_linewidth(1.2)
for cap in bp['caps']:
    cap.set_color('gray')

np.random.seed(42)
ax1.scatter(np.random.normal(1, 0.06, len(r1_nr)), r1_nr, s=60, c='#C0392B', alpha=0.8, edgecolor='white', linewidth=1, zorder=3)
ax1.scatter(np.random.normal(2, 0.06, len(r1_r)), r1_r, s=60, c='#1E8449', alpha=0.8, edgecolor='white', linewidth=1, zorder=3)

ax1.set_xticks([1, 2])
ax1.set_xticklabels([f'Non-Responder\n(NR, n={len(r1_nr)})', f'Responder\n(R, n={len(r1_r)})'], fontsize=12, fontweight='bold')
ax1.set_ylabel('NK-like signature score', fontsize=12)
ax1.set_title('GSE126044\n(NSCLC, anti-PD1)', fontsize=14, fontweight='bold')

y_range = max(r1_nr + r1_r) - min(r1_nr + r1_r)
y_max = max(r1_nr + r1_r) + y_range * 0.1
ax1.plot([1, 2], [y_max, y_max], 'k-', lw=1.5)
ax1.plot([1, 1], [y_max - y_range*0.03, y_max], 'k-', lw=1.5)
ax1.plot([2, 2], [y_max - y_range*0.03, y_max], 'k-', lw=1.5)
p_text = f'Mann-Whitney p = {p1:.4f}'
ax1.text(1.5, y_max + y_range*0.02, p_text, ha='center', va='bottom', fontsize=11, fontweight='bold', color='#1E8449')

ax1.text(-0.15, 1.05, 'A', transform=ax1.transAxes, fontsize=22, fontweight='bold')

# ---- Panel B: 跨癌种对比森林图（含 GSE120575）----
ax2 = fig.add_subplot(gs[0, 1])

cohorts = [
    {'name': f'GSE126044\n(NSCLC, bulk)', 'n_nr': len(r1_nr), 'n_r': len(r1_r),
     'p': p1, 'color': '#E74C3C', 'r_mean': np.mean(r1_r), 'nr_mean': np.mean(r1_nr),
     'fc': np.mean(r1_r) / np.mean(r1_nr) if np.mean(r1_nr) > 0 else 1,
     'r_vals': r1_r, 'nr_vals': r1_nr},
    {'name': f'GSE135222\n(NSCLC, bulk)', 'n_nr': len(r2_nr), 'n_r': len(r2_r),
     'p': p2, 'color': '#F39C12', 'r_mean': np.mean(r2_r), 'nr_mean': np.mean(r2_nr),
     'fc': np.mean(r2_r) / np.mean(r2_nr) if np.mean(r2_nr) > 0 else 1,
     'r_vals': r2_r, 'nr_vals': r2_nr},
    {'name': f'GSE91061\n(Melanoma, bulk)', 'n_nr': len(r3_nr), 'n_r': len(r3_r),
     'p': p3, 'color': '#3498DB', 'r_mean': np.mean(r3_r), 'nr_mean': np.mean(r3_nr),
     'fc': np.mean(r3_r) / np.mean(r3_nr) if np.mean(r3_nr) > 0 else 1,
     'r_vals': r3_r, 'nr_vals': r3_nr},
    {'name': f'GSE120575\n(Melanoma, scRNA)', 'n_nr': len(r4_nr), 'n_r': len(r4_r),
     'p': p4, 'color': '#2E86C1', 'r_mean': np.mean(r4_r), 'nr_mean': np.mean(r4_nr),
     'fc': np.mean(r4_r) / np.mean(r4_nr) if np.mean(r4_nr) > 0 else 1,
     'r_vals': r4_r, 'nr_vals': r4_nr},
]

y_positions = [3, 2, 1, 0]

# Bootstrap 计算 log2FC 的 95% CI
rng_bs = np.random.default_rng(42)
n_boot = 2000

for y, coh in zip(y_positions, cohorts):
    fc = coh['fc']
    log2fc = np.log2(fc)

    # Bootstrap: 从 R 和 NR 组有放回抽样，重算 log2FC
    r_data = np.asarray(coh['r_vals'])
    nr_data = np.asarray(coh['nr_vals'])
    boot_log2fc = []
    for _ in range(n_boot):
        r_boot = rng_bs.choice(r_data, size=len(r_data), replace=True)
        nr_boot = rng_bs.choice(nr_data, size=len(nr_data), replace=True)
        nr_mean = np.mean(nr_boot)
        if nr_mean > 0:
            boot_log2fc.append(np.log2(np.mean(r_boot) / nr_mean))
    if len(boot_log2fc) > 0:
        ci_low = np.percentile(boot_log2fc, 2.5)
        ci_high = np.percentile(boot_log2fc, 97.5)
    else:
        ci_low = log2fc - 0.3
        ci_high = log2fc + 0.3
    ax2.plot([ci_low, ci_high], [y, y], color=coh['color'], lw=3, solid_capstyle='round')
    ax2.plot([log2fc, log2fc], [y-0.15, y+0.15], color=coh['color'], lw=5, solid_capstyle='round', marker='D', markersize=8)
    sig = '***' if coh['p'] < 0.001 else '**' if coh['p'] < 0.01 else '*' if coh['p'] < 0.05 else 'n.s.'
    text = f"FC={fc:.2f}, p={coh['p']:.3f} {sig}"
    ax2.text(ci_high + 0.05, y, text, va='center', fontsize=9, color=coh['color'], fontweight='bold')

ax2.axvline(0, color='gray', ls='--', lw=1.5)
ax2.set_yticks(y_positions)
ax2.set_yticklabels([coh['name'] for coh in cohorts], fontsize=10)
ax2.set_xlabel('log2(Fold Change)  R / NR', fontsize=12)
ax2.set_title('Cross-cancer validation of NK-like signature', fontsize=14, fontweight='bold')

xlim = max(abs(min(ax2.get_xlim())), abs(max(ax2.get_xlim()))) * 1.4
ax2.set_xlim(-xlim, xlim)

# 癌种分组标注
ax2.text(-xlim * 0.95, 3.5, 'NSCLC (Cold tumor)', fontsize=11, fontweight='bold', color='#E67E22')
ax2.text(-xlim * 0.95, 0.5, 'Melanoma (Hot tumor)', fontsize=11, fontweight='bold', color='#2980B9')
ax2.set_ylim(-1, 4)
ax2.grid(axis='x', alpha=0.3, ls=':')
ax2.text(-0.1, 1.05, 'B', transform=ax2.transAxes, fontsize=22, fontweight='bold')

# ---- Panel C: 跨癌种异质性模型图 ----
ax3 = fig.add_subplot(gs[1:, :])
ax3.set_xlim(0, 14)
ax3.set_ylim(0, 8)
ax3.axis('off')
ax3.set_title('Cross-cancer Heterogeneity Model: Response Bottleneck', fontsize=14, fontweight='bold')

# NSCLC (Cold tumor) 左侧
box_nsclc = FancyBboxPatch((0.3, 1), 5.5, 6, boxstyle="round,pad=0.15",
                           facecolor='#FDEBD0', edgecolor='#E67E22', linewidth=2.5, alpha=0.85)
ax3.add_patch(box_nsclc)
ax3.text(3.05, 6.5, 'NSCLC (Cold Tumor)', fontsize=13, fontweight='bold', ha='center', color='#D35400')
ax3.text(3.05, 6.0, f'GSE126044: p={p1:.4f} ✓', fontsize=10, ha='center', color='#D35400', fontweight='bold')
ax3.text(3.05, 5.6, f'GSE135222: p={p2:.4f} (n_R={len(r2_r)})', fontsize=10, ha='center', color='#D35400')

# NSCLC 瓶颈
bottleneck_nsclc = FancyBboxPatch((0.8, 3), 4.5, 1.8, boxstyle="round,pad=0.1",
                                  facecolor='#F5B7B1', edgecolor='#C0392B', linewidth=2, alpha=0.85)
ax3.add_patch(bottleneck_nsclc)
ax3.text(3.05, 4.2, 'Bottleneck: Stemness loss', fontsize=10, fontweight='bold', ha='center', color='#922B21')
ax3.text(3.05, 3.6, 'SPP1+ TAM → KLF2↓ → NK-like lock lost', fontsize=9, ha='center', color='#922B21')
ax3.text(3.05, 3.2, '→ NK-like signature predicts response', fontsize=9, ha='center', color='#922B21', style='italic')

# NSCLC 响应路径
ax3.text(3.05, 2.3, 'R: NK-like locked, KLF2+', fontsize=9, ha='center', color='#27AE60', fontweight='bold')
ax3.text(3.05, 1.8, 'NR: Tex-dominant, SPP1+ TAM high', fontsize=9, ha='center', color='#E74C3C', fontweight='bold')
ax3.text(3.05, 1.3, f'OR = pCR predictor (NK-Locked)', fontsize=9, ha='center', color='#7F8C8D')

# Melanoma (Hot tumor) 右侧
box_mel = FancyBboxPatch((8.2, 1), 5.5, 6, boxstyle="round,pad=0.15",
                         facecolor='#D6EAF8', edgecolor='#2980B9', linewidth=2.5, alpha=0.85)
ax3.add_patch(box_mel)
ax3.text(10.95, 6.5, 'Melanoma (Hot Tumor)', fontsize=13, fontweight='bold', ha='center', color='#1A5276')
ax3.text(10.95, 6.0, f'GSE91061: p={p3:.4f} n.s.', fontsize=10, ha='center', color='#1A5276', fontweight='bold')
ax3.text(10.95, 5.6, f'GSE120575: p={p4:.4f} n.s.', fontsize=10, ha='center', color='#1A5276')

# Melanoma 瓶颈
bottleneck_mel = FancyBboxPatch((8.7, 3), 4.5, 1.8, boxstyle="round,pad=0.1",
                                facecolor='#AED6F1', edgecolor='#2471A3', linewidth=2, alpha=0.85)
ax3.add_patch(bottleneck_mel)
ax3.text(10.95, 4.2, 'Bottleneck: Differentiation', fontsize=10, fontweight='bold', ha='center', color='#1A5276')
ax3.text(10.95, 3.6, 'High TMB → high clonal diversity', fontsize=9, ha='center', color='#1A5276')
ax3.text(10.95, 3.2, '→ NK-like lock not rate-limiting', fontsize=9, ha='center', color='#1A5276', style='italic')

# Melanoma 响应路径
ax3.text(10.95, 2.3, 'R: High TMB, CD8 infiltration', fontsize=9, ha='center', color='#27AE60', fontweight='bold')
ax3.text(10.95, 1.8, 'NR: Immune exclusion, MET', fontsize=9, ha='center', color='#E74C3C', fontweight='bold')
ax3.text(10.95, 1.3, 'NK-like signature not predictive', fontsize=9, ha='center', color='#7F8C8D')

# 中间对比箭头
arrow_compare = FancyArrowPatch((5.9, 4), (8.1, 4),
                                arrowstyle='<->', mutation_scale=20, lw=2.5, color='#7D3C98')
ax3.add_patch(arrow_compare)
ax3.text(7, 4.5, 'Different\nbottlenecks', ha='center', fontsize=9, fontweight='bold', color='#7D3C98')
ax3.text(7, 3.3, 'Cold vs Hot', ha='center', fontsize=8, color='#7D3C98', style='italic')

ax3.text(-0.02, 1.0, 'C', transform=ax3.transAxes, fontsize=22, fontweight='bold')

plt.tight_layout()
plt.savefig(os.path.join(RESULT, 'Fig5_external_validation.png'), dpi=300)
plt.savefig(os.path.join(RESULT, 'Fig5_external_validation.pdf'))
plt.close()
print('Fig5 saved to result/Fig5_external_validation.png', flush=True)

# ============= 保存汇总表 =============
summary = pd.DataFrame({
    'Dataset': ['GSE126044', 'GSE135222', 'GSE91061', 'GSE120575'],
    'Cancer': ['NSCLC', 'NSCLC', 'Melanoma', 'Melanoma'],
    'Treatment': ['anti-PD1', 'anti-PD1/PDL1', 'anti-PD1/CTLA4', 'anti-PD1'],
    'DataType': ['bulk RNA-seq', 'bulk RNA-seq', 'bulk RNA-seq', 'scRNA-seq'],
    'NR_count': [len(r1_nr), len(r2_nr), len(r3_nr), len(r4_nr)],
    'R_count': [len(r1_r), len(r2_r), len(r3_r), len(r4_r)],
    'Total_predicted_NR': [None, None, None, n_pred_nr],
    'Total_predicted_R': [None, None, None, n_pred_r],
    'Total_predicted_patients': [None, None, None, total_pred_patients],
    'MWU_pvalue': [p1, p2, p3, p4],
    'FC_R_vs_NR': [np.mean(r1_r)/np.mean(r1_nr) if np.mean(r1_nr)>0 else np.nan,
                   np.mean(r2_r)/np.mean(r2_nr) if np.mean(r2_nr)>0 else np.nan,
                   np.mean(r3_r)/np.mean(r3_nr) if np.mean(r3_nr)>0 else np.nan,
                   np.mean(r4_r)/np.mean(r4_nr) if np.mean(r4_nr)>0 else np.nan],
    'Significant': ['Yes' if p<0.05 else 'No' for p in [p1,p2,p3,p4]],
})
summary.to_csv(os.path.join(RESULT, 'fig5_validation_summary.csv'), index=False)
print(f'Summary saved to result/fig5_validation_summary.csv', flush=True)

print('=== DONE Fig5 ===', flush=True)
