#!/usr/bin/env python3
"""
T05 Step 2 + Step 3:
  Step 2 — Bootstrap 内部验证 (1000 resamples, Model C 多变量logistic回归)
  Step 3 — 阈值敏感性分析 (克隆细胞数阈值 [3,5,10,15,20,30,50] × NK-like占比阈值 [0.3,0.4,0.5,0.6,0.7])
"""
import os, sys
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.stats as ss
import statsmodels.formula.api as smf
from sklearn.metrics import roc_auc_score
from sklearn.utils import resample
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
DATA = os.path.dirname(HERE)
ADATA = os.path.join(DATA, 'adata')
RESULT = os.path.join(DATA, 'result')
os.makedirs(RESULT, exist_ok=True)

# ============================================================
# Load h5ad & prepare per-patient frame (replicates figure1.py's logic)
# ============================================================
print('=== Loading GSE243013_T_cells.h5ad ===', flush=True)
adata = sc.read_h5ad(os.path.join(ADATA, 'GSE243013_T_cells.h5ad'))
print(f'n_cells = {adata.n_obs:,}, n_genes = {adata.n_vars:,}', flush=True)

for c in ['sub_cell_type', 'sampleID', 'pathological_response',
          'clonotype', 'clonotype_number', 'major_cell_type']:
    adata.obs[c] = adata.obs[c].astype(str)

adata.obs['is_cd8'] = adata.obs['sub_cell_type'].str.startswith('CD8T')
adata.obs['is_nklike'] = adata.obs['sub_cell_type'] == 'CD8T_NK-like_FGFBP2'
adata.obs['is_tex'] = adata.obs['sub_cell_type'].str.contains('Tex', case=False, na=False)


def compute_per_patient(ad, clone_count_thr: int, nk_frac_thr: float = 0.5):
    """按指定阈值重新计算 nk_locked 标签（clone_count_thr 控制大克隆定义，
    nk_frac_thr 控制 NK-like 占比阈值判定"锁定"）。
    返回 DataFrame: patient / nk_dominant_ratio / nk_locked / cd8t / total_clone_count / is_LUAD
    """
    rows = []
    for pid, sub in ad.obs.groupby('sampleID'):
        resp = sub['pathological_response'].iloc[0]
        if resp not in ('pCR', 'non-MPR'):
            continue
        cdm = sub['is_cd8']
        cd8n = int(cdm.sum())
        nkn = int((cdm & sub['is_nklike']).sum())
        if cd8n == 0:
            continue
        cd8p = sub[cdm]
        cc = cd8p['clonotype'].value_counts()
        big = cc[cc >= clone_count_thr]
        fracs = []
        for cl in big.index:
            clc = cd8p[cd8p['clonotype'] == cl]
            nk_in_cl = float(clc['is_nklike'].sum())
            fracs.append(nk_in_cl / len(clc))
        nk_dom = float(np.mean(fracs)) if fracs else np.nan
        nk_locked = int(nk_dom > nk_frac_thr) if not np.isnan(nk_dom) else 0
        ct = int(cc.shape[0])
        # Cancer type (供 Model C 校正)
        ct_label = sub.get('cancer_type')
        if ct_label is None:
            ct_val = 'LUAD'  # fallback
        else:
            ct_val = str(ct_label.iloc[0]) if hasattr(ct_label, 'iloc') else str(ct_label)
        rows.append(dict(
            patient=str(pid),
            response=resp,
            resp_bin=int(resp == 'pCR'),
            cd8t=cd8n,
            nklike=nkn,
            nk_pct=nkn / cd8n,
            total_clone_count=ct,
            nk_dominant_ratio=nk_dom,
            nk_locked=nk_locked,
            is_LUAD=int(ct_val == 'LUAD'),
        ))
    return pd.DataFrame(rows)


# 以发表版阈值 (clone=5, nk_frac=0.5) 作为基准
print('=== Building per-patient frame (baseline thresholds 5 / 0.5) ===', flush=True)
base = compute_per_patient(adata, clone_count_thr=5, nk_frac_thr=0.5)
print(f'n_patients = {len(base)}  pCR={int(base.resp_bin.sum())}  non-MPR={int((1-base.resp_bin).sum())}', flush=True)

# 真实 Model C 多变量 logistic 回归
print('=== Fitting Model C (real data) ===', flush=True)
m_real = smf.logit(
    'resp_bin ~ nk_locked + total_clone_count + cd8t + is_LUAD',
    data=base
).fit(disp=0)
auc_real = roc_auc_score(base['resp_bin'], m_real.predict(base))
or_real = float(np.exp(m_real.params['nk_locked']))
ci_real = np.exp(m_real.conf_int().loc['nk_locked'])
p_real = float(m_real.pvalues['nk_locked'])
print(f'Real: OR={or_real:.3f} [{ci_real[0]:.3f},{ci_real[1]:.3f}]  p={p_real:.3e}  AUC={auc_real:.3f}', flush=True)

# ============================================================
# Step 2: Bootstrap internal validation (1000 resamples)
# ============================================================
print('=== Step 2: Bootstrap 1000 resamples ===', flush=True)
N_BOOT = 1000
rng = np.random.RandomState(42)
boot_or = np.zeros(N_BOOT)
boot_auc = np.zeros(N_BOOT)
n_success = 0
for i in range(N_BOOT):
    # 放回采样 n=len(base) 个患者
    idx = rng.randint(0, len(base), size=len(base))
    bdf = base.iloc[idx].reset_index(drop=True)
    # 边界检查：若 bootstrap 样本全是 pCR 或全 non-MPR，跳过
    if bdf['resp_bin'].nunique() < 2:
        continue
    try:
        m_b = smf.logit(
            'resp_bin ~ nk_locked + total_clone_count + cd8t + is_LUAD',
            data=bdf
        ).fit(disp=0, maxiter=100)
        if np.isnan(m_b.params['nk_locked']):
            continue
        boot_or[i] = np.exp(m_b.params['nk_locked'])
        # AUC 真实值 + 抖动 (Mann-Whitney U 转换等价)
        if bdf['resp_bin'].nunique() == 2:
            boot_auc[i] = roc_auc_score(bdf['resp_bin'], m_b.predict(bdf))
        n_success += 1
    except Exception:
        continue

boot_or = boot_or[boot_or != 0]
boot_auc = boot_auc[boot_auc != 0]
# Use median + percentiles (mean is dominated by quasi-separation outliers)
or_mean = float(np.mean(boot_or[boot_or < 1e6]))  # trimmed mean (exclude extreme outliers)
or_med = float(np.median(boot_or[boot_or < 1e6]))
or_lo = float(np.percentile(boot_or[boot_or < 1e6], 2.5))
or_hi = float(np.percentile(boot_or[boot_or < 1e6], 97.5))
auc_mean = float(np.mean(boot_auc))
auc_med = float(np.median(boot_auc))
auc_lo = float(np.percentile(boot_auc, 2.5))
auc_hi = float(np.percentile(boot_auc, 97.5))
# 偏差校正 (Bootstrap bias-corrected)
auc_bc = 2 * auc_real - auc_mean  # 简单偏差校正
print(f'Bootstrap: n_success={n_success}/{N_BOOT}', flush=True)
print(f'  OR trimmed-mean={or_mean:.3f}, median={or_med:.3f}  95%CI [{or_lo:.3f}, {or_hi:.3f}]', flush=True)
print(f'  AUC mean={auc_mean:.3f}, median={auc_med:.3f}  95%CI [{auc_lo:.3f}, {auc_hi:.3f}]  bias-corrected={auc_bc:.3f}', flush=True)

# 保存
boot_csv = os.path.join(RESULT, 'T05_05_bootstrap_auc.csv')
pd.DataFrame([dict(
    metric='Model_C_nk_locked',
    n_patients=len(base),
    n_bootstrap=N_BOOT,
    n_success=int(n_success),
    or_real=or_real,
    or_lo_real=float(ci_real[0]),
    or_hi_real=float(ci_real[1]),
    p_real=p_real,
    auc_real=auc_real,
    or_boot_trimmed_mean=or_mean,
    or_boot_median=or_med,
    or_boot_lo=or_lo,
    or_boot_hi=or_hi,
    auc_boot_mean=auc_mean,
    auc_boot_median=auc_med,
    auc_boot_lo=auc_lo,
    auc_boot_hi=auc_hi,
    auc_bc=auc_bc,
)]).to_csv(boot_csv, index=False)
print(f'Saved -> {boot_csv}', flush=True)

# ============================================================
# Step 3: Threshold sensitivity analysis
# ============================================================
print('=== Step 3: Threshold sensitivity (clone × nk_frac) ===', flush=True)
clone_thrs = [3, 5, 10, 15, 20, 30, 50]
nk_frac_thrs = [0.3, 0.4, 0.5, 0.6, 0.7]
sens_rows = []
for ct_thr in clone_thrs:
    for nk_thr in nk_frac_thrs:
        df = compute_per_patient(adata, clone_count_thr=ct_thr, nk_frac_thr=nk_thr)
        n_pat = len(df)
        if df['resp_bin'].nunique() < 2 or df['nk_locked'].nunique() < 2:
            sens_rows.append(dict(
                clone_threshold=ct_thr, nk_frac_threshold=nk_thr,
                n_patients=n_pat, n_locked=int(df['nk_locked'].sum()),
                or_value=np.nan, or_lo=np.nan, or_hi=np.nan, p_value=np.nan,
                note='Insufficient variation',
            ))
            continue
        try:
            m = smf.logit(
                'resp_bin ~ nk_locked + total_clone_count + cd8t + is_LUAD',
                data=df
            ).fit(disp=0, maxiter=100)
            or_v = float(np.exp(m.params['nk_locked']))
            ci = np.exp(m.conf_int().loc['nk_locked'])
            p_v = float(m.pvalues['nk_locked'])
            sens_rows.append(dict(
                clone_threshold=ct_thr, nk_frac_threshold=nk_thr,
                n_patients=n_pat, n_locked=int(df['nk_locked'].sum()),
                or_value=or_v, or_lo=float(ci[0]), or_hi=float(ci[1]), p_value=p_v,
                note='OK',
            ))
        except Exception as e:
            sens_rows.append(dict(
                clone_threshold=ct_thr, nk_frac_threshold=nk_thr,
                n_patients=n_pat, n_locked=int(df['nk_locked'].sum()),
                or_value=np.nan, or_lo=np.nan, or_hi=np.nan, p_value=np.nan,
                note=f'Fit failed: {e}',
            ))
        last = sens_rows[-1]
        print(f'  clone={ct_thr}, nk_frac={nk_thr}: OR={last["or_value"]}, p={last["p_value"]}', flush=True)

sens_df = pd.DataFrame(sens_rows)
sens_csv = os.path.join(RESULT, 'T05_threshold_sensitivity.csv')
sens_df.to_csv(sens_csv, index=False)
print(f'Saved -> {sens_csv}', flush=True)

# Threshold curve plot: 1 row x 5 subplots (one per nk_frac)
fig, axes = plt.subplots(1, 5, figsize=(22, 4.5), sharey=True)
colors = ['#E74C3C', '#F39C12', '#27AE60', '#2E86C1', '#8E44AD']
# Cap OR values for display: clip at 1e3 for visual clarity
OR_CAP = 100.0
sens_plot = sens_df.copy()
for col in ['or_value', 'or_lo', 'or_hi']:
    sens_plot[col] = sens_plot[col].clip(upper=OR_CAP)
for i, nk_thr in enumerate(nk_frac_thrs):
    ax = axes[i]
    sub = sens_plot[sens_plot['nk_frac_threshold'] == nk_thr].sort_values('clone_threshold')
    # Use NaN mask
    valid = sub['or_value'].notna()
    ax.errorbar(sub.loc[valid, 'clone_threshold'], sub.loc[valid, 'or_value'],
                yerr=[sub.loc[valid, 'or_value'] - sub.loc[valid, 'or_lo'],
                      sub.loc[valid, 'or_hi'] - sub.loc[valid, 'or_value']],
                fmt='o-', capsize=5, capthick=2, linewidth=2, markersize=8,
                color=colors[i % len(colors)], label=f'NK-frac={nk_thr}')
    ax.axhline(1, color='gray', ls='--', lw=1, alpha=0.6)
    ax.set_xscale('log')
    ax.set_xlabel('Clone cell-number threshold', fontsize=11)
    if i == 0:
        ax.set_ylabel(f'Odds Ratio (clipped at {OR_CAP})', fontsize=11)
    ax.set_title(f'NK-frac threshold = {nk_thr}', fontsize=12, fontweight='bold')
    ax.grid(alpha=0.3, ls=':')
    ax.legend(loc='upper right', fontsize=9)
plt.suptitle('Threshold Sensitivity Analysis: OR vs Clone Threshold (colored by NK-frac)',
             fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
curve_png = os.path.join(RESULT, 't05_threshold_curve.png')
curve_pdf = os.path.join(RESULT, 't05_threshold_curve.pdf')
plt.savefig(curve_png, dpi=300, bbox_inches='tight')
plt.savefig(curve_pdf, bbox_inches='tight')
plt.close()
print(f'Saved -> {curve_png}', flush=True)

print('=== DONE T05 Step 2 + 3 ===', flush=True)
