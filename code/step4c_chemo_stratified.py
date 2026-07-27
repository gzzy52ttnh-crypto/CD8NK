"""
step4c_chemo_stratified.py
主队列 GSE243013 内部化疗方案分层分析

目的: 直接回答"不同化疗方案如何影响 NK-like CD8+ T 细胞克隆命运锁定对 pCR 的预测效力"

重大发现 (核对自 GSE243013 h5ad):
- 主队列 188 例 pCR vs non-MPR 患者中, 175 例接受 anti-PD1+化疗联合方案
- 之前 "主队列是 anti-PD1 单药" 的理解有误, OR=19.57 实际是 "化疗背景下" 的预测效力
- 主队列内化疗方案异质性可支持分层分析, 无需依赖 GSE179994 (NR 无配对的硬限制)

化疗大类分组 (n=188):
- Platinum+Taxane (Carbo/Cis + Abraxane/Paclitaxel): 123 例 (pCR=65, NR=58)
- Platinum+Pemetrexed (Carbo/Cis + Pemetrexed): 35 例 (pCR=8, NR=27) ← 与 GSE179994 Treatment A 同方案
- No chemo (单纯 anti-PD1): 13 例 (pCR=4, NR=9)
- Platinum+Gemcitabine: 6 例 (pCR=3, NR=3, excluded due to n<10)
- Other: 10 例 (pCR=4, NR=6)

分析内容:
1. 各亚组 Firth 惩罚 logistic 回归 (pCR ~ nk_dominant_ratio)
2. 跨亚组异质性检验 (Cochran's Q + I²)
3. 森林图展示各亚组 OR + 95% CI
4. 主队列 Carbo+Pemetrexed 亚组 vs GSE179994 Treatment A 跨队列对比
5. 化疗方案对 NK-dominant ratio 水平的影响 (Kruskal-Wallis)

输出:
- Fig_chemo_stratified.png/pdf (4 panels: 森林图 + NK-dominant分布 + pCR率 + 跨队列对比)
- chemo_stratified_results.csv
- chemo_stratified_patient_metrics.csv
"""
import config as cfg
from _common import firth_logistic_fit
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import chi2, kruskal, mannwhitneyu, fisher_exact
from scipy.optimize import minimize
import statsmodels.api as sm
import os
import warnings
warnings.filterwarnings('ignore')

# 复用 step2_firth_validation.py 的 Firth 回归实现
import importlib.util
firth_spec = importlib.util.spec_from_file_location(
    "step2_firth_validation",
    os.path.join(os.path.dirname(__file__), "step2_firth_validation.py"))
# 不导入整个模块(会执行分析), 只复制函数定义
# 直接内联 Firth 函数 (与 step2_firth_validation.py 保持一致)


def firth_profile_ci(X, y, beta_mle, param_idx=1, alpha=0.05, n_points=30):
    """Profile likelihood CI for Firth regression (简化版, 与 step2 一致)"""
    n, p = X.shape
    chi2_crit = chi2.ppf(1 - alpha, 1)
    # Full model log-likelihood (without Firth penalty, 用普通 MLE 近似)
    eta = X @ beta_mle
    p_i = np.clip(1.0 / (1.0 + np.exp(-eta)), 1e-10, 1 - 1e-10)
    ll_full = np.sum(y * np.log(p_i) + (1 - y) * np.log(1 - p_i))

    def neg_ll_fixed(fixed_val):
        X_reduced = X.copy()
        # 固定 param_idx 列为 fixed_val * 1, 优化截距
        # 构造: y = intercept + fixed_val * x, 只优化 intercept
        offset = fixed_val * X[:, param_idx]
        # Newton-Raphson on intercept only
        b0 = 0.0
        for _ in range(50):
            eta = b0 + offset
            p = np.clip(1.0 / (1.0 + np.exp(-eta)), 1e-10, 1 - 1e-10)
            w = p * (1 - p)
            grad = np.sum(y - p)
            hess = -np.sum(w)
            if abs(hess) < 1e-12:
                break
            b0 = b0 - grad / hess
        eta = b0 + offset
        p = np.clip(1.0 / (1.0 + np.exp(-eta)), 1e-10, 1 - 1e-10)
        ll = np.sum(y * np.log(p) + (1 - y) * np.log(1 - p))
        return -ll, b0

    # 搜索 CI 边界
    beta_hat = beta_mle[param_idx]
    se_hat = 1.0 / np.sqrt(max(np.sum(X[:, param_idx] ** 2 * 0.25), 1e-6))
    bounds = []
    for direction in [-1, 1]:
        lo, hi = beta_hat, beta_hat + direction * 5 * se_hat
        last_val = beta_hat
        for _ in range(50):
            mid = (lo + hi) / 2
            neg_ll_mid, _ = neg_ll_fixed(mid)
            lr_stat = -2 * (neg_ll_mid - (-ll_full))  # 注意符号
            # lr_stat 应该 >= 0, 当 mid 远离 MLE 时增大
            if lr_stat < chi2_crit:
                lo = mid if direction > 0 else mid
                hi = mid if direction > 0 else mid
                # 扩大搜索
                if direction > 0:
                    hi = mid + se_hat
                    lo = beta_hat
                else:
                    hi = beta_hat
                    lo = mid - se_hat
                last_val = mid
            else:
                if direction > 0:
                    hi = mid
                else:
                    lo = mid
                last_val = mid
        bounds.append(last_val)
    return min(bounds), max(bounds)


def fit_logistic_firth(y, x, use_firth=True):
    """统一接口: 拟合 logistic 回归, 返回 OR, CI, p"""
    X = np.column_stack([np.ones(len(y)), x])
    y = np.asarray(y, dtype=float)
    try:
        if use_firth:
            beta, se, ll, conv = firth_logistic_fit(X, y)
        else:
            # 标准 MLE (statsmodels)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                model = sm.Logit(y, X).fit(disp=False, maxiter=100)
            beta = model.params
            se = model.bse
            ll = model.llf
            conv = model.converged
        or_val = float(np.exp(beta[1]))
        # Wald CI (与主分析一致)
        z_crit = 1.959963984540054  # 0.975 quantile of N(0,1)
        ci_lower = float(np.exp(beta[1] - z_crit * se[1]))
        ci_upper = float(np.exp(beta[1] + z_crit * se[1]))
        # Wald p-value
        z = beta[1] / se[1] if se[1] > 0 else 0
        from scipy.stats import norm
        p_val = 2 * (1 - norm.cdf(abs(z)))
        return {
            'OR': or_val, 'CI_lower': ci_lower, 'CI_upper': ci_upper,
            'p_value': p_val, 'beta': float(beta[1]), 'se': float(se[1]),
            'log_lik': float(ll), 'converged': bool(conv),
            'method': 'Firth' if use_firth else 'MLE'
        }
    except Exception as e:
        return {
            'OR': np.nan, 'CI_lower': np.nan, 'CI_upper': np.nan,
            'p_value': np.nan, 'beta': np.nan, 'se': np.nan,
            'log_lik': np.nan, 'converged': False,
            'method': 'Firth' if use_firth else 'MLE', 'error': str(e)
        }


print("[step4c_chemo_stratified] Starting chemotherapy stratified analysis...")
os.makedirs(cfg.RESULT_DIR, exist_ok=True)

# ============================================================
# 1. 数据加载与化疗方案分组
# ============================================================
import scanpy as sc

print("Loading GSE243013 metadata (backed mode)...")
adata = sc.read_h5ad(os.path.join(cfg.ROOT_DIR, 'GSE243013_T_cells.h5ad'), backed='r')
pat_tx = adata.obs.drop_duplicates('sampleID')[['sampleID', 'chemotherapy', 'anti-PD1_therapy']].copy()
pat_tx.columns = ['patient_id', 'chemotherapy', 'anti_pd1']
del adata
print(f"  Extracted treatment info for {len(pat_tx)} patients")

ppm = pd.read_csv(cfg.result_path('per_patient_metrics.csv'))
print(f"  Loaded per_patient_metrics: {len(ppm)} patients")

merged = ppm.merge(pat_tx, on='patient_id', how='left')
print(f"  Merged: {len(merged)}, chemo missing: {merged['chemotherapy'].isna().sum()}")

# 化疗大类分组（统一引用 _common.classify_chemo）
from _common import classify_chemo

merged['chemo_class'] = merged['chemotherapy'].apply(classify_chemo)

# 只保留 pCR vs non-MPR 主分析队列
main = merged[merged['response'].isin(['pCR', 'non-MPR'])].copy()
main['response_binary'] = (main['response'] == 'pCR').astype(int)

print(f"\n=== 化疗大类分组 (pCR vs non-MPR, n={len(main)}) ===")
ct = pd.crosstab(main['chemo_class'], main['response'], margins=True)
print(ct)
print()
print("=== 各组 NK-dominant ratio 描述 ===")
print(main.groupby('chemo_class')['nk_dominant_ratio'].describe()[['count', 'mean', 'std', 'min', 'max']].round(4))

# 保存
main.to_csv(cfg.result_path('chemo_stratified_patient_metrics.csv'), index=False)
print(f"\n  Saved: chemo_stratified_patient_metrics.csv")

# ============================================================
# 2. 各亚组 Firth 惩罚 logistic 回归
# ============================================================
print("\n" + "="*60)
print("STRATIFIED LOGISTIC REGRESSION (pCR ~ nk_dominant_ratio)")
print("="*60)

# 分析亚组: 样本量 >= 10 的组
MIN_GROUP_SIZE = 10
group_counts = main['chemo_class'].value_counts()
analyzable_groups = group_counts[group_counts >= MIN_GROUP_SIZE].index.tolist()
analyzable_groups = [g for g in analyzable_groups if g != 'Unknown']
excluded_groups = [g for g in group_counts.index if g not in analyzable_groups and g != 'Unknown']
print(f"  Analyzable groups (n>={MIN_GROUP_SIZE}): {analyzable_groups}")
print(f"  Excluded groups (n<{MIN_GROUP_SIZE}): {dict(group_counts[excluded_groups])}")

stratified_results = []
# 全队列 (参考)
y_all = main['response_binary'].values
x_all = main['nk_dominant_ratio'].values
res_all = fit_logistic_firth(y_all, x_all, use_firth=True)
res_all['group'] = 'All (n=188)'
res_all['n'] = len(main)
res_all['n_pCR'] = int(y_all.sum())
res_all['n_NR'] = int((1 - y_all).sum())
res_all['excluded'] = False
res_all['exclusion_reason'] = ''
stratified_results.append(res_all)
print(f"\n  All (n={len(main)}): OR={res_all['OR']:.2f} [{res_all['CI_lower']:.2f}, {res_all['CI_upper']:.2f}], p={res_all['p_value']:.4f}")

for grp in analyzable_groups:
    sub = main[main['chemo_class'] == grp]
    y = sub['response_binary'].values
    x = sub['nk_dominant_ratio'].values
    n_pcr = int(y.sum())
    n_nr = int((1 - y).sum())
    res = fit_logistic_firth(y, x, use_firth=True)
    res['group'] = grp
    res['n'] = len(sub)
    res['n_pCR'] = n_pcr
    res['n_NR'] = n_nr
    res['excluded'] = False
    res['exclusion_reason'] = ''
    stratified_results.append(res)
    warn = ' [WARNING: n<20, results may be unstable]' if len(sub) < 20 else ''
    print(f"  {grp} (n={len(sub)}, pCR={n_pcr}, NR={n_nr}): "
          f"OR={res['OR']:.2f} [{res['CI_lower']:.2f}, {res['CI_upper']:.2f}], p={res['p_value']:.4f}{warn}")

for grp in excluded_groups:
    sub = main[main['chemo_class'] == grp]
    y = sub['response_binary'].values
    n_pcr = int(y.sum())
    n_nr = int((1 - y).sum())
    stratified_results.append({
        'group': grp,
        'n': len(sub),
        'n_pCR': n_pcr,
        'n_NR': n_nr,
        'OR': np.nan,
        'CI_lower': np.nan,
        'CI_upper': np.nan,
        'p_value': np.nan,
        'method': '',
        'converged': False,
        'excluded': True,
        'exclusion_reason': f'n<{MIN_GROUP_SIZE}'
    })
    print(f"  {grp} (n={len(sub)}, EXCLUDED): too small for analysis")

strat_df = pd.DataFrame(stratified_results)
out_cols = ['group', 'n', 'n_pCR', 'n_NR', 'OR', 'CI_lower', 'CI_upper', 'p_value', 'method', 'converged', 'excluded', 'exclusion_reason']
strat_df_out = strat_df[out_cols].copy()
strat_df_out.to_csv(cfg.result_path('chemo_stratified_results.csv'), index=False)
print(f"\n  Saved: chemo_stratified_results.csv")

# ============================================================
# 3. 跨亚组异质性检验 (Cochran's Q)
# ============================================================
print("\n" + "="*60)
print("HETEROGENEITY TEST (Cochran's Q)")
print("="*60)

# 使用各亚组 beta 和 se 计算 Cochran's Q
sub_betas = []
sub_ses = []
sub_names = []
for r in stratified_results:
    if r['group'] != 'All (n=188)' and not r.get('excluded', False) and 'beta' in r and not np.isnan(r['beta']) and r.get('se', 0) > 0:
        sub_betas.append(r['beta'])
        sub_ses.append(r['se'])
        sub_names.append(r['group'])

sub_betas = np.array(sub_betas)
sub_ses = np.array(sub_ses)
sub_vars = sub_ses ** 2

# 固定效应权重
w = 1.0 / sub_vars
beta_pooled = np.sum(w * sub_betas) / np.sum(w)
Q = np.sum(w * (sub_betas - beta_pooled) ** 2)
df_q = len(sub_betas) - 1
p_het = float(chi2.sf(Q, df_q))
I2 = max(0, (Q - df_q) / Q) * 100 if Q > 0 else 0

print(f"  Groups: {sub_names}")
print(f"  Pooled beta (fixed effect): {beta_pooled:.4f}")
print(f"  Cochran's Q = {Q:.4f}, df = {df_q}, p = {p_het:.4f}")
print(f"  I² = {I2:.1f}%")
print(f"  Interpretation: {'Significant heterogeneity' if p_het < 0.10 else 'No significant heterogeneity'}")

# ============================================================
# 4. 化疗方案对 NK-dominant ratio 水平的影响 (Kruskal-Wallis)
# ============================================================
print("\n" + "="*60)
print("NK-DOMINANT RATIO ACROSS CHEMO GROUPS (Kruskal-Wallis)")
print("="*60)

groups_for_kw = [main[main['chemo_class'] == g]['nk_dominant_ratio'].values
                 for g in analyzable_groups if len(main[main['chemo_class'] == g]) >= MIN_GROUP_SIZE]
kw_stat, kw_p = kruskal(*groups_for_kw)
print(f"  Kruskal-Wallis: H={kw_stat:.4f}, p={kw_p:.4f}")
print(f"  Interpretation: {'Significant difference' if kw_p < 0.05 else 'No significant difference'} in NK-dominant ratio across chemo groups")

# 各组 pCR 率
print("\n=== 各组 pCR 率 ===")
for grp in analyzable_groups:
    sub = main[main['chemo_class'] == grp]
    pcr_rate = sub['response_binary'].mean() * 100
    print(f"  {grp}: pCR rate = {pcr_rate:.1f}% ({int(sub['response_binary'].sum())}/{len(sub)})")

# ============================================================
# 5. 可视化 (4-panel figure)
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# --- Panel A: 森林图 ---
axA = axes[0, 0]
plot_data = stratified_results.copy()
y_pos = np.arange(len(plot_data))[::-1]  # 从上到下

colors = ['#37474F' if 'All' in r['group'] else '#1976D2' for r in plot_data]
for i, r in enumerate(plot_data):
    or_val = r['OR']
    ci_lo = r['CI_lower']
    ci_hi = r['CI_upper']
    color = colors[i]
    axA.errorbar(or_val, y_pos[i],
                 xerr=[[or_val - ci_lo], [ci_hi - or_val]],
                 fmt='s', color=color, capsize=5, markersize=10, elinewidth=2)
    label = f"{r['group']} (n={r['n']}, pCR={r['n_pCR']})"
    label += f"\nOR={or_val:.2f} [{ci_lo:.2f}, {ci_hi:.2f}], p={r['p_value']:.3f}"
    axA.text(max(ci_hi, or_val) * 1.05, y_pos[i], label, va='center', fontsize=8)

axA.axvline(x=1, color='red', linestyle='--', alpha=0.5, label='OR=1 (no effect)')
beta_pooled_exp = float(np.exp(beta_pooled))
axA.axvline(x=beta_pooled_exp, color='green', linestyle=':', alpha=0.7,
            label=f'Pooled OR={beta_pooled_exp:.2f}')
axA.set_yticks(y_pos)
axA.set_yticklabels([r['group'] for r in plot_data], fontsize=9)
axA.set_xlabel('Odds Ratio (pCR per unit NK-dominant ratio)', fontsize=10)
axA.set_xscale('log')
axA.set_title(f'Panel A: Forest Plot - NK-dominant vs pCR by Chemo Class\n'
              f'Cochran Q={Q:.2f}, p={p_het:.3f}, I²={I2:.1f}%', fontsize=11)
axA.legend(loc='lower right', fontsize=8)
axA.grid(alpha=0.3, axis='x', ls=':')
axA.spines['top'].set_visible(False); axA.spines['right'].set_visible(False)

# --- Panel B: NK-dominant ratio 分布 (按化疗方案) ---
axB = axes[0, 1]
plot_groups = [g for g in analyzable_groups if len(main[main['chemo_class'] == g]) >= MIN_GROUP_SIZE]
box_data = [main[main['chemo_class'] == g]['nk_dominant_ratio'].values for g in plot_groups]
bp = axB.boxplot(box_data, labels=[f"{g}\n(n={len(main[main['chemo_class']==g])})" for g in plot_groups],
                 patch_artist=True, widths=0.6)
box_colors = ['#FF9800', '#4CAF50', '#2196F3', '#9C27B0', '#795548']
for patch, color in zip(bp['boxes'], box_colors[:len(plot_groups)]):
    patch.set_facecolor(color); patch.set_alpha(0.6)
for i, data in enumerate(box_data):
    x_j = np.random.normal(i + 1, 0.05, size=len(data))
    axB.scatter(x_j, data, c='black', s=15, alpha=0.5, zorder=5)

axB.set_ylabel('NK-dominant ratio', fontsize=10)
axB.set_title(f'Panel B: NK-dominant Ratio Distribution by Chemo Class\n'
              f'Kruskal-Wallis: H={kw_stat:.2f}, p={kw_p:.4f}', fontsize=11)
axB.grid(alpha=0.3, axis='y', ls=':')
axB.spines['top'].set_visible(False); axB.spines['right'].set_visible(False)

# --- Panel C: 各组 pCR 率对比 ---
axC = axes[1, 0]
pcr_rates = []
ci_lo_rates = []
ci_hi_rates = []
group_labels = []
from scipy.stats import norm
z = 1.96
for grp in plot_groups:
    sub = main[main['chemo_class'] == grp]
    n = len(sub)
    k = int(sub['response_binary'].sum())
    rate = k / n
    # 正确的 Wilson score CI
    denom = 1 + z**2 / n
    center = (rate + z**2 / (2 * n)) / denom
    margin = z * np.sqrt(rate * (1 - rate) / n + z**2 / (4 * n**2)) / denom
    lo = max(0, (center - margin) * 100)
    hi = min(100, (center + margin) * 100)
    pcr_rates.append(rate * 100)
    ci_lo_rates.append(lo)
    ci_hi_rates.append(hi)
    group_labels.append(f"{grp}\n(n={n})")

x_pos = np.arange(len(plot_groups))
axC.bar(x_pos, pcr_rates, color=box_colors[:len(plot_groups)], alpha=0.7, edgecolor='black')
# yerr 使用绝对值, 确保非负
yerr_lo = np.array(pcr_rates) - np.array(ci_lo_rates)
yerr_hi = np.array(ci_hi_rates) - np.array(pcr_rates)
yerr_lo = np.maximum(yerr_lo, 0)
yerr_hi = np.maximum(yerr_hi, 0)
axC.errorbar(x_pos, pcr_rates, yerr=[yerr_lo, yerr_hi],
             fmt='none', color='black', capsize=5, lw=1.5)
for i, (rate, lo, hi) in enumerate(zip(pcr_rates, ci_lo_rates, ci_hi_rates)):
    axC.text(i, hi + 2, f'{rate:.1f}%\n[{lo:.0f}, {hi:.0f}]', ha='center', fontsize=8)

# 全队列参考线
overall_pcr = main['response_binary'].mean() * 100
axC.axhline(y=overall_pcr, color='red', linestyle='--', alpha=0.5,
            label=f'Overall pCR rate = {overall_pcr:.1f}%')
axC.set_xticks(x_pos)
axC.set_xticklabels(group_labels, fontsize=9)
axC.set_ylabel('pCR rate (%)', fontsize=10)
axC.set_ylim(0, 100)
axC.set_title('Panel C: pCR Rate by Chemo Class (Wilson CI)', fontsize=11)
axC.legend(loc='upper right', fontsize=8)
axC.grid(alpha=0.3, axis='y', ls=':')
axC.spines['top'].set_visible(False); axC.spines['right'].set_visible(False)

# --- Panel D: 主队列 Carbo+Pemetrexed vs GSE179994 Treatment A 跨队列对比 ---
axD = axes[1, 1]
# 主队列 Carbo+Pemetrexed 亚组
main_peme = main[main['chemo_class'] == 'Platinum+Pemetrexed']
main_peme_pcr = main_peme[main_peme['response'] == 'pCR']['nk_dominant_ratio'].values
main_peme_nr = main_peme[main_peme['response'] == 'non-MPR']['nk_dominant_ratio'].values

# GSE179994: 从 step4 已计算的签名得分 (post-treatment R vs NR)
g17_path = cfg.result_path('GSE179994_response_sample_scores.csv')
g17_data = None
if os.path.exists(g17_path):
    g17_df = pd.read_csv(g17_path)
    g17_data = g17_df

# 准备对比数据: 主队列 Peme 亚组 NK-dominant ratio (R vs NR)
positions = [1, 2]
data_compare = [main_peme_nr, main_peme_pcr]
labels = [f'Main cohort\nPlatinum+Pemetrexed\nNR (n={len(main_peme_nr)})',
          f'Main cohort\nPlatinum+Pemetrexed\npCR (n={len(main_peme_pcr)})']

bp2 = axD.boxplot(data_compare, positions=positions, labels=labels,
                  patch_artist=True, widths=0.5)
for patch, color in zip(bp2['boxes'], ['#F44336', '#4CAF50']):
    patch.set_facecolor(color); patch.set_alpha(0.6)
for i, data in enumerate(data_compare):
    x_j = np.random.normal(positions[i], 0.04, size=len(data))
    axD.scatter(x_j, data, c='black', s=20, alpha=0.6, zorder=5)

# 主队列 Peme 亚组 MWU 检验
if len(main_peme_pcr) >= 3 and len(main_peme_nr) >= 3:
    mwu_stat, mwu_p = mannwhitneyu(main_peme_pcr, main_peme_nr, alternative='two-sided')
else:
    mwu_p = np.nan

# 主队列 Peme 亚组 OR
res_peme = [r for r in stratified_results if r['group'] == 'Platinum+Pemetrexed'][0]
axD.set_title(f'Panel D: Main Cohort Platinum+Pemetrexed Subgroup\n'
              f'(Same regimen as GSE179994 Treatment A)\n'
              f'MWU p={mwu_p:.4f} | Firth OR={res_peme["OR"]:.2f} '
              f'[{res_peme["CI_lower"]:.2f}, {res_peme["CI_upper"]:.2f}], p={res_peme["p_value"]:.3f}',
              fontsize=10)
axD.set_ylabel('NK-dominant ratio', fontsize=10)
axD.grid(alpha=0.3, axis='y', ls=':')
axD.spines['top'].set_visible(False); axD.spines['right'].set_visible(False)

# 添加 GSE179994 参考文本
axD.text(0.02, 0.98,
         f'GSE179994 Treatment A (same regimen):\n'
         f'AUC=0.4722, p=0.9399 (n=14, post-treatment)\n'
         f'(Different metric: signature score, not TCR clonotype)',
         transform=axD.transAxes, fontsize=8, va='top',
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

fig.suptitle('Chemotherapy Regimen Stratified Analysis (GSE243013 Main Cohort, n=188)\n'
             'Impact of chemo class on NK-dominant clonotype ratio predictive power for pCR',
             fontsize=13, y=1.00)

plt.tight_layout()
fig.savefig(cfg.result_path('Fig_chemo_stratified.pdf'), dpi=150, bbox_inches='tight')
fig.savefig(cfg.result_path('Fig_chemo_stratified.png'), dpi=150, bbox_inches='tight')
plt.close(fig)
print(f"\n  Saved: Fig_chemo_stratified.png/pdf")

# ============================================================
# 6. 总结
# ============================================================
print("\n" + "="*70)
print("CHEMO STRATIFIED ANALYSIS SUMMARY")
print("="*70)
print(f"""
Dataset: GSE243013 main cohort (n=188, pCR vs non-MPR)
- 175/188 patients received anti-PD1 + chemotherapy (NOT monotherapy as previously thought)
- 13 patients received anti-PD1 monotherapy (No chemo)

Chemo class distribution:
{ct.to_string()}
""")

print("Stratified OR results:")
for r in stratified_results:
    print(f"  {r['group']:30s} (n={r['n']:3d}, pCR={r['n_pCR']:3d}): "
          f"OR={r['OR']:7.2f} [{r['CI_lower']:6.2f}, {r['CI_upper']:7.2f}], p={r['p_value']:.4f}")

print(f"""
Heterogeneity: Cochran Q={Q:.2f}, p={p_het:.4f}, I²={I2:.1f}%
NK-dominant ratio across groups: Kruskal-Wallis H={kw_stat:.2f}, p={kw_p:.4f}

Key findings:
1. Overall OR (n=188) = {res_all['OR']:.2f} reflects MIXED chemo regimens (175/188 had chemo)
2. Stratified analysis reveals whether NK-dominant predictive power is uniform across regimens
3. Platinum+Pemetrexed subgroup (n=35) is the SAME regimen as GSE179994 Treatment A
   → Direct cross-cohort validation of chemo effect on NK-like predictive power
4. {'Significant heterogeneity detected' if p_het < 0.10 else 'No significant heterogeneity'} across chemo classes
""")

print("[step4c_chemo_stratified] Done.")
