"""
step2_firth_validation.py
Step 2 补充：Firth's Penalized Logistic Regression 稳健性验证
目的：排除因样本量或数据结构导致的"完全分离"偏差，确认 OR 估计的稳健性
自实现 Firth 惩罚逻辑回归（最小化惩罚负对数似然 + 轮廓似然 CI）
"""
import config as cfg
from _common import firth_logistic_fit
import pandas as pd
import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2
import statsmodels.api as sm
import os
import warnings
warnings.filterwarnings('ignore')

print("[step2_firth_validation] Starting...")
os.makedirs(cfg.RESULT_DIR, exist_ok=True)


def firth_profile_likelihood_ci(X, y, beta, se, param_idx=1, alpha=0.05, n_points=25):
    """
    轮廓似然置信区间（profile likelihood CI）
    对参数 param_idx 固定其值，优化其余参数，计算惩罚似然比
    使用 chisq(1, 1-alpha) 作为阈值
    """
    n, p = X.shape
    chi2_crit = chi2.ppf(1 - alpha, 1)
    
    # Full model penalized log-likelihood (at MLE)
    eta_full = X @ beta
    p_full = 1.0 / (1.0 + np.exp(-eta_full))
    p_full = np.clip(p_full, 1e-10, 1 - 1e-10)
    w_full = p_full * (1 - p_full)
    I_full = X.T @ (X * w_full[:, None])
    sign_full, log_det_full = np.linalg.slogdet(I_full)
    ll_full = np.sum(y * np.log(p_full) + (1 - y) * np.log(1 - p_full))
    ll_pen_full = ll_full + 0.5 * log_det_full if sign_full > 0 else ll_full
    
    def constrained_ll_pen(beta_val):
        """固定 beta[param_idx] = beta_val，优化其余参数的惩罚似然"""
        beta_c = beta.copy()
        beta_c[param_idx] = beta_val
        
        # Other indices
        other_idx = [i for i in range(p) if i != param_idx]
        X_other = X[:, other_idx]
        
        # 用迭代优化其余参数（简化：Newton 几步）
        for _ in range(30):
            eta = X @ beta_c
            p_i = 1.0 / (1.0 + np.exp(-eta))
            p_i = np.clip(p_i, 1e-10, 1 - 1e-10)
            w = p_i * (1 - p_i)
            
            WX = X * w[:, None]
            I = X.T @ WX
            try:
                I_inv = np.linalg.inv(I)
            except:
                I_inv = np.linalg.pinv(I)
            
            h = np.sum((X @ I_inv) * WX, axis=1)
            h = np.clip(h, 1e-10, 1 - 1e-10)
            
            # 只更新非固定参数
            score = X.T @ (y - p_i + h * (0.5 - p_i))
            score_other = score[other_idx]
            I_other = I[np.ix_(other_idx, other_idx)]
            
            try:
                step_other = np.linalg.solve(I_other, score_other)
            except:
                step_other = np.linalg.lstsq(I_other, score_other, rcond=None)[0]
            
            # 限幅
            max_step = 0.5
            step_norm = np.max(np.abs(step_other))
            if step_norm > max_step:
                step_other = step_other * max_step / step_norm
            
            beta_c[other_idx] += step_other
            
            if np.max(np.abs(step_other)) < 1e-7:
                break
        
        # Compute penalized LL at constrained solution
        eta_c = X @ beta_c
        p_c = 1.0 / (1.0 + np.exp(-eta_c))
        p_c = np.clip(p_c, 1e-10, 1 - 1e-10)
        w_c = p_c * (1 - p_c)
        I_c = X.T @ (X * w_c[:, None])
        sign_c, log_det_c = np.linalg.slogdet(I_c)
        ll_c = np.sum(y * np.log(p_c) + (1 - y) * np.log(1 - p_c))
        ll_pen_c = ll_c + 0.5 * log_det_c if sign_c > 0 else ll_c
        
        return ll_pen_c
    
    # 搜索区间：从 MLE 向两侧各扩展 4 个 SE
    beta_hat = beta[param_idx]
    se_hat = se[param_idx]
    low = beta_hat - 4 * se_hat
    high = beta_hat + 4 * se_hat
    
    # 似然比统计量: 2 * (ll_full - ll_constrained)
    # CI 定义: 2 * (ll_full - ll_constrained) < chi2_crit
    # 即 ll_constrained > ll_full - chi2_crit / 2
    target = ll_pen_full - chi2_crit / 2.0
    
    # 搜索下界
    ci_lower = low
    for b in np.linspace(beta_hat, low, n_points):
        ll = constrained_ll_pen(b)
        if ll < target:
            ci_lower = b
            break
        ci_lower = b
    
    # 搜索上界
    ci_upper = high
    for b in np.linspace(beta_hat, high, n_points):
        ll = constrained_ll_pen(b)
        if ll < target:
            ci_upper = b
            break
        ci_upper = b
    
    return ci_lower, ci_upper, ll_pen_full


# ── 1. 数据准备 ──
df = pd.read_csv(cfg.result_path('per_patient_metrics.csv'))
df['y'] = (df['response'] == 'pCR').astype(int)
df = df.dropna(subset=['y', 'nk_dominant_ratio'])

X_simple = sm.add_constant(df['nk_dominant_ratio'].values)
y = df['y'].values

print(f"样本数: n={len(y)}, pCR={int(y.sum())}, non-MPR={int(len(y)-y.sum())}")
print(f"NK-dominant ratio range: [{df['nk_dominant_ratio'].min():.4f}, {df['nk_dominant_ratio'].max():.4f}]")

# ── 2. 普通 Logistic 回归（statsmodels MLE） ──
print("\n=== 普通 Logistic 回归 (statsmodels MLE) ===")
model_sm = sm.Logit(y, X_simple)
result_sm = model_sm.fit(disp=0, maxiter=5000)
print(result_sm.summary())

or_sm = np.exp(result_sm.params[1])
ci_sm = np.exp(result_sm.conf_int()[1])
p_sm = result_sm.pvalues[1]
print(f"\nOR (nk_dominant_ratio): {or_sm:.4f}")
print(f"95% CI (Wald): [{ci_sm[0]:.4f}, {ci_sm[1]:.4f}]")
print(f"p-value (Wald): {p_sm:.6e}")

# 检查是否存在完全分离
pred = result_sm.predict(X_simple)
sep_pos = np.all(pred[y == 1] > 0.99) or np.all(pred[y == 0] < 0.01)
print(f"疑似完全分离: {'是' if sep_pos else '否'}")

# ── 3. Firth 惩罚逻辑回归 ──
print("\n=== Firth's Penalized Logistic Regression ===")
beta_firth, se_firth, ll_pen_firth, converged_f = firth_logistic_fit(X_simple, y)
print(f"Converged: {converged_f}")
print(f"beta_firth (intercept, nk_dominant): {beta_firth}")
print(f"se_firth: {se_firth}")

or_firth = np.exp(beta_firth[1])
ci_firth_wald_lower = np.exp(beta_firth[1] - 1.96 * se_firth[1])
ci_firth_wald_upper = np.exp(beta_firth[1] + 1.96 * se_firth[1])

# Wald p (基于惩罚 SE)
from scipy.stats import norm
z_firth = beta_firth[1] / se_firth[1] if se_firth[1] > 0 else 0
p_firth_wald = 2 * (1 - norm.cdf(abs(z_firth)))

print(f"\nOR (Firth): {or_firth:.4f}")
print(f"95% CI (Wald-based): [{ci_firth_wald_lower:.4f}, {ci_firth_wald_upper:.4f}]")
print(f"p-value (Wald-based): {p_firth_wald:.6e}")

# 轮廓似然 CI（更准确，优先使用）
print("\n计算轮廓似然置信区间 (Profile Likelihood CI)...")
ci_pl_lower, ci_pl_upper, ll_pen_full = firth_profile_likelihood_ci(
    X_simple, y, beta_firth, se_firth, param_idx=1, alpha=0.05
)
or_pl_lower = np.exp(ci_pl_lower)
or_pl_upper = np.exp(ci_pl_upper)
print(f"95% CI (Profile Likelihood): [{or_pl_lower:.4f}, {or_pl_upper:.4f}]")

# ── 4. 多变量模型（含临床协变量）也做 Firth ──
print("\n=== 多变量模型 Firth 验证（Clinical + NK-dominant）===")

df_mv = df.copy()
clinical_feats = []
if 'age' in df_mv.columns:
    df_mv['age_num'] = pd.to_numeric(df_mv['age'], errors='coerce')
    clinical_feats.append('age_num')
if 'stage' in df_mv.columns:
    stage_series = df_mv['stage'].astype(str).str.strip()
    stage_series = stage_series.replace({'': 'Other', 'unknown': 'Other', 'nan': 'Other'})
    def cs(s):
        s_up = str(s).upper()
        if s_up.startswith('IV'): return 'IV'
        elif s_up.startswith('III'): return 'III'
        elif s_up.startswith('II'): return 'II'
        elif s_up.startswith('I') and not s_up.startswith('II'): return 'I'
        else: return 'Other'
    df_mv['stage_group'] = stage_series.apply(cs)
    stage_d = pd.get_dummies(df_mv['stage_group'], prefix='stage', drop_first=True).astype(int)
    for c in stage_d.columns:
        df_mv[c] = stage_d[c].values
        clinical_feats.append(c)
if 'histology' in df_mv.columns:
    df_mv['hist_luad'] = (df_mv['histology'].astype(str).str.upper() == 'LUAD').astype(int)
    clinical_feats.append('hist_luad')
if 'sex' in df_mv.columns or 'gender' in df_mv.columns:
    sex_c = 'sex' if 'sex' in df_mv.columns else 'gender'
    df_mv['sex_male'] = (df_mv[sex_c].astype(str).str.upper().isin(['M', 'MALE', '男'])).astype(int)
    clinical_feats.append('sex_male')

df_mv = df_mv.dropna(subset=['y'] + clinical_feats + ['nk_dominant_ratio'])
y_mv = df_mv['y'].values

# 特征矩阵: const + clinical + nk
X_cols = clinical_feats + ['nk_dominant_ratio']
X_mv = df_mv[X_cols].values
X_mv = sm.add_constant(X_mv)

print(f"多变量模型 n={len(y_mv)}, features={X_cols}")

# 普通 MLE
model_mv_sm = sm.Logit(y_mv, X_mv)
result_mv_sm = model_mv_sm.fit(disp=0, maxiter=5000)
or_mv_sm = np.exp(result_mv_sm.params[-1])
ci_mv_sm = np.exp(result_mv_sm.conf_int()[-1])
p_mv_sm = result_mv_sm.pvalues[-1]
print(f"普通 Logistic OR (NK-dominant): {or_mv_sm:.4f}, 95%CI=[{ci_mv_sm[0]:.4f}, {ci_mv_sm[1]:.4f}], p={p_mv_sm:.6e}")

# Firth
beta_mv_f, se_mv_f, ll_mv_f, conv_mv = firth_logistic_fit(X_mv, y_mv)
print(f"Firth converged: {conv_mv}")
or_mv_f = np.exp(beta_mv_f[-1])
ci_mv_f_wald_l = np.exp(beta_mv_f[-1] - 1.96 * se_mv_f[-1])
ci_mv_f_wald_u = np.exp(beta_mv_f[-1] + 1.96 * se_mv_f[-1])
z_mv_f = beta_mv_f[-1] / se_mv_f[-1] if se_mv_f[-1] > 0 else 0
p_mv_f = 2 * (1 - norm.cdf(abs(z_mv_f)))
print(f"Firth OR (NK-dominant): {or_mv_f:.4f}, 95%CI(Wald)=[{ci_mv_f_wald_l:.4f}, {ci_mv_f_wald_u:.4f}], p={p_mv_f:.6e}")

# Profile likelihood CI for multivariate
ci_mv_pl_l, ci_mv_pl_u, _ = firth_profile_likelihood_ci(
    X_mv, y_mv, beta_mv_f, se_mv_f, param_idx=X_mv.shape[1]-1, alpha=0.05
)
print(f"Firth 95%CI (Profile Likelihood): [{np.exp(ci_mv_pl_l):.4f}, {np.exp(ci_mv_pl_u):.4f}]")

# ── 5. 保存对比表 ──
comparison_data = []

# 单变量对比
comparison_data.append({
    'model': 'Univariate (NK-dominant only)',
    'method': 'Standard Logistic (MLE)',
    'OR': round(or_sm, 4),
    'CI_lower': round(ci_sm[0], 4),
    'CI_upper': round(ci_sm[1], 4),
    'CI_method': 'Wald',
    'p_value': f'{p_sm:.6e}',
    'n_patients': len(y),
    'converged': True,
})

comparison_data.append({
    'model': 'Univariate (NK-dominant only)',
    'method': 'Firth Penalized',
    'OR': round(or_firth, 4),
    'CI_lower': round(or_pl_lower, 4),
    'CI_upper': round(or_pl_upper, 4),
    'CI_method': 'Profile Likelihood',
    'p_value': f'{p_firth_wald:.6e}',
    'n_patients': len(y),
    'converged': converged_f,
})

comparison_data.append({
    'model': 'Univariate (NK-dominant only)',
    'method': 'Firth Penalized (Wald)',
    'OR': round(or_firth, 4),
    'CI_lower': round(ci_firth_wald_lower, 4),
    'CI_upper': round(ci_firth_wald_upper, 4),
    'CI_method': 'Wald (penalized SE)',
    'p_value': f'{p_firth_wald:.6e}',
    'n_patients': len(y),
    'converged': converged_f,
})

# 多变量对比
comparison_data.append({
    'model': 'Multivariate (Clinical + NK-dominant)',
    'method': 'Standard Logistic (MLE)',
    'OR': round(or_mv_sm, 4),
    'CI_lower': round(ci_mv_sm[0], 4),
    'CI_upper': round(ci_mv_sm[1], 4),
    'CI_method': 'Wald',
    'p_value': f'{p_mv_sm:.6e}',
    'n_patients': len(y_mv),
    'converged': True,
})

comparison_data.append({
    'model': 'Multivariate (Clinical + NK-dominant)',
    'method': 'Firth Penalized',
    'OR': round(or_mv_f, 4),
    'CI_lower': round(np.exp(ci_mv_pl_l), 4),
    'CI_upper': round(np.exp(ci_mv_pl_u), 4),
    'CI_method': 'Profile Likelihood',
    'p_value': f'{p_mv_f:.6e}',
    'n_patients': len(y_mv),
    'converged': conv_mv,
})

df_comp = pd.DataFrame(comparison_data)
df_comp.to_csv(cfg.result_path('firth_validation_comparison.csv'), index=False)
print("\nfirth_validation_comparison.csv saved")

print("\n=== 结果对比总结 ===")
print(df_comp.to_string())

print("\n[step2_firth_validation] Done.")
