"""
step5_clinical_translation.py
复现 Fig6 (clinical model): Nomogram + Bootstrap calibration + DCA + AUC comparison.
"""
import config as cfg
import _common
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from scipy.stats import norm
import os
import warnings
warnings.filterwarnings('ignore')

print("[step5_clinical_translation] Starting...")

os.makedirs(cfg.RESULT_DIR, exist_ok=True)

# ── Load data ──
df = pd.read_csv(cfg.result_path('per_patient_metrics.csv'))
print(f"Loaded patient metrics: {df.shape}")

# Prepare features
df_model = df.copy()
# Binary response: 1 = pCR, 0 = non-MPR
df_model['y'] = (df_model['response'] == 'pCR').astype(int)
df_model['nk_dominant_ratio_scaled'] = df_model['nk_dominant_ratio'] * 100  # scale to 0-100

# Handle missing values
df_model = df_model.dropna(subset=['y'])

# Clinical features (统一引用 _common.build_clinical_features)
df_model, clinical_features = _common.build_clinical_features(df_model)

# Drop rows with NaN in features
feature_cols = clinical_features + ['nk_dominant_ratio_scaled']
df_model = df_model.dropna(subset=feature_cols)

# Merge SPP1 data for consistent patient subset with step5_5
# Ensure all models (clinical, NK, IRS) use the same patient subset (n=179)
df_mye = pd.read_csv(cfg.result_path('myeloid_per_patient.csv'))
df_model = df_model.merge(df_mye[['patient_id', 'SPP1_TAM_ratio']],
                          left_on='patient_id', right_on='patient_id', how='inner')
df_model['SPP1_TAM_ratio'] = pd.to_numeric(df_model['SPP1_TAM_ratio'], errors='coerce')
df_model = df_model.dropna(subset=['SPP1_TAM_ratio']).copy()

X_clinical = df_model[clinical_features].values if len(clinical_features) > 0 else np.zeros((len(df_model), 0))
X_full = df_model[feature_cols].values
y = df_model['y'].values

print(f"Model n: {len(y)}, pCR: {y.sum()}, non-MPR: {(1-y).sum()}")
print(f"Clinical features ({len(clinical_features)}): {clinical_features}")

# ── Fit models ──
has_clinical = len(clinical_features) > 0
if has_clinical:
    model_clinical = LogisticRegression(penalty=None, max_iter=5000, solver='lbfgs')
    model_clinical.fit(X_clinical, y)
    pred_clinical = model_clinical.predict_proba(X_clinical)[:, 1]
    auc_clinical = roc_auc_score(y, pred_clinical)
    print(f"AUC clinical: {auc_clinical:.4f}")
else:
    print("[WARNING] No clinical features available, skipping clinical-only model")
    model_clinical = None
    pred_clinical = np.zeros_like(y, dtype=float)
    auc_clinical = np.nan

model_full = LogisticRegression(penalty=None, max_iter=5000, solver='lbfgs')
model_full.fit(X_full, y)
pred_full = model_full.predict_proba(X_full)[:, 1]
auc_full = roc_auc_score(y, pred_full)

print(f"AUC full: {auc_full:.4f}")

# ── IRS 模型（审计 P3 修复：整合 IRS Nomogram 对比） ──
# IRS = nk_dominant_ratio × (1 - SPP1_TAM_ratio_normalized)
# 使用已合并的 df_model（n=179，与 step5_5 一致）
df_model_irs = df_model.copy()

spp1_min = df_model_irs['SPP1_TAM_ratio'].min()
spp1_max = df_model_irs['SPP1_TAM_ratio'].max()
df_model_irs['SPP1_norm'] = (df_model_irs['SPP1_TAM_ratio'] - spp1_min) / (spp1_max - spp1_min) if spp1_max > spp1_min else 0.0
df_model_irs['IRS'] = df_model_irs['nk_dominant_ratio'] * (1.0 - df_model_irs['SPP1_norm'])
df_model_irs['irs_scaled'] = df_model_irs['IRS'] * 100

feature_cols_irs = clinical_features + ['irs_scaled']
df_model_irs = df_model_irs.dropna(subset=feature_cols_irs)
X_irs = df_model_irs[feature_cols_irs].values
y_irs = df_model_irs['y'].values

model_irs = LogisticRegression(penalty=None, max_iter=5000, solver='lbfgs')
model_irs.fit(X_irs, y_irs)
pred_irs = model_irs.predict_proba(X_irs)[:, 1]
auc_irs = roc_auc_score(y_irs, pred_irs)
print(f"AUC Clinical+IRS: {auc_irs:.4f} (n={len(y_irs)}, IRS 整合)")
print(f"  ΔAUC (IRS - NK): {auc_irs - auc_full:.4f}")

# 5-fold CV out-of-fold predictions (for calibration and DCA)
from sklearn.model_selection import StratifiedKFold
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
oof_full = np.zeros(len(y))
oof_clinical = np.zeros(len(y))
for train_idx, val_idx in skf.split(X_full, y):
    m_full_fold = LogisticRegression(penalty=None, max_iter=5000, solver='lbfgs')
    m_full_fold.fit(X_full[train_idx], y[train_idx])
    oof_full[val_idx] = m_full_fold.predict_proba(X_full[val_idx])[:, 1]
    
    if has_clinical:
        m_clin_fold = LogisticRegression(penalty=None, max_iter=5000, solver='lbfgs')
        m_clin_fold.fit(X_clinical[train_idx], y[train_idx])
        oof_clinical[val_idx] = m_clin_fold.predict_proba(X_clinical[val_idx])[:, 1]
    else:
        oof_clinical[val_idx] = y.mean()

auc_oof_full = roc_auc_score(y, oof_full)
auc_oof_clinical = roc_auc_score(y, oof_clinical) if has_clinical else np.nan
print(f"OOF AUC clinical: {auc_oof_clinical:.4f}" if has_clinical else "OOF AUC clinical: N/A")
print(f"OOF AUC full: {auc_oof_full:.4f}")

# Coefficients
coef_clinical = model_clinical.coef_[0] if has_clinical else np.array([])
coef_full = model_full.coef_[0]
intercept_clinical = model_clinical.intercept_[0] if has_clinical else 0
intercept_full = model_full.intercept_[0]

# Nomogram: 统一使用 |β × range| 归一化（与 step5_5 一致，适用于跨尺度变量）
lp_contrib = np.abs(coef_full) * (X_full.max(axis=0) - X_full.min(axis=0))
max_lp = lp_contrib.max() if lp_contrib.max() > 0 else 1
var_max_points = (lp_contrib / max_lp * 100).astype(int)

# ============================================================
# Save nomogram data (both Clinical only and Clinical+NK models)
# ============================================================
nomogram_data = []

# Clinical only model coefficients
if has_clinical:
    lp_contrib_cli = np.abs(coef_clinical) * (X_clinical.max(axis=0) - X_clinical.min(axis=0))
    max_lp_cli = lp_contrib_cli.max() if lp_contrib_cli.max() > 0 else 1
    var_max_pts_cli = (lp_contrib_cli / max_lp_cli * 100).astype(int)
    for i, feat in enumerate(clinical_features):
        nomogram_data.append({
            'variable': feat,
            'coefficient': round(coef_clinical[i], 4) if i < len(coef_clinical) else '',
            'max_points': int(var_max_pts_cli[i]) if i < len(var_max_pts_cli) else '',
            'model': 'Clinical_only'
        })

# Clinical+NK model coefficients
for i, feat in enumerate(clinical_features):
    nomogram_data.append({
        'variable': feat,
        'coefficient': round(coef_full[i], 4) if i < len(coef_full) else '',
        'max_points': int(var_max_points[i]) if i < len(var_max_points) else '',
        'model': 'Clinical+NK'
    })
nomogram_data.append({
    'variable': 'nk_dominant_ratio_scaled',
    'coefficient': round(coef_full[-1], 4),
    'max_points': int(var_max_points[-1]),
    'model': 'Clinical+NK'
})

df_nomogram = pd.DataFrame(nomogram_data)
df_nomogram.to_csv(cfg.result_path('fig6_nomogram_data.csv'), index=False)
print("fig6_nomogram_data.csv written")

# ============================================================
# FIGURE 6: Clinical Model (4 Panels)
# ============================================================
fig6, axes6 = plt.subplots(2, 2, figsize=(14, 12))

# Panel A: Nomogram (hand-drawn with matplotlib)
ax6A = axes6[0, 0]
total_max_points = int(var_max_points.sum())
n_vars = len(feature_cols)

ax6A.set_xlim(0, total_max_points + 60)
ax6A.set_ylim(-0.5, n_vars + 1.5)

# Points axis
ax6A.axhline(y=n_vars + 0.5, xmin=0, xmax=total_max_points / (total_max_points + 60), color='black', linewidth=1)
ax6A.text(total_max_points // 2, n_vars + 0.8, 'Points', ha='center', fontsize=10)

points_step = 20
for i in range(0, total_max_points + 1, points_step):
    x_frac = i / (total_max_points + 60)
    ax6A.axvline(x=i, ymin=(n_vars + 0.3) / (n_vars + 2), ymax=(n_vars + 0.7) / (n_vars + 2), color='black', linewidth=0.8)
    ax6A.text(i, n_vars + 0.6, str(i), ha='center', fontsize=7)

# Variables
for v_idx, feat in enumerate(feature_cols):
    y_pos = n_vars - v_idx
    ax6A.axhline(y=y_pos, xmin=0, xmax=var_max_points[v_idx] / (total_max_points + 60), color='black', linewidth=1)
    ax6A.text(-5, y_pos + 0.3, feat, ha='right', fontsize=8, va='center')
    
    # Get coefficient-based points mapping
    coef = coef_full[v_idx] if v_idx < len(coef_full) else 0
    vmax_pts = var_max_points[v_idx]
    # Map variable range to points (min→0, max→vmax_pts)
    var_min = df_model[feat].min()
    var_max = df_model[feat].max()
    var_range = var_max - var_min if var_max > var_min else 1
    
    n_ticks = min(5, vmax_pts // 20 + 1)
    for ti in range(n_ticks):
        tick_pt = int(ti * vmax_pts / max(n_ticks - 1, 1))
        val = var_min + (tick_pt / max(vmax_pts, 1)) * var_range
        ax6A.text(tick_pt, y_pos - 0.3, f'{val:.1f}', ha='center', fontsize=6)

# Total points → linear predictor → probability
ax6A.axhline(y=0, xmin=0, xmax=total_max_points / (total_max_points + 60), color='black', linewidth=1)
ax6A.text(total_max_points // 2, 0.3, 'Total Points', ha='center', fontsize=9)
ax6A.text(total_max_points // 2, -0.5, 'Linear Predictor', ha='center', fontsize=9)

# Probability scale (using actual model intercept + coefficient range)
lp_at_zero = intercept_full  # when all variables = 0 points (minimum values)
# Total points = sum of points → LP = intercept + sum(coef * x)
# Each point corresponds to max_lp / 100 units of LP change（max_lp 来自 |β×range| 归一化）
lp_per_point = max_lp / 100.0
lp_at_max = intercept_full + total_max_points * lp_per_point
lp_min = min(lp_at_zero, lp_at_max)
lp_max = max(lp_at_zero, lp_at_max)
prob_ticks = [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95]
ax6A.text(total_max_points // 2, -1.3, 'pCR Probability', ha='center', fontsize=9)
for pt in prob_ticks:
    lp = np.log(pt / (1 - pt))
    if lp_min <= lp <= lp_max:
        x_pos = (lp - lp_min) / (lp_max - lp_min) * total_max_points
        ax6A.text(x_pos, -1.5, f'{pt:.2f}', ha='center', fontsize=7)

ax6A.set_title('Panel A: Nomogram (pCR ~ NK-dominant + Clinical)')
ax6A.axis('off')

# Panel B: Calibration curve (5-fold CV out-of-fold predictions)
ax6B = axes6[0, 1]
bins = np.linspace(0, 1, 11)

from sklearn.calibration import calibration_curve
prob_true, prob_pred = calibration_curve(y, oof_full, n_bins=10, strategy='uniform')
ax6B.plot(prob_pred, prob_true, 'o-', color='steelblue', linewidth=2, label='Model (5-fold CV OOF)')
ax6B.plot([0, 1], [0, 1], 'k--', label='Perfect calibration')
ax6B.fill_between(prob_pred, 
                  np.maximum(prob_true - 0.1, 0),
                  np.minimum(prob_true + 0.1, 1),
                  alpha=0.1, color='steelblue')

ax6B.legend(fontsize=8)
ax6B.set_title('Panel B: Calibration Curve (5-fold CV OOF)')
ax6B.set_xlabel('Predicted Probability'); ax6B.set_ylabel('Observed Proportion')
ax6B.set_xlim(0, 1); ax6B.set_ylim(0, 1)

# Panel C: DCA Decision Curve (5-fold CV OOF predictions)
ax6C = axes6[1, 0]
thresholds = np.linspace(0.01, 0.99, 99)

nb_treat_all = y.mean() - (1 - y.mean()) * (thresholds / (1 - thresholds))
nb_treat_none = np.zeros_like(thresholds)

nb_model = []
nb_clinical = []
for t in thresholds:
    treat = oof_full >= t
    tp = ((treat == 1) & (y == 1)).sum()
    fp = ((treat == 1) & (y == 0)).sum()
    n = len(y)
    nb = (tp / n) - (fp / n) * (t / (1 - t))
    nb_model.append(nb)
    
    treat_c = oof_clinical >= t
    tp_c = ((treat_c == 1) & (y == 1)).sum()
    fp_c = ((treat_c == 1) & (y == 0)).sum()
    nb_c = (tp_c / n) - (fp_c / n) * (t / (1 - t))
    nb_clinical.append(nb_c)

ax6C.plot(thresholds, nb_treat_all, 'b-', label='Treat All', linewidth=1.5)
ax6C.plot(thresholds, nb_treat_none, 'k-', label='Treat None', linewidth=1.5)
ax6C.plot(thresholds, nb_clinical, 'g-', label='Clinical only (OOF)', linewidth=1.5)
ax6C.plot(thresholds, nb_model, 'r-', label='Clinical + NK (OOF)', linewidth=2)
ax6C.set_xlabel('Threshold Probability'); ax6C.set_ylabel('Net Benefit')
ax6C.set_title('Panel C: Decision Curve Analysis (5-fold CV OOF)')
ax6C.legend(fontsize=8)
ax6C.set_xlim(0, 1)

# Panel D: AUC comparison bar (含 IRS 模型，审计 P3 修复)
ax6D = axes6[1, 1]
bars = [auc_clinical, auc_full, auc_irs]
bar_labels = ['Clinical only', 'Clinical + NK-Locked', 'Clinical + IRS']
colors_bar = ['#95E1D3', '#FF6B6B', '#9B59B6']

ax6D.bar(bar_labels, bars, color=colors_bar, width=0.5)
for i, (v, l) in enumerate(zip(bars, bar_labels)):
    ax6D.text(i, v + 0.01, f'AUC={v:.4f}', ha='center', fontsize=10, fontweight='bold')

# Paired bootstrap test for AUC difference (NK vs Clinical) - refit-based
np.random.seed(42)
n_boot_auc = 2000
boot_diff = []
boot_diff_irs = []
boot_diff_irs_clin = []
ci_diff_lower = np.nan  # 显式初始化，避免 in dir() 反模式
ci_diff_upper = np.nan
p_boot = np.nan
p_boot_irs = np.nan
p_boot_irs_clin = np.nan

for _ in range(n_boot_auc):
    idx = np.random.choice(len(y), size=len(y), replace=True)
    yb = y[idx]
    if len(set(yb)) < 2:
        continue

    try:
        m_c = LogisticRegression(penalty=None, max_iter=5000, solver='lbfgs')
        m_c.fit(X_clinical[idx], yb)
        auc_c_b = roc_auc_score(yb, m_c.predict_proba(X_clinical[idx])[:, 1])
    except:
        continue

    try:
        m_n = LogisticRegression(penalty=None, max_iter=5000, solver='lbfgs')
        m_n.fit(X_full[idx], yb)
        auc_n_b = roc_auc_score(yb, m_n.predict_proba(X_full[idx])[:, 1])
    except:
        continue

    try:
        m_i = LogisticRegression(penalty=None, max_iter=5000, solver='lbfgs')
        m_i.fit(X_irs[idx], yb)
        auc_i_b = roc_auc_score(yb, m_i.predict_proba(X_irs[idx])[:, 1])
    except:
        continue

    boot_diff.append(auc_n_b - auc_c_b)
    boot_diff_irs.append(auc_i_b - auc_n_b)
    boot_diff_irs_clin.append(auc_i_b - auc_c_b)

boot_diff = np.array(boot_diff)
boot_diff_irs = np.array(boot_diff_irs)
boot_diff_irs_clin = np.array(boot_diff_irs_clin)


def _bootstrap_p_twosided(d):
    """标准双尾 bootstrap p 值: 2 * min(P(≤0), P(≥0))"""
    if len(d) < 100:
        return np.nan
    p_le = np.mean(d <= 0)
    p_ge = np.mean(d >= 0)
    p = 2.0 * min(p_le, p_ge)
    return min(max(p, 1.0 / len(d)), 1.0)


if len(boot_diff) > 100:
    p_boot = _bootstrap_p_twosided(boot_diff)
    ci_diff_lower = np.percentile(boot_diff, 2.5)
    ci_diff_upper = np.percentile(boot_diff, 97.5)

if len(boot_diff_irs) > 100:
    p_boot_irs = _bootstrap_p_twosided(boot_diff_irs)

if len(boot_diff_irs_clin) > 100:
    p_boot_irs_clin = _bootstrap_p_twosided(boot_diff_irs_clin)

print(f"  (bootstrap successful iterations: {len(boot_diff)}/{n_boot_auc})")

print(f"\nAUC comparison (paired bootstrap):")
print(f"  Clinical only AUC = {auc_clinical:.4f}")
print(f"  Clinical + NK-Locked AUC = {auc_full:.4f} (n={len(y)})")
print(f"  Clinical + IRS AUC = {auc_irs:.4f} (n={len(y)})")
print(f"  NK vs Clinical: ΔAUC = {auc_full - auc_clinical:.4f}, p = {p_boot:.4f}")
print(f"  IRS vs NK: ΔAUC = {auc_irs - auc_full:.4f}, p = {p_boot_irs:.4f}")

ax6D.set_title(f'Panel D: AUC Comparison\n'
               f'NK vs Clinical: p={p_boot:.4f}\n'
               f'IRS vs NK: p={p_boot_irs:.4f}')
ax6D.set_ylabel('AUC')
ax6D.set_ylim(0, 1)

plt.tight_layout()
fig6.savefig(cfg.result_path('Fig6_clinical_model.pdf'), dpi=150, bbox_inches='tight')
fig6.savefig(cfg.result_path('Fig6_clinical_model.png'), dpi=150, bbox_inches='tight')
plt.close(fig6)
print("Fig6_clinical_model saved")

# ── Save clinical model statistics to CSV (含 IRS 对比，审计 P3) ──
fig6_stats = {
    'auc_clinical': auc_clinical if has_clinical else np.nan,
    'auc_clinical_nk': auc_full,
    'auc_clinical_irs': auc_irs,
    'auc_oof_clinical': auc_oof_clinical if has_clinical else np.nan,
    'auc_oof_clinical_nk': auc_oof_full,
    'delta_auc_nk': auc_full - (auc_clinical if has_clinical else 0),
    'delta_auc_oof_nk_vs_clinical': auc_oof_full - (auc_oof_clinical if has_clinical else 0),
    'delta_auc_irs_vs_nk': auc_irs - auc_full,
    'bootstrap_p_nk_vs_clinical': p_boot,
    'bootstrap_p_irs_vs_clinical': p_boot_irs_clin,
    'bootstrap_p_irs_vs_nk': p_boot_irs,
    'bootstrap_ci_lower': ci_diff_lower,
    'bootstrap_ci_upper': ci_diff_upper,
    'n_patients': len(y),
    'n_patients_irs': len(y_irs),
    'n_pcr': int(y.sum()),
    'n_non_mpr': int((1-y).sum()),
}
fig6_stats_df = pd.DataFrame({'metric': list(fig6_stats.keys()), 'value': list(fig6_stats.values())})
fig6_stats_df.to_csv(cfg.result_path('fig6_model_stats.csv'), index=False)
print("fig6_model_stats.csv saved (含 IRS 对比)")

print("[step5_clinical_translation] Done.")
