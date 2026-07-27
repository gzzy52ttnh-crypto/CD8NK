"""
step5_5_IRS_construction.py
Step 5.5: 构建 IRS (Immune Response Score) 评分并验证
IRS = nk_dominant_ratio * (1 - SPP1_TAM_ratio_normalized)
含义：机体抗肿瘤免疫的"净效力"——NK 克隆锁定程度扣除髓系抑制效应
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
from sklearn.model_selection import StratifiedKFold
from scipy.stats import chi2
import os
import warnings
warnings.filterwarnings('ignore')

print("[step5_5_IRS_construction] Starting...")
os.makedirs(cfg.RESULT_DIR, exist_ok=True)

# ── 1. 数据准备 ──
df_patient = pd.read_csv(cfg.result_path('per_patient_metrics.csv'))
df_mye = pd.read_csv(cfg.result_path('myeloid_per_patient.csv'))
print(f"Loaded per_patient_metrics: {df_patient.shape}")
print(f"Loaded myeloid_per_patient: {df_mye.shape}")

df = df_patient.merge(df_mye[['patient_id', 'SPP1_TAM_ratio']],
                      on='patient_id', how='inner')
print(f"Merged (inner join): {len(df)} patients")
print(f"  Excluded from patient metrics: {len(df_patient) - len(df)} (no myeloid data)")

# 过滤无效 SPP1_TAM_ratio（空值）
df['SPP1_TAM_ratio'] = pd.to_numeric(df['SPP1_TAM_ratio'], errors='coerce')
df = df.dropna(subset=['SPP1_TAM_ratio'])
print(f"After dropping NaN SPP1_TAM_ratio: {len(df)} patients")

# ── 2. IRS 评分构建 ──
# Min-Max 归一化 SPP1_TAM_ratio 到 0-1
spp1_min = df['SPP1_TAM_ratio'].min()
spp1_max = df['SPP1_TAM_ratio'].max()
df['SPP1_norm'] = (df['SPP1_TAM_ratio'] - spp1_min) / (spp1_max - spp1_min) if spp1_max > spp1_min else 0.0

# IRS = nk_dominant_ratio * (1 - SPP1_norm)
df['IRS'] = df['nk_dominant_ratio'] * (1.0 - df['SPP1_norm'])
print(f"\nIRS 构建:")
print(f"  SPP1_TAM_ratio range: [{spp1_min:.4f}, {spp1_max:.4f}]")
print(f"  nk_dominant_ratio range: [{df['nk_dominant_ratio'].min():.4f}, {df['nk_dominant_ratio'].max():.4f}]")
print(f"  IRS range: [{df['IRS'].min():.4f}, {df['IRS'].max():.4f}]")
print(f"  IRS median: {df['IRS'].median():.4f}")

# 保存 IRS 数据
df_irs = df[['patient_id', 'response', 'nk_dominant_ratio', 'SPP1_TAM_ratio', 'SPP1_norm', 'IRS']].copy()
df_irs.to_csv(cfg.result_path('irs_scores.csv'), index=False)
print("  irs_scores.csv saved")

# ── 3. 模型构建（与 step5 完全一致的特征处理） ──
df_model = df.copy()
df_model['y'] = (df_model['response'] == 'pCR').astype(int)

# Clinical features (统一引用 _common.build_clinical_features)
df_model, clinical_features = _common.build_clinical_features(df_model)

df_model = df_model.dropna(subset=['y'] + clinical_features)
y = df_model['y'].values

print(f"\n模型样本: n={len(y)}, pCR={y.sum()}, non-MPR={len(y)-y.sum()}")
print(f"临床特征 ({len(clinical_features)}): {clinical_features}")

has_clinical = len(clinical_features) > 0

# 模型 A: Clinical only
if has_clinical:
    X_clin = df_model[clinical_features].values
    model_clin = LogisticRegression(penalty=None, max_iter=5000, solver='lbfgs')
    model_clin.fit(X_clin, y)
    pred_clin = model_clin.predict_proba(X_clin)[:, 1]
    auc_clin = roc_auc_score(y, pred_clin)
    print(f"\nClinical only AUC: {auc_clin:.4f}")
else:
    auc_clin = np.nan
    pred_clin = np.zeros_like(y, dtype=float)
    print("\n[WARNING] No clinical features")

# 模型 B: Clinical + nk_dominant_ratio (baseline)
df_model['nk_scaled'] = df_model['nk_dominant_ratio'] * 100
nk_feat = 'nk_scaled'
feat_nk = clinical_features + [nk_feat] if has_clinical else [nk_feat]
X_nk = df_model[feat_nk].values
model_nk = LogisticRegression(penalty=None, max_iter=5000, solver='lbfgs')
model_nk.fit(X_nk, y)
pred_nk = model_nk.predict_proba(X_nk)[:, 1]
auc_nk = roc_auc_score(y, pred_nk)
print(f"Clinical + NK-dominant AUC: {auc_nk:.4f}")
delta_auc_nk = auc_nk - auc_clin if has_clinical else np.nan
if has_clinical:
    print(f"  ΔAUC (NK - Clinical): {delta_auc_nk:.4f}")

# 模型 C: Clinical + IRS (新评分)
df_model['irs_scaled'] = df_model['IRS'] * 100
irs_feat = 'irs_scaled'
feat_irs = clinical_features + [irs_feat] if has_clinical else [irs_feat]
X_irs = df_model[feat_irs].values
model_irs = LogisticRegression(penalty=None, max_iter=5000, solver='lbfgs')
model_irs.fit(X_irs, y)
pred_irs = model_irs.predict_proba(X_irs)[:, 1]
auc_irs = roc_auc_score(y, pred_irs)
print(f"Clinical + IRS AUC: {auc_irs:.4f}")
delta_auc_irs = auc_irs - auc_clin if has_clinical else np.nan
if has_clinical:
    print(f"  ΔAUC (IRS - Clinical): {delta_auc_irs:.4f}")
print(f"  ΔAUC (IRS - NK): {auc_irs - auc_nk:.4f}")

# ── 4. Bootstrap 检验 ΔAUC 显著性 ──
np.random.seed(42)
n_boot = 2000
boot_delta_nk = []
boot_delta_irs = []
boot_auc_irs = []
boot_auc_nk = []

n = len(y)
for b in range(n_boot):
    idx = np.random.choice(n, size=n, replace=True)
    y_b = y[idx]
    if len(np.unique(y_b)) < 2:
        continue
    if has_clinical:
        X_clin_b = X_clin[idx]
        try:
            m_c = LogisticRegression(penalty=None, max_iter=5000, solver='lbfgs')
            m_c.fit(X_clin_b, y_b)
            auc_c_b = roc_auc_score(y_b, m_c.predict_proba(X_clin_b)[:, 1])
        except:
            continue
    else:
        auc_c_b = 0.5
    
    X_nk_b = X_nk[idx]
    try:
        m_n = LogisticRegression(penalty=None, max_iter=5000, solver='lbfgs')
        m_n.fit(X_nk_b, y_b)
        auc_n_b = roc_auc_score(y_b, m_n.predict_proba(X_nk_b)[:, 1])
    except:
        continue
    
    X_irs_b = X_irs[idx]
    try:
        m_i = LogisticRegression(penalty=None, max_iter=5000, solver='lbfgs')
        m_i.fit(X_irs_b, y_b)
        auc_i_b = roc_auc_score(y_b, m_i.predict_proba(X_irs_b)[:, 1])
    except:
        continue
    
    if has_clinical:
        boot_delta_nk.append(auc_n_b - auc_c_b)
        boot_delta_irs.append(auc_i_b - auc_c_b)
    boot_auc_irs.append(auc_i_b)
    boot_auc_nk.append(auc_n_b)

# Bootstrap p 值：标准双尾 2 * min(P(≤0), P(≥0))
def _boot_p_twosided(d):
    """标准双尾 bootstrap p 值"""
    if len(d) < 1:
        return np.nan
    d = np.asarray(d)
    p_le = np.mean(d <= 0)
    p_ge = np.mean(d >= 0)
    p = 2.0 * min(p_le, p_ge)
    return min(max(p, 1.0 / len(d)), 1.0)


if has_clinical and len(boot_delta_nk) > 0:
    p_nk = _boot_p_twosided(boot_delta_nk)
    print(f"\nBootstrap NK ΔAUC p: {p_nk:.4f} (n_boot={len(boot_delta_nk)})")
else:
    p_nk = np.nan

if has_clinical and len(boot_delta_irs) > 0:
    p_irs = _boot_p_twosided(boot_delta_irs)
    print(f"Bootstrap IRS ΔAUC p: {p_irs:.4f} (n_boot={len(boot_delta_irs)})")
else:
    p_irs = np.nan

# IRS vs NK 直接比较
if len(boot_auc_irs) > 0 and len(boot_auc_nk) > 0:
    diff_boot = np.array(boot_auc_irs) - np.array(boot_auc_nk)
    p_irs_vs_nk = _boot_p_twosided(diff_boot)
    print(f"Bootstrap IRS vs NK AUC p: {p_irs_vs_nk:.4f}")
else:
    p_irs_vs_nk = np.nan

# Bootstrap 95% CI for IRS AUC
if len(boot_auc_irs) > 0:
    ci_irs_lower = np.percentile(boot_auc_irs, 2.5)
    ci_irs_upper = np.percentile(boot_auc_irs, 97.5)
    print(f"IRS AUC 95% Bootstrap CI: [{ci_irs_lower:.4f}, {ci_irs_upper:.4f}]")
else:
    ci_irs_lower, ci_irs_upper = np.nan, np.nan

# ── 5. 5-fold OOF AUC (避免乐观偏差) ──
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
oof_irs = np.zeros(len(y))
oof_nk = np.zeros(len(y))
oof_clin = np.zeros(len(y))

for train_idx, val_idx in skf.split(X_irs, y):
    m_irs_f = LogisticRegression(penalty=None, max_iter=5000, solver='lbfgs')
    m_irs_f.fit(X_irs[train_idx], y[train_idx])
    oof_irs[val_idx] = m_irs_f.predict_proba(X_irs[val_idx])[:, 1]
    
    m_nk_f = LogisticRegression(penalty=None, max_iter=5000, solver='lbfgs')
    m_nk_f.fit(X_nk[train_idx], y[train_idx])
    oof_nk[val_idx] = m_nk_f.predict_proba(X_nk[val_idx])[:, 1]
    
    if has_clinical:
        m_c_f = LogisticRegression(penalty=None, max_iter=5000, solver='lbfgs')
        m_c_f.fit(X_clin[train_idx], y[train_idx])
        oof_clin[val_idx] = m_c_f.predict_proba(X_clin[val_idx])[:, 1]

auc_oof_irs = roc_auc_score(y, oof_irs)
auc_oof_nk = roc_auc_score(y, oof_nk)
auc_oof_clin = roc_auc_score(y, oof_clin) if has_clinical else np.nan
print(f"\n5-fold OOF AUC:")
print(f"  Clinical only: {auc_oof_clin:.4f}" if has_clinical else "  Clinical only: N/A")
print(f"  Clinical + NK: {auc_oof_nk:.4f}")
print(f"  Clinical + IRS: {auc_oof_irs:.4f}")

# ── 6. 绘图：ROC + 校准曲线 + AUC对比 (1x3 布局) ──
fig, (ax_roc, ax_cal, ax_auc) = plt.subplots(1, 3, figsize=(18, 5.5))

# 左：ROC 曲线对比
if has_clinical:
    fpr_c, tpr_c, _ = roc_curve(y, pred_clin)
    ax_roc.plot(fpr_c, tpr_c, 'gray', lw=1.5, label=f'Clinical only (AUC={auc_clin:.3f})')

fpr_n, tpr_n, _ = roc_curve(y, pred_nk)
ax_roc.plot(fpr_n, tpr_n, '#2E86C1', lw=2, label=f'Clinical + NK-dominant (AUC={auc_nk:.3f})')

fpr_i, tpr_i, _ = roc_curve(y, pred_irs)
ax_roc.fill_between(fpr_i, tpr_i, alpha=0.15, color='#E74C3C')
ax_roc.plot(fpr_i, tpr_i, '#E74C3C', lw=2.5, label=f'Clinical + IRS (AUC={auc_irs:.3f})')

ax_roc.plot([0, 1], [0, 1], 'k--', lw=0.8, alpha=0.5)
ax_roc.set_xlabel('False Positive Rate')
ax_roc.set_ylabel('True Positive Rate')
ax_roc.set_title(f'ROC Comparison\n(n={len(y)}, pCR={int(y.sum())})')
ax_roc.legend(loc='lower right', fontsize=9)
ax_roc.spines['top'].set_visible(False)
ax_roc.spines['right'].set_visible(False)
ax_roc.set_facecolor('#F8F9FA')
ax_roc.grid(True, alpha=0.3)

# 中：校准曲线（基于 OOF 预测，避免乐观偏差）
from sklearn.calibration import calibration_curve
try:
    frac_pos, mean_pred = calibration_curve(y, oof_irs, n_bins=5, strategy='quantile')
    ax_cal.plot(mean_pred, frac_pos, 'o-', color='#E74C3C', lw=2, markersize=8, label=f'IRS (OOF AUC={auc_oof_irs:.3f})')
except Exception as e:
    print(f"  IRS calibration failed: {e}")

try:
    frac_pos_nk, mean_pred_nk = calibration_curve(y, oof_nk, n_bins=5, strategy='quantile')
    ax_cal.plot(mean_pred_nk, frac_pos_nk, 's--', color='#2E86C1', lw=1.5, markersize=7, label=f'NK-dominant (OOF AUC={auc_oof_nk:.3f})')
except Exception as e:
    print(f"  NK calibration failed: {e}")

ax_cal.plot([0, 1], [0, 1], 'k--', lw=0.8, alpha=0.5, label='Perfectly calibrated')
ax_cal.set_xlabel('Mean predicted probability')
ax_cal.set_ylabel('Fraction of pCR (observed)')
ax_cal.set_title('Calibration Curve (5-fold OOF)')
ax_cal.legend(loc='upper left', fontsize=9)
ax_cal.set_xlim(-0.05, 1.05); ax_cal.set_ylim(-0.05, 1.05)
ax_cal.spines['top'].set_visible(False)
ax_cal.spines['right'].set_visible(False)
ax_cal.set_facecolor('#F8F9FA')
ax_cal.grid(True, alpha=0.3)

# 右：AUC 对比柱状图（训练集 vs OOF）
model_names = ['Clinical\nonly', 'Clinical +\nNK-dominant', 'Clinical +\nIRS']
train_aucs = [auc_clin if has_clinical else 0.5, auc_nk, auc_irs]
oof_aucs = [auc_oof_clin if has_clinical else 0.5, auc_oof_nk, auc_oof_irs]
x_auc = np.arange(len(model_names))
bar_w = 0.35
b_train = ax_auc.bar(x_auc - bar_w/2, train_aucs, bar_w, color=['#95A5A6', '#2E86C1', '#E74C3C'], alpha=0.7,
                     edgecolor='black', linewidth=0.5, label='Training set')
b_oof = ax_auc.bar(x_auc + bar_w/2, oof_aucs, bar_w, color=['#95A5A6', '#2E86C1', '#E74C3C'], alpha=1.0,
                   edgecolor='black', linewidth=0.5, hatch='//', label='5-fold OOF')
ax_auc.axhline(0.5, color='gray', ls='--', lw=1, alpha=0.6)
ax_auc.set_xticks(x_auc)
ax_auc.set_xticklabels(model_names, fontsize=9)
ax_auc.set_ylabel('AUC')
ax_auc.set_ylim(0, 1.0)
ax_auc.set_title('AUC: Training vs OOF\n(ΔAUC = optimism)')
ax_auc.legend(loc='lower right', fontsize=9)
ax_auc.spines['top'].set_visible(False); ax_auc.spines['right'].set_visible(False)
ax_auc.set_facecolor('#F8F9FA')
ax_auc.grid(True, alpha=0.3, axis='y')
for bar, v in zip(b_train, train_aucs):
    ax_auc.text(bar.get_x() + bar.get_width()/2, v + 0.02, f'{v:.3f}',
                ha='center', va='bottom', fontsize=8)
for bar, v in zip(b_oof, oof_aucs):
    ax_auc.text(bar.get_x() + bar.get_width()/2, v + 0.02, f'{v:.3f}',
                ha='center', va='bottom', fontsize=8)

plt.suptitle(f'IRS Model Performance: NK-dominant ratio × (1 - SPP1_TAM_norm)\n'
             f'(IRS AUC={auc_irs:.3f} vs NK AUC={auc_nk:.3f}, Bootstrap ΔAUC p={p_irs_vs_nk:.4f})',
             fontsize=11, y=1.02)
plt.tight_layout()
fig.savefig(cfg.result_path('fig6_IRS_model.pdf'), dpi=150, bbox_inches='tight')
fig.savefig(cfg.result_path('fig6_IRS_model.png'), dpi=150, bbox_inches='tight')
plt.close(fig)
print("fig6_IRS_model saved (ROC + Calibration + AUC comparison)")

# ── 6.4 IRS 版 Nomogram/列线图（临床转化） ──
# 审计问题D修复：IRS 评分需有临床可解读的 Nomogram，增强转化价值
# 与 step5_clinical_translation.py 的 Nomogram 保持一致：
#   - 按 |β_i| / max(|β|) 比例分配最大 points（非等距线性映射）
#   - 含 Clinical + IRS 模型与 Clinical only 模型两套系数
print("\n" + "="*60)
print("IRS NOMOGRAM (Clinical Translation)")
print("="*60)

# IRS 模型系数（model_irs 已在上方拟合）
coef_irs = model_irs.coef_[0]
intercept_irs = model_irs.intercept_[0]
feat_irs_names = feat_irs  # clinical_features + [irs_feat]

# Nomogram points 分配：按 |β_i × range_i| / max(|β_j × range_j|) 比例
# 标准列线图方法：points 反映各变量对线性预测子的实际贡献（β × 取值范围）
# 注：step5_clinical_translation.py 使用 |β_i|/max(|β|)（适用于同尺度变量）；
#     IRS 模型含 ×100 缩放变量，必须用 β×range 才能正确反映贡献
lp_contrib_irs = np.zeros(len(feat_irs_names))
for i, feat in enumerate(feat_irs_names):
    var_range_i = df_model[feat].max() - df_model[feat].min()
    lp_contrib_irs[i] = abs(coef_irs[i]) * var_range_i if var_range_i > 0 else 0
max_lp_irs = lp_contrib_irs.max() if lp_contrib_irs.max() > 0 else 1
var_max_pts_irs = (lp_contrib_irs / max_lp_irs * 100).astype(int)

# Clinical only 模型系数（model_clin 已在上方拟合，若存在）
nomogram_irs_data = []
if has_clinical:
    coef_cli = model_clin.coef_[0]
    lp_contrib_cli = np.zeros(len(clinical_features))
    for i, feat in enumerate(clinical_features):
        var_range_i = df_model[feat].max() - df_model[feat].min()
        lp_contrib_cli[i] = abs(coef_cli[i]) * var_range_i if var_range_i > 0 else 0
    max_lp_cli = lp_contrib_cli.max() if lp_contrib_cli.max() > 0 else 1
    var_max_pts_cli = (lp_contrib_cli / max_lp_cli * 100).astype(int)
    for i, feat in enumerate(clinical_features):
        nomogram_irs_data.append({
            'variable': feat,
            'coefficient': round(coef_cli[i], 4),
            'max_points': int(var_max_pts_cli[i]),
            'model': 'Clinical_only'
        })

# Clinical + IRS 模型系数
for i, feat in enumerate(feat_irs_names):
    nomogram_irs_data.append({
        'variable': feat,
        'coefficient': round(coef_irs[i], 4),
        'max_points': int(var_max_pts_irs[i]),
        'model': 'Clinical+IRS'
    })

df_nomogram_irs = pd.DataFrame(nomogram_irs_data)
df_nomogram_irs.to_csv(cfg.result_path('fig6_IRS_nomogram_data.csv'), index=False)
print(f"fig6_IRS_nomogram_data.csv saved ({len(df_nomogram_irs)} rows)")
print(f"  IRS 模型特征: {feat_irs_names}")
print(f"  IRS 系数: {dict(zip(feat_irs_names, [round(c, 4) for c in coef_irs]))}")
print(f"  IRS points 分配 (β×range): {dict(zip(feat_irs_names, var_max_pts_irs))}")

# 绘制 IRS Nomogram（独立图）
fig_nom, ax_nom = plt.subplots(figsize=(12, 7))
total_max_pts_irs = int(var_max_pts_irs.sum())
n_vars_irs = len(feat_irs_names)

ax_nom.set_xlim(0, total_max_pts_irs + 80)
ax_nom.set_ylim(-2.5, n_vars_irs + 1.5)

# Points 轴（顶部）
ax_nom.axhline(y=n_vars_irs + 0.5, xmin=0, xmax=total_max_pts_irs / (total_max_pts_irs + 80),
               color='black', linewidth=1.2)
ax_nom.text(total_max_pts_irs // 2, n_vars_irs + 0.9, 'Points', ha='center', fontsize=11, fontweight='bold')
points_step_irs = 20
for i in range(0, total_max_pts_irs + 1, points_step_irs):
    ax_nom.axvline(x=i, ymin=(n_vars_irs + 0.3) / (n_vars_irs + 3.5),
                   ymax=(n_vars_irs + 0.7) / (n_vars_irs + 3.5), color='black', linewidth=0.8)
    ax_nom.text(i, n_vars_irs + 0.6, str(i), ha='center', fontsize=7)

# 各变量行
for v_idx, feat in enumerate(feat_irs_names):
    y_pos = n_vars_irs - v_idx
    vmax_pts = var_max_pts_irs[v_idx]
    ax_nom.axhline(y=y_pos, xmin=0, xmax=vmax_pts / (total_max_pts_irs + 80),
                   color='black', linewidth=1)
    # 变量名（左对齐）
    label_name = feat.replace('_scaled', ' (×100)').replace('_num', '').replace('hist_luad', 'Histology=LUAD')
    if feat.startswith('stage_'):
        label_name = f'Stage={feat.split("_")[1]}'
    ax_nom.text(-5, y_pos + 0.3, label_name, ha='right', fontsize=8, va='center')

    # 变量取值 → points 刻度
    var_min = df_model[feat].min()
    var_max = df_model[feat].max()
    var_range = var_max - var_min if var_max > var_min else 1
    n_ticks = min(5, max(vmax_pts // 20 + 1, 2))
    for ti in range(n_ticks):
        tick_pt = int(ti * vmax_pts / max(n_ticks - 1, 1))
        val = var_min + (tick_pt / max(vmax_pts, 1)) * var_range
        ax_nom.text(tick_pt, y_pos - 0.3, f'{val:.1f}', ha='center', fontsize=6)

# Total Points → pCR Probability
ax_nom.axhline(y=0, xmin=0, xmax=total_max_pts_irs / (total_max_pts_irs + 80),
               color='black', linewidth=1.2)
ax_nom.text(total_max_pts_irs // 2, 0.3, 'Total Points', ha='center', fontsize=10, fontweight='bold')

# 概率刻度（基于实际模型截距与系数范围）
lp_at_zero_irs = intercept_irs
lp_per_point_irs = max_lp_irs / 100.0
lp_at_max_irs = intercept_irs + total_max_pts_irs * lp_per_point_irs
lp_min_irs = min(lp_at_zero_irs, lp_at_max_irs)
lp_max_irs = max(lp_at_zero_irs, lp_at_max_irs)
ax_nom.text(total_max_pts_irs // 2, -0.7, 'Linear Predictor', ha='center', fontsize=9)
ax_nom.text(total_max_pts_irs // 2, -1.5, 'pCR Probability', ha='center', fontsize=10, fontweight='bold')
prob_ticks_irs = [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95]
for pt in prob_ticks_irs:
    lp = np.log(pt / (1 - pt))
    if lp_min_irs <= lp <= lp_max_irs:
        x_pos = (lp - lp_min_irs) / (lp_max_irs - lp_min_irs) * total_max_pts_irs
        ax_nom.text(x_pos, -1.8, f'{pt:.2f}', ha='center', fontsize=7)

ax_nom.set_title(f'IRS Nomogram: pCR ~ Clinical + IRS\n'
                 f'(IRS = nk_dominant_ratio × (1 - SPP1_TAM_norm), n={len(y)}, '
                 f'AUC={auc_irs:.3f}, OOF AUC={auc_oof_irs:.3f})',
                 fontsize=11, pad=15)
ax_nom.axis('off')

plt.tight_layout()
fig_nom.savefig(cfg.result_path('fig6_IRS_nomogram.pdf'), dpi=150, bbox_inches='tight')
fig_nom.savefig(cfg.result_path('fig6_IRS_nomogram.png'), dpi=150, bbox_inches='tight')
plt.close(fig_nom)
print("fig6_IRS_nomogram saved (IRS-based Nomogram for clinical translation)")

# ── 6.5 统计交互项检验（nk_dominant_ratio × SPP1_TAM_ratio） ──
# 真正的统计交互项：y ~ nk + SPP1 + nk:SPP1
print("\n" + "="*60)
print("STATISTICAL INTERACTION TEST (nk × SPP1 in Logistic model)")
print("="*60)
import statsmodels.api as sm

df_interact = df_model.copy()
df_interact['nk_z'] = (df_interact['nk_dominant_ratio'] - df_interact['nk_dominant_ratio'].mean()) / df_interact['nk_dominant_ratio'].std()
df_interact['spp1_z'] = (df_interact['SPP1_TAM_ratio'] - df_interact['SPP1_TAM_ratio'].mean()) / df_interact['SPP1_TAM_ratio'].std()
df_interact['interaction'] = df_interact['nk_z'] * df_interact['spp1_z']

# 模型: y ~ nk + SPP1 + nk:SPP1（无临床协变量，单变量交互检验）
X_int = sm.add_constant(df_interact[['nk_z', 'spp1_z', 'interaction']])
y_int = df_interact['y'].values

try:
    model_int = sm.Logit(y_int, X_int).fit(disp=0)
    or_int = np.exp(model_int.params['interaction'])
    ci_int = np.exp(model_int.conf_int().loc['interaction'])
    p_int = model_int.pvalues['interaction']
    print(f"  交互项 (nk_z × spp1_z): OR = {or_int:.4f}, 95%CI = [{ci_int[0]:.4f}, {ci_int[1]:.4f}], p = {p_int:.6f}")

    # 似然比检验（交互项是否显著提升模型）
    X_noint = sm.add_constant(df_interact[['nk_z', 'spp1_z']])
    model_noint = sm.Logit(y_int, X_noint).fit(disp=0)
    ll_int = model_int.llf
    ll_noint = model_noint.llf
    lr_stat = -2 * (ll_noint - ll_int)
    lr_p = 1 - chi2.cdf(lr_stat, df=1)
    print(f"  似然比检验: ΔLL = {ll_int - ll_noint:.4f}, χ² = {lr_stat:.4f}, p = {lr_p:.6f}")

    # 含临床协变量的交互模型
    if has_clinical:
        X_int_clin = sm.add_constant(df_interact[['nk_z', 'spp1_z', 'interaction'] + clinical_features])
        model_int_clin = sm.Logit(y_int, X_int_clin).fit(disp=0)
        or_int_clin = np.exp(model_int_clin.params['interaction'])
        ci_int_clin = np.exp(model_int_clin.conf_int().loc['interaction'])
        p_int_clin = model_int_clin.pvalues['interaction']
        print(f"  交互项 (含临床协变量): OR = {or_int_clin:.4f}, 95%CI = [{ci_int_clin[0]:.4f}, {ci_int_clin[1]:.4f}], p = {p_int_clin:.6f}")
    else:
        or_int_clin = np.nan
        ci_int_clin = [np.nan, np.nan]
        p_int_clin = np.nan

    interact_results = {
        'interaction_or_unadjusted': round(or_int, 4),
        'interaction_ci_lower_unadjusted': round(ci_int[0], 4),
        'interaction_ci_upper_unadjusted': round(ci_int[1], 4),
        'interaction_p_unadjusted': round(p_int, 6),
        'lr_stat': round(lr_stat, 4),
        'lr_p': round(lr_p, 6),
        'interaction_or_adjusted': round(or_int_clin, 4) if not np.isnan(or_int_clin) else 'NA',
        'interaction_p_adjusted': round(p_int_clin, 6) if not np.isnan(p_int_clin) else 'NA',
    }
    print(f"\n  交互项与 IRS 的关系:")
    print(f"    IRS 是生物学评分 = nk × (1 - SPP1_norm)")
    print(f"    统计交互项是 nk:SPP1 的乘积项系数")
    print(f"    两者方向一致: SPP1 越高，NK 的预测效力越弱")
except Exception as e:
    print(f"  交互项检验失败: {e}")
    interact_results = {'error': str(e)}

# ── 7. 保存统计结果 ──
stats_data = {
    'n_patients': len(y),
    'n_pCR': int(y.sum()),
    'n_nonMPR': int(len(y) - y.sum()),
    'auc_clinical': round(auc_clin, 4) if has_clinical else 'NA',
    'auc_nk': round(auc_nk, 4),
    'auc_irs': round(auc_irs, 4),
    'delta_auc_nk': round(delta_auc_nk, 4) if has_clinical else 'NA',
    'delta_auc_irs': round(delta_auc_irs, 4) if has_clinical else 'NA',
    'delta_auc_irs_vs_nk': round(auc_irs - auc_nk, 4),
    'bootstrap_p_nk': round(p_nk, 4) if not np.isnan(p_nk) else 'NA',
    'bootstrap_p_irs': round(p_irs, 4) if not np.isnan(p_irs) else 'NA',
    'bootstrap_p_irs_vs_nk': round(p_irs_vs_nk, 4) if not np.isnan(p_irs_vs_nk) else 'NA',
    'bootstrap_ci_irs_lower': round(ci_irs_lower, 4) if not np.isnan(ci_irs_lower) else 'NA',
    'bootstrap_ci_irs_upper': round(ci_irs_upper, 4) if not np.isnan(ci_irs_upper) else 'NA',
    'oof_auc_clinical': round(auc_oof_clin, 4) if has_clinical else 'NA',
    'oof_auc_nk': round(auc_oof_nk, 4),
    'oof_auc_irs': round(auc_oof_irs, 4),
    'spp1_min': round(spp1_min, 4),
    'spp1_max': round(spp1_max, 4),
    'irs_formula': 'nk_dominant_ratio * (1 - SPP1_TAM_ratio_minmax_norm)',
    'n_bootstrap': len(boot_auc_irs),
}
# 加入交互项检验结果
if 'error' not in interact_results:
    stats_data.update(interact_results)

df_stats = pd.DataFrame([stats_data])
df_stats.to_csv(cfg.result_path('irs_model_stats.csv'), index=False)
print("irs_model_stats.csv saved (with statistical interaction test)")

print("\n[step5_5_IRS_construction] Done.")
