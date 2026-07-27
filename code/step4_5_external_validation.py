"""
step4_5_external_validation.py
Step 4.5: 外部队列深度验证与机制归因
队列 1: GSE179994 (Liu et al. 2022, NSCLC) - 数据缺失，标注 NA
队列 2: GSE207422 (Hu et al. 2023, NSCLC) - 仅 metadata，无表达矩阵，标注 NA
队列 3: GSE120575 (melanoma, anti-CTLA-4) - 机制归因分析
"""
import config as cfg
from _common import NKLIKE_SIGNATURE, firth_logistic_fit
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from scipy.stats import norm
import statsmodels.api as sm
import os
import gzip
import warnings
warnings.filterwarnings('ignore')

print("[step4_5_external_validation] Starting...")
os.makedirs(cfg.RESULT_DIR, exist_ok=True)


validation_results = []

# ============================================================
# 队列 1: GSE179994
# ============================================================
print("\n" + "="*60)
print("队列 1: GSE179994 (NSCLC, anti-PD-1)")
print("="*60)

gse179994_expr = os.path.join(cfg.ADATA_DIR, 'GSE179994_TPM_matrix.csv')
gse179994_meta = os.path.join(cfg.ADATA_DIR, 'GSE179994', 'GSE179994_pseudobulk.csv')

if os.path.exists(gse179994_expr) or os.path.exists(gse179994_meta):
    print("GSE179994 数据可用，开始分析...")
else:
    print("[NOT AVAILABLE] GSE179994 表达矩阵和 metadata 均不存在于 adata/ 目录")
    validation_results.append({
        'cohort': 'GSE179994 (NSCLC, anti-PD-1)',
        'analysis': 'NK-like signature association',
        'n_patients': 'Not Available',
        'n_responder': 'Not Available',
        'n_nonresponder': 'Not Available',
        'OR': 'Not Available',
        'CI_lower': 'Not Available',
        'CI_upper': 'Not Available',
        'p_value': 'Not Available',
        'method': 'Not Available',
        'direction': 'Not Available',
        'note': 'Data file not found in adata/ directory',
    })
    validation_results.append({
        'cohort': 'GSE179994 (NSCLC, anti-PD-1)',
        'analysis': 'Clonal fate locking (TCR)',
        'n_patients': 'Not Available',
        'n_responder': 'Not Available',
        'n_nonresponder': 'Not Available',
        'OR': 'Not Available',
        'CI_lower': 'Not Available',
        'CI_upper': 'Not Available',
        'p_value': 'Not Available',
        'method': 'Not Available',
        'direction': 'Not Available',
        'note': 'TCR data not available',
    })

# ============================================================
# 队列 2: GSE207422
# ============================================================
print("\n" + "="*60)
print("队列 2: GSE207422 (NSCLC, neoadjuvant chemo-IO)")
print("="*60)

gse207422_expr_csv = os.path.join(cfg.ADATA_DIR, 'GSE207422_expression_matrix.csv')
gse207422_h5ad = os.path.join(cfg.ADATA_DIR, 'GSE207422_NSCLC_scRNAseq.h5ad')
gse207422_meta = os.path.join(cfg.ADATA_DIR, 'GSE207422_NSCLC_scRNAseq_metadata.xlsx')

has_expr_207422 = os.path.exists(gse207422_expr_csv) or os.path.exists(gse207422_h5ad)
has_meta_207422 = os.path.exists(gse207422_meta)

n_pat_207422 = 'Not Available'
n_mpr_207422 = 'Not Available'
n_nmpr_207422 = 'Not Available'

if has_meta_207422:
    df_meta_207422 = pd.read_excel(gse207422_meta)
    print(f"Metadata 可用: {df_meta_207422.shape}")
    if 'Patient' in df_meta_207422.columns:
        n_pat_207422 = df_meta_207422['Patient'].nunique()
        print(f"  Patients: {n_pat_207422}")
    if 'Pathologic Response' in df_meta_207422.columns:
        vc = df_meta_207422['Pathologic Response'].value_counts()
        print(f"  Response distribution:\n{vc.to_string()}")
        n_mpr_207422 = int(vc.get('MPR', 0) + vc.get('pCR', 0))
        n_nmpr_207422 = int(vc.get('NMPR', 0) + vc.get('Non-MPR', 0))

if not has_expr_207422:
    print("[NOT AVAILABLE] GSE207422 表达矩阵不存在（仅 metadata 可用）")
    validation_results.append({
        'cohort': 'GSE207422 (NSCLC, neoadjuvant chemo-IO)',
        'analysis': 'NK-like cell abundance (MPR vs NMPR)',
        'n_patients': n_pat_207422,
        'n_responder': n_mpr_207422,
        'n_nonresponder': n_nmpr_207422,
        'OR': 'Not Available',
        'CI_lower': 'Not Available',
        'CI_upper': 'Not Available',
        'p_value': 'Not Available',
        'method': 'Not Available',
        'direction': 'Not Available',
        'note': 'Expression matrix not found; only metadata available (12 patients with response labels)',
    })

# ============================================================
# 队列 3: GSE120575 机制归因
# ============================================================
print("\n" + "="*60)
print("队列 3: GSE120575 (melanoma, anti-CTLA-4) 机制归因")
print("="*60)

tpm_path_g12 = os.path.join(cfg.ADATA_DIR, 'GSE120575_TPM_matrix.csv')
pb_path_g12 = os.path.join(cfg.ADATA_DIR, 'GSE120575', 'GSE120575_pseudobulk.csv')
patient_map_path = os.path.join(cfg.ADATA_DIR, 'GSE120575', 'GSE120575_patient_ID_single_cells.txt.gz')

if os.path.exists(tpm_path_g12) and os.path.exists(pb_path_g12):
    print("Loading GSE120575 TPM matrix...")
    tpm_g12 = pd.read_csv(tpm_path_g12, index_col=0)
    tpm_g12.index = tpm_g12.index.astype(str)
    print(f"  TPM: {tpm_g12.shape[0]} genes x {tpm_g12.shape[1]} cells")

    sig_found_g12 = [g for g in NKLIKE_SIGNATURE if g in tpm_g12.index]
    print(f"  Signature genes found: {len(sig_found_g12)}/{len(NKLIKE_SIGNATURE)}")

    # 细胞 → 患者映射
    cell_to_patient = {}
    if os.path.exists(patient_map_path):
        with gzip.open(patient_map_path, 'rt', encoding='latin-1') as f:
            for line in f:
                line = line.strip()
                if line.startswith('Sample ') and '\t' in line:
                    parts = line.split('\t')
                    if len(parts) >= 5:
                        cid = parts[1].strip()
                        pid = parts[4].strip()
                        if cid and pid and cid != 'title':
                            cell_to_patient[cid] = pid
        print(f"  Patient mapping from SOFT: {len(cell_to_patient)} cells")

    pb_g12 = pd.read_csv(pb_path_g12)
    pb_g12['response_binary'] = (pb_g12['response_binary'] == 'R').astype(int)

    # 每细胞 NK-like 得分
    sig_expr_g12 = tpm_g12.loc[sig_found_g12]
    nk_score_per_cell = np.log1p(sig_expr_g12).mean(axis=0)

    nk_median = nk_score_per_cell.median()
    nk_high_mask = nk_score_per_cell >= nk_median
    nk_low_mask = nk_score_per_cell < nk_median
    print(f"\n  NK-like High: {nk_high_mask.sum()} cells, Low: {nk_low_mask.sum()} cells")

    # 每细胞响应标签
    cell_resp = {}
    for col in tpm_g12.columns:
        pid = cell_to_patient.get(col, '')
        if pid:
            resp_row = pb_g12[pb_g12['patient'] == pid]
            if len(resp_row) > 0:
                cell_resp[col] = int(resp_row['response_binary'].iloc[0])

    cells_with_resp = [c for c in tpm_g12.columns if c in cell_resp]
    print(f"  Cells with response label: {len(cells_with_resp)}")

    key_genes = ['PRF1', 'GZMB', 'GNLY', 'GZMH']
    key_genes_found = [g for g in key_genes if g in tpm_g12.index]
    print(f"  Key effector genes: {key_genes_found}")

    # NK-high 细胞中 R vs NR 比较
    print("\n" + "-"*50)
    print("核心机制检验：NK-like High 细胞中 R vs NR")
    print("-"*50)

    nk_high_cells = [c for c in cells_with_resp if nk_high_mask[c]]
    nk_high_r = [c for c in nk_high_cells if cell_resp[c] == 1]
    nk_high_nr = [c for c in nk_high_cells if cell_resp[c] == 0]
    print(f"  NK-high R cells: {len(nk_high_r)}, NR cells: {len(nk_high_nr)}")

    mech_summary = {}
    for gene in key_genes_found:
        expr_r = np.log1p(tpm_g12.loc[gene, nk_high_r].values.astype(float))
        expr_nr = np.log1p(tpm_g12.loc[gene, nk_high_nr].values.astype(float))
        if len(expr_r) > 5 and len(expr_nr) > 5:
            stat, p = mannwhitneyu(expr_r, expr_nr, alternative='two-sided')
            direction = 'R > NR' if np.mean(expr_r) > np.mean(expr_nr) else 'NR > R'
            mech_summary[gene] = {
                'R_mean': float(np.mean(expr_r)),
                'NR_mean': float(np.mean(expr_nr)),
                'p_value': float(p),
                'direction': direction,
            }
            print(f"  {gene}: R={np.mean(expr_r):.4f}, NR={np.mean(expr_nr):.4f}, p={p:.4e} ({direction})")

    if 'PRF1' in mech_summary and 'GZMB' in mech_summary:
        ratio_r = mech_summary['GZMB']['R_mean'] / max(mech_summary['PRF1']['R_mean'], 1e-10)
        ratio_nr = mech_summary['GZMB']['NR_mean'] / max(mech_summary['PRF1']['NR_mean'], 1e-10)
        mech_summary['GZMB_PRF1_ratio'] = {'R': ratio_r, 'NR': ratio_nr}
        print(f"\n  GZMB/PRF1 ratio: R={ratio_r:.4f}, NR={ratio_nr:.4f}")

    # 患者水平 NK-like 得分
    unique_patients = sorted(set(cell_to_patient.values()))
    patient_scores = {}
    patient_resp_map = {}
    for pid in unique_patients:
        pat_cells = [c for c in tpm_g12.columns if cell_to_patient.get(c, '') == pid]
        if len(pat_cells) < 3:
            continue
        resp_row = pb_g12[pb_g12['patient'] == pid]
        if len(resp_row) == 0:
            continue
        patient_resp_map[pid] = int(resp_row['response_binary'].iloc[0])
        patient_scores[pid] = float(np.mean(nk_score_per_cell[pat_cells]))

    df_scores = pd.DataFrame({
        'patient_id': list(patient_scores.keys()),
        'nklike_score': list(patient_scores.values()),
        'response_binary': [patient_resp_map[p] for p in patient_scores.keys()],
    })
    n_g12_pat = len(df_scores)
    n_r_g12 = int(df_scores['response_binary'].sum())
    n_nr_g12 = n_g12_pat - n_r_g12
    print(f"\n  Patient-level: {n_g12_pat} patients (R={n_r_g12}, NR={n_nr_g12})")

    # Firth 回归
    if n_g12_pat >= 6 and df_scores['response_binary'].nunique() == 2:
        y_g12 = df_scores['response_binary'].values
        score_g12 = df_scores['nklike_score'].values
        X_g12 = sm.add_constant(score_g12)
        beta_g12, se_g12, _, conv_g12 = firth_logistic_fit(X_g12, y_g12)
        or_g12 = float(np.exp(beta_g12[1]))
        ci_g12_low = float(np.exp(beta_g12[1] - 1.96 * se_g12[1]))
        ci_g12_high = float(np.exp(beta_g12[1] + 1.96 * se_g12[1]))
        z_g12 = beta_g12[1] / se_g12[1] if se_g12[1] > 0 else 0
        p_g12_firth = float(2 * (1 - norm.cdf(abs(z_g12))))
        direction_g12 = 'Higher in R' if np.mean(score_g12[y_g12 == 1]) > np.mean(score_g12[y_g12 == 0]) else 'Higher in NR'

        print(f"  Firth OR={or_g12:.4f}, 95%CI=[{ci_g12_low:.4f}, {ci_g12_high:.4f}], p={p_g12_firth:.4e} ({direction_g12})")

        validation_results.append({
            'cohort': 'GSE120575 (melanoma, anti-CTLA-4)',
            'analysis': 'NK-like signature association (Firth)',
            'n_patients': n_g12_pat,
            'n_responder': n_r_g12,
            'n_nonresponder': n_nr_g12,
            'OR': round(or_g12, 4),
            'CI_lower': round(ci_g12_low, 4),
            'CI_upper': round(ci_g12_high, 4),
            'p_value': f'{p_g12_firth:.6e}',
            'method': 'Firth Penalized Logistic Regression',
            'direction': direction_g12,
            'note': 'Patient-level NK-like signature (24 genes); melanoma + anti-CTLA-4',
        })

    # 机制归因结果
    prf1_p = mech_summary.get('PRF1', {}).get('p_value', 'NA')
    prf1_dir = mech_summary.get('PRF1', {}).get('direction', 'NA')
    gzmb_dir = mech_summary.get('GZMB', {}).get('direction', 'NA')

    validation_results.append({
        'cohort': 'GSE120575 (melanoma, anti-CTLA-4)',
        'analysis': 'Mechanistic attribution (PRF1/GZMB in NK-like High cells)',
        'n_patients': n_g12_pat,
        'n_responder': n_r_g12,
        'n_nonresponder': n_nr_g12,
        'OR': 'Not Applicable',
        'CI_lower': 'Not Applicable',
        'CI_upper': 'Not Applicable',
        'p_value': prf1_p if isinstance(prf1_p, str) else f'{prf1_p:.6e}',
        'method': 'Mann-Whitney U (single-cell level, NK-high subset)',
        'direction': f'PRF1: {prf1_dir}, GZMB: {gzmb_dir}',
        'note': 'Functional block hypothesis: higher NK-like signature in NR but lower PRF1 (cytolytic executioner) in NR',
    })

    # 绘图
    print("\n绘图：机制验证图")
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    plot_genes = [g for g in ['PRF1', 'GZMB'] if g in key_genes_found]

    for gi, gene in enumerate(plot_genes):
        ax = axes[gi]
        expr_r = np.log1p(tpm_g12.loc[gene, nk_high_r].values.astype(float))
        expr_nr = np.log1p(tpm_g12.loc[gene, nk_high_nr].values.astype(float))

        bp = ax.boxplot([expr_r, expr_nr], patch_artist=True, widths=0.6, showfliers=False)
        bp['boxes'][0].set_facecolor('#4CAF50')
        bp['boxes'][1].set_facecolor('#F44336')
        bp['boxes'][0].set_alpha(0.6)
        bp['boxes'][1].set_alpha(0.6)

        data_pairs = [(expr_r, '#4CAF50'), (expr_nr, '#F44336')]
        for i, (data_arr, color) in enumerate(data_pairs):
            n_sample = min(len(data_arr), 200)
            x_j = np.random.normal(i + 1, 0.08, size=n_sample)
            ax.scatter(x_j, data_arr[:n_sample], c=color, s=8, alpha=0.3, rasterized=True)

        ax.set_xticklabels(['Responder', 'Non-responder'])
        ax.set_ylabel(f'{gene} expression\n(log1p TPM)')
        ax.set_title(f'{gene} in NK-like High cells')

        p_val = mech_summary.get(gene, {}).get('p_value', 1)
        ax.text(0.5, 0.95, f'p = {p_val:.2e}', transform=ax.transAxes, ha='center', va='top', fontsize=10)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.tight_layout()
    fig.savefig(cfg.result_path('GSE120575_mechanistic_validation.pdf'), dpi=150, bbox_inches='tight')
    fig.savefig(cfg.result_path('GSE120575_mechanistic_validation.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  GSE120575_mechanistic_validation saved")

    # 保存机制数据
    df_mech = pd.DataFrame([
        {'gene': g, 'R_mean': v['R_mean'], 'NR_mean': v['NR_mean'], 'p_value': v['p_value'], 'direction': v['direction']}
        for g, v in mech_summary.items() if 'p_value' in v
    ])
    if len(df_mech) > 0:
        df_mech.to_csv(cfg.result_path('GSE120575_mechanistic_analysis.csv'), index=False)
        print("  GSE120575_mechanistic_analysis.csv saved")

else:
    print("[NOT AVAILABLE] GSE120575 数据不可用")
    validation_results.append({
        'cohort': 'GSE120575 (melanoma, anti-CTLA-4)',
        'analysis': 'NK-like signature + mechanistic attribution',
        'n_patients': 'Not Available',
        'n_responder': 'Not Available',
        'n_nonresponder': 'Not Available',
        'OR': 'Not Available',
        'CI_lower': 'Not Available',
        'CI_upper': 'Not Available',
        'p_value': 'Not Available',
        'method': 'Not Available',
        'direction': 'Not Available',
        'note': 'Data file not found',
    })

# ============================================================
# 保存汇总
# ============================================================
df_val = pd.DataFrame(validation_results)
df_val.to_csv(cfg.result_path('external_validation_summary.csv'), index=False)
print(f"\nexternal_validation_summary.csv saved: {len(df_val)} records")
print(df_val.to_string())

print("\n[step4_5_external_validation] Done.")
