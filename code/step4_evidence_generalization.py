"""
step4_evidence_generalization.py
复现 Fig5 (external validation) + Fig7 (spatial interaction) + FigS6/S7.
"""
import config as cfg
from _common import NKLIKE_SIGNATURE
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, spearmanr
import os
import warnings
warnings.filterwarnings('ignore')

print("[step4_evidence_generalization] Starting...")

os.makedirs(cfg.RESULT_DIR, exist_ok=True)

# ============================================================
# FIGURE 5: External Validation (6 Panels)
# 从原始 TPM/UMI 矩阵重新计算 NK-like 24 基因签名得分
# Panel A: GSE120575 boxplot (melanoma, anti-CTLA-4)
# Panel B: GSE120575 ROC
# Panel C: GSE179994 boxplot (NSCLC, pre/post anti-PD-1)
# Panel D: GSE207422 boxplot (NSCLC, anti-PD-1)
# Panel E: Cross-cohort AUC comparison
# Panel F: Top signature genes heatmap
# ============================================================
fig5, axes5 = plt.subplots(2, 3, figsize=(20, 12))

validation_results = []

# --- GSE120575: 从 TPM 矩阵重新计算签名得分 ---
tpm_path = os.path.join(cfg.ROOT_DIR, 'adata', 'GSE120575_TPM_matrix.csv')
pb_path = os.path.join(cfg.ADATA_DIR, 'GSE120575', 'GSE120575_pseudobulk.csv')

if os.path.exists(tpm_path) and os.path.exists(pb_path):
    print("Loading GSE120575 TPM matrix (may take a moment)...")
    tpm = pd.read_csv(tpm_path, index_col=0)
    tpm.index = tpm.index.astype(str)
    print(f"  TPM matrix: {tpm.shape[0]} genes × {tpm.shape[1]} samples")

    sig_found = [g for g in NKLIKE_SIGNATURE if g in tpm.index]
    print(f"  Signature genes found: {len(sig_found)}/{len(NKLIKE_SIGNATURE)}")

    # 从 GEO metadata 文件获取 cell → patient 映射
    # 文件格式：SOFT 格式，每行以 Sample <name> 开头，后面跟 characteristics
    patient_map_path = os.path.join(cfg.ADATA_DIR, 'GSE120575', 'GSE120575_patient_ID_single_cells.txt.gz')
    sample_cols = tpm.columns.tolist()
    patient_map = {}
    
    if os.path.exists(patient_map_path):
        import gzip
        cell_to_patient = {}
        with gzip.open(patient_map_path, 'rt', encoding='latin-1') as f:
            for line in f:
                line = line.strip()
                if line.startswith('Sample ') and '\t' in line:
                    parts = line.split('\t')
                    if len(parts) >= 5:
                        cell_id = parts[1].strip()
                        pat_id = parts[4].strip()
                        if cell_id and pat_id and cell_id != 'title':
                            cell_to_patient[cell_id] = pat_id
        
        for col in sample_cols:
            if col in cell_to_patient:
                patient_map[col] = cell_to_patient[col]
        
        n_mapped = len(patient_map)
        print(f"  Patient mapping from SOFT file: {n_mapped}/{len(sample_cols)} cells mapped")
    
    if len(patient_map) == 0:
        print("  [WARNING] Could not parse patient IDs, using regex fallback (may be incorrect)")
        import re
        pat_re = re.compile(r'(?:^|_)(P\d+)(?:_|$)')
        for col in sample_cols:
            col_str = str(col)
            m = pat_re.search(col_str)
            if m:
                patient_map[col] = m.group(1)
            else:
                parts = col_str.split('_')
                if len(parts) >= 2 and parts[1].startswith('P'):
                    patient_map[col] = parts[1]
                elif len(parts) >= 1:
                    patient_map[col] = parts[0]

    # 患者水平 pseudobulk（均值）
    sig_expr = tpm.loc[sig_found]
    patient_scores = {}
    for pat in set(patient_map.values()):
        pat_cols = [c for c, p in patient_map.items() if p == pat]
        if len(pat_cols) == 0:
            continue
        pat_expr = sig_expr[pat_cols].mean(axis=1)
        # 签名得分 = 基因表达均值（log1p 后）
        score = np.log1p(pat_expr).mean()
        patient_scores[pat] = score

    df_scores = pd.DataFrame({'patient_id': list(patient_scores.keys()),
                              'nklike_score': list(patient_scores.values())})

    # 合并响应标签
    pb = pd.read_csv(pb_path)
    pb['response_binary'] = (pb['response_binary'] == 'R').astype(int)
    df_merged = df_scores.merge(pb[['patient', 'response_binary', 'response']],
                                 left_on='patient_id', right_on='patient', how='inner')
    df_merged = df_merged.dropna(subset=['nklike_score', 'response_binary'])
    print(f"  Merged patients: {len(df_merged)} (R={int(df_merged['response_binary'].sum())}, NR={int((1-df_merged['response_binary']).sum())})")

    if len(df_merged) >= 6 and df_merged['response_binary'].nunique() == 2:
        from sklearn.metrics import roc_auc_score
        from scipy.stats import mannwhitneyu

        y_g12 = df_merged['response_binary'].values
        score_g12 = df_merged['nklike_score'].values
        auc_g12 = roc_auc_score(y_g12, score_g12)

        r_scores = df_merged.loc[df_merged['response_binary']==1, 'nklike_score'].values
        nr_scores = df_merged.loc[df_merged['response_binary']==0, 'nklike_score'].values
        stat_g12, p_g12 = mannwhitneyu(r_scores, nr_scores, alternative='two-sided')

        validation_results.append({
            'cohort': 'GSE120575 (melanoma, anti-CTLA-4)',
            'n_patients': len(df_merged),
            'n_responder': int(df_merged['response_binary'].sum()),
            'n_nonresponder': int((1-df_merged['response_binary']).sum()),
            'n_signature_genes': len(sig_found),
            'metric': 'AUC',
            'value': round(auc_g12, 4),
            'p_value': round(p_g12, 6),
            'direction': 'Higher in R' if np.mean(r_scores) > np.mean(nr_scores) else 'Higher in NR'
        })
        print(f"  GSE120575: AUC={auc_g12:.4f}, MWU p={p_g12:.4f}")

        ax5A = axes5[0, 0]
        bp_data = [nr_scores, r_scores]
        bp = ax5A.boxplot(bp_data, labels=['Non-responder', 'Responder'], patch_artist=True, widths=0.5)
        for patch, color in zip(bp['boxes'], ['#F44336', '#4CAF50']):
            patch.set_facecolor(color); patch.set_alpha(0.6)
        for i, (data, color) in enumerate(zip([nr_scores, r_scores], ['#F44336', '#4CAF50'])):
            x_j = np.random.normal(i+1, 0.05, size=len(data))
            ax5A.scatter(x_j, data, c=color, s=20, zorder=5, alpha=0.8, edgecolors='black', linewidth=0.3)
        ax5A.set_ylabel('NK-like signature score (24 genes)')
        ax5A.set_title(f'Panel A: GSE120575 Melanoma (anti-CTLA-4)\nAUC={auc_g12:.3f}, MWU p={p_g12:.4f}')

        # ── Panel B: GSE120575 ROC curve ──
        from sklearn.metrics import roc_curve
        ax5B = axes5[0, 1]
        fpr_v, tpr_v, _ = roc_curve(y_g12, score_g12)
        ax5B.plot(fpr_v, tpr_v, '#E74C3C', lw=2.5, label=f'NK-like signature (AUC={auc_g12:.3f})')
        ax5B.plot([0, 1], [0, 1], 'k--', lw=0.8, alpha=0.5, label='Random (AUC=0.5)')
        ax5B.set_xlabel('False Positive Rate')
        ax5B.set_ylabel('True Positive Rate')
        ax5B.set_title(f'Panel B: GSE120575 ROC Curve\n(n={len(df_merged)})')
        ax5B.legend(loc='lower right', fontsize=9)
        ax5B.spines['top'].set_visible(False); ax5B.spines['right'].set_visible(False)

        # ── Panel F: Top signature genes expression heatmap (GSE120575 R vs NR) ──
        ax5F = axes5[1, 2]
        r_pat_ids = df_merged.loc[df_merged['response_binary']==1, 'patient_id'].values
        nr_pat_ids = df_merged.loc[df_merged['response_binary']==0, 'patient_id'].values
        r_pat_cols = [c for c, p in patient_map.items() if p in r_pat_ids]
        nr_pat_cols = [c for c, p in patient_map.items() if p in nr_pat_ids]

        gene_diff = {}
        for g in sig_found:
            r_mean = sig_expr.loc[g, r_pat_cols].mean() if len(r_pat_cols) > 0 else 0
            nr_mean = sig_expr.loc[g, nr_pat_cols].mean() if len(nr_pat_cols) > 0 else 0
            gene_diff[g] = r_mean - nr_mean

        top_genes = sorted(gene_diff, key=lambda x: -abs(gene_diff[x]))[:12]
        top_genes_sorted = sorted(top_genes, key=lambda x: gene_diff[x])

        heatmap_data = np.zeros((len(top_genes_sorted), 2))
        for i, g in enumerate(top_genes_sorted):
            heatmap_data[i, 0] = sig_expr.loc[g, nr_pat_cols].mean() if len(nr_pat_cols) > 0 else 0
            heatmap_data[i, 1] = sig_expr.loc[g, r_pat_cols].mean() if len(r_pat_cols) > 0 else 0

        row_mean = heatmap_data.mean(axis=1, keepdims=True)
        row_std = heatmap_data.std(axis=1, keepdims=True)
        row_std[row_std == 0] = 1
        heatmap_z = (heatmap_data - row_mean) / row_std

        im = ax5F.imshow(heatmap_z, cmap='RdBu_r', aspect='auto', vmin=-1.5, vmax=1.5)
        ax5F.set_xticks([0, 1]); ax5F.set_xticklabels(['Non-responder', 'Responder'])
        ax5F.set_yticks(range(len(top_genes_sorted))); ax5F.set_yticklabels(top_genes_sorted, fontsize=8)
        ax5F.set_title('Panel F: Top 12 Signature Genes\n(R vs NR expression, z-scored)')
        plt.colorbar(im, ax=ax5F, label='z-score')
    else:
        axes5[0, 0].text(0.5, 0.5, 'GSE120575: insufficient data',
                          ha='center', va='center', transform=axes5[0, 0].transAxes, fontsize=12, color='grey')
        axes5[0, 0].set_title('Panel A: GSE120575 (insufficient)')
        for ax, panel in zip([axes5[0,1], axes5[1,2]], ['B', 'F']):
            ax.text(0.5, 0.5, f'Panel {panel}: Data unavailable',
                    ha='center', va='center', transform=ax.transAxes, fontsize=11, color='grey')
else:
    axes5[0, 0].text(0.5, 0.5, 'GSE120575 TPM data not found',
                      ha='center', va='center', transform=axes5[0, 0].transAxes, fontsize=12, color='grey')
    axes5[0, 0].set_title('Panel A: (data missing)')
    for ax, panel in zip([axes5[0,1], axes5[1,2]], ['B', 'F']):
        ax.text(0.5, 0.5, f'Panel {panel}: Data unavailable',
                ha='center', va='center', transform=ax.transAxes, fontsize=11, color='grey')

# --- GSE179994: NSCLC tumor scRNA-seq, anti-PD-1 + chemo, 36 patients ---
# Liu et al. 2022, Nature Cancer. Complete T cell rawCounts (150,849 cells)
# Response labels from Supplementary Table 1 (MOESM3)
print("\nProcessing GSE179994 (Liu 2022 Nature Cancer)...")
g17_meta_path = os.path.join(cfg.ADATA_DIR, 'GSE179994_Tcell.metadata.tsv.gz')
g17_scores_path = os.path.join(cfg.ADATA_DIR, 'GSE179994', 'GSE179994_all_signature_scores.csv')

g17_auc = np.nan
g17_p = np.nan

# Response labels from paper Supplementary Table 1 (sample-level, accounts for mixed responses)
G17_RESPONSE_MAP = {
    'P1.post.1': 'Responder',      # P001 LN metastasis (responsive)
    'P1.post.2': 'Responder',      # P001 LN metastasis (responsive, 2nd biopsy)
    'P1.post.3': 'Non-responder',  # P001 lung primary (refractory)
    'P10.post.1': 'Responder',     # P010 LN metastasis
    'P13.post.1': 'Responder',     # P013 liver metastasis (responsive)
    'P13.post.3': 'Non-responder', # P013 new LN metastasis (nonresponsive)
    'P19.post.1': 'Responder',     # P019 LN metastasis
    'P29.post.1': 'Responder',     # P029 left lung
    'P30.post.1': 'Responder',     # P030 right lung
    'P33.post.1': 'Responder',     # P033 right lung
    'P35.post.1': 'Responder',     # P035 right lung
    'P36.post.1': 'Non-responder', # P036 right lung
    'P37.post.1': 'Non-responder', # P037 left lung
    'P38.post.1': 'Non-responder', # P038 left lung
}

if os.path.exists(g17_meta_path) and os.path.exists(g17_scores_path):
    df_g17_meta = pd.read_csv(g17_meta_path, sep='\t')
    df_g17_scores = pd.read_csv(g17_scores_path)
    df_g17 = df_g17_meta.merge(df_g17_scores, left_on='cellid', right_on='cell', how='inner')
    print(f"  Total cells: {len(df_g17)}, patients: {df_g17['patient'].nunique()}")

    # Timepoint: pre (ut) vs post (tr)
    df_g17['timepoint'] = df_g17['sample'].apply(
        lambda x: 'pre' if '.pre' in str(x) or '.ut' in str(x)
        else ('post' if '.post' in str(x) or '.tr' in str(x) else 'unknown'))

    # Add response labels to post-treatment samples
    df_g17['response'] = df_g17['sample'].map(G17_RESPONSE_MAP)

    # CD8 T cells only
    df_g17_cd8 = df_g17[df_g17['celltype'] == 'CD8'].copy()
    df_g17_post = df_g17_cd8[(df_g17_cd8['timepoint'] == 'post') & (df_g17_cd8['response'].notna())].copy()
    print(f"  CD8 post-treatment cells with response: {len(df_g17_post)}")
    print(f"  Responder: {(df_g17_post['response']=='Responder').sum()}, Non-responder: {(df_g17_post['response']=='Non-responder').sum()}")

    # Sample-level scores (pseudobulk)
    sample_scores_g17 = df_g17_post.groupby(['sample', 'response']).agg(
        mean_score=('score', 'mean'),
        n_cells=('score', 'count')
    ).reset_index()

    r_scores_g17 = sample_scores_g17[sample_scores_g17['response'] == 'Responder']['mean_score'].values
    nr_scores_g17 = sample_scores_g17[sample_scores_g17['response'] == 'Non-responder']['mean_score'].values

    print(f"  Responder samples: n={len(r_scores_g17)}, mean={np.mean(r_scores_g17):.4f}")
    print(f"  Non-responder samples: n={len(nr_scores_g17)}, mean={np.mean(nr_scores_g17):.4f}")

    if len(r_scores_g17) >= 3 and len(nr_scores_g17) >= 3:
        stat_g17, p_g17 = mannwhitneyu(r_scores_g17, nr_scores_g17, alternative='two-sided')
        y_g17 = np.array([1] * len(r_scores_g17) + [0] * len(nr_scores_g17))
        score_g17 = np.concatenate([r_scores_g17, nr_scores_g17])
        auc_g17 = roc_auc_score(y_g17, score_g17)
        g17_auc = auc_g17
        g17_p = p_g17
        direction_g17 = 'Higher in R' if np.mean(r_scores_g17) > np.mean(nr_scores_g17) else 'Higher in NR'
        print(f"  AUC={auc_g17:.4f}, MWU p={p_g17:.4f}, direction={direction_g17}")

        validation_results.append({
            'cohort': 'GSE179994 (NSCLC, tumor, anti-PD-1+chemo)',
            'n_patients': 11,  # 8 R + 3 NR (post-treatment)
            'n_responder': len(r_scores_g17),
            'n_nonresponder': len(nr_scores_g17),
            'n_signature_genes': 24,
            'metric': 'AUC',
            'value': round(auc_g17, 4),
            'p_value': round(p_g17, 6),
            'direction': direction_g17,
            'reason': 'Sample-level (pseudobulk) NK-like signature score; response from Supp Table 1'
        })

    # Also compute pre vs post dynamics for paired patients
    cd8_paired = df_g17_cd8[df_g17_cd8['patient'].isin(
        set(df_g17_cd8[df_g17_cd8['timepoint']=='pre']['patient']) &
        set(df_g17_cd8[df_g17_cd8['timepoint']=='post']['patient'])
    )].copy()
    pat_tp_g17 = cd8_paired.groupby(['patient', 'timepoint'])['score'].mean().reset_index()

    # Panel C: GSE179994 R vs NR boxplot + pre/post dynamics
    ax5C = axes5[0, 2]

    if len(r_scores_g17) >= 3 and len(nr_scores_g17) >= 3:
        bp_data_g17 = [nr_scores_g17, r_scores_g17]
        bp_g17 = ax5C.boxplot(bp_data_g17, labels=['Non-responder', 'Responder'],
                              patch_artist=True, widths=0.5, positions=[1, 1.8])
        for patch, color in zip(bp_g17['boxes'], ['#F44336', '#4CAF50']):
            patch.set_facecolor(color); patch.set_alpha(0.6)
        for i, (data, color) in enumerate(zip([nr_scores_g17, r_scores_g17], ['#F44336', '#4CAF50'])):
            x_j = np.random.normal([1, 1.8][i], 0.05, size=len(data))
            ax5C.scatter(x_j, data, c=color, s=30, zorder=5, alpha=0.8, edgecolors='black', linewidth=0.3)

        # Add pre/post paired lines on the right side
        ax5C2 = ax5C.twinx()
        x_offset = 3.0
        pat_ids = pat_tp_g17['patient'].unique()
        for idx, pat in enumerate(sorted(pat_ids)):
            pre_val = pat_tp_g17[(pat_tp_g17['patient']==pat) & (pat_tp_g17['timepoint']=='pre')]['score'].values
            post_val = pat_tp_g17[(pat_tp_g17['patient']==pat) & (pat_tp_g17['timepoint']=='post')]['score'].values
            if len(pre_val) > 0 and len(post_val) > 0:
                x_pre = x_offset + idx * 0.3
                x_post = x_pre + 0.15
                ax5C2.plot([x_pre, x_post], [pre_val[0], post_val[0]], 'k-', alpha=0.4, lw=0.8)
                ax5C2.scatter([x_pre, x_post], [pre_val[0], post_val[0]], c=['#2196F3', '#FF9800'], s=15, zorder=5)
        ax5C2.set_ylabel('Pre/Post paired score', fontsize=8, color='grey')
        ax5C2.tick_params(axis='y', labelcolor='grey', labelsize=7)
        ax5C2.set_xlim(2.5, x_offset + len(pat_ids) * 0.3 + 0.5)

        ax5C.set_ylabel('NK-like signature score (24 genes)')
        ax5C.set_xlim(0.5, 2.5)
        ax5C.set_title(f'Panel C: GSE179994 NSCLC (anti-PD-1+chemo)\n'
                       f'R vs NR: AUC={auc_g17:.3f}, p={p_g17:.3f}\n'
                       f'Right: pre→post dynamics (n={len(pat_ids)} paired)')
        print(f"  GSE179994 Panel C: R vs NR boxplot + dynamics drawn")
    else:
        ax5C.text(0.5, 0.5, 'GSE179994: insufficient response groups',
                  ha='center', va='center', transform=ax5C.transAxes, fontsize=11, color='grey')
        ax5C.set_title('Panel C: GSE179994 (insufficient)')
else:
    ax5C = axes5[0, 2]
    ax5C.text(0.5, 0.5, 'GSE179994 data not found',
              ha='center', va='center', transform=ax5C.transAxes, fontsize=12, color='grey')
    ax5C.set_title('Panel C: GSE179994 (unavailable)')

# --- GSE207422: NSCLC scRNA-seq with RECIST response ---
print("\nProcessing GSE207422 (Hu 2023 Genome Medicine)...")
g20_meta_path = os.path.join(cfg.ADATA_DIR, 'GSE207422_NSCLC_scRNAseq_metadata.xlsx')
g20_umi_path = os.path.join(cfg.ADATA_DIR, 'GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz')

g20_auc = np.nan
g20_p = np.nan

if os.path.exists(g20_meta_path) and os.path.exists(g20_umi_path):
    df_g20_meta = pd.read_excel(g20_meta_path)
    print(f"  Metadata: {df_g20_meta.shape[0]} samples")
    print(f"  RECIST values: {df_g20_meta['RECIST'].value_counts().to_dict()}")
    
    print("  Reading UMI matrix (signature genes only)...")
    import gzip
    sig_set = set(NKLIKE_SIGNATURE)
    sig_expr_g20 = {}
    cell_cols = None
    
    with gzip.open(g20_umi_path, 'rt') as f:
        header = f.readline().strip()
        cell_cols = header.split('\t')
        cell_cols = cell_cols[1:]
        
        for line in f:
            parts = line.strip().split('\t')
            gene = parts[0]
            if gene in sig_set:
                expr_vals = np.array([int(float(v)) for v in parts[1:]])
                sig_expr_g20[gene] = expr_vals
    
    print(f"  Signature genes found: {len(sig_expr_g20)}/{len(NKLIKE_SIGNATURE)}")
    print(f"  Number of cells: {len(cell_cols)}")
    
    if len(sig_expr_g20) > 0 and cell_cols:
        sig_found_g20 = list(sig_expr_g20.keys())
        
        cell_to_patient = {}
        for col in cell_cols:
            parts = str(col).split('_')
            if len(parts) >= 3 and parts[0] == 'BD' and parts[1].startswith('immune'):
                pat_num = parts[1].replace('immune', '')
                patient_id = f'P{pat_num}'
                cell_to_patient[col] = patient_id
        
        print(f"  Cell-to-patient mapping: {len(cell_to_patient)}/{len(cell_cols)}")
        
        patient_scores_g20 = {}
        for pat in set(cell_to_patient.values()):
            pat_idx = [i for i, c in enumerate(cell_cols) if c in cell_to_patient and cell_to_patient[c] == pat]
            if len(pat_idx) == 0:
                continue
            
            pat_expr_sum = 0
            n_genes = 0
            for gene in sig_found_g20:
                gene_vals = sig_expr_g20[gene][pat_idx]
                pat_expr_sum += np.log1p(gene_vals).mean()
                n_genes += 1
            
            if n_genes > 0:
                patient_scores_g20[pat] = pat_expr_sum / n_genes
        
        df_scores_g20 = pd.DataFrame({'patient_id': list(patient_scores_g20.keys()),
                                       'nklike_score': list(patient_scores_g20.values())})
        
        df_g20_meta['patient_id'] = df_g20_meta['Patient']
        df_merged_g20 = df_scores_g20.merge(df_g20_meta[['patient_id', 'RECIST']],
                                             on='patient_id', how='inner')
        
        df_merged_g20['response_binary'] = df_merged_g20['RECIST'].apply(
            lambda x: 1 if str(x).strip() in ['PR'] else 0)
        
        r_mask_g20 = df_merged_g20['response_binary'] == 1
        nr_mask_g20 = df_merged_g20['response_binary'] == 0
        
        print(f"  Merged: {len(df_merged_g20)} patients (PR={r_mask_g20.sum()}, non-PR={nr_mask_g20.sum()})")
        
        if len(df_merged_g20) >= 6 and df_merged_g20['response_binary'].nunique() == 2:
            r_scores_g20 = df_merged_g20.loc[r_mask_g20, 'nklike_score'].values
            nr_scores_g20 = df_merged_g20.loc[nr_mask_g20, 'nklike_score'].values
            
            auc_g20 = roc_auc_score(df_merged_g20['response_binary'].values, df_merged_g20['nklike_score'].values)
            stat_g20, p_g20 = mannwhitneyu(r_scores_g20, nr_scores_g20, alternative='two-sided')
            
            validation_results.append({
                'cohort': 'GSE207422 (NSCLC, tumor, anti-PD-1)',
                'n_patients': len(df_merged_g20),
                'n_responder': int(r_mask_g20.sum()),
                'n_nonresponder': int(nr_mask_g20.sum()),
                'n_signature_genes': len(sig_found_g20),
                'metric': 'AUC',
                'value': round(auc_g20, 4),
                'p_value': round(p_g20, 6),
                'direction': 'Higher in R' if np.mean(r_scores_g20) > np.mean(nr_scores_g20) else 'Higher in NR'
            })
            print(f"  GSE207422: AUC={auc_g20:.4f}, MWU p={p_g20:.4f}")
            
            ax5D = axes5[1, 0]
            bp_data_g20 = [nr_scores_g20, r_scores_g20]
            bp_g20 = ax5D.boxplot(bp_data_g20, labels=['Non-responder (SD)', 'Responder (PR)'], 
                                  patch_artist=True, widths=0.5)
            for patch, color in zip(bp_g20['boxes'], ['#F44336', '#4CAF50']):
                patch.set_facecolor(color); patch.set_alpha(0.6)
            for i, (data, color) in enumerate(zip([nr_scores_g20, r_scores_g20], ['#F44336', '#4CAF50'])):
                x_j = np.random.normal(i+1, 0.05, size=len(data))
                ax5D.scatter(x_j, data, c=color, s=20, zorder=5, alpha=0.8, edgecolors='black', linewidth=0.3)
            ax5D.set_ylabel('NK-like signature score (24 genes)')
            ax5D.set_title(f'Panel D: GSE207422 NSCLC (anti-PD-1)\nAUC={auc_g20:.3f}, MWU p={p_g20:.4f}')
            
            g20_auc = auc_g20
            g20_p = p_g20
        else:
            ax5D = axes5[1, 0]
            ax5D.text(0.5, 0.5, f'GSE207422: insufficient groups\n({r_mask_g20.sum()} PR, {nr_mask_g20.sum()} non-PR)',
                      ha='center', va='center', transform=ax5D.transAxes, fontsize=11, color='grey')
            ax5D.set_title('Panel D: GSE207422 (insufficient)')
    else:
        ax5D = axes5[1, 0]
        ax5D.text(0.5, 0.5, 'GSE207422: no signature genes found',
                  ha='center', va='center', transform=ax5D.transAxes, fontsize=12, color='grey')
        ax5D.set_title('Panel D: GSE207422 (no signature genes)')
else:
    ax5D = axes5[1, 0]
    ax5D.text(0.5, 0.5, 'GSE207422 data not found',
              ha='center', va='center', transform=ax5D.transAxes, fontsize=12, color='grey')
    ax5D.set_title('Panel D: GSE207422 (unavailable)')

# ── Panel E: Cross-cohort AUC comparison bar chart ──
ax5E = axes5[1, 1]
main_auc = np.nan
main_auc_path = os.path.join(cfg.RESULT_DIR, 'irs_model_stats.csv')
if os.path.exists(main_auc_path):
    try:
        df_irs_stats = pd.read_csv(main_auc_path)
        if 'auc_irs' in df_irs_stats.columns:
            main_auc = float(df_irs_stats['auc_irs'].iloc[0])
    except Exception:
        pass

cohort_names = []
auc_values = []
colors_bar = []

if not np.isnan(main_auc):
    cohort_names.append('GSE243013\n(NSCLC, anti-PD1)')
    auc_values.append(main_auc)
    colors_bar.append('#2E86C1')
if not np.isnan(g17_auc):
    cohort_names.append('GSE179994\n(NSCLC, anti-PD1+chemo)')
    auc_values.append(g17_auc)
    colors_bar.append('#8E44AD')
if 'auc_g12' in locals() and not np.isnan(auc_g12):
    cohort_names.append('GSE120575\n(Melanoma, anti-CTLA-4)')
    auc_values.append(auc_g12)
    colors_bar.append('#E74C3C')
if not np.isnan(g20_auc):
    cohort_names.append('GSE207422\n(NSCLC, anti-PD1)')
    auc_values.append(g20_auc)
    colors_bar.append('#27AE60')

if cohort_names:
    bars = ax5E.bar(range(len(cohort_names)), auc_values, color=colors_bar, alpha=0.7, edgecolor='black', linewidth=0.5)
    ax5E.axhline(0.5, color='gray', ls='--', lw=1, alpha=0.6, label='Random (0.5)')
    ax5E.set_xticks(range(len(cohort_names)))
    ax5E.set_xticklabels(cohort_names, fontsize=9)
    ax5E.set_ylabel('AUC')
    ax5E.set_ylim(0, 1.0)
    ax5E.set_title('Panel E: Cross-cohort AUC Comparison\n(NK-like signature predicts response)')
    for i, (bar, v) in enumerate(zip(bars, auc_values)):
        ax5E.text(bar.get_x() + bar.get_width()/2, v + 0.02, f'{v:.3f}',
                  ha='center', va='bottom', fontsize=10, fontweight='bold')
    ax5E.legend(loc='upper right', fontsize=8)
    ax5E.spines['top'].set_visible(False); ax5E.spines['right'].set_visible(False)
else:
    ax5E.text(0.5, 0.5, 'No AUC data available',
              ha='center', va='center', transform=ax5E.transAxes, fontsize=12, color='grey')
    ax5E.set_title('Panel E: AUC Comparison (unavailable)')

# 外部队列数据可用性说明（保留记录但不再占用图表面板）
for cohort_name in ['GSE126044 (NSCLC, anti-PD-1)', 'GSE135222 (NSCLC, anti-PD-1/L1)', 'GSE91061 (melanoma, anti-CTLA-4/PD-1)']:
    validation_results.append({
        'cohort': cohort_name.split('(')[0].strip(),
        'n_patients': 'Not Available',
        'n_responder': 'Not Available',
        'n_nonresponder': 'Not Available',
        'n_signature_genes': 'Not Available',
        'metric': 'AUC',
        'value': 'Not Available',
        'p_value': 'Not Available',
        'direction': 'Not Available',
        'reason': 'No matching scRNA-seq data with response labels in local adata directory',
    })

if validation_results:
    df_val = pd.DataFrame(validation_results)
    df_val.to_csv(cfg.result_path('fig5_validation_summary.csv'), index=False)
    print(f"fig5_validation_summary.csv saved: {len(df_val)} cohorts")

plt.tight_layout()
fig5.savefig(cfg.result_path('Fig5_external_validation.pdf'), dpi=150, bbox_inches='tight')
fig5.savefig(cfg.result_path('Fig5_external_validation.png'), dpi=150, bbox_inches='tight')
plt.close(fig5)
print("Fig5_external_validation saved")

# ============================================================
# FIGURE 7: Spatial Interaction (4 Panels)
# Panel A: Stemness-Effector Decoupling vs SPP1 (KEY RESULT — promoted from FigS)
# Panel B: SPP1 vs NK Cyto (supporting)
# Panel C: SPP1 vs Stemness (supporting)
# Panel D: NK Cyto by response (clinical relevance)
# ============================================================
spatial_path = cfg.result_path('spatial_roi_scores.csv')
spatial_archive = os.path.join(cfg.ROOT_DIR, 'archive_v_original', 'result', 'spatial_roi_scores.csv')

fig7, axes7 = plt.subplots(2, 2, figsize=(14, 12))

spatial_available = False
if os.path.exists(spatial_path):
    df_sp = pd.read_csv(spatial_path)
    spatial_available = True
elif os.path.exists(spatial_archive):
    df_sp = pd.read_csv(spatial_archive)
    spatial_available = True

if spatial_available:
    print(f"Loaded spatial data: {df_sp.shape}")
    print(f"Spatial columns: {list(df_sp.columns)}")

    # 明确指定数值列（避免模糊匹配选错列）
    # spp1_tam 是 SPP1/TAM 比值（数值），spp1_rec_interact 是交互得分（数值）
    # 排除 'spp1_grp', 'rec_grp', 'lr_group' 等分类列
    spp1_numeric_cols = ['spp1_tam', 'spp1_rec_interact', 'SPP1']
    cyto_numeric_cols = ['nk_cyto_tcell', 'nk_cyto_score']
    stem_numeric_cols = ['stem_tcell', 'stem_score']

    spp1_cols = [c for c in spp1_numeric_cols if c in df_sp.columns]
    cyto_cols = [c for c in cyto_numeric_cols if c in df_sp.columns]
    stem_cols = [c for c in stem_numeric_cols if c in df_sp.columns]

    print(f"  SPP1 numeric cols: {spp1_cols}")
    print(f"  Cyto numeric cols: {cyto_cols}")
    print(f"  Stem numeric cols: {stem_cols}")

    # ── Panel A: Stemness-Effector Decoupling vs SPP1 (promoted from FigS_spatial_decoupling) ──
    # Decoupling_Index = z(stem_tcell) - z(nk_cyto_tcell); >0 = stemness-high/cyto-low (arrest)
    ax7A = axes7[0, 0]
    if 'spp1_rec_interact' in df_sp.columns and 'decoupling_index' in df_sp.columns:
        from scipy import stats as ss
        from matplotlib.lines import Line2D

        x_dec = df_sp['spp1_rec_interact'].values
        y_dec = df_sp['decoupling_index'].values
        valid_dec = ~(np.isnan(x_dec) | np.isnan(y_dec))

        if valid_dec.sum() > 3:
            # Spearman correlation
            r_dec, p_dec = spearmanr(x_dec[valid_dec], y_dec[valid_dec])

            # Color by SPP1 group (values: "High SPP1", "Low SPP1", possibly others)
            spp1_grp_col = 'spp1_grp' if 'spp1_grp' in df_sp.columns else None
            colors_dec = []
            if spp1_grp_col:
                for g in df_sp.loc[valid_dec, spp1_grp_col]:
                    if str(g).strip() in ['High SPP1', 'High']:
                        colors_dec.append('#E74C3C')
                    elif str(g).strip() in ['Low SPP1', 'Low']:
                        colors_dec.append('#2E86C1')
                    else:
                        colors_dec.append('#95A5A6')
            else:
                colors_dec = ['#95A5A6'] * valid_dec.sum()

            ax7A.scatter(x_dec[valid_dec], y_dec[valid_dec], c=colors_dec,
                        s=60, alpha=0.7, edgecolors='black', linewidth=0.5)

            # Linear regression line
            slope_d, intercept_d, r_val_d, p_val_d, _ = ss.linregress(x_dec[valid_dec], y_dec[valid_dec])
            x_line = np.linspace(x_dec[valid_dec].min(), x_dec[valid_dec].max(), 100)
            ax7A.plot(x_line, slope_d * x_line + intercept_d, 'k--', lw=1.5, alpha=0.8)

            # Zero reference line
            ax7A.axhline(0, color='gray', ls=':', alpha=0.5)

            # Stats annotation
            ax7A.text(0.05, 0.95,
                     f'Spearman r={r_dec:.3f}, p={p_dec:.4f}\nLinear R²={r_val_d**2:.3f}',
                     transform=ax7A.transAxes, fontsize=10, verticalalignment='top',
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

            # Legend
            legend_elements = [
                Line2D([0], [0], marker='o', color='w', markerfacecolor='#E74C3C', markersize=10, label='High SPP1 ROI'),
                Line2D([0], [0], marker='o', color='w', markerfacecolor='#2E86C1', markersize=10, label='Low SPP1 ROI'),
                Line2D([0], [0], marker='o', color='w', markerfacecolor='#95A5A6', markersize=10, label='Mid SPP1 ROI'),
            ]
            ax7A.legend(handles=legend_elements, loc='lower right', fontsize=8)
            ax7A.set_xlabel('SPP1 Interaction Score (spp1_rec_interact)')
            ax7A.set_ylabel('Decoupling Index (z[stem] - z[cytotoxicity])')
            ax7A.set_title(f'Panel A: Stemness-Effector Decoupling vs SPP1\n'
                          f'(DI>0: stemness-high/cyto-low = arrest state)')
            print(f"  [Panel A] Decoupling: Spearman r={r_dec:.4f}, p={p_dec:.4f}, n={valid_dec.sum()} ROI")
        else:
            ax7A.text(0.5, 0.5, 'Insufficient valid decoupling data', ha='center', va='center',
                     transform=ax7A.transAxes, color='grey')
            ax7A.set_title('Panel A: Decoupling vs SPP1')
    else:
        # Fallback: compute decoupling_index on the fly if not present
        if 'stem_tcell' in df_sp.columns and 'nk_cyto_tcell' in df_sp.columns and 'spp1_rec_interact' in df_sp.columns:
            stem_z = (df_sp['stem_tcell'] - df_sp['stem_tcell'].mean()) / df_sp['stem_tcell'].std()
            cyto_z = (df_sp['nk_cyto_tcell'] - df_sp['nk_cyto_tcell'].mean()) / df_sp['nk_cyto_tcell'].std()
            df_sp['decoupling_index'] = stem_z - cyto_z
            r_dec, p_dec = spearmanr(df_sp['spp1_rec_interact'], df_sp['decoupling_index'])
            ax7A.scatter(df_sp['spp1_rec_interact'], df_sp['decoupling_index'], alpha=0.6, s=40, c='steelblue')
            ax7A.axhline(0, color='gray', ls=':', alpha=0.5)
            ax7A.text(0.05, 0.95, f'Spearman r={r_dec:.3f}, p={p_dec:.4f}',
                     transform=ax7A.transAxes, fontsize=10, verticalalignment='top',
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            ax7A.set_xlabel('SPP1 Interaction Score')
            ax7A.set_ylabel('Decoupling Index (z[stem] - z[cyto])')
            ax7A.set_title(f'Panel A: Decoupling vs SPP1 (computed on-the-fly)\nr={r_dec:.3f}, p={p_dec:.4f}')
            print(f"  [Panel A] Decoupling (on-the-fly): r={r_dec:.4f}, p={p_dec:.4f}")
        else:
            ax7A.text(0.5, 0.5, 'decoupling_index / spp1_rec_interact not found',
                     ha='center', va='center', transform=ax7A.transAxes, color='grey')
            ax7A.set_title('Panel A: Decoupling vs SPP1 (unavailable)')

    # ── Panel B: SPP1 signal vs NK cyto score scatter (was Panel A) ──
    ax7B = axes7[0, 1]
    if spp1_cols and cyto_cols:
        spp1_vals = df_sp[spp1_cols[0]].values
        cyto_vals = df_sp[cyto_cols[0]].values
        valid = ~(np.isnan(spp1_vals) | np.isnan(cyto_vals))
        if valid.sum() > 2:
            sr7a, sp7a = spearmanr(spp1_vals[valid], cyto_vals[valid])
            ax7B.scatter(spp1_vals[valid], cyto_vals[valid], alpha=0.7, s=60, c='steelblue',
                        edgecolors='white', linewidth=0.5)
            # 回归线
            z7b = np.polyfit(spp1_vals[valid], cyto_vals[valid], 1)
            xs7b = np.linspace(spp1_vals[valid].min(), spp1_vals[valid].max(), 50)
            ax7B.plot(xs7b, np.polyval(z7b, xs7b), 'k--', lw=1.5, alpha=0.7)
            ax7B.set_xlabel(f'{spp1_cols[0]} (SPP1/TAM ratio, log2)')
            ax7B.set_ylabel(f'{cyto_cols[0]} (NK-cyto/Tcell, log2)')
            ax7B.set_title(f'Panel B: SPP1 vs NK Cytotoxicity\nSpearman r={sr7a:.4f}, p={sp7a:.4f}')
            ax7B.grid(alpha=0.3, ls=':')
            print(f"  [Panel B] SPP1 vs Cyto: r={sr7a:.4f}, p={sp7a:.4f}, n={valid.sum()}")
        else:
            ax7B.text(0.5, 0.5, f'Insufficient valid data ({valid.sum()})', ha='center', va='center',
                     transform=ax7B.transAxes, color='grey')
            ax7B.set_title('Panel B: SPP1 vs NK Cyto')
    else:
        ax7B.text(0.5, 0.5, 'SPP1/cyto columns not found', ha='center', va='center',
                 transform=ax7B.transAxes, color='grey')
        ax7B.set_title('Panel B: SPP1 vs NK Cyto')

    # ── Panel C: SPP1 vs Stemness scatter ──
    ax7C = axes7[1, 0]
    if spp1_cols and stem_cols:
        spp1_v = df_sp[spp1_cols[0]].values
        stem_v = df_sp[stem_cols[0]].values
        valid_c = ~(np.isnan(spp1_v) | np.isnan(stem_v))
        if valid_c.sum() > 2:
            sr7c, sp7c = spearmanr(spp1_v[valid_c], stem_v[valid_c])
            ax7C.scatter(spp1_v[valid_c], stem_v[valid_c], alpha=0.7, s=60, c='darkorange',
                        edgecolors='white', linewidth=0.5)
            # 回归线
            z7c = np.polyfit(spp1_v[valid_c], stem_v[valid_c], 1)
            xs7c = np.linspace(spp1_v[valid_c].min(), spp1_v[valid_c].max(), 50)
            ax7C.plot(xs7c, np.polyval(z7c, xs7c), 'k--', lw=1.5, alpha=0.7)
            ax7C.set_xlabel(f'{spp1_cols[0]} (SPP1/TAM ratio, log2)')
            ax7C.set_ylabel(f'{stem_cols[0]} (Stem/Tcell, log2)')
            ax7C.set_title(f'Panel C: SPP1 vs Stemness\nSpearman r={sr7c:.4f}, p={sp7c:.4f}')
            ax7C.grid(alpha=0.3, ls=':')
            print(f"  [Panel C] SPP1 vs Stem: r={sr7c:.4f}, p={sp7c:.4f}, n={valid_c.sum()}")
        else:
            ax7C.text(0.5, 0.5, f'Insufficient valid data ({valid_c.sum()})', ha='center', va='center',
                     transform=ax7C.transAxes, color='grey')
            ax7C.set_title('Panel C: SPP1 vs Stemness')
    else:
        ax7C.text(0.5, 0.5, 'SPP1/stem columns not found', ha='center', va='center',
                 transform=ax7C.transAxes, color='grey')
        ax7C.set_title('Panel C: SPP1 vs Stemness')

    # ── Panel D: NK Cyto by response (clinical relevance, was Panel B) ──
    ax7D = axes7[1, 1]
    resp_col = [c for c in df_sp.columns if 'response' in str(c).lower() or 'group' in str(c).lower() or 'R/NR' in str(c)]
    # Avoid using spp1_grp as response column
    resp_col = [c for c in resp_col if c != 'spp1_grp' and c != 'rec_grp' and c != 'lr_group']
    if resp_col and cyto_cols:
        resp_vals = df_sp[resp_col[0]].values
        resp_mask = np.isin(resp_vals, ['Responder', 'R', 'pCR', 'Response', 'responder'])
        nonresp_mask = np.isin(resp_vals, ['Non-responder', 'NR', 'non-MPR', 'Non-response', 'nonresponder', 'Progressive'])
        bp7b_data = []
        bp7b_labels = []
        if resp_mask.sum() > 0:
            bp7b_data.append(df_sp.loc[resp_mask, cyto_cols[0]].dropna().values)
            bp7b_labels.append('Responder')
        if nonresp_mask.sum() > 0:
            bp7b_data.append(df_sp.loc[nonresp_mask, cyto_cols[0]].dropna().values)
            bp7b_labels.append('Non-responder')
        if bp7b_data and len(bp7b_data) == 2:
            bp7b = ax7D.boxplot(bp7b_data, patch_artist=True)
            ax7D.set_xticklabels(bp7b_labels)
            for b, c in zip(bp7b['boxes'], ['#4CAF50', '#F44336']):
                b.set_facecolor(c)
            ax7D.set_title(f'Panel D: {cyto_cols[0]} by response')
            ax7D.set_ylabel(cyto_cols[0])
        else:
            ax7D.text(0.5, 0.5, f'Only {len(bp7b_data)} groups found', ha='center', va='center',
                     transform=ax7D.transAxes, color='grey')
            ax7D.set_title('Panel D: by response')
    else:
        ax7D.text(0.5, 0.5, 'Response/cyto column not found', ha='center', va='center',
                 transform=ax7D.transAxes, color='grey')
        ax7D.set_title('Panel D: NK Cyto by response')
else:
    for i, ax in enumerate(axes7.flatten()):
        ax.text(0.5, 0.5, f'Panel {chr(65+i)}: Spatial data not available',
               ha='center', va='center', transform=ax.transAxes, fontsize=10, color='grey')
        ax.set_title(f'Panel {chr(65+i)}: (unavailable)')

plt.tight_layout()
fig7.savefig(cfg.result_path('Fig7_spatial_interaction.pdf'), dpi=150, bbox_inches='tight')
fig7.savefig(cfg.result_path('Fig7_spatial_interaction.png'), dpi=150, bbox_inches='tight')
plt.close(fig7)
print("Fig7_spatial_interaction saved")

# ============================================================
# FIGURE S6 & S7: Spatial overview and validation
# ============================================================
for fig_name, title in [('FigS6_spatial_overview', 'Spatial Overview'),
                         ('FigS7_spatial_validation', 'Spatial Validation')]:
    fig, ax = plt.subplots(figsize=(8, 6))
    if spatial_available:
        # Plot spatial data distribution overview
        numeric_cols = df_sp.select_dtypes(include=[np.number]).columns[:8]
        if len(numeric_cols) > 0:
            df_sp[numeric_cols].boxplot(ax=ax, rot=45)
            ax.set_title(f'{title}: Score distributions')
    else:
        ax.text(0.5, 0.5, f'{title}\n(spatial data not available)',
               ha='center', va='center', transform=ax.transAxes, fontsize=12, color='grey')
        ax.set_title(f'{title}: (unavailable)')
    plt.tight_layout()
    fig.savefig(cfg.result_path(f'{fig_name}.pdf'), dpi=150, bbox_inches='tight')
    fig.savefig(cfg.result_path(f'{fig_name}.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"{fig_name} saved")

print("[step4_evidence_generalization] Done.")
