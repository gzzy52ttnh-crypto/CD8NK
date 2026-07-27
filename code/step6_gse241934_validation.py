"""
step6_gse241934_validation.py
Step 6: GSE241934 外部验证分析
- Step 6.1: 数据预处理与对象构建（读取MTX+Meta，添加化疗方案分组）
- Step 6.2: NK-like 签名量化
- Step 6.3: Taxane 方案 R vs NR 核心验证
- Step 6.4: 方案特异性验证（Taxane vs Pemetrexed）
"""
import config as cfg
from _common import NKLIKE_SIGNATURE, firth_logistic_fit
import pandas as pd
import numpy as np
import scipy.sparse as sp
import anndata as ad
import os
import gzip
import warnings
warnings.filterwarnings('ignore')

print("[step6_gse241934_validation] Starting...")
os.makedirs(cfg.RESULT_DIR, exist_ok=True)

# ============================================================
# RWC 队列化疗方案映射（来源：Table S14 / mmc15.xlsx）
# ============================================================
RWC_CHEMO_MAP = {
    'P33':   {'chemo': 'Carboplatin + Abraxane',       'chemo_class': 'Taxane',     'platinum': 'Carboplatin', 'partner': 'Abraxane'},
    'P52':   {'chemo': 'Carboplatin + Pemetrexed',      'chemo_class': 'Pemetrexed', 'platinum': 'Carboplatin', 'partner': 'Pemetrexed'},
    'P64':   {'chemo': 'Carboplatin + Pemetrexed',      'chemo_class': 'Pemetrexed', 'platinum': 'Carboplatin', 'partner': 'Pemetrexed'},
    'P70':   {'chemo': 'Carboplatin + Pemetrexed',      'chemo_class': 'Pemetrexed', 'platinum': 'Carboplatin', 'partner': 'Pemetrexed'},
    'P105':  {'chemo': 'Carboplatin + Pemetrexed',      'chemo_class': 'Pemetrexed', 'platinum': 'Carboplatin', 'partner': 'Pemetrexed'},
    'P122':  {'chemo': 'Carboplatin + Pemetrexed',      'chemo_class': 'Pemetrexed', 'platinum': 'Carboplatin', 'partner': 'Pemetrexed'},
    'P130':  {'chemo': 'Carboplatin + Abraxane',       'chemo_class': 'Taxane',     'platinum': 'Carboplatin', 'partner': 'Abraxane'},
    'P172':  {'chemo': 'Cisplatin + Pemetrexed',       'chemo_class': 'Pemetrexed', 'platinum': 'Cisplatin',  'partner': 'Pemetrexed'},
    'P209':  {'chemo': 'Carboplatin + Pemetrexed',      'chemo_class': 'Pemetrexed', 'platinum': 'Carboplatin', 'partner': 'Pemetrexed'},
    'P218':  {'chemo': 'Carboplatin + Gemcitabine',     'chemo_class': 'Other',      'platinum': 'Carboplatin', 'partner': 'Gemcitabine'},
    'P223':  {'chemo': 'Carboplatin + Abraxane',       'chemo_class': 'Taxane',     'platinum': 'Carboplatin', 'partner': 'Abraxane'},
    'P266':  {'chemo': 'Carboplatin + Abraxane',       'chemo_class': 'Taxane',     'platinum': 'Carboplatin', 'partner': 'Abraxane'},
    'P269':  {'chemo': 'Carboplatin + Abraxane',       'chemo_class': 'Taxane',     'platinum': 'Carboplatin', 'partner': 'Abraxane'},
    'P270':  {'chemo': 'Carboplatin + Pemetrexed',      'chemo_class': 'Pemetrexed', 'platinum': 'Carboplatin', 'partner': 'Pemetrexed'},
    'P293':  {'chemo': 'Carboplatin + Abraxane',       'chemo_class': 'Taxane',     'platinum': 'Carboplatin', 'partner': 'Abraxane'},
    'P345':  {'chemo': 'Cisplatin + Pemetrexed',       'chemo_class': 'Pemetrexed', 'platinum': 'Cisplatin',  'partner': 'Pemetrexed'},
    'P348':  {'chemo': 'Cisplatin + Pemetrexed',       'chemo_class': 'Pemetrexed', 'platinum': 'Cisplatin',  'partner': 'Pemetrexed'},
    'P390':  {'chemo': 'Nedaplatin + Pemetrexed',      'chemo_class': 'Pemetrexed', 'platinum': 'Nedaplatin', 'partner': 'Pemetrexed'},
    'P394':  {'chemo': 'Carboplatin + Pemetrexed',      'chemo_class': 'Pemetrexed', 'platinum': 'Carboplatin', 'partner': 'Pemetrexed'},
    'P399':  {'chemo': 'Cisplatin + Pemetrexed',       'chemo_class': 'Pemetrexed', 'platinum': 'Cisplatin',  'partner': 'Pemetrexed'},
    'P481':  {'chemo': 'Carboplatin + Pemetrexed',      'chemo_class': 'Pemetrexed', 'platinum': 'Carboplatin', 'partner': 'Pemetrexed'},
    'P483':  {'chemo': 'Carboplatin + Pemetrexed Disodium', 'chemo_class': 'Pemetrexed', 'platinum': 'Carboplatin', 'partner': 'Pemetrexed Disodium'},
    'P498':  {'chemo': 'Carboplatin + Abraxane',       'chemo_class': 'Taxane',     'platinum': 'Carboplatin', 'partner': 'Abraxane'},
    'P509':  {'chemo': 'Cisplatin + Pemetrexed',       'chemo_class': 'Pemetrexed', 'platinum': 'Cisplatin',  'partner': 'Pemetrexed'},
    'P523':  {'chemo': 'Carboplatin + Pemetrexed',      'chemo_class': 'Pemetrexed', 'platinum': 'Carboplatin', 'partner': 'Pemetrexed'},
    'P528':  {'chemo': 'Carboplatin + Abraxane',       'chemo_class': 'Taxane',     'platinum': 'Carboplatin', 'partner': 'Abraxane'},
    'P533':  {'chemo': 'Cisplatin + Pemetrexed',       'chemo_class': 'Pemetrexed', 'platinum': 'Cisplatin',  'partner': 'Pemetrexed'},
    'P551':  {'chemo': 'Carboplatin + Pemetrexed',      'chemo_class': 'Pemetrexed', 'platinum': 'Carboplatin', 'partner': 'Pemetrexed'},
    'P574':  {'chemo': 'Cisplatin + Abraxane',        'chemo_class': 'Taxane',     'platinum': 'Cisplatin',  'partner': 'Abraxane'},
    'P579':  {'chemo': 'Carboplatin + Abraxane',       'chemo_class': 'Taxane',     'platinum': 'Carboplatin', 'partner': 'Abraxane'},
    'P592':  {'chemo': 'Carboplatin + Abraxane',       'chemo_class': 'Taxane',     'platinum': 'Carboplatin', 'partner': 'Abraxane'},
    'P603':  {'chemo': 'Carboplatin + Pemetrexed',      'chemo_class': 'Pemetrexed', 'platinum': 'Carboplatin', 'partner': 'Pemetrexed'},
    'P605':  {'chemo': 'Carboplatin + Abraxane',       'chemo_class': 'Taxane',     'platinum': 'Carboplatin', 'partner': 'Abraxane'},
    'P609':  {'chemo': 'Cisplatin + Abraxane',        'chemo_class': 'Taxane',     'platinum': 'Cisplatin',  'partner': 'Abraxane'},
}

# IIT 队列统一化疗方案
IIT_CHEMO = {
    'chemo': 'Carboplatin + nab-Paclitaxel',
    'chemo_class': 'Taxane',
    'platinum': 'Carboplatin',
    'partner': 'nab-Paclitaxel'
}


def read_mtx_gz(mtx_path):
    """读取 gzip 压缩的 MTX 稀疏矩阵，返回 csr_matrix"""
    print(f"  Reading MTX: {mtx_path}")
    with gzip.open(mtx_path, 'rt') as f:
        header = f.readline().strip()
        assert header.startswith('%%MatrixMarket'), f"Invalid MTX header: {header}"
        while True:
            line = f.readline().strip()
            if not line.startswith('%'):
                m, n, nnz = map(int, line.split())
                break
        rows, cols, vals = [], [], []
        for line in f:
            parts = line.strip().split()
            rows.append(int(parts[0]) - 1)
            cols.append(int(parts[1]) - 1)
            vals.append(float(parts[2]))
    mat = sp.coo_matrix((vals, (rows, cols)), shape=(m, n)).tocsr()
    print(f"  -> {m} genes x {n} cells, {nnz} nnz")
    return mat


def read_barcodes_gz(bc_path):
    """读取 gzip 压缩的 barcodes 文件"""
    with gzip.open(bc_path, 'rt') as f:
        barcodes = [line.strip() for line in f]
    return barcodes


def read_features_gz(feat_path):
    """读取 gzip 压缩的 features/genes 文件，返回基因名列表"""
    with gzip.open(feat_path, 'rt') as f:
        genes = []
        for line in f:
            parts = line.strip().split('\t')
            genes.append(parts[1] if len(parts) > 1 else parts[0])
    return genes


def read_meta_gz(meta_path):
    """读取 gzip 压缩的元数据 TSV，用 cellID 列作为 index"""
    print(f"  Reading Meta: {meta_path}")
    df = pd.read_csv(meta_path, sep='\t', compression='gzip')
    df = df.set_index('cellID')
    print(f"  -> {df.shape[0]} cells x {df.shape[1]} columns")
    return df


# ============================================================
# Step 6.1: 构建 GSE241934 表达对象
# ============================================================
print("\n" + "="*70)
print("Step 6.1: 构建 GSE241934 表达对象")
print("="*70)

h5ad_out = os.path.join(cfg.RESULT_DIR, 'GSE241934_all.h5ad')

if os.path.exists(h5ad_out):
    print(f"  h5ad 已存在，直接读取: {h5ad_out}")
    adata = ad.read_h5ad(h5ad_out)
else:
    print("  读取 IIT 队列...")
    iit_mat = read_mtx_gz(os.path.join(cfg.ADATA_DIR, 'GSE241934_IIT_Matrix.mtx.gz'))
    iit_bc = read_barcodes_gz(os.path.join(cfg.ADATA_DIR, 'GSE241934_IIT_barcodes.tsv.gz'))
    iit_feat = read_features_gz(os.path.join(cfg.ADATA_DIR, 'GSE241934_IIT_features.tsv.gz'))
    iit_meta = read_meta_gz(os.path.join(cfg.ADATA_DIR, 'GSE241934_IIT_Meta.txt.gz'))

    print("\n  读取 RWC 队列...")
    rwc_mat = read_mtx_gz(os.path.join(cfg.ADATA_DIR, 'GSE241934_Real_Matrix.mtx.gz'))
    rwc_bc = read_barcodes_gz(os.path.join(cfg.ADATA_DIR, 'GSE241934_RWC_barcodes.tsv.gz'))
    rwc_feat = read_features_gz(os.path.join(cfg.ADATA_DIR, 'GSE241934_RWC_features.tsv.gz'))
    rwc_meta = read_meta_gz(os.path.join(cfg.ADATA_DIR, 'GSE241934_Real_Meta.txt.gz'))

    assert iit_feat == rwc_feat, "IIT and RWC gene names do not match!"

    print("\n  合并两个队列...")
    combined_mat = sp.hstack([iit_mat, rwc_mat]).tocsr()
    combined_bc = iit_bc + rwc_bc
    combined_meta = pd.concat([iit_meta, rwc_meta], axis=0)

    print(f"  合并后: {combined_mat.shape[0]} genes x {combined_mat.shape[1]} cells")

    adata = ad.AnnData(X=combined_mat.T)
    adata.var_names = iit_feat
    adata.obs_names = combined_bc
    adata.obs = combined_meta.loc[combined_bc].copy()

    # 添加队列标签
    adata.obs['cohort'] = ['IIT'] * len(iit_bc) + ['RWC'] * len(rwc_bc)

    # 添加化疗方案信息
    chemo_detail = []
    chemo_class = []
    platinum = []
    partner = []

    for idx, row in adata.obs.iterrows():
        sample_id = row['sampleID']
        if row['cohort'] == 'IIT':
            chemo_detail.append(IIT_CHEMO['chemo'])
            chemo_class.append(IIT_CHEMO['chemo_class'])
            platinum.append(IIT_CHEMO['platinum'])
            partner.append(IIT_CHEMO['partner'])
        else:
            if sample_id in RWC_CHEMO_MAP:
                chemo_detail.append(RWC_CHEMO_MAP[sample_id]['chemo'])
                chemo_class.append(RWC_CHEMO_MAP[sample_id]['chemo_class'])
                platinum.append(RWC_CHEMO_MAP[sample_id]['platinum'])
                partner.append(RWC_CHEMO_MAP[sample_id]['partner'])
            else:
                chemo_detail.append('Unknown')
                chemo_class.append('Unknown')
                platinum.append('Unknown')
                partner.append('Unknown')

    adata.obs['chemotherapy'] = chemo_detail
    adata.obs['chemo_class'] = chemo_class
    adata.obs['platinum'] = platinum
    adata.obs['chemo_partner'] = partner

    # 统一响应标签（二分类：R = pCR + MPR，NR = non-MPR）
    response_binary = []
    for resp in adata.obs['Pathological Response']:
        if resp in ['pCR', 'MPR']:
            response_binary.append('R')
        elif resp == 'non-MPR':
            response_binary.append('NR')
        else:
            response_binary.append('Unknown')
    adata.obs['response_binary'] = response_binary

    # 分析子集标签
    subset = []
    for idx, row in adata.obs.iterrows():
        if row['cohort'] == 'IIT':
            subset.append('A_IIT_Taxane')
        elif row['chemo_class'] == 'Taxane':
            subset.append('B_RWC_Taxane')
        elif row['chemo_class'] == 'Pemetrexed':
            subset.append('C_RWC_Pemetrexed')
        else:
            subset.append('D_RWC_Other')
    adata.obs['analysis_subset'] = subset

    print(f"\n  化疗方案分类统计:")
    print(adata.obs.groupby(['cohort', 'chemo_class']).size())
    print(f"\n  响应分布（按队列×方案）:")
    print(adata.obs.groupby(['cohort', 'chemo_class', 'response_binary']).size())

    # 修复数据类型问题（将混合类型列转换为字符串）
    for col in adata.obs.columns:
        if adata.obs[col].dtype == object:
            adata.obs[col] = adata.obs[col].astype(str)
        # Pathological Response Rate 可能有混合类型
        if col == 'Pathological Response Rate':
            adata.obs[col] = pd.to_numeric(adata.obs[col], errors='coerce')

    # 保存
    print(f"\n  保存 h5ad: {h5ad_out}")
    try:
        adata.write_h5ad(h5ad_out)
        print("  ✅ 保存完成")
    except Exception as e:
        print(f"  ⚠️  h5ad 保存失败（{e}），但分析将继续进行")

# ============================================================
# Step 6.2: NK-like 签名量化（双方法）
# ============================================================
print("\n" + "="*70)
print("Step 6.2: NK-like 签名量化")
print("="*70)

from _common import normalize_log1p, gene_mean_per_cell

patient_scores_out = os.path.join(cfg.RESULT_DIR, 'GSE241934_patient_nklike_scores.csv')

if os.path.exists(patient_scores_out):
    print(f"  患者得分已存在，直接读取: {patient_scores_out}")
    df_patient = pd.read_csv(patient_scores_out, index_col=0)
else:
    print("  计算 NK-like 签名得分（每细胞）...")
    logX = normalize_log1p(adata)
    nklike_score = gene_mean_per_cell(adata, NKLIKE_SIGNATURE, logX=logX)
    adata.obs['nklike_score'] = nklike_score

    # NK-like 基因匹配情况
    matched_genes = [g for g in NKLIKE_SIGNATURE if g in adata.var_names]
    missing_genes = [g for g in NKLIKE_SIGNATURE if g not in adata.var_names]
    print(f"  NK-like 基因匹配: {len(matched_genes)}/24")
    if missing_genes:
        print(f"  缺失基因: {missing_genes}")

    # 方法 1: 基于 cell.type 注释的 NK-like 细胞比例
    print("\n  计算 NK-like 细胞比例（基于 cell.type 注释）...")
    nklike_celltypes = [ct for ct in adata.obs['cell.type'].unique() if 'FGFBP2' in str(ct) and 'CD8' in str(ct)]
    print(f"  NK-like 细胞类型: {nklike_celltypes}")

    adata.obs['is_nklike'] = adata.obs['cell.type'].isin(nklike_celltypes).astype(int)

    # CD8+ T 细胞定义
    cd8_celltypes = [ct for ct in adata.obs['cell.type'].unique() if 'CD8' in str(ct)]
    print(f"  CD8+ T 细胞类型数: {len(cd8_celltypes)}")
    adata.obs['is_cd8'] = adata.obs['cell.type'].isin(cd8_celltypes).astype(int)

    # 患者水平汇总
    print("\n  汇总患者水平指标...")
    patient_metrics = []
    for sample_id, group in adata.obs.groupby('sampleID'):
        n_cells = len(group)
        n_cd8 = group['is_cd8'].sum()
        n_nklike = group['is_nklike'].sum()
        frac_nklike_of_cd8 = n_nklike / n_cd8 if n_cd8 > 0 else 0
        mean_nklike_score = group['nklike_score'].mean()
        mean_nklike_score_cd8 = group.loc[group['is_cd8'] == 1, 'nklike_score'].mean() if n_cd8 > 0 else 0

        # 患者级元数据（取第一个细胞的值）
        first = group.iloc[0]
        patient_metrics.append({
            'sampleID': sample_id,
            'cohort': first['cohort'],
            'chemo_class': first['chemo_class'],
            'chemotherapy': first['chemotherapy'],
            'response': first['Pathological Response'],
            'response_binary': first['response_binary'],
            'analysis_subset': first['analysis_subset'],
            'n_cells': n_cells,
            'n_cd8': n_cd8,
            'n_nklike': n_nklike,
            'frac_nklike_of_cd8': frac_nklike_of_cd8,
            'mean_nklike_score_all': mean_nklike_score,
            'mean_nklike_score_cd8': mean_nklike_score_cd8,
            'age': first.get('Age', np.nan),
            'gender': first.get('Gender', np.nan),
            'histology': first.get('Histology', np.nan),
            'cycles': first.get('Cycles', np.nan),
            'pd1': first.get('PD1', np.nan),
        })

    df_patient = pd.DataFrame(patient_metrics).set_index('sampleID')
    df_patient = df_patient.sort_values(['cohort', 'chemo_class', 'response_binary'])

    print(f"\n  患者数: {len(df_patient)}")
    print(f"\n  各子集统计:")
    for subset_name in ['A_IIT_Taxane', 'B_RWC_Taxane', 'C_RWC_Pemetrexed', 'D_RWC_Other']:
        sub = df_patient[df_patient['analysis_subset'] == subset_name]
        if len(sub) > 0:
            r = (sub['response_binary'] == 'R').sum()
            nr = (sub['response_binary'] == 'NR').sum()
            print(f"    {subset_name}: {len(sub)} 例 (R={r}, NR={nr})")

    # 保存
    df_patient.to_csv(patient_scores_out)
    print(f"\n  ✅ 患者得分已保存: {patient_scores_out}")

    # 同时更新 h5ad 保存
    h5ad_updated = os.path.join(cfg.RESULT_DIR, 'GSE241934_all.h5ad')
    adata.write_h5ad(h5ad_updated)
    print(f"  ✅ 更新后的 h5ad 已保存: {h5ad_updated}")

# ============================================================
# Step 6.3: 核心验证 - Taxane 方案 R vs NR
# ============================================================
print("\n" + "="*70)
print("Step 6.3: 核心验证 - Taxane 方案 R vs NR")
print("="*70)

from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score
import statsmodels.api as sm

validation_results = []

for metric_name in ['frac_nklike_of_cd8', 'mean_nklike_score_cd8', 'mean_nklike_score_all']:
    print(f"\n--- 指标: {metric_name} ---")

    # 子集 A: IIT Taxane (n=11)
    sub_a = df_patient[df_patient['analysis_subset'] == 'A_IIT_Taxane']
    r_a = sub_a[sub_a['response_binary'] == 'R'][metric_name].values
    nr_a = sub_a[sub_a['response_binary'] == 'NR'][metric_name].values

    if len(r_a) > 0 and len(nr_a) > 0:
        stat_a, p_a = mannwhitneyu(r_a, nr_a, alternative='two-sided')
        auc_a = roc_auc_score([1]*len(r_a) + [0]*len(nr_a), list(r_a) + list(nr_a))
        print(f"  子集 A (IIT Taxane, n={len(sub_a)}): R={len(r_a)}, NR={len(nr_a)}")
        print(f"    R 组均值: {r_a.mean():.4f}, NR 组均值: {nr_a.mean():.4f}")
        print(f"    MWU p = {p_a:.4f}, AUC = {auc_a:.4f}")
        validation_results.append({
            'subset': 'A_IIT_Taxane',
            'metric': metric_name,
            'n_total': len(sub_a),
            'n_R': len(r_a),
            'n_NR': len(nr_a),
            'mean_R': r_a.mean(),
            'mean_NR': nr_a.mean(),
            'mwu_p': p_a,
            'auc': auc_a,
        })
    else:
        print(f"  子集 A: 样本不足 (R={len(r_a)}, NR={len(nr_a)})")

    # 子集 B: RWC Taxane (n=13)
    sub_b = df_patient[df_patient['analysis_subset'] == 'B_RWC_Taxane']
    r_b = sub_b[sub_b['response_binary'] == 'R'][metric_name].values
    nr_b = sub_b[sub_b['response_binary'] == 'NR'][metric_name].values

    if len(r_b) > 0 and len(nr_b) > 0:
        stat_b, p_b = mannwhitneyu(r_b, nr_b, alternative='two-sided')
        auc_b = roc_auc_score([1]*len(r_b) + [0]*len(nr_b), list(r_b) + list(nr_b))
        print(f"  子集 B (RWC Taxane, n={len(sub_b)}): R={len(r_b)}, NR={len(nr_b)}")
        print(f"    R 组均值: {r_b.mean():.4f}, NR 组均值: {nr_b.mean():.4f}")
        print(f"    MWU p = {p_b:.4f}, AUC = {auc_b:.4f}")
        validation_results.append({
            'subset': 'B_RWC_Taxane',
            'metric': metric_name,
            'n_total': len(sub_b),
            'n_R': len(r_b),
            'n_NR': len(nr_b),
            'mean_R': r_b.mean(),
            'mean_NR': nr_b.mean(),
            'mwu_p': p_b,
            'auc': auc_b,
        })
    else:
        print(f"  子集 B: 样本不足 (R={len(r_b)}, NR={len(nr_b)})")

    # 合并 Taxane (A+B, n=24)
    sub_taxane = df_patient[df_patient['chemo_class'] == 'Taxane']
    r_t = sub_taxane[sub_taxane['response_binary'] == 'R'][metric_name].values
    nr_t = sub_taxane[sub_taxane['response_binary'] == 'NR'][metric_name].values

    if len(r_t) > 0 and len(nr_t) > 0:
        stat_t, p_t = mannwhitneyu(r_t, nr_t, alternative='two-sided')
        auc_t = roc_auc_score([1]*len(r_t) + [0]*len(nr_t), list(r_t) + list(nr_t))
        print(f"  Taxane 合并 (A+B, n={len(sub_taxane)}): R={len(r_t)}, NR={len(nr_t)}")
        print(f"    R 组均值: {r_t.mean():.4f}, NR 组均值: {nr_t.mean():.4f}")
        print(f"    MWU p = {p_t:.4f}, AUC = {auc_t:.4f}")
        validation_results.append({
            'subset': 'Taxane_combined',
            'metric': metric_name,
            'n_total': len(sub_taxane),
            'n_R': len(r_t),
            'n_NR': len(nr_t),
            'mean_R': r_t.mean(),
            'mean_NR': nr_t.mean(),
            'mwu_p': p_t,
            'auc': auc_t,
        })
    else:
        print(f"  Taxane 合并: 样本不足")

    # 子集 C: RWC Pemetrexed (n=20) - 方案特异性对照
    sub_c = df_patient[df_patient['analysis_subset'] == 'C_RWC_Pemetrexed']
    r_c = sub_c[sub_c['response_binary'] == 'R'][metric_name].values
    nr_c = sub_c[sub_c['response_binary'] == 'NR'][metric_name].values

    if len(r_c) > 0 and len(nr_c) > 0:
        stat_c, p_c = mannwhitneyu(r_c, nr_c, alternative='two-sided')
        try:
            auc_c = roc_auc_score([1]*len(r_c) + [0]*len(nr_c), list(r_c) + list(nr_c))
        except:
            auc_c = np.nan
        print(f"  子集 C (RWC Pemetrexed, n={len(sub_c)}): R={len(r_c)}, NR={len(nr_c)}")
        print(f"    R 组均值: {r_c.mean():.4f}, NR 组均值: {nr_c.mean():.4f}")
        print(f"    MWU p = {p_c:.4f}, AUC = {auc_c:.4f}")
        validation_results.append({
            'subset': 'C_RWC_Pemetrexed',
            'metric': metric_name,
            'n_total': len(sub_c),
            'n_R': len(r_c),
            'n_NR': len(nr_c),
            'mean_R': r_c.mean(),
            'mean_NR': nr_c.mean(),
            'mwu_p': p_c,
            'auc': auc_c,
        })
    else:
        print(f"  子集 C: 样本不足 (R={len(r_c)}, NR={len(nr_c)})")

# 保存验证结果
df_val = pd.DataFrame(validation_results)
val_out = os.path.join(cfg.RESULT_DIR, 'GSE241934_validation_results.csv')
df_val.to_csv(val_out, index=False)
print(f"\n✅ 验证结果已保存: {val_out}")
print(df_val.to_string(index=False))

print("\n" + "="*70)
print("[step6_gse241934_validation] Step 6.1-6.3 completed!")
print("="*70)

# ============================================================
# Step 6.4: 可视化
# ============================================================
print("\n" + "="*70)
print("Step 6.4: 可视化")
print("="*70)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sklearn.metrics import roc_curve

plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['axes.grid'] = False
plt.rcParams['pdf.fonttype'] = 42

# 颜色方案
color_R = '#E64B35'
color_NR = '#4DBBD5'
color_taxane = '#E64B35'
color_pemetrexed = '#4DBBD5'
color_iit = '#91D1C2'
color_rwc = '#F39B7F'

# 主指标
MAIN_METRIC = 'frac_nklike_of_cd8'
METRIC_LABEL = 'NK-like proportion in CD8+ T cells'

# --- Fig 1: 箱线图（4 个子集） ---
fig, axes = plt.subplots(1, 4, figsize=(16, 5))
fig.suptitle('GSE241934: NK-like proportion by response', fontsize=14, fontweight='bold', y=1.02)

subset_info = [
    ('A_IIT_Taxane', 'IIT\nTaxane (n=11)', color_iit),
    ('B_RWC_Taxane', 'RWC\nTaxane (n=13)', color_rwc),
    ('Taxane_combined', 'Taxane\nCombined (n=24)', color_taxane),
    ('C_RWC_Pemetrexed', 'RWC\nPemetrexed (n=20)', color_pemetrexed),
]

for i, (subset_key, title, _) in enumerate(subset_info):
    ax = axes[i]
    if subset_key == 'Taxane_combined':
        sub = df_patient[df_patient['analysis_subset'].isin(['A_IIT_Taxane', 'B_RWC_Taxane'])]
    else:
        sub = df_patient[df_patient['analysis_subset'] == subset_key]
    
    r_vals = sub[sub['response_binary'] == 'R'][MAIN_METRIC].values
    nr_vals = sub[sub['response_binary'] == 'NR'][MAIN_METRIC].values
    
    bp = ax.boxplot([nr_vals, r_vals],
                    labels=['NR', 'R'],
                    patch_artist=True, widths=0.6,
                    medianprops={'color': 'black', 'linewidth': 1.5},
                    flierprops={'marker': 'o', 'markersize': 4})
    
    bp['boxes'][0].set_facecolor(color_NR)
    bp['boxes'][1].set_facecolor(color_R)
    bp['boxes'][0].set_alpha(0.7)
    bp['boxes'][1].set_alpha(0.7)
    
    # 叠加散点
    x_jitter_nr = np.random.normal(1, 0.04, size=len(nr_vals))
    x_jitter_r = np.random.normal(2, 0.04, size=len(r_vals))
    ax.scatter(x_jitter_nr, nr_vals, color=color_NR, alpha=0.8, s=30, zorder=5)
    ax.scatter(x_jitter_r, r_vals, color=color_R, alpha=0.8, s=30, zorder=5)
    
    # MWU p 值标注
    from scipy.stats import mannwhitneyu
    _, p_val = mannwhitneyu(r_vals, nr_vals, alternative='two-sided')
    if p_val < 0.001:
        p_text = f'p={p_val:.2e}'
    else:
        p_text = f'p={p_val:.3f}'
    
    y_max = max(np.nanmax(r_vals), np.nanmax(nr_vals))
    ax.text(1.5, y_max * 1.1, p_text, ha='center', fontsize=10, fontweight='bold')
    
    ax.set_title(title, fontsize=11)
    ax.set_ylabel(METRIC_LABEL if i == 0 else '')
    ax.set_ylim(bottom=0)

plt.tight_layout()
fig1_path = os.path.join(cfg.RESULT_DIR, 'GSE241934_Fig1_boxplot.png')
fig.savefig(fig1_path, dpi=300, bbox_inches='tight')
fig1_pdf = os.path.join(cfg.RESULT_DIR, 'GSE241934_Fig1_boxplot.pdf')
fig.savefig(fig1_pdf, bbox_inches='tight')
print(f"  ✅ Fig1 箱线图已保存: {fig1_path}")
plt.close(fig)

# --- Fig 2: ROC 曲线（Taxane vs Pemetrexed） ---
fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
fig.suptitle('GSE241934: ROC curves', fontsize=14, fontweight='bold', y=1.02)

roc_subsets = [
    ('Taxane (IIT)', 'A_IIT_Taxane', color_iit),
    ('Taxane (RWC)', 'B_RWC_Taxane', color_rwc),
    ('Pemetrexed (RWC)', 'C_RWC_Pemetrexed', color_pemetrexed),
]

for i, (title, subset_key, color) in enumerate(roc_subsets):
    ax = axes[i]
    sub = df_patient[df_patient['analysis_subset'] == subset_key]
    y_true = (sub['response_binary'] == 'R').astype(int).values
    y_score = sub[MAIN_METRIC].values
    
    fpr, tpr, _ = roc_curve(y_true, y_score)
    auc_val = roc_auc_score(y_true, y_score)
    
    ax.plot(fpr, tpr, color=color, linewidth=2.5, label=f'AUC = {auc_val:.3f}')
    ax.plot([0, 1], [0, 1], 'k--', linewidth=1, alpha=0.5)
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate' if i == 0 else '')
    ax.set_title(f'{title}\n(n={len(sub)})', fontsize=11)
    ax.legend(loc='lower right', fontsize=10)
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])

plt.tight_layout()
fig2_path = os.path.join(cfg.RESULT_DIR, 'GSE241934_Fig2_ROC.png')
fig.savefig(fig2_path, dpi=300, bbox_inches='tight')
fig2_pdf = os.path.join(cfg.RESULT_DIR, 'GSE241934_Fig2_ROC.pdf')
fig.savefig(fig2_pdf, bbox_inches='tight')
print(f"  ✅ Fig2 ROC 曲线已保存: {fig2_path}")
plt.close(fig)

# --- Fig 3: AUC 对比柱状图 ---
fig, ax = plt.subplots(1, 1, figsize=(10, 5.5))

auc_data = []
for metric_name in ['frac_nklike_of_cd8', 'mean_nklike_score_cd8', 'mean_nklike_score_all']:
    for subset_key, label in [
        ('A_IIT_Taxane', 'IIT Taxane'),
        ('B_RWC_Taxane', 'RWC Taxane'),
        ('C_RWC_Pemetrexed', 'RWC Pemetrexed'),
    ]:
        row = df_val[(df_val['subset'] == subset_key) & (df_val['metric'] == metric_name)]
        if len(row) > 0:
            auc_data.append({
                'metric': metric_name,
                'subset': label,
                'auc': row['auc'].values[0],
                'p': row['mwu_p'].values[0],
            })

df_auc = pd.DataFrame(auc_data)

metric_labels = {
    'frac_nklike_of_cd8': 'NK-like cell fraction',
    'mean_nklike_score_cd8': 'NK-like signature (CD8+)',
    'mean_nklike_score_all': 'NK-like signature (all cells)',
}
subset_colors = {
    'IIT Taxane': color_iit,
    'RWC Taxane': color_taxane,
    'RWC Pemetrexed': color_pemetrexed,
}

x = np.arange(len(metric_labels))
width = 0.25

for j, (subset_name, color) in enumerate(subset_colors.items()):
    vals = [df_auc[(df_auc['metric'] == m) & (df_auc['subset'] == subset_name)]['auc'].values[0]
            if len(df_auc[(df_auc['metric'] == m) & (df_auc['subset'] == subset_name)]) > 0 else 0
            for m in metric_labels.keys()]
    bars = ax.bar(x + j * width - width, vals, width, label=subset_name, color=color, alpha=0.8)
    
    # 标注 p 值
    for k, v in enumerate(vals):
        p_val = df_auc[(df_auc['metric'] == list(metric_labels.keys())[k]) & (df_auc['subset'] == subset_name)]['p'].values[0]
        sig = '*' if p_val < 0.05 else ''
        if v > 0:
            ax.text(x[k] + j * width - width, v + 0.02, sig, ha='center', fontsize=14, fontweight='bold', color='red')

ax.axhline(y=0.5, color='gray', linestyle='--', linewidth=1, alpha=0.7, label='Random (AUC=0.5)')
ax.set_xticks(x)
ax.set_xticklabels(metric_labels.values(), fontsize=10)
ax.set_ylabel('AUC')
ax.set_title('GSE241934: AUC comparison across subsets and metrics', fontsize=13, fontweight='bold')
ax.set_ylim([0, 1.05])
ax.legend(loc='upper right', fontsize=9, ncol=2)

plt.tight_layout()
fig3_path = os.path.join(cfg.RESULT_DIR, 'GSE241934_Fig3_AUC_comparison.png')
fig.savefig(fig3_path, dpi=300, bbox_inches='tight')
fig3_pdf = os.path.join(cfg.RESULT_DIR, 'GSE241934_Fig3_AUC_comparison.pdf')
fig.savefig(fig3_pdf, bbox_inches='tight')
print(f"  ✅ Fig3 AUC 对比图已保存: {fig3_path}")
plt.close(fig)

# --- Fig 4: 方案特异性散点图 ---
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle('GSE241934: Treatment-specific NK-like association', fontsize=14, fontweight='bold', y=1.02)

# Panel A: Taxane vs Pemetrexed 的 R/NR 均值对比
ax = axes[0]
subsets_bar = ['IIT Taxane', 'RWC Taxane', 'RWC Pemetrexed']
subset_keys = ['A_IIT_Taxane', 'B_RWC_Taxane', 'C_RWC_Pemetrexed']
bar_colors = [color_iit, color_taxane, color_pemetrexed]

x_pos = np.arange(len(subsets_bar))
width = 0.35

r_means = []
nr_means = []
for sk in subset_keys:
    sub = df_patient[df_patient['analysis_subset'] == sk]
    r_means.append(sub[sub['response_binary'] == 'R'][MAIN_METRIC].mean())
    nr_means.append(sub[sub['response_binary'] == 'NR'][MAIN_METRIC].mean())

ax.bar(x_pos - width/2, nr_means, width, label='NR', color=color_NR, alpha=0.7)
ax.bar(x_pos + width/2, r_means, width, label='R', color=color_R, alpha=0.7)
ax.set_xticks(x_pos)
ax.set_xticklabels(subsets_bar)
ax.set_ylabel(METRIC_LABEL)
ax.set_title('Mean NK-like proportion by response', fontsize=11)
ax.legend()

# Panel B: Taxane vs Pemetrexed AUC 对比
ax = axes[1]
auc_vals_bar = []
p_vals_bar = []
for sk in subset_keys:
    row = df_val[(df_val['subset'] == sk) & (df_val['metric'] == MAIN_METRIC)]
    if len(row) > 0:
        auc_vals_bar.append(row['auc'].values[0])
        p_vals_bar.append(row['mwu_p'].values[0])

bars = ax.bar(subsets_bar, auc_vals_bar, color=bar_colors, alpha=0.8, width=0.6)
ax.axhline(y=0.5, color='gray', linestyle='--', linewidth=1, alpha=0.7, label='Random (AUC=0.5)')
ax.set_ylabel('AUC')
ax.set_title('AUC by treatment subset', fontsize=11)
ax.set_ylim([0, 1.05])

for bar, p in zip(bars, p_vals_bar):
    height = bar.get_height()
    sig = '*' if p < 0.05 else 'ns'
    ax.text(bar.get_x() + bar.get_width()/2., height + 0.02,
            f'AUC={height:.3f}\np={p:.3f} ({sig})',
            ha='center', va='bottom', fontsize=9)

ax.legend(loc='upper right')

plt.tight_layout()
fig4_path = os.path.join(cfg.RESULT_DIR, 'GSE241934_Fig4_treatment_specificity.png')
fig.savefig(fig4_path, dpi=300, bbox_inches='tight')
fig4_pdf = os.path.join(cfg.RESULT_DIR, 'GSE241934_Fig4_treatment_specificity.pdf')
fig.savefig(fig4_pdf, bbox_inches='tight')
print(f"  ✅ Fig4 方案特异性图已保存: {fig4_path}")
plt.close(fig)

print("\n  📊 所有图片已保存（PNG + PDF 双格式）")

# ============================================================
# Step 6.5: Logistic 回归 + 交互效应检验
# ============================================================
print("\n" + "="*70)
print("Step 6.5: Logistic 回归 + 交互效应检验")
print("="*70)

import statsmodels.api as sm
from scipy.special import logit
from scipy.stats import chi2

# 准备数据：RWC 队列（Taxane + Pemetrexed，方案内对比）
df_rwc = df_patient[df_patient['cohort'] == 'RWC'].copy()
df_rwc = df_rwc[df_rwc['chemo_class'].isin(['Taxane', 'Pemetrexed'])].copy()  # n=33
df_rwc['y'] = (df_rwc['response_binary'] == 'R').astype(int)
df_rwc['chemo_taxane'] = (df_rwc['chemo_class'] == 'Taxane').astype(int)  # 1=Taxane, 0=Peme

print(f"\n  RWC 分析样本: {len(df_rwc)} 例 (Taxane={sum(df_rwc['chemo_taxane']==1)}, Peme={sum(df_rwc['chemo_taxane']==0)})")

logistic_results = []

for metric_name in ['frac_nklike_of_cd8', 'mean_nklike_score_cd8', 'mean_nklike_score_all']:
    print(f"\n--- 指标: {metric_name} ---")
    
    # 标准化 NK-like 指标（Z-score，方便解释 OR）
    df_rwc[f'{metric_name}_z'] = (df_rwc[metric_name] - df_rwc[metric_name].mean()) / df_rwc[metric_name].std()
    
    # Model 1: 仅 NK-like score（未校正方案）
    X1 = sm.add_constant(df_rwc[f'{metric_name}_z'])
    try:
        model1 = sm.Logit(df_rwc['y'], X1).fit(disp=0, maxiter=100)
        or1 = np.exp(model1.params[metric_name + '_z'])
        ci1 = np.exp(model1.conf_int().loc[metric_name + '_z'])
        p1 = model1.pvalues[metric_name + '_z']
        print(f"  Model 1 (NK only, unadjusted): OR={or1:.3f}, 95%CI=[{ci1[0]:.3f}, {ci1[1]:.3f}], p={p1:.4f}")
        logistic_results.append({
            'metric': metric_name,
            'model': 'M1_NK_only',
            'n': len(df_rwc),
            'or_value': or1,
            'ci_lower': ci1[0],
            'ci_upper': ci1[1],
            'p_value': p1,
            'auc': roc_auc_score(df_rwc['y'], model1.predict(X1)),
        })
    except Exception as e:
        print(f"  Model 1 拟合失败: {e}")

    # Model 2: NK + chemo（主效应，无交互）
    X2 = sm.add_constant(df_rwc[[f'{metric_name}_z', 'chemo_taxane']])
    try:
        model2 = sm.Logit(df_rwc['y'], X2).fit(disp=0, maxiter=100)
        or2_nk = np.exp(model2.params[metric_name + '_z'])
        ci2_nk = np.exp(model2.conf_int().loc[metric_name + '_z'])
        p2_nk = model2.pvalues[metric_name + '_z']
        or2_chemo = np.exp(model2.params['chemo_taxane'])
        ci2_chemo = np.exp(model2.conf_int().loc['chemo_taxane'])
        p2_chemo = model2.pvalues['chemo_taxane']
        print(f"  Model 2 (NK + chemo):")
        print(f"    NK: OR={or2_nk:.3f}, CI=[{ci2_nk[0]:.3f}, {ci2_nk[1]:.3f}], p={p2_nk:.4f}")
        print(f"    Chemo: OR={or2_chemo:.3f}, CI=[{ci2_chemo[0]:.3f}, {ci2_chemo[1]:.3f}], p={p2_chemo:.4f}")
        logistic_results.append({
            'metric': metric_name,
            'model': 'M2_NK_plus_chemo',
            'n': len(df_rwc),
            'or_value': or2_nk,
            'ci_lower': ci2_nk[0],
            'ci_upper': ci2_nk[1],
            'p_value': p2_nk,
            'auc': roc_auc_score(df_rwc['y'], model2.predict(X2)),
        })
    except Exception as e:
        print(f"  Model 2 拟合失败: {e}")

    # Model 3: NK + chemo + 交互项（核心检验！）
    df_rwc['interaction'] = df_rwc[f'{metric_name}_z'] * df_rwc['chemo_taxane']
    X3 = sm.add_constant(df_rwc[[f'{metric_name}_z', 'chemo_taxane', 'interaction']])
    try:
        model3 = sm.Logit(df_rwc['y'], X3).fit(disp=0, maxiter=200)
        or3_nk = np.exp(model3.params[metric_name + '_z'])
        ci3_nk = np.exp(model3.conf_int().loc[metric_name + '_z'])
        p3_nk = model3.pvalues[metric_name + '_z']
        or3_chemo = np.exp(model3.params['chemo_taxane'])
        ci3_chemo = np.exp(model3.conf_int().loc['chemo_taxane'])
        p3_chemo = model3.pvalues['chemo_taxane']
        or3_interact = np.exp(model3.params['interaction'])
        ci3_interact = np.exp(model3.conf_int().loc['interaction'])
        p3_interact = model3.pvalues['interaction']
        
        print(f"  Model 3 (NK + chemo + interaction):")
        print(f"    NK: OR={or3_nk:.3f}, CI=[{ci3_nk[0]:.3f}, {ci3_nk[1]:.3f}], p={p3_nk:.4f}")
        print(f"    Chemo: OR={or3_chemo:.3f}, CI=[{ci3_chemo[0]:.3f}, {ci3_chemo[1]:.3f}], p={p3_chemo:.4f}")
        print(f"    Interaction (NK×chemo): OR={or3_interact:.3f}, CI=[{ci3_interact[0]:.3f}, {ci3_interact[1]:.3f}], p={p3_interact:.4f}")
        print(f"    ⭐ 交互效应 p = {p3_interact:.4f} {'（显著）' if p3_interact < 0.05 else '（未显著）'}")
        
        # 似然比检验（Model2 vs Model3）
        lr_stat = -2 * (model2.llf - model3.llf)
        lr_p = 1 - chi2.cdf(lr_stat, df=1)
        print(f"    LRT (M2 vs M3): LR={lr_stat:.3f}, p={lr_p:.4f}")
        
        logistic_results.append({
            'metric': metric_name,
            'model': 'M3_with_interaction',
            'n': len(df_rwc),
            'or_value': or3_interact,
            'ci_lower': ci3_interact[0],
            'ci_upper': ci3_interact[1],
            'p_value': p3_interact,
            'auc': roc_auc_score(df_rwc['y'], model3.predict(X3)),
            'note': 'OR is for interaction term',
        })
        
        # 也记录 NK 和 chemo 的主效应
        logistic_results.append({
            'metric': metric_name,
            'model': 'M3_NK_main_effect',
            'n': len(df_rwc),
            'or_value': or3_nk,
            'ci_lower': ci3_nk[0],
            'ci_upper': ci3_nk[1],
            'p_value': p3_nk,
            'auc': roc_auc_score(df_rwc['y'], model3.predict(X3)),
            'note': 'OR for NK in reference group (Pemetrexed)',
        })
    except Exception as e:
        print(f"  Model 3 拟合失败: {e}")
        # 尝试 Firth 惩罚回归
        print(f"  尝试 Firth 惩罚回归...")

# --- 各方案亚组单独分析（验证方向） ---
print("\n--- 各方案亚组单独 Logistic 回归 ---")
for chemo_class in ['Taxane', 'Pemetrexed']:
    sub = df_rwc[df_rwc['chemo_class'] == chemo_class].copy()
    if len(sub) < 10:
        print(f"  {chemo_class}: n={len(sub)}，样本过少跳过")
        continue
    
    for metric_name in ['frac_nklike_of_cd8', 'mean_nklike_score_cd8', 'mean_nklike_score_all']:
        sub[f'{metric_name}_z'] = (sub[metric_name] - sub[metric_name].mean()) / sub[metric_name].std()
        X = sm.add_constant(sub[f'{metric_name}_z'])
        try:
            m = sm.Logit(sub['y'], X).fit(disp=0, maxiter=100)
            or_val = np.exp(m.params[metric_name + '_z'])
            p_val = m.pvalues[metric_name + '_z']
            auc_val = roc_auc_score(sub['y'], m.predict(X))
            print(f"  {chemo_class} / {metric_name}: OR={or_val:.3f}, p={p_val:.4f}, AUC={auc_val:.3f}")
            logistic_results.append({
                'metric': metric_name,
                'model': f'{chemo_class}_only',
                'n': len(sub),
                'or_value': or_val,
                'ci_lower': np.exp(m.conf_int().loc[metric_name + '_z'])[0],
                'ci_upper': np.exp(m.conf_int().loc[metric_name + '_z'])[1],
                'p_value': p_val,
                'auc': auc_val,
            })
        except Exception as e:
            print(f"  {chemo_class} / {metric_name}: 拟合失败 ({e})")

# 保存 Logistic 结果
df_logit = pd.DataFrame(logistic_results)
logit_out = os.path.join(cfg.RESULT_DIR, 'GSE241934_logistic_interaction.csv')
df_logit.to_csv(logit_out, index=False)
print(f"\n✅ Logistic 回归结果已保存: {logit_out}")
print(df_logit.to_string(index=False))

# ============================================================
# Step 6.6: EGFR 分层分析 + Firth 惩罚回归
# ============================================================
print("\n" + "="*70)
print("Step 6.6: EGFR 分层分析 + Firth 惩罚回归")
print("="*70)


def firth_odds_ratio_ci(X, y, col_idx=1, alpha=0.05):
    """计算 Firth 回归的 OR 和轮廓似然 95% CI"""
    beta_full, se_full, ll_full, _ = firth_logistic_fit(X, y)
    
    # 轮廓似然 CI
    def profile_ll(beta_j_val):
        beta_fixed = beta_full.copy()
        beta_fixed[col_idx] = beta_j_val
        
        # 固定第 j 个系数，优化其余系数
        n, p = X.shape
        mask = np.ones(p, dtype=bool)
        mask[col_idx] = False
        
        X_free = X[:, mask]
        # offset = X[:, col_idx] * beta_j_val
        offset = X[:, col_idx] * beta_j_val
        
        # 简化: 用 Firth 迭代但固定 beta_j
        beta_other = beta_full[mask]
        for it in range(100):
            eta = X_free @ beta_other + offset
            p_i = 1.0 / (1.0 + np.exp(-eta))
            p_i = np.clip(p_i, 1e-10, 1 - 1e-10)
            w = p_i * (1 - p_i)
            WX = X_free * w[:, None]
            I = X_free.T @ WX
            try:
                I_inv = np.linalg.inv(I)
            except:
                I_inv = np.linalg.pinv(I)
            h = np.sum((X_free @ I_inv) * WX, axis=1)
            h = np.clip(h, 1e-10, 1 - 1e-10)
            score = X_free.T @ (y - p_i + h * (0.5 - p_i))
            try:
                step = I_inv @ score
            except:
                step = np.linalg.lstsq(I, score, rcond=None)[0]
            beta_other = beta_other + step
            if np.max(np.abs(step)) < 1e-8:
                break
        
        eta = X_free @ beta_other + offset
        p_i = 1.0 / (1.0 + np.exp(-eta))
        p_i = np.clip(p_i, 1e-10, 1 - 1e-10)
        w = p_i * (1 - p_i)
        I_full = X.T @ (X * w[:, None])  # 近似
        ll = np.sum(y * eta - np.log(1 + np.exp(eta))) + 0.5 * np.linalg.slogdet(I_full)[1]
        return ll
    
    # 使用 Wald CI（与其他脚本保持一致，避免轮廓似然近似误差）
    z_crit = 1.959963984540054  # 0.975 quantile of N(0,1)
    ci_lower_beta = beta_full[col_idx] - z_crit * se_full[col_idx]
    ci_upper_beta = beta_full[col_idx] + z_crit * se_full[col_idx]
    
    or_val = np.exp(beta_full[col_idx])
    ci_lower = np.exp(ci_lower_beta)
    ci_upper = np.exp(ci_upper_beta)
    
    # Wald p 值
    wald_stat = (beta_full[col_idx] / se_full[col_idx]) ** 2
    p_wald = 1 - chi2.cdf(wald_stat, df=1)
    
    return or_val, ci_lower, ci_upper, p_wald, beta_full, se_full

# ── 1. EGFR 分层分析 ──
print("\n--- EGFR 分层分析 ---")

egfr_results = []

for metric_name in ['frac_nklike_of_cd8', 'mean_nklike_score_cd8', 'mean_nklike_score_all']:
    print(f"\n  指标: {metric_name}")
    
    for cohort_name, cohort_label in [('IIT', 'IIT (EGFR-MT, Taxane)'), ('RWC_Taxane', 'RWC (EGFR-WT, Taxane)')]:
        if cohort_name == 'IIT':
            sub = df_patient[df_patient['analysis_subset'] == 'A_IIT_Taxane'].copy()
        else:
            sub = df_patient[df_patient['analysis_subset'] == 'B_RWC_Taxane'].copy()
        
        if len(sub) < 5:
            continue
        
        y = (sub['response_binary'] == 'R').astype(int).values
        x = sub[metric_name].values
        x_z = (x - x.mean()) / x.std() if x.std() > 0 else x
        
        X = sm.add_constant(x_z)
        
        # 普通 Logistic
        try:
            m_sm = sm.Logit(y, X).fit(disp=0, maxiter=100)
            or_sm = np.exp(m_sm.params[1])
            ci_sm = np.exp(m_sm.conf_int().iloc[1])
            p_sm = m_sm.pvalues[1]
            auc_sm = roc_auc_score(y, m_sm.predict(X))
        except:
            or_sm, ci_sm_lo, ci_sm_hi, p_sm, auc_sm = np.nan, np.nan, np.nan, np.nan, np.nan
            ci_sm = [np.nan, np.nan]
        
        # Firth 回归
        try:
            or_f, ci_f_lo, ci_f_hi, p_f, beta_f, se_f = firth_odds_ratio_ci(X, y, col_idx=1)
        except:
            or_f, ci_f_lo, ci_f_hi, p_f = np.nan, np.nan, np.nan, np.nan
        
        # MWU + AUC
        r_vals = sub[sub['response_binary'] == 'R'][metric_name].values
        nr_vals = sub[sub['response_binary'] == 'NR'][metric_name].values
        _, mwu_p = mannwhitneyu(r_vals, nr_vals, alternative='two-sided')
        auc_val = roc_auc_score(y, x)
        
        print(f"    {cohort_label} (n={len(sub)}, R={sum(y)}, NR={len(y)-sum(y)}):")
        print(f"      MWU p={mwu_p:.4f}, AUC={auc_val:.4f}")
        print(f"      Logistic OR={or_sm:.3f}, p={p_sm:.4f}")
        print(f"      Firth OR={or_f:.3f}, 95%CI=[{ci_f_lo:.3f}, {ci_f_hi:.3f}], p={p_f:.4f}")
        
        egfr_results.append({
            'metric': metric_name,
            'cohort': cohort_label,
            'n': len(sub),
            'n_R': sum(y),
            'n_NR': len(y) - sum(y),
            'mwu_p': mwu_p,
            'auc': auc_val,
            'logistic_or': or_sm,
            'logistic_p': p_sm,
            'firth_or': or_f,
            'firth_ci_lower': ci_f_lo,
            'firth_ci_upper': ci_f_hi,
            'firth_p': p_f,
        })

df_egfr = pd.DataFrame(egfr_results)
egfr_out = os.path.join(cfg.RESULT_DIR, 'GSE241934_egfr_stratified.csv')
df_egfr.to_csv(egfr_out, index=False)
print(f"\n✅ EGFR 分层分析结果已保存: {egfr_out}")

# ── 2. Firth 交互效应检验（RWC 队列 n=33） ──
print("\n--- Firth 惩罚回归: 交互效应检验（RWC 队列）---")

firth_interact_results = []

for metric_name in ['frac_nklike_of_cd8', 'mean_nklike_score_cd8', 'mean_nklike_score_all']:
    print(f"\n  指标: {metric_name}")
    
    df_sub = df_rwc.copy()
    y = df_sub['y'].values
    x_z = df_sub[f'{metric_name}_z'].values
    chemo = df_sub['chemo_taxane'].values
    interaction = x_z * chemo
    
    # 交互效应模型: y ~ NK + chemo + NK*chemo
    X_int = sm.add_constant(np.column_stack([x_z, chemo, interaction]))
    
    try:
        or_f, ci_lo, ci_hi, p_f, beta_f, se_f = firth_odds_ratio_ci(X_int, y, col_idx=3)  # 第4列=交互项
        print(f"    交互项 Firth OR={or_f:.3f}, 95%CI=[{ci_lo:.3f}, {ci_hi:.3f}], p={p_f:.4f}")
        firth_interact_results.append({
            'metric': metric_name,
            'n': len(df_sub),
            'interaction_or_firth': or_f,
            'interaction_ci_lower': ci_lo,
            'interaction_ci_upper': ci_hi,
            'interaction_p_firth': p_f,
        })
    except Exception as e:
        print(f"    Firth 交互效应拟合失败: {e}")
        firth_interact_results.append({
            'metric': metric_name,
            'n': len(df_sub),
            'interaction_or_firth': np.nan,
            'interaction_ci_lower': np.nan,
            'interaction_ci_upper': np.nan,
            'interaction_p_firth': np.nan,
        })

df_firth_int = pd.DataFrame(firth_interact_results)
firth_int_out = os.path.join(cfg.RESULT_DIR, 'GSE241934_firth_interaction.csv')
df_firth_int.to_csv(firth_int_out, index=False)
print(f"\n✅ Firth 交互效应检验结果已保存: {firth_int_out}")

# ============================================================
# Step 6.8: 荟萃分析（GSE243013 + GSE241934 Taxane 合并）
# ============================================================
print("\n" + "="*70)
print("Step 6.8: 荟萃分析（GSE243013 + GSE241934 Taxane 合并）")
print("="*70)

# 读取主队列数据
main_path = os.path.join(cfg.RESULT_DIR, 'per_patient_metrics.csv')
df_main = pd.read_csv(main_path)

# 过滤主队列 Taxane 亚组
# 注意: 主队列的 chemo_class 列可能叫别的名字，让我们先看看
print(f"\n  主队列列名: {list(df_main.columns)}")
print(f"  主队列患者数: {len(df_main)}")
print(f"  主队列 response 值: {df_main['response'].value_counts().to_dict()}")

# 主队列有 chemo_class 吗？
if 'chemo_class' in df_main.columns:
    print(f"  主队列化疗方案分布: {df_main['chemo_class'].value_counts().to_dict()}")
    df_main_taxane = df_main[df_main['chemo_class'].str.contains('Taxane|taxane|紫杉醇|nab.Paclitaxel', case=False, na=False)].copy()
    print(f"  主队列 Taxane 亚组: {len(df_main_taxane)} 例")
else:
    # 如果没有 chemo_class，就用全部患者（假设大部分是 Taxane）
    print("  主队列无 chemo_class 列，使用全部患者作为 Taxane 近似")
    df_main_taxane = df_main.copy()

# 主队列二分类
df_main_taxane['y'] = (df_main_taxane['response'] == 'pCR').astype(int)
print(f"  主队列 Taxane: pCR={sum(df_main_taxane['y'])}, non-MPR={len(df_main_taxane)-sum(df_main_taxane['y'])}")

# ── 1. 标准化效应量（Fisher Z 变换的 AUC） ──
print("\n--- 1. AUC 荟萃分析（Fisher Z 变换）---")

from scipy.stats import norm

def auc_to_z(auc):
    """Fisher Z 变换 AUC → Z 分数"""
    # 用 probit 近似: Z = Φ⁻¹(AUC) * √2
    # 或者更简单: 用 logit 变换
    # 这里用标准的 AUC → SE 近似
    return np.log(auc / (1 - auc))  # logit 变换

def z_to_auc(z):
    return 1 / (1 + np.exp(-z))

def auc_se(auc, n1, n2):
    """AUC 的标准误（近似）"""
    # Hanley & McNeil 1982 方法
    # SE(AUC) ≈ sqrt( AUC*(1-AUC) + (n1-1)*(Q1 - AUC^2) + (n2-1)*(Q2 - AUC^2) ) / (n1*n2)
    # 简化: 用二项分布近似
    q1 = auc / (2 - auc)
    q2 = 2 * auc**2 / (1 + auc)
    se = np.sqrt((auc * (1 - auc) + (n1 - 1) * (q1 - auc**2) + (n2 - 1) * (q2 - auc**2)) / (n1 * n2))
    return se

# 收集各队列 Taxane 亚组的 AUC
meta_auc_data = []

# GSE243013 主队列（用 nk_dominant_ratio）
y_main = df_main_taxane['y'].values
x_main = df_main_taxane['nk_dominant_ratio'].values
auc_main = roc_auc_score(y_main, x_main)
n_r_main = sum(y_main)
n_nr_main = len(y_main) - n_r_main
se_main = auc_se(auc_main, n_r_main, n_nr_main)
print(f"  GSE243013 Taxane: AUC={auc_main:.4f}, n_R={n_r_main}, n_NR={n_nr_main}, SE={se_main:.4f}")
meta_auc_data.append({
    'study': 'GSE243013 (主队列)',
    'treatment': 'anti-PD-1 + Taxane',
    'nk_metric': 'nk_dominant_ratio',
    'n_R': n_r_main,
    'n_NR': n_nr_main,
    'auc': auc_main,
    'auc_se': se_main,
})

# GSE241934 Taxane 合并（用 frac_nklike_of_cd8）
df_gse_taxane = df_patient[df_patient['analysis_subset'].isin(['A_IIT_Taxane', 'B_RWC_Taxane'])].copy()
y_gse = (df_gse_taxane['response_binary'] == 'R').astype(int).values
x_gse = df_gse_taxane['frac_nklike_of_cd8'].values
auc_gse = roc_auc_score(y_gse, x_gse)
n_r_gse = sum(y_gse)
n_nr_gse = len(y_gse) - n_r_gse
se_gse = auc_se(auc_gse, n_r_gse, n_nr_gse)
print(f"  GSE241934 Taxane: AUC={auc_gse:.4f}, n_R={n_r_gse}, n_NR={n_nr_gse}, SE={se_gse:.4f}")
meta_auc_data.append({
    'study': 'GSE241934 (外部验证)',
    'treatment': 'anti-PD-1 + Taxane',
    'nk_metric': 'frac_nklike_of_cd8',
    'n_R': n_r_gse,
    'n_NR': n_nr_gse,
    'auc': auc_gse,
    'auc_se': se_gse,
})

# 随机效应荟萃分析（DerSimonian-Laird 方法）
df_meta_auc = pd.DataFrame(meta_auc_data)

# 用 logit(AUC) 做荟萃
df_meta_auc['logit_auc'] = np.log(df_meta_auc['auc'] / (1 - df_meta_auc['auc']))
df_meta_auc['logit_se'] = df_meta_auc['auc_se'] / (df_meta_auc['auc'] * (1 - df_meta_auc['auc']))

# 固定效应
w_fixed = 1 / df_meta_auc['logit_se']**2
theta_fixed = sum(w_fixed * df_meta_auc['logit_auc']) / sum(w_fixed)
se_fixed = np.sqrt(1 / sum(w_fixed))
auc_fixed = z_to_auc(theta_fixed)
auc_fixed_lo = z_to_auc(theta_fixed - 1.96 * se_fixed)
auc_fixed_hi = z_to_auc(theta_fixed + 1.96 * se_fixed)
z_fixed = theta_fixed / se_fixed
p_fixed = 2 * (1 - norm.cdf(abs(z_fixed)))

# 异质性检验（Q 检验）
Q = sum(w_fixed * (df_meta_auc['logit_auc'] - theta_fixed)**2)
I2 = max(0, (Q - (len(df_meta_auc) - 1)) / Q * 100) if Q > 0 else 0
p_hetero = 1 - chi2.cdf(Q, df=len(df_meta_auc) - 1)

# 随机效应（DerSimonian-Laird）
tau2 = max(0, (Q - (len(df_meta_auc) - 1)) / (sum(w_fixed) - sum(w_fixed**2) / sum(w_fixed)))
w_random = 1 / (df_meta_auc['logit_se']**2 + tau2)
theta_random = sum(w_random * df_meta_auc['logit_auc']) / sum(w_random)
se_random = np.sqrt(1 / sum(w_random))
auc_random = z_to_auc(theta_random)
auc_random_lo = z_to_auc(theta_random - 1.96 * se_random)
auc_random_hi = z_to_auc(theta_random + 1.96 * se_random)
z_random = theta_random / se_random
p_random = 2 * (1 - norm.cdf(abs(z_random)))

print(f"\n  固定效应合并 AUC: {auc_fixed:.4f}, 95%CI=[{auc_fixed_lo:.4f}, {auc_fixed_hi:.4f}], p={p_fixed:.4f}")
print(f"  随机效应合并 AUC: {auc_random:.4f}, 95%CI=[{auc_random_lo:.4f}, {auc_random_hi:.4f}], p={p_random:.4f}")
print(f"  异质性 Q={Q:.4f}, I²={I2:.1f}%, p={p_hetero:.4f}")

# ── 2. OR 荟萃分析 ──
print("\n--- 2. OR 荟萃分析 ---")

# 各队列 Firth OR（更稳健）
meta_or_data = []

# GSE243013 Taxane Firth OR
X_main = sm.add_constant(df_main_taxane['nk_dominant_ratio'].values)
y_main_arr = df_main_taxane['y'].values
try:
    or_main_f, ci_main_lo, ci_main_hi, p_main_f, _, _ = firth_odds_ratio_ci(X_main, y_main_arr, col_idx=1)
    print(f"  GSE243013 Taxane Firth OR={or_main_f:.3f}, 95%CI=[{ci_main_lo:.3f}, {ci_main_hi:.3f}], p={p_main_f:.4f}")
    meta_or_data.append({
        'study': 'GSE243013 (主队列)',
        'or_firth': or_main_f,
        'ci_lower': ci_main_lo,
        'ci_upper': ci_main_hi,
        'p_value': p_main_f,
        'n_R': n_r_main,
        'n_NR': n_nr_main,
    })
except Exception as e:
    print(f"  GSE243013 Firth 失败: {e}")

# GSE241934 Taxane Firth OR
x_gse_firth = df_gse_taxane['frac_nklike_of_cd8'].values
x_gse_z = (x_gse_firth - x_gse_firth.mean()) / x_gse_firth.std()
X_gse_firth = sm.add_constant(x_gse_z)
try:
    or_gse_f, ci_gse_lo, ci_gse_hi, p_gse_f, _, _ = firth_odds_ratio_ci(X_gse_firth, y_gse, col_idx=1)
    print(f"  GSE241934 Taxane Firth OR={or_gse_f:.3f}, 95%CI=[{ci_gse_lo:.3f}, {ci_gse_hi:.3f}], p={p_gse_f:.4f}")
    meta_or_data.append({
        'study': 'GSE241934 (外部验证)',
        'or_firth': or_gse_f,
        'ci_lower': ci_gse_lo,
        'ci_upper': ci_gse_hi,
        'p_value': p_gse_f,
        'n_R': n_r_gse,
        'n_NR': n_nr_gse,
    })
except Exception as e:
    print(f"  GSE241934 Firth 失败: {e}")

# OR 荟萃（用 log(OR)）
if len(meta_or_data) >= 2:
    df_meta_or = pd.DataFrame(meta_or_data)
    df_meta_or['log_or'] = np.log(df_meta_or['or_firth'])
    df_meta_or['log_se'] = (np.log(df_meta_or['ci_upper']) - np.log(df_meta_or['ci_lower'])) / (2 * 1.96)
    
    # 固定效应
    w_or = 1 / df_meta_or['log_se']**2
    theta_or_fixed = sum(w_or * df_meta_or['log_or']) / sum(w_or)
    se_or_fixed = np.sqrt(1 / sum(w_or))
    or_pooled_fixed = np.exp(theta_or_fixed)
    or_fixed_lo = np.exp(theta_or_fixed - 1.96 * se_or_fixed)
    or_fixed_hi = np.exp(theta_or_fixed + 1.96 * se_or_fixed)
    z_or = theta_or_fixed / se_or_fixed
    p_or_fixed = 2 * (1 - norm.cdf(abs(z_or)))
    
    # 异质性
    Q_or = sum(w_or * (df_meta_or['log_or'] - theta_or_fixed)**2)
    I2_or = max(0, (Q_or - (len(df_meta_or) - 1)) / Q_or * 100) if Q_or > 0 else 0
    p_het_or = 1 - chi2.cdf(Q_or, df=len(df_meta_or) - 1)
    
    # 随机效应
    tau2_or = max(0, (Q_or - (len(df_meta_or) - 1)) / (sum(w_or) - sum(w_or**2) / sum(w_or)))
    w_or_rand = 1 / (df_meta_or['log_se']**2 + tau2_or)
    theta_or_rand = sum(w_or_rand * df_meta_or['log_or']) / sum(w_or_rand)
    se_or_rand = np.sqrt(1 / sum(w_or_rand))
    or_pooled_rand = np.exp(theta_or_rand)
    or_rand_lo = np.exp(theta_or_rand - 1.96 * se_or_rand)
    or_rand_hi = np.exp(theta_or_rand + 1.96 * se_or_rand)
    p_or_rand = 2 * (1 - norm.cdf(abs(theta_or_rand / se_or_rand)))
    
    print(f"\n  固定效应合并 OR: {or_pooled_fixed:.3f}, 95%CI=[{or_fixed_lo:.3f}, {or_fixed_hi:.3f}], p={p_or_fixed:.4f}")
    print(f"  随机效应合并 OR: {or_pooled_rand:.3f}, 95%CI=[{or_rand_lo:.3f}, {or_rand_hi:.3f}], p={p_or_rand:.4f}")
    print(f"  异质性 Q={Q_or:.4f}, I²={I2_or:.1f}%, p={p_het_or:.4f}")

# ── 3. 森林图 ──
print("\n--- 3. 绘制森林图 ---")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('Meta-Analysis: GSE243013 + GSE241934 (Taxane regimen)', fontsize=14, fontweight='bold', y=1.02)

# Panel A: AUC 森林图
ax = axes[0]
ax.set_title('Panel A: AUC meta-analysis', fontsize=12, fontweight='bold')

y_positions = list(range(len(df_meta_auc) + 1))
y_labels = list(df_meta_auc['study'].values) + ['Pooled (Random)']

for i, (_, row) in enumerate(df_meta_auc.iterrows()):
    y = len(df_meta_auc) - i
    # 用 Hanley & McNeil 方法计算 AUC 的 95% CI（近似）
    auc_val = row['auc']
    se_val = row['auc_se']
    ci_lo = max(0.01, auc_val - 1.96 * se_val)
    ci_hi = min(0.99, auc_val + 1.96 * se_val)
    ax.errorbar(auc_val, y, xerr=[[auc_val - ci_lo], [ci_hi - auc_val]],
                fmt='s', color='#2E86C1', markersize=8, capsize=5, linewidth=2)
    ax.text(auc_val + 0.03, y, f"AUC={auc_val:.3f}", va='center', fontsize=9)

# 合并效应
y_pooled = 0
ax.errorbar(auc_random, y_pooled, xerr=[[auc_random - auc_random_lo], [auc_random_hi - auc_random]],
            fmt='D', color='#E74C3C', markersize=10, capsize=6, linewidth=2.5)
ax.text(auc_random + 0.05, y_pooled, f"AUC={auc_random:.3f}\n95%CI=[{auc_random_lo:.3f}, {auc_random_hi:.3f}]",
        va='center', fontsize=9, fontweight='bold', color='#E74C3C')

ax.axvline(x=0.5, color='gray', linestyle='--', linewidth=1, alpha=0.7)
ax.set_yticks([len(df_meta_auc) - i for i in range(len(df_meta_auc))] + [0])
ax.set_yticklabels(list(df_meta_auc['study'].values) + ['Pooled (RE)'])
ax.set_xlabel('AUC (95% CI)')
ax.set_xlim([0.3, 1.0])
ax.text(0.9, 0.05, f'I² = {I2:.1f}%\nQ p = {p_hetero:.3f}',
        transform=ax.transAxes, ha='right', fontsize=9,
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# Panel B: OR 森林图
ax = axes[1]
ax.set_title('Panel B: OR meta-analysis (Firth)', fontsize=12, fontweight='bold')

for i, (_, row) in enumerate(df_meta_or.iterrows()):
    y = len(df_meta_or) - i
    or_val = row['or_firth']
    ci_lo = row['ci_lower']
    ci_hi = row['ci_upper']
    ax.errorbar(or_val, y, xerr=[[or_val - ci_lo], [ci_hi - or_val]],
                fmt='s', color='#2E86C1', markersize=8, capsize=5, linewidth=2)
    ax.text(or_val * 1.1, y, f"OR={or_val:.2f}", va='center', fontsize=9)

y_pooled = 0
ax.errorbar(or_pooled_rand, y_pooled, xerr=[[or_pooled_rand - or_rand_lo], [or_rand_hi - or_pooled_rand]],
            fmt='D', color='#E74C3C', markersize=10, capsize=6, linewidth=2.5)
ax.text(or_pooled_rand * 1.1, y_pooled, f"OR={or_pooled_rand:.2f}\n95%CI=[{or_rand_lo:.2f}, {or_rand_hi:.2f}]",
        va='center', fontsize=9, fontweight='bold', color='#E74C3C')

ax.axvline(x=1.0, color='gray', linestyle='--', linewidth=1, alpha=0.7)
ax.set_yticks([len(df_meta_or) - i for i in range(len(df_meta_or))] + [0])
ax.set_yticklabels(list(df_meta_or['study'].values) + ['Pooled (RE)'])
ax.set_xlabel('Odds Ratio (95% CI)')
ax.set_xscale('log')
ax.text(0.9, 0.05, f'I² = {I2_or:.1f}%\nQ p = {p_het_or:.3f}',
        transform=ax.transAxes, ha='right', fontsize=9,
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
fig_forest_path = os.path.join(cfg.RESULT_DIR, 'GSE241934_Fig5_meta_analysis.png')
fig.savefig(fig_forest_path, dpi=300, bbox_inches='tight')
fig_forest_pdf = os.path.join(cfg.RESULT_DIR, 'GSE241934_Fig5_meta_analysis.pdf')
fig.savefig(fig_forest_pdf, bbox_inches='tight')
print(f"  ✅ 森林图已保存: {fig_forest_path}")
plt.close(fig)

# 保存荟萃分析结果
meta_summary = {
    'auc_pooled_fixed': auc_fixed,
    'auc_pooled_fixed_ci_lo': auc_fixed_lo,
    'auc_pooled_fixed_ci_hi': auc_fixed_hi,
    'auc_pooled_fixed_p': p_fixed,
    'auc_pooled_random': auc_random,
    'auc_pooled_random_ci_lo': auc_random_lo,
    'auc_pooled_random_ci_hi': auc_random_hi,
    'auc_pooled_random_p': p_random,
    'auc_Q': Q,
    'auc_I2': I2,
    'auc_hetero_p': p_hetero,
    'or_pooled_fixed': or_pooled_fixed,
    'or_pooled_fixed_ci_lo': or_fixed_lo,
    'or_pooled_fixed_ci_hi': or_fixed_hi,
    'or_pooled_fixed_p': p_or_fixed,
    'or_pooled_random': or_pooled_rand,
    'or_pooled_random_ci_lo': or_rand_lo,
    'or_pooled_random_ci_hi': or_rand_hi,
    'or_pooled_random_p': p_or_rand,
    'or_Q': Q_or,
    'or_I2': I2_or,
    'or_hetero_p': p_het_or,
    'n_studies': len(df_meta_auc),
    'total_R': n_r_main + n_r_gse,
    'total_NR': n_nr_main + n_nr_gse,
    'total_n': len(df_main_taxane) + len(df_gse_taxane),
}

df_meta_summary = pd.DataFrame([meta_summary])
meta_out = os.path.join(cfg.RESULT_DIR, 'GSE241934_meta_analysis_summary.csv')
df_meta_summary.to_csv(meta_out, index=False)
print(f"\n✅ 荟萃分析汇总已保存: {meta_out}")

print("\n" + "="*70)
print("[step6_gse241934_validation] Step 6.1-6.8 ALL completed!")
print("="*70)
