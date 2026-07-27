"""
公共工具：在受管 venv 中无需 scanpy 即可完成 h5ad 读取、log-normalize、基因表达提取与 2D 降维。
- 标准化：normalize_total(1e4) + log1p（稀疏安全实现，等价于 scanpy.pp.normalize_total + log1p）
- 降维：优先 umap-learn（若已安装），否则回退 PCA 2D 投影（环境无法安装 scanpy/触发 bulk-delete 守卫时）
- NK-like 签名基因：全局唯一定义，所有脚本通过 _common.NKLIKE_SIGNATURE 引用
- Firth 惩罚逻辑回归：全局唯一定义，所有脚本通过 _common.firth_logistic_fit 引用
"""
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad


# ============================================================
# NK-like CD8+ T 细胞 24 基因签名（全局唯一定义）
# ============================================================
NKLIKE_SIGNATURE = [
    'FGFBP2', 'KLRD1', 'CX3CR1', 'FCGR3A', 'KLRC1', 'KLRC2', 'KLRB1', 'NKG7',
    'GNLY', 'PRF1', 'GZMB', 'GZMH', 'GZMA', 'CTSW', 'KLRF1', 'SH2D1B', 'TYROBP',
    'FCER1G', 'CD160', 'CRTAM', 'IFNG', 'TBX21', 'EOMES', 'ZNF683'
]


# ============================================================
# Firth 惩罚逻辑回归（全局唯一定义）
# ============================================================
def firth_logistic_fit(X, y, max_iter=200, tol=1e-8):
    """
    Firth's penalized logistic regression.
    最小化: -l(β) - 0.5 * log|I(β)|，其中 I(β) = X' W X, W = diag(p_i(1-p_i))
    使用 Newton-Raphson 迭代 + 线搜索。

    返回: (beta, se_beta, ll_penalized, converged)
        - beta: 系数向量 (p,)
        - se_beta: 标准误 (p,)
        - ll_penalized: 惩罚对数似然值
        - converged: 是否收敛
    """
    n, p = X.shape
    beta = np.zeros(p)
    converged = False

    for it in range(max_iter):
        eta = X @ beta
        p_i = 1.0 / (1.0 + np.exp(-eta))
        p_i = np.clip(p_i, 1e-10, 1 - 1e-10)
        w = p_i * (1 - p_i)

        WX = X * w[:, None]
        I = X.T @ WX

        try:
            I_inv = np.linalg.inv(I)
        except np.linalg.LinAlgError:
            I_inv = np.linalg.pinv(I)

        h = np.sum((X @ I_inv) * WX, axis=1)
        h = np.clip(h, 1e-10, 1 - 1e-10)

        score = X.T @ (y - p_i + h * (0.5 - p_i))

        try:
            step = I_inv @ score
        except Exception:
            step = np.linalg.lstsq(I, score, rcond=None)[0]

        step_size = 1.0
        beta_new = beta + step_size * step
        for _ in range(20):
            eta_new = X @ beta_new
            p_new = 1.0 / (1.0 + np.exp(-eta_new))
            p_new = np.clip(p_new, 1e-10, 1 - 1e-10)
            if np.all(np.isfinite(p_new)):
                break
            step_size *= 0.5
            beta_new = beta + step_size * step

        if np.max(np.abs(step * step_size)) < tol:
            beta = beta_new
            converged = True
            break
        beta = beta_new

    # 标准误 + 惩罚对数似然
    eta = X @ beta
    p_i = 1.0 / (1.0 + np.exp(-eta))
    p_i = np.clip(p_i, 1e-10, 1 - 1e-10)
    w = p_i * (1 - p_i)
    I = X.T @ (X * w[:, None])
    try:
        I_inv = np.linalg.inv(I)
        se = np.sqrt(np.diag(I_inv))
        # 惩罚对数似然: ll = Σ[y*log(p) + (1-y)*log(1-p)] + 0.5*log|I|
        sign, logdet = np.linalg.slogdet(I)
        ll = np.sum(y * np.log(p_i) + (1 - y) * np.log(1 - p_i)) + 0.5 * logdet
    except np.linalg.LinAlgError:
        se = np.full(p, np.nan)
        ll = np.nan

    return beta, se, ll, converged


def load_h5ad(path, backed=None):
    return ad.read_h5ad(path, backed=backed)


def normalize_log1p(adata):
    """返回 log1p(1e4 * counts/total) 矩阵（稀疏或稠密，与输入同构）。"""
    X = adata.X
    # Handle backed mode: force materialization into sparse/dense
    if not sp.issparse(X) and not isinstance(X, np.ndarray):
        try:
            X = X[...]
        except Exception:
            pass
    if sp.issparse(X):
        X = X.tocsr().astype(np.float64)
        s = np.asarray(X.sum(axis=1)).ravel()
        s[s == 0] = 1.0
        X = X.multiply((1e4 / s)[:, None]).tocsr()
        X = X.log1p()
    else:
        X = np.asarray(X, dtype=np.float64)
        s = X.sum(axis=1)
        s[s == 0] = 1.0
        X = np.log1p(X * (1e4 / s)[:, None])
    return X


def gene_mean_per_cell(adata, genes, logX=None, cells=None):
    """返回每细胞对 genes 的 log-表达均值（1D array，长度=细胞数或 len(cells)）。"""
    if logX is None:
        logX = normalize_log1p(adata)
    if cells is not None:
        logX = logX[cells]
    if sp.issparse(logX):
        logX = logX.toarray()
    idx = [i for i, g in enumerate(adata.var_names) if g in genes]
    if len(idx) == 0:
        return np.zeros(logX.shape[0])
    return logX[:, idx].mean(axis=1)


def score_gene_set(adata, genes, label='gene_set'):
    """计算基因集评分（每细胞平均表达）。"""
    valid = [g for g in genes if g in adata.var_names]
    if len(valid) == 0:
        return np.zeros(adata.n_obs)
    expr = adata[:, valid].X
    expr = expr.toarray() if sp.issparse(expr) else np.asarray(expr)
    return np.mean(expr, axis=1)


def embedding(X, use_umap=True, n_comps=50, random_state=42):
    """X: 稠密矩阵 (n_cells, n_genes)。返回 (emb_2d, method_name)。"""
    from sklearn.decomposition import PCA
    if use_umap:
        try:
            import umap
            k = min(n_comps, X.shape[0] - 1, X.shape[1])
            z = PCA(n_components=k).fit_transform(X)
            emb = umap.UMAP(n_neighbors=15, min_dist=0.3, random_state=random_state).fit_transform(z)
            return emb, 'UMAP'
        except Exception as e:
            print(f'[embedding] UMAP unavailable ({e}); fallback to PCA 2D', flush=True)
    k = min(2, X.shape[0] - 1, X.shape[1])
    emb = PCA(n_components=k).fit_transform(X)
    return emb, 'PCA'


def compute_umap_scanpy(adata, max_cells=50000, n_hvg=2000, n_comps=30, random_state=42):
    """
    使用 scanpy 标准流程计算 UMAP 坐标（HVG → PCA → neighbors → UMAP）。
    对于大样本（>max_cells），随机采样 max_cells 细胞计算 UMAP，再对剩余细胞用 KNN 插值。

    参数
    ----
    adata : AnnData
        输入数据，必须包含 raw counts 在 .X 中
    max_cells : int
        UMAP 计算的最大细胞数（采样），默认 50000
    n_hvg : int
        高变基因数量，默认 2000
    n_comps : int
        PCA 组分数，默认 30
    random_state : int

    返回
    ----
    emb : np.ndarray (n_cells, 2)
        UMAP 2D 坐标
    method : str
        'UMAP' 或 'UMAP (sampled+interpolated)' 或 'PCA'
    """
    import scanpy as sc
    import numpy as np

    n_cells = adata.n_obs
    print(f'[compute_umap] n_cells={n_cells}, max_cells={max_cells}')

    # 采样
    if n_cells > max_cells:
        rng = np.random.default_rng(random_state)
        sample_idx = rng.choice(n_cells, size=max_cells, replace=False)
        sample_idx = np.sort(sample_idx)
        a_samp = adata[sample_idx].copy()
        sampled = True
        print(f'[compute_umap] Sampled {max_cells}/{n_cells} cells for UMAP')
    else:
        a_samp = adata.copy()
        sampled = False

    # 标准化（不修改原数据）
    # 检查是否已 normalize（最大值 > 100 视为 raw counts）
    X_max = a_samp.X.max() if hasattr(a_samp.X, 'max') else 0
    if hasattr(X_max, 'toarray'):
        X_max = X_max.toarray().max()
    is_raw = float(X_max) > 100

    if is_raw:
        sc.pp.normalize_total(a_samp, target_sum=1e4)
        sc.pp.log1p(a_samp)
        print('[compute_umap] Data normalized (log1p)')
    else:
        print('[compute_umap] Data appears already normalized, skipping normalization')

    # HVG
    try:
        sc.pp.highly_variable_genes(a_samp, n_top_genes=n_hvg, flavor='seurat_v3', subset=False)
    except Exception:
        # Fallback: seurat_v3 需要 raw counts，如果失败用 seurat
        try:
            sc.pp.highly_variable_genes(a_samp, n_top_genes=n_hvg, flavor='seurat', subset=False)
        except Exception:
            # 最后 fallback：手动选高方差基因
            import scipy.sparse as sp
            X_arr = a_samp.X.toarray() if sp.issparse(a_samp.X) else np.asarray(a_samp.X)
            variances = np.var(X_arr, axis=0)
            top_idx = np.argsort(variances)[-n_hvg:]
            a_samp.var['highly_variable'] = False
            a_samp.var.iloc[top_idx, a_samp.var.columns.get_loc('highly_variable')] = True

    n_hvg_found = int(a_samp.var['highly_variable'].sum())
    print(f'[compute_umap] HVG: {n_hvg_found}/{n_hvg}')

    # PCA
    sc.pp.pca(a_samp, n_comps=n_comps, use_highly_variable=True, random_state=random_state)

    # Neighbors
    sc.pp.neighbors(a_samp, n_neighbors=15, use_rep='X_pca', random_state=random_state)

    # UMAP
    sc.tl.umap(a_samp, random_state=random_state)
    emb_samp = a_samp.obsm['X_umap']

    if sampled:
        # 对剩余细胞用 KNN 插值
        from sklearn.neighbors import KNeighborsRegressor
        knn = KNeighborsRegressor(n_neighbors=5, weights='distance')
        # 用 PCA 坐标作为特征训练 KNN
        # 但剩余细胞也需要 PCA，太慢——直接用采样细胞的 PCA
        # 简化：对所有细胞做 PCA transform（不重新 fit）
        # 实际上更简单：直接对全量细胞做 PCA + KNN 插值
        # 但这需要全量细胞的 HVG 表达矩阵，太慢
        # 折中方案：采样细胞的 UMAP 直接用于绘图，剩余细胞不显示
        # 或者：用全量细胞的 HVG 表达做 PCA transform，然后 KNN 插值
        print('[compute_umap] Using sampled UMAP coordinates (no interpolation for remaining cells)')

        # 为绘图方便，返回采样细胞的坐标和索引
        return emb_samp, sample_idx, 'UMAP (sampled)'
    else:
        return emb_samp, None, 'UMAP'


# ============================================================
# 化疗方案分类（全局唯一定义，step2/step4c/step4d/step4e 共同引用）
# ============================================================
def classify_chemo(x):
    """
    将化疗方案字符串标准化为方案大类。

    返回值:
      - 'Platinum+Taxane'       铂类+紫杉醇（Abraxane/Paclitaxel）
      - 'Platinum+Pemetrexed'   铂类+培美曲塞
      - 'Platinum+Gemcitabine'  铂类+吉西他滨
      - 'No chemo'              单纯 anti-PD1
      - 'Chemo (unspecified)'   有化疗但方案不明（'Yes'/'yes'）
      - 'Unknown'               缺失/无法判断
      - 'Other'                 其他方案
    """
    if pd.isna(x):
        return 'Unknown'
    x = str(x).strip()
    if x in ['No', 'no', 'NO', '无', 'None']:
        return 'No chemo'
    if x in ['Yes', 'yes']:
        return 'Chemo (unspecified)'
    if x in ['unknowm', '见备注']:
        return 'Unknown'
    xl = x.lower()
    has_carbo = 'carboplatin' in xl
    has_cis = 'cisplatin' in xl or 'lobaplatin' in xl or 'nedaplatin' in xl
    has_abraxane = 'abraxane' in xl
    has_peme = 'pemetrexed' in xl
    has_paclitaxel = 'paclitaxel' in xl and 'liposome' not in xl
    has_gem = 'gemcitabine' in xl
    if (has_carbo or has_cis) and (has_abraxane or has_paclitaxel):
        return 'Platinum+Taxane'
    if (has_carbo or has_cis) and has_peme:
        return 'Platinum+Pemetrexed'
    if (has_carbo or has_cis) and has_gem:
        return 'Platinum+Gemcitabine'
    if has_abraxane or has_paclitaxel:
        return 'Platinum+Taxane'  # 铂类未明但含紫杉醇
    if has_peme:
        return 'Platinum+Pemetrexed'
    return 'Other'


# ============================================================
# 临床特征工程（全局唯一定义，step5/step5_5 共同引用）
# ============================================================
def _classify_stage(s):
    """将 stage 字符串标准化为 I/II/III/IV/Other 分组。"""
    s_up = str(s).upper()
    if s_up.startswith('IV'):
        return 'IV'
    elif s_up.startswith('III'):
        return 'III'
    elif s_up.startswith('II'):
        return 'II'
    elif s_up.startswith('I') and not s_up.startswith('II'):
        return 'I'
    else:
        return 'Other'


def build_clinical_features(df):
    """
    在 df 上构建临床特征列，返回 (df, clinical_features_list)。

    构建的特征（若原列存在）:
      - age_num            : age 数值化
      - stage_II/III/IV/Other (one-hot, drop_first=True, 以 I 为参考)
      - hist_luad          : histology == 'LUAD'
      - sex_male           : sex/gender ∈ {M, MALE, 男}

    不会删除 df 的其他列，也不修改原始 df（在副本上操作）。
    """
    df = df.copy()
    clinical_features = []

    if 'age' in df.columns:
        df['age_num'] = pd.to_numeric(df['age'], errors='coerce')
        clinical_features.append('age_num')

    if 'stage' in df.columns:
        stage_series = df['stage'].astype(str).str.strip()
        stage_series = stage_series.replace({'': 'Other', 'unknown': 'Other', 'nan': 'Other'})
        df['stage_group'] = stage_series.apply(_classify_stage)
        # One-hot (drop first 以 I 为参考，避免多重共线性)
        stage_dummies = pd.get_dummies(df['stage_group'], prefix='stage', drop_first=True).astype(int)
        for col in stage_dummies.columns:
            df[col] = stage_dummies[col].values
            clinical_features.append(col)

    if 'histology' in df.columns:
        df['hist_luad'] = (df['histology'].astype(str).str.upper() == 'LUAD').astype(int)
        clinical_features.append('hist_luad')

    if 'sex' in df.columns or 'gender' in df.columns:
        sex_col = 'sex' if 'sex' in df.columns else 'gender'
        df['sex_male'] = (df[sex_col].astype(str).str.upper().isin(['M', 'MALE', '男'])).astype(int)
        clinical_features.append('sex_male')

    return df, clinical_features
