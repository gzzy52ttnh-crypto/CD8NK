"""
公共工具：在受管 venv 中无需 scanpy 即可完成 h5ad 读取、log-normalize、基因表达提取与 2D 降维。
- 标准化：normalize_total(1e4) + log1p（稀疏安全实现，等价于 scanpy.pp.normalize_total + log1p）
- 降维：优先 umap-learn（若已安装），否则回退 PCA 2D 投影（环境无法安装 scanpy/触发 bulk-delete 守卫时）
"""
import os
import numpy as np
import scipy.sparse as sp
import anndata as ad


def load_h5ad(path, backed=None):
    return ad.read_h5ad(path, backed=backed)


def normalize_log1p(adata):
    """返回 log1p(1e4 * counts/total) 矩阵（稀疏或稠密，与输入同构）。"""
    X = adata.X
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
