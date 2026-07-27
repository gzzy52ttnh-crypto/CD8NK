import os

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULT_DIR = os.path.join(ROOT_DIR, 'result')
ADATA_DIR = os.path.join(ROOT_DIR, 'adata')
# 原始输入数据（上游 pipeline 产物，通过符号链接指向 adata/）
H5AD_T = os.path.join(ROOT_DIR, 'GSE243013_T_cells.h5ad')
H5AD_IMMUNE = os.path.join(ROOT_DIR, 'GSE243013_immune.h5ad')
# 过程数据（step3 运行时从 immune.h5ad 拆分生成，属于中间产物，放 result 目录）
H5AD_NK = os.path.join(RESULT_DIR, 'GSE243013_NK_like_CD8.h5ad')
H5AD_MYELOID = os.path.join(RESULT_DIR, 'GSE243013_myeloid.h5ad')


def result_path(filename):
    return os.path.join(RESULT_DIR, filename)


# Obs column name mapping (actual names in h5ad)
COL_SAMPLE = 'sampleID'
COL_CELL_TYPE = 'sub_cell_type'
COL_RESPONSE = 'pathological_response'
COL_RESPONSE_RATE = 'pathological_response_rate'
COL_CLONOTYPE = 'clonotype'
COL_CLONOTYPE_NUM = 'clonotype_number'
COL_EXPANSION = 'expansion'
COL_AGE = 'age'
COL_SEX = 'gender'
COL_STAGE = 'pre_treatment_staging'
COL_HISTOLOGY = 'cancer_type'
COL_SMOKING = 'smoking_history'


# ============================================================
# 分析阈值（全局唯一定义，step2/step4d/step6c/step6e 等共同引用）
# ============================================================
# 大克隆阈值：clone_size >= BIG_CLONE_THRESHOLD 视为"大克隆"
# 选择依据见 step2_core_index_construction.py 阈值敏感性分析（threshold=5 平衡稳定性与统计效力）
BIG_CLONE_THRESHOLD = 5

# NK-dominant 克隆判定阈值：nk_ratio >= NK_DOMINANT_RATIO_THRESHOLD 视为 NK-dominant
# 同时用于 Tex-dominant (tex_ratio >= 0.5) 和 Mixed (nk+tex >= 0.5) 判定
NK_DOMINANT_RATIO_THRESHOLD = 0.5
