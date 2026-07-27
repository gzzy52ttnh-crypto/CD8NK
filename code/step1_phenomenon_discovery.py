"""
step1_phenomenon_discovery.py
复现 Table1 (baseline) + Fig1 (overview) + FigS1 (supplement).
"""
import config as cfg
import os
import sys
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import mannwhitneyu, fisher_exact
from scipy.stats import chi2_contingency

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import compute_umap_scanpy

print("[step1_phenomenon_discovery] Starting...")

os.makedirs(cfg.RESULT_DIR, exist_ok=True)

# ── Load data ──
adata = sc.read_h5ad(cfg.H5AD_T)
print(f"Loaded T cells: {adata.shape}")

# ── Compute UMAP if missing ──
if 'X_umap' not in adata.obsm or adata.obsm['X_umap'].shape[0] != adata.n_obs:
    print("[step1] X_umap missing in obsm, computing via scanpy pipeline (sampled)...")
    try:
        emb_samp, sample_idx, method = compute_umap_scanpy(adata, max_cells=50000, n_hvg=2000)
        # 全量细胞 UMAP（采样细胞有坐标，其他用 NaN）
        full_emb = np.full((adata.n_obs, 2), np.nan)
        full_emb[sample_idx] = emb_samp
        adata.obsm['X_umap'] = full_emb
        print(f"[step1] UMAP computed ({method}): {emb_samp.shape[0]} sampled cells")
    except Exception as e:
        print(f"[step1] UMAP computation failed: {e}")
        adata.obsm['X_umap'] = np.full((adata.n_obs, 2), np.nan)

# ── Pre-processing ──
# Filter to pCR vs non-MPR for baseline
obs = adata.obs.copy()
obs['patient'] = obs[cfg.COL_SAMPLE]
# Fix typo in original data: 'unknowm' → 'unknown'
if cfg.COL_SMOKING in obs.columns:
    obs[cfg.COL_SMOKING] = obs[cfg.COL_SMOKING].astype(str).replace({'unknowm': 'unknown'})
valid_resp = obs[cfg.COL_RESPONSE].isin(['pCR', 'non-MPR'])
obs_table = obs[valid_resp].drop_duplicates(subset='patient').copy()

# ============================================================
# TABLE 1: Baseline characteristics
# ============================================================
table1_rows = []

# Demographics and clinical vars
cat_vars = {
    'sex': cfg.COL_SEX,
    'stage': cfg.COL_STAGE,
    'histology': cfg.COL_HISTOLOGY,
    'smoking': cfg.COL_SMOKING,
}
cont_vars = {
    'age': cfg.COL_AGE,
}

pcr_mask = obs_table[cfg.COL_RESPONSE] == 'pCR'
nmp_mask = obs_table[cfg.COL_RESPONSE] == 'non-MPR'
n_pcr = pcr_mask.sum()
n_nmp = nmp_mask.sum()

for display_name, col in cat_vars.items():
    if col not in obs_table.columns:
        continue
    all_cats = obs_table[col].astype(str).value_counts()
    for cat_val in all_cats.index:
        row = {'Variable': f'{display_name} = {cat_val}'}
        n_pcr_cat = (obs_table.loc[pcr_mask, col].astype(str) == cat_val).sum()
        n_nmp_cat = (obs_table.loc[nmp_mask, col].astype(str) == cat_val).sum()
        row['pCR (n=%d)' % n_pcr] = f'{n_pcr_cat} ({100*n_pcr_cat/n_pcr:.1f}%)' if n_pcr > 0 else '0'
        row['non-MPR (n=%d)' % n_nmp] = f'{n_nmp_cat} ({100*n_nmp_cat/n_nmp:.1f}%)' if n_nmp > 0 else '0'
        # Fisher exact for 2x2
        if len(all_cats) == 2:
            a = n_pcr_cat; b = n_pcr - n_pcr_cat
            c = n_nmp_cat; d = n_nmp - n_nmp_cat
            _, fp = fisher_exact([[a,b],[c,d]])
            row['p_value'] = round(fp, 4)
        else:
            row['p_value'] = ''
        row['test'] = 'Fisher' if len(all_cats) == 2 else ''
        table1_rows.append(row)

for display_name, col in cont_vars.items():
    if col not in obs_table.columns:
        continue
    row = {'Variable': display_name}
    v_pcr = obs_table.loc[pcr_mask, col].dropna().astype(float)
    v_nmp = obs_table.loc[nmp_mask, col].dropna().astype(float)
    row['pCR (n=%d)' % n_pcr] = f'{np.median(v_pcr):.1f} [{np.percentile(v_pcr,25):.1f}-{np.percentile(v_pcr,75):.1f}]'
    row['non-MPR (n=%d)' % n_nmp] = f'{np.median(v_nmp):.1f} [{np.percentile(v_nmp,25):.1f}-{np.percentile(v_nmp,75):.1f}]'
    _, mwp = mannwhitneyu(v_pcr, v_nmp, alternative='two-sided')
    row['p_value'] = round(mwp, 4)
    row['test'] = 'MWU'
    table1_rows.append(row)

df_table1 = pd.DataFrame(table1_rows)

# BH FDR 多重比较校正
from statsmodels.stats.multitest import multipletests
pvals_raw = pd.to_numeric(df_table1['p_value'], errors='coerce').values
mask_valid = ~np.isnan(pvals_raw)
pvals_adj = np.full(len(pvals_raw), np.nan)
if mask_valid.sum() > 0:
    _, pvals_bh, _, _ = multipletests(pvals_raw[mask_valid], method='fdr_bh')
    pvals_adj[mask_valid] = pvals_bh
df_table1['p_value_BH'] = [round(x, 4) if not np.isnan(x) else '' for x in pvals_adj]

df_table1.to_csv(cfg.result_path('Table1_baseline.csv'), index=False)
print("Table1_baseline.csv written (with BH FDR correction)")

# ============================================================
# FIGURE 1: Overview (4 Panels)
# ============================================================
fig1, axes = plt.subplots(2, 2, figsize=(14, 12))
axA, axB, axC, axD = axes[0,0], axes[0,1], axes[1,0], axes[1,1]

# Panel A: UMAP by cell_type (sub_cell_type)
cell_types = adata.obs[cfg.COL_CELL_TYPE]
# Focus on major CD8 subtypes
cd8_types = [ct for ct in cell_types.unique() if 'CD8T' in str(ct)]
colors_cd8 = plt.cm.tab20(np.linspace(0, 1, len(cd8_types)))
# Show only CD8T cells in UMAP, grey others
is_cd8 = adata.obs[cfg.COL_CELL_TYPE].isin(cd8_types)
others = adata[~is_cd8].copy()
cd8 = adata[is_cd8].copy()
umap_arr = adata.obsm.get('X_umap', None)
if umap_arr is not None and not np.all(np.isnan(umap_arr)):
    # 采样细胞有 UMAP 坐标，其他为 NaN——只画有效坐标的细胞
    valid_others = ~np.isnan(umap_arr[~is_cd8.values, 0])
    others_umap = umap_arr[~is_cd8.values][valid_others]
    if len(others_umap) > 0:
        axA.scatter(others_umap[:, 0], others_umap[:, 1],
                    c='lightgrey', s=0.5, alpha=0.3, rasterized=True)
    cd8_umap = umap_arr[is_cd8.values]
    cd8_valid = ~np.isnan(cd8_umap[:, 0])
    for i, ct in enumerate(cd8_types):
        ct_mask = (cd8.obs[cfg.COL_CELL_TYPE] == ct).values & cd8_valid
        if ct_mask.sum() > 0:
            axA.scatter(cd8_umap[ct_mask, 0], cd8_umap[ct_mask, 1],
                        c=[colors_cd8[i]], s=0.5, label=ct, rasterized=True)
    axA.legend(fontsize=6, markerscale=4, loc='upper right', bbox_to_anchor=(1.5, 1))
    n_plotted = int((~np.isnan(umap_arr[:, 0])).sum())
    axA.set_title(f'Panel A: CD8+ T cell UMAP by sub_cell_type\n(n={n_plotted} sampled cells)')
else:
    axA.text(0.5, 0.5, 'UMAP coordinates not available',
             ha='center', va='center', transform=axA.transAxes, fontsize=12, color='grey')
    axA.set_title('Panel A: CD8+ T cell UMAP (unavailable)')
axA.set_xlabel('UMAP1'); axA.set_ylabel('UMAP2')

# Panel B: Pie chart of CD8+ subtypes proportion
cd8_counts = cd8.obs[cfg.COL_CELL_TYPE].value_counts()
pie_colors = [colors_cd8[cd8_types.index(ct)] for ct in cd8_counts.index]
axB.pie(cd8_counts.values, labels=cd8_counts.index, autopct='%1.1f%%',
        colors=pie_colors, textprops={'fontsize': 7})
axB.set_title('Panel B: CD8+ subtype proportions')

# Panel C: NK-like proportion pCR vs non-MPR
nk_like_mask = adata.obs[cfg.COL_CELL_TYPE] == 'CD8T_NK-like_FGFBP2'
patients = adata.obs[cfg.COL_SAMPLE].unique()
patient_nk_ratio = {}
patient_has_data = set()
for p in patients:
    p_mask = adata.obs[cfg.COL_SAMPLE] == p
    total = p_mask.sum()
    nk_count = (p_mask & nk_like_mask).sum()
    if total > 0:
        patient_nk_ratio[p] = nk_count / total
        patient_has_data.add(p)

pcr_patients = obs_table.loc[pcr_mask, 'patient'].values
nmp_patients = obs_table.loc[nmp_mask, 'patient'].values

nk_pcr = [patient_nk_ratio[p] for p in pcr_patients if p in patient_nk_ratio]
nk_nmp = [patient_nk_ratio[p] for p in nmp_patients if p in patient_nk_ratio]
n_pcr_excl = len(pcr_patients) - len(nk_pcr)
n_nmp_excl = len(nmp_patients) - len(nk_nmp)
if n_pcr_excl > 0 or n_nmp_excl > 0:
    print(f"  [Panel C] Excluded patients with no CD8T cells: pCR={n_pcr_excl}, non-MPR={n_nmp_excl}")

bp = axC.boxplot([nk_pcr, nk_nmp], patch_artist=True)
axC.set_xticklabels(['pCR', 'non-MPR'])
bp['boxes'][0].set_facecolor('#4CAF50'); bp['boxes'][1].set_facecolor('#F44336')
_, mwp_nk = mannwhitneyu(nk_pcr, nk_nmp, alternative='two-sided')
axC.set_title(f'Panel C: NK-like ratio (pCR vs non-MPR)\nMWU p={mwp_nk:.4f}')
axC.set_ylabel('NK-like CD8+ proportion')

# Panel D: Cell count per patient histogram
patient_cell_counts = adata.obs[cfg.COL_SAMPLE].value_counts().values
axD.hist(patient_cell_counts, bins=30, color='steelblue', edgecolor='white')
axD.set_title('Panel D: Cell count per patient')
axD.set_xlabel('Number of cells'); axD.set_ylabel('Number of patients')

plt.tight_layout()
fig1.savefig(cfg.result_path('Fig1_overview.pdf'), dpi=150, bbox_inches='tight')
fig1.savefig(cfg.result_path('Fig1_overview.png'), dpi=150, bbox_inches='tight')
plt.close(fig1)
print("Fig1_overview saved")

# ============================================================
# FIGURE S1: Supplement (3-4 Panels)
# ============================================================
figS1, axesS = plt.subplots(2, 2, figsize=(16, 14))

# Panel A: CD8 subtypes boxplot matrix (pCR vs non-MPR)
axSA = axesS[0,0]
subtype_data = {}
for ct in cd8_types:
    ct_mask = adata.obs[cfg.COL_CELL_TYPE] == ct
    ratios_pcr = []
    ratios_nmp = []
    for p in pcr_patients:
        p_mask = adata.obs[cfg.COL_SAMPLE] == p
        total = p_mask.sum()
        ct_count = (p_mask & ct_mask).sum()
        ratios_pcr.append(ct_count / total if total > 0 else 0)
    for p in nmp_patients:
        p_mask = adata.obs[cfg.COL_SAMPLE] == p
        total = p_mask.sum()
        ct_count = (p_mask & ct_mask).sum()
        ratios_nmp.append(ct_count / total if total > 0 else 0)
    subtype_data[ct] = (ratios_pcr, ratios_nmp)

positions = []
labels_sa = []
all_data_sa = []
for i, ct in enumerate(cd8_types):
    rp, rn = subtype_data[ct]
    positions.extend([i*3 + 1, i*3 + 2])
    all_data_sa.extend([rp, rn])
    labels_sa.extend([f'{ct}\npCR', f'{ct}\nnon'])

bp_sa = axSA.boxplot(all_data_sa, positions=positions, patch_artist=True, widths=0.7)
for i, ct in enumerate(cd8_types):
    bp_sa['boxes'][i*2].set_facecolor('#4CAF50')
    bp_sa['boxes'][i*2+1].set_facecolor('#F44336')
axSA.set_xticks([i*3 + 1.5 for i in range(len(cd8_types))])
axSA.set_xticklabels(cd8_types, rotation=45, ha='right', fontsize=7)
axSA.set_title('Panel A: CD8 subtype proportions pCR vs non-MPR')
axSA.set_ylabel('Proportion')

# Panel B: QC metrics (UMI count, gene count, mito%)
axSB1 = axesS[0,1]
umi_counts = adata.obs['total_counts'].values
axSB1.hist(np.log10(umi_counts + 1), bins=50, color='steelblue', edgecolor='white', alpha=0.7)
axSB1.set_title('Panel B: log10(UMI count) distribution')
axSB1.set_xlabel('log10(UMI+1)')

axSB2 = axesS[1,0]
gene_counts = adata.obs['n_genes_by_counts'].values
axSB2.hist(np.log10(gene_counts + 1), bins=50, color='darkorange', edgecolor='white', alpha=0.7)
axSB2.set_title('Panel B: log10(gene count) distribution')
axSB2.set_xlabel('log10(genes+1)')

axSB3 = axesS[1,1]
mito_pct = adata.obs['pct_counts_mt'].values
axSB3.hist(mito_pct, bins=50, color='darkred', edgecolor='white', alpha=0.7)
axSB3.set_title('Panel B: Mito% distribution')
axSB3.set_xlabel('Mito%')

# NOTE: Panel C (marker gene feature plots) and Panel D (Shannon entropy)
# require further data not currently accessible in single-pass;
# documented as pending in log.

figS1.tight_layout()
figS1.savefig(cfg.result_path('FigS1_supplement.pdf'), dpi=150, bbox_inches='tight')
figS1.savefig(cfg.result_path('FigS1_supplement.png'), dpi=150, bbox_inches='tight')
plt.close(figS1)
print("FigS1_supplement saved")

print("[step1_phenomenon_discovery] Done.")
