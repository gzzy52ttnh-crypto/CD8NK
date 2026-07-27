"""
step3_mechanism_exploration.py
复现 Fig3 (myeloid) + Fig4 (mechanism) + FigS5 (interaction network).
Extracts NK-like CD8 and myeloid subsets from immune.h5ad on first run.
"""
import config as cfg
from _common import normalize_log1p, load_h5ad, compute_umap_scanpy
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, spearmanr
import os
import warnings
warnings.filterwarnings('ignore')

print("[step3_mechanism_exploration] Starting...")

os.makedirs(cfg.RESULT_DIR, exist_ok=True)

# ── Extract / load NK-like CD8 and Myeloid h5ad ──
immune = sc.read_h5ad(cfg.H5AD_IMMUNE)
print(f"Loaded immune: {immune.shape}")

# Extract myeloid
myeloid = immune[immune.obs['major_cell_type'] == 'Myeloid cell'].copy()
print(f"Myeloid subset: {myeloid.shape}")
myeloid.write(cfg.H5AD_MYELOID)
print(f"Written {cfg.H5AD_MYELOID}")

# Extract NK-like CD8 (all CD8T subtypes)
cd8_subtypes = [s for s in immune.obs[cfg.COL_CELL_TYPE].unique() if 'CD8T' in str(s)]
nk_cd8 = immune[immune.obs[cfg.COL_CELL_TYPE].isin(cd8_subtypes)].copy()
print(f"NK-like CD8 subset: {nk_cd8.shape}")
nk_cd8.write(cfg.H5AD_NK)
print(f"Written {cfg.H5AD_NK}")

# ── Compute UMAP for myeloid if missing ──
if 'X_umap' not in myeloid.obsm or myeloid.obsm['X_umap'].shape[0] != myeloid.n_obs:
    print("[step3] X_umap missing in myeloid.obsm, computing via scanpy pipeline (sampled)...")
    try:
        emb_samp, sample_idx, method = compute_umap_scanpy(myeloid, max_cells=50000, n_hvg=2000)
        full_emb = np.full((myeloid.n_obs, 2), np.nan)
        full_emb[sample_idx] = emb_samp
        myeloid.obsm['X_umap'] = full_emb
        print(f"[step3] Myeloid UMAP computed ({method}): {emb_samp.shape[0]} sampled cells")
    except Exception as e:
        print(f"[step3] Myeloid UMAP computation failed: {e}")
        myeloid.obsm['X_umap'] = np.full((myeloid.n_obs, 2), np.nan)

# ── Load per-patient metrics ──
df_patient = pd.read_csv(cfg.result_path('per_patient_metrics.csv'))

# ============================================================
# FIGURE 3: Myeloid (4 Panels)
# ============================================================
fig3, axes3 = plt.subplots(2, 2, figsize=(14, 12))

# Panel A: Myeloid UMAP
umap_mye = myeloid.obsm.get('X_umap', None)
if umap_mye is not None and not np.all(np.isnan(umap_mye)):
    subtypes_mye = myeloid.obs[cfg.COL_CELL_TYPE].unique()
    colors_mye = plt.cm.tab20(np.linspace(0, 1, len(subtypes_mye)))
    axA3 = axes3[0, 0]
    valid_global = ~np.isnan(umap_mye[:, 0])
    for i, st in enumerate(subtypes_mye):
        mask = (myeloid.obs[cfg.COL_CELL_TYPE] == st).values & valid_global
        if mask.sum() > 0:
            axA3.scatter(umap_mye[mask, 0], umap_mye[mask, 1],
                         c=[colors_mye[i]], s=0.5, label=st, rasterized=True)
    axA3.legend(fontsize=5, markerscale=3, loc='upper right', bbox_to_anchor=(1.5, 1))
    n_plotted_mye = int(valid_global.sum())
    axA3.set_title(f'Panel A: Myeloid UMAP\n(n={n_plotted_mye} sampled cells)')
    axA3.set_xlabel('UMAP1'); axA3.set_ylabel('UMAP2')
else:
    axes3[0, 0].text(0.5, 0.5, 'No UMAP coordinates', ha='center', va='center', transform=axes3[0,0].transAxes)
    axes3[0, 0].set_title('Panel A: Myeloid UMAP (unavailable)')

# Panel B: Stacked bar of myeloid proportions per patient
mye_obs = myeloid.obs
patients_m = mye_obs[cfg.COL_SAMPLE].unique()
mye_subtypes_all = mye_obs[cfg.COL_CELL_TYPE].value_counts().index[:8]  # top 8

patient_mye_props = []
patient_mye_resp = []
for p in patients_m:
    p_data = mye_obs[mye_obs[cfg.COL_SAMPLE] == p]
    total = len(p_data)
    props = {}
    for st in mye_subtypes_all:
        props[st] = (p_data[cfg.COL_CELL_TYPE] == st).sum() / total if total > 0 else 0
    patient_mye_props.append(props)
    patient_mye_resp.append(p_data[cfg.COL_RESPONSE].iloc[0] if len(p_data) > 0 else 'unknown')

# Sort patients by response
sorted_idx = np.argsort(patient_mye_resp)
sorted_props = [patient_mye_props[i] for i in sorted_idx]
sorted_resp = [patient_mye_resp[i] for i in sorted_idx]

axB3 = axes3[0, 1]
bottom = np.zeros(len(sorted_props))
for i, st in enumerate(mye_subtypes_all):
    vals = [p[st] for p in sorted_props]
    axB3.bar(range(len(vals)), vals, bottom=bottom, label=st, width=0.8)
    bottom += vals
axB3.set_title('Panel B: Myeloid proportions per patient')
axB3.set_xlabel('Patients (sorted by response)')
axB3.legend(fontsize=5, loc='upper right', bbox_to_anchor=(1.45, 1))

# Panel C: SPP1+ TAM ratio pCR vs non-MPR
# ── STRICT GENE-EXPRESSION DEFINITION: SPP1 > 0 in Mac/Mono cells, denominator = Mac/Mono ──
spp1_tam_ratios = {}
mac_mono_counts = {}

# Use myeloid h5ad as primary source for both SPP1 expression and cell type annotation
adata_mac = myeloid
mac_mono_mask = adata_mac.obs[cfg.COL_CELL_TYPE].str.contains('Mφ|Mac|Mono', regex=True, na=False).values

if 'SPP1' not in adata_mac.var_names:
    raise ValueError("SPP1 gene not found in myeloid h5ad var_names — cannot compute SPP1+ TAM ratio")

spp1_idx = list(adata_mac.var_names).index('SPP1')
spp1_expr_all = adata_mac.X[:, spp1_idx]
spp1_expr_all = spp1_expr_all.toarray().flatten() if hasattr(spp1_expr_all, 'toarray') else np.asarray(spp1_expr_all).flatten()
is_spp1_pos = spp1_expr_all > 0  # raw count > 0 = SPP1+

for p in patients_m:
    p_mask = (adata_mac.obs[cfg.COL_SAMPLE] == p).values
    p_mac_mono_mask = p_mask & mac_mono_mask
    total_mac = p_mac_mono_mask.sum()
    mac_mono_counts[p] = total_mac
    spp1_count = is_spp1_pos[p_mac_mono_mask].sum()
    spp1_tam_ratios[p] = spp1_count / total_mac if total_mac > 0 else np.nan

spp1_vals_list = [v for v in spp1_tam_ratios.values() if not np.isnan(v)]
print(f"SPP1+ TAM ratio (SPP1>0 in Mac/Mono, denom=Mac/Mono): range=[{min(spp1_vals_list):.4f}, {max(spp1_vals_list):.4f}], median={np.median(spp1_vals_list):.4f}")

obs_table_m = mye_obs.drop_duplicates(subset=cfg.COL_SAMPLE)
pcr_p_m = obs_table_m.loc[obs_table_m[cfg.COL_RESPONSE]=='pCR', cfg.COL_SAMPLE].values
nmp_p_m = obs_table_m.loc[obs_table_m[cfg.COL_RESPONSE]=='non-MPR', cfg.COL_SAMPLE].values

# Filter out patients with no myeloid data (NaN ratio) — same fix as step1 N10
spp1_pcr = [spp1_tam_ratios[p] for p in pcr_p_m if p in spp1_tam_ratios and not np.isnan(spp1_tam_ratios[p])]
spp1_nmp = [spp1_tam_ratios[p] for p in nmp_p_m if p in spp1_tam_ratios and not np.isnan(spp1_tam_ratios[p])]
n_pcr_excl_m = len(pcr_p_m) - len(spp1_pcr)
n_nmp_excl_m = len(nmp_p_m) - len(spp1_nmp)
if n_pcr_excl_m > 0 or n_nmp_excl_m > 0:
    print(f"  [Panel C] Excluded patients with no Mac/Mono data: pCR={n_pcr_excl_m}, non-MPR={n_nmp_excl_m}")

axC3 = axes3[1, 0]
bp_c3 = axC3.boxplot([spp1_pcr, spp1_nmp], patch_artist=True)
axC3.set_xticklabels(['pCR', 'non-MPR'])
bp_c3['boxes'][0].set_facecolor('#4CAF50'); bp_c3['boxes'][1].set_facecolor('#F44336')
_, mwp_spp1 = mannwhitneyu(spp1_pcr, spp1_nmp, alternative='two-sided')
axC3.set_title(f'Panel C: SPP1+ TAM ratio\nMWU p={mwp_spp1:.4e}')
axC3.set_ylabel('SPP1+ TAM proportion')

# Panel D: SPP1+ TAM ratio vs NK-dominant ratio scatter
patient_spp1_map = {p: v for p, v in spp1_tam_ratios.items() if not np.isnan(v)}
nk_ratios_d3 = []
spp1_vals_d3 = []
for _, row in df_patient.iterrows():
    pid = row['patient_id']
    if pid in patient_spp1_map:
        nk_ratios_d3.append(row['nk_dominant_ratio'])
        spp1_vals_d3.append(patient_spp1_map[pid])

axD3 = axes3[1, 1]
if len(nk_ratios_d3) > 2:
    sr3, sp3 = spearmanr(nk_ratios_d3, spp1_vals_d3)
else:
    sr3, sp3 = 0, 1
axD3.scatter(nk_ratios_d3, spp1_vals_d3, alpha=0.6, c='steelblue', s=20)
axD3.set_xlabel('NK-dominant ratio'); axD3.set_ylabel('SPP1+ TAM ratio')
axD3.set_title(f'Panel D: SPP1+ TAM vs NK-dominant\nSpearman r={sr3:.4f}, p={sp3:.4f}')

plt.tight_layout()
fig3.savefig(cfg.result_path('Fig3_myeloid.pdf'), dpi=150, bbox_inches='tight')
fig3.savefig(cfg.result_path('Fig3_myeloid.png'), dpi=150, bbox_inches='tight')
plt.close(fig3)
print("Fig3_myeloid saved")

# ── Save myeloid per-patient ──
mye_per_patient = []
for p in patients_m:
    p_data = mye_obs[mye_obs[cfg.COL_SAMPLE] == p]
    total = len(p_data)
    row_data = {'patient_id': p, 'response': p_data[cfg.COL_RESPONSE].iloc[0]}
    for st in mye_subtypes_all:
        row_data[st] = round((p_data[cfg.COL_CELL_TYPE] == st).sum() / total, 4) if total > 0 else 0
    row_data['mac_mono_count'] = int(mac_mono_counts.get(p, 0))
    row_data['SPP1_TAM_ratio'] = round(spp1_tam_ratios.get(p, np.nan), 4) if not np.isnan(spp1_tam_ratios.get(p, np.nan)) else ''
    mye_per_patient.append(row_data)

df_mye = pd.DataFrame(mye_per_patient)
df_mye.to_csv(cfg.result_path('myeloid_per_patient.csv'), index=False)
print("myeloid_per_patient.csv written")

# ============================================================
# FIGURE 4: Mechanism (NK-like vs Tex) (5 Panels, 2×3)
# Panel A: Pathway heatmap | Panel B: KLF2          | Panel C: Volcano
# Panel D: SPP1 receptors   | Panel E: TCF7 dysfunction | Panel F: (hidden)
# Panel E showcases "environment disrupts cell function" (TCF7+ cytotoxicity
# suppressed in High SPP1 env — integrates step3b core mechanism deep-dive)
# ============================================================
fig4, axes4 = plt.subplots(2, 3, figsize=(21, 12))

nk_obs = nk_cd8.obs
nk_like_mask = nk_obs[cfg.COL_CELL_TYPE] == 'CD8T_NK-like_FGFBP2'
tex_mask = nk_obs[cfg.COL_CELL_TYPE].isin(['CD8T_Tex_CXCL13', 'CD8T_terminal_Tex_LAYN'])

# Panel A: ssGSEA pathway heatmap placeholder
ax4A = axes4[0, 0]
# Build simple score: cytotoxicity = GZMB + GZMH + PRF1, stemness = TCF7 + LEF1, exhaustion = LAG3 + TIGIT + PDCD1
gene_sets = {
    'Cytotoxicity': ['GZMB', 'GZMH', 'PRF1', 'NKG7', 'GNLY'],
    'Stemness': ['TCF7', 'LEF1', 'SELL', 'IL7R'],
    'Exhaustion': ['LAG3', 'TIGIT', 'PDCD1', 'HAVCR2', 'CTLA4'],
    'Tissue_Residence': ['ITGAE', 'CD69', 'ZNF683', 'CXCR6'],
}

pathway_scores = {}
for name, genes in gene_sets.items():
    valid_genes = [g for g in genes if g in nk_cd8.var_names]
    if valid_genes:
        expr = nk_cd8[:, valid_genes].X.toarray() if hasattr(nk_cd8.X, 'toarray') else nk_cd8[:, valid_genes].X
        pathway_scores[name] = np.mean(expr, axis=1)
    else:
        pathway_scores[name] = np.zeros(nk_cd8.shape[0])

score_data = []
for name in pathway_scores:
    nk_like_score = np.mean(pathway_scores[name][nk_like_mask.values])
    tex_score = np.mean(pathway_scores[name][tex_mask.values])
    score_data.append({'pathway': name, 'NK-like_mean': round(nk_like_score, 4), 'Tex_mean': round(tex_score, 4)})

# Heatmap
score_matrix = np.array([[d['NK-like_mean'], d['Tex_mean']] for d in score_data])
im4a = ax4A.imshow(score_matrix, aspect='auto', cmap='RdYlBu_r')
ax4A.set_xticks([0, 1]); ax4A.set_xticklabels(['NK-like', 'Tex'])
ax4A.set_yticks(range(len(score_data))); ax4A.set_yticklabels([d['pathway'] for d in score_data])
for i in range(len(score_data)):
    for j in range(2):
        ax4A.text(j, i, f'{score_matrix[i,j]:.3f}', ha='center', va='center', fontsize=8)
ax4A.set_title('Panel A: Pathway scores (mean expression)')
plt.colorbar(im4a, ax=ax4A)

# Panel B: KLF2 & stemness genes in NK-like cells: High vs Low SPP1+ TAM (per-patient)
ax4B = axes4[0, 1]

stem_genes = ['KLF2', 'TCF7', 'LEF1', 'SELL', 'IL7R']
stem_genes_in_data = [g for g in stem_genes if g in nk_cd8.var_names]

high_thresh = np.percentile([v for v in spp1_tam_ratios.values() if not np.isnan(v)], 60)
low_thresh = np.percentile([v for v in spp1_tam_ratios.values() if not np.isnan(v)], 40)
high_spp1_pats = [p for p, v in spp1_tam_ratios.items() if not np.isnan(v) and v >= high_thresh]
low_spp1_pats = [p for p, v in spp1_tam_ratios.items() if not np.isnan(v) and v <= low_thresh]

per_pat_kfl2 = []
groups_list = []
n_cells_list = []
klf2_in_data = 'KLF2' in nk_cd8.var_names
if klf2_in_data:
    klf2_idx = list(nk_cd8.var_names).index('KLF2')
    for p in df_patient['patient_id']:
        p_nk_mask = nk_like_mask & (nk_cd8.obs[cfg.COL_SAMPLE] == p).values
        n_nk = p_nk_mask.sum()
        if n_nk == 0:
            continue
        p_klf2_expr = nk_cd8.X[p_nk_mask, klf2_idx]
        p_klf2_expr = p_klf2_expr.toarray().flatten() if hasattr(p_klf2_expr, 'toarray') else np.asarray(p_klf2_expr).flatten()
        mean_klf2 = np.mean(p_klf2_expr)
        per_pat_kfl2.append(mean_klf2)
        if p in high_spp1_pats:
            groups_list.append('High SPP1+ TAM')
        elif p in low_spp1_pats:
            groups_list.append('Low SPP1+ TAM')
        else:
            groups_list.append('Mid')
        n_cells_list.append(n_nk)

if klf2_in_data and len([g for g in groups_list if g != 'Mid']) >= 4:
    high_vals = [per_pat_kfl2[i] for i, g in enumerate(groups_list) if g == 'High SPP1+ TAM']
    low_vals = [per_pat_kfl2[i] for i, g in enumerate(groups_list) if g == 'Low SPP1+ TAM']
    stat_mw, p_mw = mannwhitneyu(high_vals, low_vals, alternative='two-sided')
    bp_data = [low_vals, high_vals]
    bp = ax4B.boxplot(bp_data, labels=['Low SPP1+ TAM', 'High SPP1+ TAM'], patch_artist=True, widths=0.5)
    colors = ['#2E86C1', '#E74C3C']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    for i, (data_list, color) in enumerate(zip([low_vals, high_vals], colors)):
        x_jitter = np.random.normal(i+1, 0.05, size=len(data_list))
        ax4B.scatter(x_jitter, data_list, c=color, s=20, zorder=5, alpha=0.8, edgecolors='black', linewidth=0.5)
    ax4B.set_ylabel('KLF2 expression (NK-like cells)')
    ax4B.set_title(f'Panel B: KLF2 in NK-like cells\nHigh vs Low SPP1+ TAM\nMWU p={p_mw:.4f}')
    print(f"\nKLF2 per-patient (NK-like only): High(n={len(high_vals)})={np.mean(high_vals):.4f}, Low(n={len(low_vals)})={np.mean(low_vals):.4f}, p={p_mw:.4f}")
else:
    ax4B.text(0.5, 0.5, 'KLF2 data not available', ha='center', va='center', transform=ax4B.transAxes)
    ax4B.set_title('Panel B: KLF2 in NK-like cells')

# Panel C: Volcano plot (NK-like vs Tex DEG)
nk_like_raw = nk_cd8[nk_like_mask].X
tex_raw = nk_cd8[tex_mask].X

nk_like_dense = nk_like_raw.toarray() if hasattr(nk_like_raw, 'toarray') else np.asarray(nk_like_raw)
tex_dense = tex_raw.toarray() if hasattr(tex_raw, 'toarray') else np.asarray(tex_raw)

nk_mean = np.mean(nk_like_dense, axis=0)
tex_mean = np.mean(tex_dense, axis=0)

log2fc = np.log2(nk_mean + 1) - np.log2(tex_mean + 1)

from scipy.stats import mannwhitneyu as mwu
from statsmodels.stats.multitest import multipletests

min_expr = 0.1
expressed = (nk_mean > min_expr) | (tex_mean > min_expr)
gene_indices = np.where(expressed)[0]
ngenes_test = len(gene_indices)
print(f"  Volcano: testing {ngenes_test}/{nk_cd8.shape[1]} genes (mean expr > {min_expr} in at least one group)")

pvals_all = np.ones(nk_cd8.shape[1])
for idx_i, gi in enumerate(gene_indices):
    if idx_i % 5000 == 0 and idx_i > 0:
        print(f"    ... {idx_i}/{ngenes_test} genes tested")
    _, pv = mwu(nk_like_dense[:, gi], tex_dense[:, gi], alternative='two-sided')
    pvals_all[gi] = pv

_, fdr_all, _, _ = multipletests(pvals_all, method='fdr_bh')
neg_log10fdr_all = -np.log10(fdr_all + 1e-300)

top_genes = ['FGFBP2', 'CXCL13', 'TCF7', 'NKG7', 'GZMB']

ax4C = axes4[0, 2]
plot_mask = expressed
sig_mask = (fdr_all < 0.05) & (np.abs(log2fc) > 0.5) & plot_mask
nonsig_mask = (~sig_mask) & plot_mask
ax4C.scatter(log2fc[nonsig_mask], neg_log10fdr_all[nonsig_mask], alpha=0.3, s=3, c='grey', rasterized=True)
ax4C.scatter(log2fc[sig_mask], neg_log10fdr_all[sig_mask], alpha=0.6, s=5, c='red', rasterized=True)
var_names = nk_cd8.var_names
for g in top_genes:
    if g in var_names:
        idx = list(var_names).index(g)
        ax4C.scatter(log2fc[idx], neg_log10fdr_all[idx], s=30, c='red', edgecolors='black', zorder=5)
        ax4C.annotate(g, (log2fc[idx], neg_log10fdr_all[idx]), fontsize=8, fontweight='bold')
ax4C.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=0.5, label='FDR=0.05')
ax4C.set_xlabel('log2FC (NK-like / Tex)'); ax4C.set_ylabel('-log10(FDR)')
ax4C.set_title(f'Panel C: Volcano plot (NK-like vs Tex, n={ngenes_test} genes)')

# Panel D: SPP1 receptors expression
spp1_receptors = ['CD44', 'ITGB1', 'ITGAV']
nk_like_obs = nk_cd8[nk_like_mask]
tex_obs = nk_cd8[tex_mask]

rec_data = []
rec_labels = []
for rec in spp1_receptors:
    if rec in nk_cd8.var_names:
        idx = list(nk_cd8.var_names).index(rec)
        nk_expr = nk_like_dense[:, idx]
        tex_expr_r = tex_dense[:, idx]
        rec_data.append(nk_expr)
        rec_data.append(tex_expr_r)
        rec_labels.append(f'{rec}\nNK-like')
        rec_labels.append(f'{rec}\nTex')

ax4D = axes4[1, 0]
if rec_data:
    bp4d = ax4D.boxplot(rec_data, patch_artist=True, widths=0.6)
    for i in range(0, len(rec_data), 2):
        bp4d['boxes'][i].set_facecolor('#FF6B6B')
        bp4d['boxes'][i+1].set_facecolor('#4ECDC4')
    ax4D.set_xticklabels(rec_labels, rotation=45, ha='right', fontsize=8)
    ax4D.set_title('Panel D: SPP1 receptor expression')
    ax4D.set_ylabel('Expression')

# Panel E: TCF7 dysfunction — cytotoxic gene expression in TCF7+ NK-like cells
# under High vs Low SPP1 environment (median split, consistent with step3b)
# Hypothesis: High SPP1 environment disrupts TCF7→cytotoxic differentiation program
ax4E = axes4[1, 1]

# SPP1 environment split (median, consistent with step3b_tcf7_dysfunction.py)
spp1_vals_for_median = [v for v in spp1_tam_ratios.values() if not np.isnan(v)]
spp1_median_val = np.median(spp1_vals_for_median)
high_spp1_med_pats = [p for p, v in spp1_tam_ratios.items() if not np.isnan(v) and v >= spp1_median_val]
low_spp1_med_pats = [p for p, v in spp1_tam_ratios.items() if not np.isnan(v) and v < spp1_median_val]
print(f"\n[Panel E] SPP1 median split: High(≥{spp1_median_val:.4f})={len(high_spp1_med_pats)}pats, Low(<{spp1_median_val:.4f})={len(low_spp1_med_pats)}pats")

# TCF7 expression in NK-like cells (sparse: ~84% cells = 0)
if 'TCF7' in nk_cd8.var_names:
    tcf7_gidx = list(nk_cd8.var_names).index('TCF7')
    tcf7_expr_all = nk_cd8.X[:, tcf7_gidx]
    tcf7_expr_all = tcf7_expr_all.toarray().flatten() if hasattr(tcf7_expr_all, 'toarray') else np.asarray(tcf7_expr_all).flatten()
    # TCF7_High = TCF7 > 0 (nonzero cells, matching step3b)
    tcf7_hi_mask_arr = tcf7_expr_all > 0
    print(f"  TCF7+ NK-like cells: {tcf7_hi_mask_arr.sum()}/{len(tcf7_hi_mask_arr)} ({100*tcf7_hi_mask_arr.mean():.1f}%)")

    # SPP1 environment for each cell
    pat_arr_e = nk_cd8.obs[cfg.COL_SAMPLE].values
    high_env_mask_e = np.isin(pat_arr_e, high_spp1_med_pats) & tcf7_hi_mask_arr
    low_env_mask_e = np.isin(pat_arr_e, low_spp1_med_pats) & tcf7_hi_mask_arr
    print(f"  TCF7+ in High SPP1 env: {high_env_mask_e.sum()} cells; Low SPP1 env: {low_env_mask_e.sum()} cells")

    cyto_genes_e = ['GZMB', 'PRF1', 'GNLY']
    cyto_data_e = {}
    for gene in cyto_genes_e:
        if gene in nk_cd8.var_names:
            gidx_e = list(nk_cd8.var_names).index(gene)
            gexpr_e = nk_cd8.X[:, gidx_e]
            gexpr_e = gexpr_e.toarray().flatten() if hasattr(gexpr_e, 'toarray') else np.asarray(gexpr_e).flatten()
            cyto_data_e[gene] = {
                'High_SPP1': gexpr_e[high_env_mask_e],
                'Low_SPP1': gexpr_e[low_env_mask_e],
            }

    # Violin plot: for each gene, Low vs High SPP1 side by side
    np.random.seed(42)
    plot_data_e = []
    positions_e = []
    labels_e = []
    colors_e = []
    p_vals_e = []
    for i, gene in enumerate(cyto_genes_e):
        if gene in cyto_data_e:
            lo_vals = cyto_data_e[gene]['Low_SPP1']
            hi_vals = cyto_data_e[gene]['High_SPP1']
            # subsample for plotting speed
            lo_sample = np.random.choice(lo_vals, min(5000, len(lo_vals)), replace=False) if len(lo_vals) > 5000 else lo_vals
            hi_sample = np.random.choice(hi_vals, min(5000, len(hi_vals)), replace=False) if len(hi_vals) > 5000 else hi_vals
            plot_data_e.extend([lo_sample, hi_sample])
            positions_e.extend([i*3 + 1, i*3 + 2])
            labels_e.extend([f'{gene}\nLow SPP1', f'{gene}\nHigh SPP1'])
            colors_e.extend(['#2E86C1', '#E74C3C'])
            if len(hi_vals) > 10 and len(lo_vals) > 10:
                _, mwu_p_e = mannwhitneyu(hi_vals, lo_vals, alternative='two-sided')
                p_vals_e.append((i, mwu_p_e, max(np.percentile(lo_vals, 95), np.percentile(hi_vals, 95))))
                print(f"  {gene} (TCF7+): High mean={np.mean(hi_vals):.4f}, Low mean={np.mean(lo_vals):.4f}, MWU p={mwu_p_e:.4e}")

    if plot_data_e:
        parts_e = ax4E.violinplot(plot_data_e, positions=positions_e, showmeans=True, showmedians=True, widths=0.8)
        for i, pc in enumerate(parts_e['bodies']):
            pc.set_facecolor(colors_e[i])
            pc.set_alpha(0.6)
        # p-value annotations
        for gene_i, pv, y_ref in p_vals_e:
            ax4E.text(gene_i*3 + 1.5, y_ref * 1.15, f'p={pv:.1e}',
                     ha='center', fontsize=7, fontweight='bold',
                     bbox=dict(boxstyle='round,pad=0.2', facecolor='wheat', alpha=0.7))
        ax4E.set_xticks(positions_e)
        ax4E.set_xticklabels(labels_e, fontsize=7)
        ax4E.set_ylabel('Expression in TCF7$^{+}$ NK-like cells')
        ax4E.set_title('Panel E: TCF7$^{+}$ cytotoxicity vs SPP1 env\n(High SPP1 → differentiation arrest)')
        ax4E.spines['top'].set_visible(False)
        ax4E.spines['right'].set_visible(False)
    else:
        ax4E.text(0.5, 0.5, 'Cytotoxic genes not available', ha='center', va='center', transform=ax4E.transAxes)
        ax4E.set_title('Panel E: TCF7 dysfunction')
else:
    ax4E.text(0.5, 0.5, 'TCF7 gene not found', ha='center', va='center', transform=ax4E.transAxes)
    ax4E.set_title('Panel E: TCF7 dysfunction')

# Panel F: (reserved — hide unused axis)
axes4[1, 2].axis('off')

plt.tight_layout()
fig4.savefig(cfg.result_path('Fig4_mechanism.pdf'), dpi=150, bbox_inches='tight')
fig4.savefig(cfg.result_path('Fig4_mechanism.png'), dpi=150, bbox_inches='tight')
plt.close(fig4)
print("Fig4_mechanism saved")

# ============================================================
# Ligand-Receptor Scoring (migrated from figure3.py)
# ============================================================
print("Generating ligand-receptor scores...")
lr_pairs = [
    ('SPP1', 'CD44'), ('SPP1', 'ITGAV'), ('SPP1', 'ITGB1'),
    ('CXCL8', 'CXCR2'), ('CCL2', 'CCR2'), ('IL1B', 'IL1R1'),
    ('TNF', 'TNFRSF1A'), ('VEGFA', 'FLT1'),
]

# SPP1+ TAM in myeloid (same definition as Panel C: SPP1 raw count > 0, in Mac/Mono)
mac_mask = immune.obs['major_cell_type'] == 'Myeloid cell'
m_tmp = immune[mac_mask].copy()
mac_sub_mask = m_tmp.obs[cfg.COL_CELL_TYPE].str.contains('Mφ|Mac|Mono', regex=True, na=False)
if 'SPP1' in m_tmp.var_names:
    gi = list(m_tmp.var_names).index('SPP1')
    spp1_raw_col = m_tmp.X[:, gi]
    spp1_raw = spp1_raw_col.toarray().ravel() if hasattr(spp1_raw_col, 'toarray') else np.asarray(spp1_raw_col).ravel()
else:
    spp1_raw = np.zeros(m_tmp.n_obs)
mac_spp1_raw = spp1_raw[mac_sub_mask.values]
ist_mask = np.zeros(m_tmp.n_obs, dtype=bool)
ist_mask[mac_sub_mask.values] = (mac_spp1_raw > 0)
logX_mye = normalize_log1p(m_tmp)

# Load T cells for receptor expression
t_adata = load_h5ad(cfg.H5AD_T)
t_adata.var_names = t_adata.var_names.astype(str)
t_logX = normalize_log1p(t_adata)
t_adata.obs['is_nklike'] = t_adata.obs[cfg.COL_CELL_TYPE] == 'CD8T_NK-like_FGFBP2'
t_adata.obs['is_tex'] = t_adata.obs[cfg.COL_CELL_TYPE].isin(['CD8T_Tex_CXCL13', 'CD8T_terminal_Tex_LAYN'])

ligand_expr = {}
receptor_expr_nk = {}
receptor_expr_tex = {}
interaction_nk = {}
interaction_tex = {}

for ligand, receptor in lr_pairs:
    pair_key = f'{ligand}-{receptor}'
    if ligand in m_tmp.var_names:
        gi = list(m_tmp.var_names).index(ligand)
        lig_col = logX_mye[:, gi]
        lig_vals = lig_col.toarray().ravel() if hasattr(lig_col, 'toarray') else np.asarray(lig_col).ravel()
        lig_mean = float(lig_vals[ist_mask].mean())
    else:
        lig_mean = 0
    ligand_expr[pair_key] = lig_mean

    if receptor in t_adata.var_names:
        ri = list(t_adata.var_names).index(receptor)
        t_arr = t_logX.toarray() if hasattr(t_logX, 'toarray') else np.asarray(t_logX)
        nk_mask = t_adata.obs['is_nklike'].values
        tex_mask = t_adata.obs['is_tex'].values
        rec_nk = float(t_arr[nk_mask, ri].mean()) if nk_mask.sum() > 0 else 0
        rec_tex = float(t_arr[tex_mask, ri].mean()) if tex_mask.sum() > 0 else 0
    else:
        rec_nk = rec_tex = 0
    receptor_expr_nk[pair_key] = rec_nk
    receptor_expr_tex[pair_key] = rec_tex
    interaction_nk[pair_key] = lig_mean * rec_nk
    interaction_tex[pair_key] = lig_mean * rec_tex

lr_df_out = pd.DataFrame({
    'pair': list(ligand_expr.keys()),
    'ligand_expr_SPP1_TAM': list(ligand_expr.values()),
    'receptor_expr_NKlike': list(receptor_expr_nk.values()),
    'receptor_expr_Tex': list(receptor_expr_tex.values()),
    'interaction_NKlike': list(interaction_nk.values()),
    'interaction_Tex': list(interaction_tex.values()),
})
lr_df_out.to_csv(cfg.result_path('ligand_receptor_scores.csv'), index=False)
print(f"ligand_receptor_scores.csv saved: {len(lr_df_out)} pairs")
print("  NK-like side (假说对应：SPP1+ TAM → NK-like):")
for k, v in sorted(interaction_nk.items(), key=lambda x: -x[1]):
    print(f"    {k}: NK_score={v:.4f} (L={ligand_expr[k]:.2f}, R_NK={receptor_expr_nk[k]:.2f})")
print("  Tex side (对照：SPP1+ TAM → Tex):")
for k, v in sorted(interaction_tex.items(), key=lambda x: -x[1]):
    print(f"    {k}: Tex_score={v:.4f} (L={ligand_expr[k]:.2f}, R_Tex={receptor_expr_tex[k]:.2f})")

# ============================================================
# FIGURE S5: Interaction Network (2 Panels: network + top interactions bar)
# ============================================================
lr_path = cfg.result_path('ligand_receptor_scores.csv')
if os.path.exists(lr_path):
    df_lr = pd.read_csv(lr_path)
    print(f"Loaded ligand_receptor_scores: {df_lr.shape}")
    
    # 双面板布局：左网络图，右 top interactions 条形图
    figS5, (axS5_net, axS5_bar) = plt.subplots(1, 2, figsize=(16, 7),
                                                gridspec_kw={'width_ratios': [1.2, 1]})
    import networkx as nx
    G = nx.Graph()
    node_colors_s5 = {}
    node_sizes_s5 = {}
    edge_data_s5 = []
    for _, row in df_lr.iterrows():
        pair_str = str(row.get('pair', ''))
        if '-' in pair_str:
            source, target = pair_str.split('-', 1)
        else:
            source = row.get('ligand', 'Unknown')
            target = row.get('receptor', 'Unknown')
        score = float(row.get('interaction_NKlike', 1.0))
        score_tex = float(row.get('interaction_Tex', 1.0))
        G.add_edge(source, target, weight=score, weight_tex=score_tex)
        node_colors_s5[source] = '#FF6B6B'   # ligand (TAM side)
        node_colors_s5[target] = '#4ECDC4'   # receptor (T/NK side)
        # node size proportional to degree
        node_sizes_s5[source] = node_sizes_s5.get(source, 300) + 200
        node_sizes_s5[target] = node_sizes_s5.get(target, 300) + 200
        edge_data_s5.append({'pair': pair_str, 'source': source, 'target': target,
                             'NKlike': score, 'Tex': score_tex})
    
    # 使用 kamada_kawai 布局（比 spring 更紧凑）
    try:
        pos = nx.kamada_kawai_layout(G)
    except Exception:
        pos = nx.spring_layout(G, seed=42, k=1.5)
    
    edges = G.edges()
    weights = [G[u][v]['weight'] * 1.5 for u, v in edges]
    node_colors = [node_colors_s5.get(n, '#95E1D3') for n in G.nodes()]
    node_sizes = [node_sizes_s5.get(n, 400) for n in G.nodes()]
    
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, ax=axS5_net, alpha=0.85, edgecolors='black', linewidths=0.8)
    nx.draw_networkx_edges(G, pos, width=weights, alpha=0.6, ax=axS5_net, edge_color='gray')
    nx.draw_networkx_labels(G, pos, font_size=9, font_weight='bold', ax=axS5_net)
    
    # 边权重标签（只标 top 3）
    edge_list = list(G.edges(data=True))
    edge_list.sort(key=lambda x: -x[2]['weight'])
    for u, v, d in edge_list[:3]:
        mid_x = (pos[u][0] + pos[v][0]) / 2
        mid_y = (pos[u][1] + pos[v][1]) / 2
        axS5_net.text(mid_x, mid_y, f"{d['weight']:.2f}", fontsize=8, color='red',
                      ha='center', va='center', bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='red'))
    
    # Highlight key edges
    for u, v, d in G.edges(data=True):
        if ('SPP1' in str(u) and 'ITGB1' in str(v)) or ('SPP1' in str(v) and 'ITGB1' in str(u)):
            nx.draw_networkx_edges(G, pos, edgelist=[(u,v)], width=4, edge_color='red', ax=axS5_net)
        if ('SPP1' in str(u) and 'CD44' in str(v)) or ('SPP1' in str(v) and 'CD44' in str(u)):
            nx.draw_networkx_edges(G, pos, edgelist=[(u,v)], width=4, edge_color='orange', ax=axS5_net)
    
    axS5_net.set_title('Ligand-Receptor Network\n(node size ∝ degree, edge width ∝ NK-like interaction score)', fontsize=11)
    axS5_net.axis('off')
    
    # 右侧条形图：top interactions NK-like vs Tex
    df_edges = pd.DataFrame(edge_data_s5).sort_values('NKlike', ascending=True)
    y_pos = np.arange(len(df_edges))
    bar_h = 0.35
    axS5_bar.barh(y_pos + bar_h/2, df_edges['NKlike'], bar_h, color='#4ECDC4', label='NK-like', alpha=0.8, edgecolor='black', linewidth=0.5)
    axS5_bar.barh(y_pos - bar_h/2, df_edges['Tex'], bar_h, color='#FF6B6B', label='Tex', alpha=0.8, edgecolor='black', linewidth=0.5)
    axS5_bar.set_yticks(y_pos)
    axS5_bar.set_yticklabels(df_edges['pair'], fontsize=9)
    axS5_bar.set_xlabel('Interaction Score (ligand_expr × receptor_expr)')
    axS5_bar.set_title('L-R Interaction Scores\n(NK-like vs Tex side-by-side)', fontsize=11)
    axS5_bar.legend(loc='lower right', fontsize=9)
    axS5_bar.spines['top'].set_visible(False); axS5_bar.spines['right'].set_visible(False)
    # 数值标签
    for i, (nk, tex) in enumerate(zip(df_edges['NKlike'], df_edges['Tex'])):
        if nk > 0.05:
            axS5_bar.text(nk + 0.05, i + bar_h/2, f'{nk:.2f}', va='center', fontsize=7, color='#1A5276')
        if tex > 0.05:
            axS5_bar.text(tex + 0.05, i - bar_h/2, f'{tex:.2f}', va='center', fontsize=7, color='#922B21')
    
    plt.tight_layout()
    figS5.savefig(cfg.result_path('FigS5_interaction_network.pdf'), dpi=150, bbox_inches='tight')
    figS5.savefig(cfg.result_path('FigS5_interaction_network.png'), dpi=150, bbox_inches='tight')
    plt.close(figS5)
    print("FigS5_interaction_network saved")
else:
    print("[WARNING] ligand_receptor_scores.csv not found, skipping FigS5")

# ── Save step3 statistics to CSV ──
step3_stats = {
    'spp1_tam_mwu_p': mwp_spp1,
    'spp1_tam_pcr_median': float(np.median(spp1_pcr)),
    'spp1_tam_nmp_median': float(np.median(spp1_nmp)),
    'spp1_tam_pcr_n': len(spp1_pcr),
    'spp1_tam_nmp_n': len(spp1_nmp),
    'spp1_tam_total_n': len(spp1_pcr) + len(spp1_nmp),
    'spp1_tam_missing_n': 188 - (len(spp1_pcr) + len(spp1_nmp)),
    'spp1_tam_spearman_r': sr3,
    'spp1_tam_spearman_p': sp3,
    'klf2_high_spp1_mean': float(np.mean(high_vals)) if 'high_vals' in dir() and len(high_vals) > 0 else np.nan,
    'klf2_low_spp1_mean': float(np.mean(low_vals)) if 'low_vals' in dir() and len(low_vals) > 0 else np.nan,
    'klf2_mwu_p': p_mw if 'p_mw' in dir() else np.nan,
}
step3_stats_df = pd.DataFrame({'metric': list(step3_stats.keys()), 'value': list(step3_stats.values())})
step3_stats_df.to_csv(cfg.result_path('step3_stats.csv'), index=False)
print("step3_stats.csv saved")

# ============================================================
# 补充分析: Taxane vs Pemetrexed 化疗方案分层下 NK-like 功能基因表达对比
# 直接回答"方案特异性不是源于 NK-like 本身功能差异"
# ============================================================
print("\n" + "="*60)
print("Taxane vs Pemetrexed: NK-like functional gene comparison")
print("="*60)

from _common import classify_chemo

# 从 immune 中提取化疗方案映射
chemo_map = immune.obs.groupby(cfg.COL_SAMPLE)['chemotherapy'].first().to_dict()
nk_cd8.obs['chemo_class'] = nk_cd8.obs[cfg.COL_SAMPLE].map(chemo_map).apply(classify_chemo)

# 功能基因集
FUNCTIONAL_GENE_SETS = {
    'Effector': ['GNLY', 'PRF1', 'GZMB', 'GZMH', 'GZMA', 'IFNG'],
    'Chemokine': ['CX3CR1', 'CCL5', 'XCL1', 'XCL2', 'NKG7'],
    'Receptor': ['KLRD1', 'KLRC1', 'KLRC2', 'KLRB1', 'KLRF1', 'CD160', 'CRTAM', 'CTSW', 'SH2D1B', 'FCGR3A'],
    'Transcription': ['TBX21', 'EOMES', 'KLF2', 'TCF7', 'ZNF683', 'TYROBP'],
}

# 标准化表达矩阵
from _common import normalize_log1p
X_nk = normalize_log1p(nk_cd8)
var_names = list(nk_cd8.var_names)

# 分组
taxane_mask = (nk_cd8.obs['chemo_class'] == 'Platinum+Taxane').values
peme_mask = (nk_cd8.obs['chemo_class'] == 'Platinum+Pemetrexed').values
print(f"  Taxane NK-like cells: {taxane_mask.sum()}")
print(f"  Pemetrexed NK-like cells: {peme_mask.sum()}")
print(f"  Taxane patients: {nk_cd8.obs.loc[taxane_mask, cfg.COL_SAMPLE].nunique()}")
print(f"  Pemetrexed patients: {nk_cd8.obs.loc[peme_mask, cfg.COL_SAMPLE].nunique()}")

func_results = []
fig_s6, axes_s6 = plt.subplots(2, 2, figsize=(12, 10))
axes_flat = axes_s6.flatten()

for idx, (gs_name, genes) in enumerate(FUNCTIONAL_GENE_SETS.items()):
    avail = [g for g in genes if g in var_names]
    if len(avail) == 0:
        print(f"  {gs_name}: no genes found, skipping")
        continue
    gene_idx = [var_names.index(g) for g in avail]
    X_subset = X_nk[:, gene_idx]
    if hasattr(X_subset, 'toarray'):
        X_subset = X_subset.toarray()
    # 细胞水平均值
    cell_mean = X_subset.mean(axis=1)

    taxane_vals = cell_mean[taxane_mask]
    peme_vals = cell_mean[peme_mask]

    # 细胞水平 MWU
    try:
        u_stat, p_cell = mannwhitneyu(taxane_vals, peme_vals, alternative='two-sided')
    except:
        p_cell = np.nan

    # 患者水平（避免伪重复）
    nk_cd8_obs = nk_cd8.obs.copy()
    nk_cd8_obs['cell_mean'] = cell_mean
    patient_mean = nk_cd8_obs.groupby(cfg.COL_SAMPLE)['cell_mean'].mean()
    patient_chemo = nk_cd8_obs.groupby(cfg.COL_SAMPLE)['chemo_class'].first()
    taxane_pat = patient_mean[patient_chemo == 'Platinum+Taxane'].values
    peme_pat = patient_mean[patient_chemo == 'Platinum+Pemetrexed'].values
    try:
        _, p_patient = mannwhitneyu(taxane_pat, peme_pat, alternative='two-sided')
    except:
        p_patient = np.nan

    log2fc = np.log2((taxane_vals.mean() + 1e-6) / (peme_vals.mean() + 1e-6))
    print(f"  {gs_name}: Taxane={taxane_vals.mean():.4f}, Peme={peme_vals.mean():.4f}, log2FC={log2fc:.4f}, p_cell={p_cell:.4f}, p_patient={p_patient:.4f}")

    func_results.append({
        'gene_set': gs_name,
        'n_genes': len(avail),
        'taxane_mean': round(taxane_vals.mean(), 4),
        'peme_mean': round(peme_vals.mean(), 4),
        'log2fc': round(log2fc, 4),
        'mwu_p_cell_level': round(p_cell, 6) if not np.isnan(p_cell) else 'NA',
        'mwu_p_patient_level': round(p_patient, 6) if not np.isnan(p_patient) else 'NA',
        'n_taxane_cells': int(taxane_mask.sum()),
        'n_peme_cells': int(peme_mask.sum()),
        'n_taxane_patients': int(nk_cd8.obs.loc[taxane_mask, cfg.COL_SAMPLE].nunique()),
        'n_peme_patients': int(nk_cd8.obs.loc[peme_mask, cfg.COL_SAMPLE].nunique()),
    })

    # 箱线图
    ax = axes_flat[idx]
    df_plot = pd.DataFrame({
        'value': np.concatenate([taxane_vals, peme_vals]),
        'group': ['Taxane'] * len(taxane_vals) + ['Pemetrexed'] * len(peme_vals)
    })
    df_plot.boxplot(column='value', by='group', ax=ax, grid=False)
    ax.set_title(f'{gs_name} (n_genes={len(avail)})\npatient-level p={p_patient:.4f}' if not np.isnan(p_patient) else f'{gs_name}')
    ax.set_ylabel('Mean log-norm expression')

# 保存
df_func = pd.DataFrame(func_results)
df_func.to_csv(cfg.result_path('chemo_functional_comparison.csv'), index=False)
print(f"\n  chemo_functional_comparison.csv saved")

fig_s6.suptitle('Taxane vs Pemetrexed: NK-like functional gene expression', fontsize=14, y=1.02)
plt.tight_layout()
fig_s6.savefig(cfg.result_path('FigS6_chemo_functional.png'), dpi=150, bbox_inches='tight')
fig_s6.savefig(cfg.result_path('FigS6_chemo_functional.pdf'), dpi=150, bbox_inches='tight')
plt.close(fig_s6)
print("  FigS6_chemo_functional saved")

# 更新 step3_stats
if func_results:
    effector_p = [r['mwu_p_patient_level'] for r in func_results if r['gene_set'] == 'Effector']
    if effector_p and effector_p[0] != 'NA':
        step3_stats['chemo_functional_effector_p_patient'] = effector_p[0]
    step3_stats_df = pd.DataFrame({'metric': list(step3_stats.keys()), 'value': list(step3_stats.values())})
    step3_stats_df.to_csv(cfg.result_path('step3_stats.csv'), index=False)

print("[step3_mechanism_exploration] Done.")
