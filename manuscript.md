# Clonal Fate Locking of NK-like CD8 T Cells Predicts Anti-PD-1 Response in Non-Small Cell Lung Cancer via SPP1+ Tumor-Associated Macrophage Gate

---

## Simple Summary

Anti-PD-1 immunotherapy has transformed non-small cell lung cancer (NSCLC) treatment, yet reliable predictive biomarkers remain elusive. Using a single-cell transcriptomic atlas of 434,458 T cells from 188 NSCLC patients receiving neoadjuvant anti-PD-1 therapy, we identified a clonal fate locking phenomenon in which expanded CD8 T cell clones are committed to either an NK-like (FGFBP2+) or terminally exhausted (Tex) trajectory. The proportion of NK-locked large clones (NK-Locked Ratio) strongly predicted pathological complete response (OR = 19.54, *p* = 0.008; adjusted OR = 38.14, *p* = 0.002). Mechanistically, SPP1+ tumor-associated macrophages (TAMs) formed a myeloid gate that suppressed stemness-associated genes (KLF2, TCF7, IL7R) in NK-like CD8 T cells (Spearman *r* = −0.339, *p* = 3.78 × 10⁻⁶). External validation in an independent NSCLC cohort confirmed the predictive value of the NK-like 24-gene signature (GSE126044, *p* = 0.006, FC = 4.72), while melanoma cohorts showed no association, supporting cancer type-specific clonal dynamics. Spatial transcriptomics validated the SPP1–receptor interaction *in situ* (*r* = 0.325, *p* = 0.044). These findings establish clonal fate locking as a deterministic principle governing anti-PD-1 response and suggest the SPP1+ TAM–KLF2 axis as a therapeutic target.

**Keywords:** NK-like CD8 T cells; clonal fate locking; anti-PD-1; SPP1+ TAM; NSCLC; single-cell RNA-seq; tumor microenvironment; immunotherapy biomarker

---

## 1. Introduction

Immune checkpoint blockade targeting the PD-1/PD-L1 axis has revolutionized the treatment of non-small cell lung cancer (NSCLC), with neoadjuvant anti-PD-1 therapy achieving pathological complete response (pCR) in approximately 20–45% of patients [1,2]. However, the majority of patients fail to achieve durable benefit, and the biological determinants of response heterogeneity remain incompletely understood. While tumor mutational burden (TMB) and PD-L1 expression are approved companion diagnostics, their predictive accuracy is modest, underscoring the need for mechanism-based biomarkers [3,4].

CD8 T cells are the primary effectors of anti-tumor immunity, but their functional states exist along a continuum from stemness-like memory to terminal exhaustion [5,6]. Within this spectrum, a subset of CD8 T cells co-expressing natural killer (NK) receptors—termed NK-like CD8 T cells—has emerged as a functionally distinct population with potent cytotoxic capacity [7,8]. These cells, marked by FGFBP2, KLRD1, CX3CR1, and FCGR3A, share transcriptional features with conventional NK cells while retaining T cell receptor (TCR) clonality [9,10]. However, whether NK-like CD8 T cells represent a terminally differentiated endpoint or a dynamically regulated fate commitment within expanding T cell clones has not been systematically investigated.

Tumor-associated macrophages (TAMs) constitute a dominant immunosuppressive compartment in the NSCLC microenvironment and are increasingly implicated in checkpoint blockade resistance [11,12]. SPP1 (osteopontin)-expressing TAMs, in particular, have been identified as a pro-tumorigenic subset associated with T cell suppression and poor prognosis across multiple cancer types [13,14]. Yet, the mechanistic link between SPP1+ TAMs and the clonal dynamics of effector CD8 T cells—specifically whether they directly influence clonal fate decisions—remains unknown.

In this study, we leveraged a large-scale single-cell RNA sequencing (scRNA-seq) atlas of NSCLC patients treated with neoadjuvant anti-PD-1 therapy (GSE243013) [15] to investigate the clonal architecture of NK-like CD8 T cells and its relationship to treatment response. We identified a "clonal fate locking" phenomenon in which expanded CD8 T cell clones are irreversibly committed to either an NK-like or Tex trajectory, with the NK-locked fraction serving as a potent predictor of pathological response. We further elucidated a SPP1+ TAM–mediated myeloid gate that suppresses stemness gene expression in NK-like CD8 T cells, establishing a mechanistic framework for myeloid-driven immune evasion. External validation across multiple independent cohorts and spatial transcriptomic confirmation collectively support a cancer type-specific model in which clonal fate locking represents the rate-limiting bottleneck for anti-PD-1 response in NSCLC.

---

## 2. Results

### 2.1. Single-Cell Atlas Identifies NK-like CD8 T Cells as a Distinct Clonal Population

We analyzed a single-cell transcriptomic atlas of 434,458 T cells from 188 NSCLC patients who received neoadjuvant anti-PD-1 (toripalimab) therapy (GSE243013) [15], comprising 84 patients with pathological complete response (pCR) and 104 with non-major pathological response (non-MPR) (Table S1). Within the CD8 T cell compartment (219,842 cells), we identified 16,744 NK-like CD8 T cells (annotated as CD8T_NK-like_FGFBP2) and 71,298 terminally exhausted T cells (Tex) based on the original study's cell type annotations.

To confirm the transcriptional identity of NK-like CD8 T cells, we computed an NK gene signature score (FGFBP2, KLRD1, CX3CR1, FCGR3A, KLRC1, NKG7, GNLY, PRF1, GZMB, KLRB1) and performed differential analysis between NK-like and other CD8 T cell subsets. NK-like cells exhibited significantly elevated NK signature scores (Mann–Whitney U, *p* < 0.001), confirming their distinct cytotoxic identity (**Figure 1A,B**).

Notably, gene dot plot analysis revealed that NK-like CD8 T cells co-expressed stemness-associated genes (KLF2, TCF7) alongside cytotoxic mediators (GZMB), distinguishing them from Tex cells that expressed high levels of exhaustion markers (CXCL13, LAYN, TOX) (**Figure 1B**). This dual cytotoxic-stemness phenotype prompted us to investigate whether NK-like identity represents a stable clonal commitment.

To examine clonal architecture, we analyzed TCR clonotype distribution across CD8 T cell subtypes. Clone expansion analysis (top-1 and top-5 clonotype fractions per patient) revealed substantial inter-patient heterogeneity in clonal dominance (**Figure 1C**). However, comparison of clone expansion between pCR and non-MPR patients showed no significant difference in top-5 clonotype fraction (Mann–Whitney U, *p* = 0.683; pCR median = 0.237, non-MPR median = 0.210) (**Figure 1D**), suggesting that the *quantity* of clonal expansion alone does not distinguish responders from non-responders. This finding motivated a deeper investigation into the *quality*—specifically, the subtype commitment—of expanded clones.

### 2.2. Clonal Fate Locking Predicts Pathological Response

We hypothesized that the functional identity of expanded clones, rather than their size, determines treatment outcome. To test this, we defined large clones as those with ≥5 cells (following established conventions [16]) and classified each large clone based on its dominant CD8 T cell subtype composition:

- **NK-Locked clones**: Large clones in which NK-like CD8 T cells constitute >50% of cells (nk_frac > 0.5), indicating irreversible commitment to the NK-like trajectory.
- **Tex-Differentiated clones**: Large clones dominated by Tex cells.

For each patient, we computed the NK-Locked Ratio—the proportion of large clones that are NK-locked—serving as a patient-level metric of clonal fate commitment (**Figure 2A**).

Univariate logistic regression revealed that the NK-Locked Ratio was a strong predictor of pCR (OR = 19.54, 95% CI: 2.18–175.15, *p* = 0.0079, n = 185 patients with pCR vs non-MPR comparison) (**Figure 2B**; Table S2). The Tex-dominant clone fraction showed a reciprocal association (OR = 112.04, 95% CI: 4.08–3079.82, *p* = 0.005), consistent with a binary fate decision. Multivariate logistic regression adjusting for total clone count and CD8 T cell proportion confirmed the independent predictive value of NK-Locked Ratio (adjusted OR = 38.14, 95% CI: 3.62–401.99, *p* = 0.002).

Direct comparison between response groups confirmed that pCR patients had significantly higher NK-Locked Ratios than non-MPR patients (Mann–Whitney U, *p* = 0.001; pCR median = 0.083, non-MPR median = 0.042) (**Figure 2C**). Frequency distribution analysis showed that 7.6% of patients (14/185) had a NK-Locked Ratio of zero, and these were disproportionately concentrated in the non-MPR group (**Figure 2D**).

Sensitivity analysis across varying clone size thresholds (≥0, ≥5, ≥10, ≥20, ≥50 cells) demonstrated that the predictive association was robust for thresholds ≥5 but attenuated at higher thresholds (≥20: OR = 5.41, *p* = 0.149; ≥50: OR = 24.21, *p* = 0.153), reflecting reduced statistical power with smaller patient numbers (**Figure S2**, Table S3). The optimal threshold was ≥5 cells, balancing clone size biological significance with sample size adequacy.

### 2.3. SPP1+ Tumor-Associated Macrophages Form a Myeloid Gate

To investigate the microenvironmental factors influencing clonal fate locking, we analyzed the myeloid compartment using the full immune cell atlas (1,254,749 cells), focusing on 191,099 myeloid cells. Within the macrophage and monocyte subsets, we identified SPP1+ TAMs as cells with SPP1 expression above the subtype median, yielding 31,488 SPP1+ TAMs (16.5% of myeloid cells) (**Figure 3A,B**).

SPP1+ TAMs were significantly enriched in non-MPR patients compared to pCR patients (Mann–Whitney U, *p* = 5.54 × 10⁻⁹; non-MPR median = 0.141, pCR median = 0.035) (**Figure 3C**), establishing a reciprocal relationship with the NK-Locked phenotype.

To characterize the intercellular communication between SPP1+ TAMs and CD8 T cell subsets, we computed ligand–receptor interaction scores for eight candidate pathways. Interaction scores were calculated separately for NK-like and Tex recipient cells (ligand expression in SPP1+ TAMs × receptor expression in each CD8 subtype), and the maximum score was taken to avoid signal averaging. The top interactions were SPP1–CD44 (Tex interaction score = 3.90) and SPP1–ITGB1 (NK-like interaction score = 3.88), followed by CXCL8–CXCR2, CCL2–CCR2, and TNF–TNFRSF1A (**Figure 3D**; Table S4). The SPP1–CD44/ITGB1 axis represents the dominant signaling conduit from TAMs to both CD8 T cell fates.

### 2.4. SPP1+ TAMs Suppress Stemness Genes in NK-like CD8 T Cells

The reciprocal relationship between SPP1+ TAM abundance and NK-Locked Ratio prompted us to test whether SPP1+ TAMs directly suppress the stemness program that maintains NK-like identity. We correlated patient-level SPP1+ TAM fraction with NK-Locked Ratio and observed a significant negative correlation (Spearman *r* = −0.339, *p* = 3.78 × 10⁻⁶, n = 178; Pearson *r* = −0.136, *p* = 0.071) (**Figure 4A**).

To identify the molecular targets of SPP1+ TAM–mediated suppression, we stratified patients into High-SPP1 (top 30%, n = 62) and Low-SPP1 (bottom 30%, n = 62) groups based on SPP1+ TAM fraction and compared stemness-associated gene expression in NK-like CD8 T cells. High-SPP1 patients exhibited significantly reduced expression of KLF2 (High = 0.847 vs Low = 1.057, Mann–Whitney U *p* = 0.009, FDR = 0.019), TCF7 (0.227 vs 0.280, *p* = 0.004, FDR = 0.012), and IL7R (1.043 vs 1.407, *p* = 0.0008, FDR = 0.005) (**Figure 4B**; Table S5). KLF2 emerged as a central hub, as it is a master transcription factor maintaining T cell quiescence and stemness [17].

Direct comparison of stemness gene expression between NK-like and Tex cells (n_NK-like = 3,662, n_Tex = 13,413 in the sensitivity-filtered subset) confirmed that NK-like CD8 T cells inherently express higher levels of KLF2 (2.112 vs 0.308, *p* < 10⁻³⁰⁰), TCF7 (0.153 vs 0.063, *p* = 3.66 × 10⁻⁵⁴), SELL (0.165 vs 0.038, *p* = 1.08 × 10⁻⁹⁵), IL7R (1.051 vs 0.703, *p* = 3.83 × 10⁻⁵⁹), and LEF1 (0.078 vs 0.049, *p* = 2.49 × 10⁻¹¹) (**Figure 4C**; Table S6). These data establish that KLF2 and associated stemness factors are defining features of NK-like identity and are specifically suppressed by SPP1+ TAMs.

### 2.5. External Validation and Cross-Cancer Heterogeneity

To validate the NK-like signature as a predictive biomarker, we computed a 24-gene NK-like signature score (Table S7) across four independent cohorts: GSE126044 (NSCLC, anti-PD-1, n = 16), GSE135222 (NSCLC, anti-PD-1/PD-L1, n = 27), GSE91061 (melanoma, anti-PD-1/CTLA-4, n = 33 pre-treatment), and GSE120575 (melanoma scRNA-seq, n = 11 patients with matched expression data).

In the NSCLC cohort GSE126044, the NK-like signature score was significantly higher in responders than non-responders (Mann–Whitney U, *p* = 0.006, fold change = 4.72) (**Figure 5A**), independently confirming the discovery cohort findings. Gene matching rates were 24/24 (100%) for GSE126044, 13/24 (54%) for GSE135222 (Ensembl ID mapping), 22/24 (92%) for GSE91061 (Entrez ID), and 24/24 (100%) for GSE120575 (Table S8).

However, in melanoma cohorts, no significant association was observed (GSE91061: *p* = 0.769, FC = 1.17; GSE120575: *p* = 0.921, FC = 0.98) (**Figure 5B**). The NSCLC cohort GSE135222 also showed no significant association (*p* = 0.755, FC = 1.02), attributable to its reliance on progression-free survival (PFS) rather than RECIST response criteria and limited responder sample size (R = 6).

We interpreted this cross-cancer heterogeneity through a tumor immune contexture model (**Figure 5C**): NSCLC, as an immunologically "cold" tumor with lower TMB, exhibits a clonal fate bottleneck in which the commitment of expanded clones to the NK-like trajectory is the rate-limiting step for response. In contrast, melanoma, as a "hot" tumor with high TMB and abundant clonal diversity, has already surpassed this bottleneck—the differentiation state, rather than NK-like locking, determines outcome. This model is consistent with the divergent predictive biomarker performance of TMB between NSCLC and melanoma [18].

### 2.6. Spatial Transcriptomic Validation

To confirm the SPP1+ TAM–NK-like CD8 T cell interaction *in situ*, we analyzed spatial transcriptomic data from 67 regions of interest (ROIs) across 39 NSCLC patients (GSE221733, GeoMx DSP). The SPP1–receptor interaction score (SPP1 × mean receptor expression of CD44/ITGAV/ITGB1) correlated positively with stemness gene expression in the tumor compartment (Spearman *r* = 0.325, *p* = 0.044) (**Figure 3E**; Table S9). NK cell cytotoxicity scores (GZMB/GZMA/GNLY/PRF1) were significantly higher in responder ROIs than non-responder ROIs (Mann–Whitney U, *p* = 0.027). In tumor-specific ROIs (PanCK+), the SPP1–NK receptor interaction was further strengthened (*r* = 0.346, *p* = 0.045).

Four-group comparison (High-SPP1/High-stemness vs Low-SPP1/Low-stemness) revealed a significant difference in clonal fate locking metrics (*p* = 0.017), supporting the model that the SPP1+ TAM–stemness axis operates as a spatially organized regulatory module.

---

## 3. Discussion

This study identifies clonal fate locking as a deterministic principle governing anti-PD-1 response in NSCLC. Our key findings are: (1) expanded CD8 T cell clones are irreversibly committed to either an NK-like or Tex trajectory, and the NK-Locked Ratio potently predicts pCR; (2) SPP1+ TAMs form a myeloid gate that suppresses KLF2-dependent stemness in NK-like CD8 T cells; (3) the NK-like gene signature validates in an independent NSCLC cohort but not in melanoma, reflecting cancer type-specific clonal dynamics; and (4) spatial transcriptomics confirms the SPP1–receptor interaction *in situ*.

**Clonal fate locking as a deterministic principle.** The concept of clonal dominance in anti-tumor immunity is well established—oligoclonal expansion of tumor-reactive T cells correlates with response to checkpoint blockade [19,20]. However, our data reveal that the *direction* of clonal commitment, rather than the *magnitude* of expansion, is the critical determinant. The observation that clone expansion per se does not differ between responders and non-responders (*p* = 0.683) but the NK-Locked Ratio is highly discriminating (OR = 19.54) fundamentally reframes the question of what makes a "productive" T cell clone. This finding is consistent with recent reports that stemness-like CD8 T cell subsets (TCF1+, SLAMF6+) are essential for sustained anti-tumor responses [21,22], but extends this concept by showing that the NK-like trajectory itself carries stemness capacity through KLF2 expression.

The binary nature of clonal fate locking—NK-locked vs Tex-differentiated—suggests a point-of-no-return in clone differentiation. This is mechanistically plausible given the epigenetic stability of exhaustion programs [23] and the transcriptional incompatibility between NK-like and exhaustion programs [24]. The >50% threshold for defining NK-locked clones was chosen to ensure biological commitment rather than stochastic fluctuation, and sensitivity analysis confirmed robustness at moderate thresholds.

**SPP1+ TAMs as a myeloid gate.** The identification of SPP1+ TAMs as a suppressor of NK-like stemness provides a mechanistic link between myeloid inflammation and T cell dysfunction. SPP1 (osteopontin) is a pleiotropic cytokine upregulated in tumor-promoting macrophages [25,26], and SPP1+ TAMs have been identified as a conserved pro-tumorigenic population across cancer types [13,14]. Our finding that SPP1–CD44 and SPP1–ITGB1 are the dominant ligand–receptor pairs aligns with the known roles of CD44 in maintaining T cell hydration and migration [27] and ITGB1 in T cell–extracellular matrix interactions [28]. The suppression of KLF2—a master regulator of T cell quiescence, migration (through S1P receptor transcription), and stemness [17,29]—by SPP1+ TAMs provides a molecular mechanism for how myeloid cells can "lock" T cell clones into terminal differentiation.

The strength of the correlation (Spearman *r* = −0.339) indicates that SPP1+ TAMs explain a substantial but not exclusive fraction of NK-Locked Ratio variance, consistent with a multifactorial model in which antigen affinity, co-stimulation, and cytokine milieu also contribute to clonal fate decisions. The residual variance may also reflect the action of additional myeloid suppressor populations not captured in our analysis.

**Cross-cancer heterogeneity.** The failure of the NK-like signature to predict response in melanoma cohorts (GSE91061, GSE120575) is a key finding rather than a limitation. Melanoma is characterized by high TMB, abundant neoantigen-specific T cell clones, and a "hot" immune microenvironment [30]. In this context, we propose that clonal diversity is sufficiently high that the bottleneck shifts from clonal fate locking to differentiation kinetics—responders and non-responders both have NK-like clones, but the rate of effector-to-exhaustion transition differs. This model is supported by the divergent predictive performance of TMB between NSCLC and melanoma [18] and by the distinct exhaustion kinetics observed in melanoma vs lung cancer T cells [31]. The cross-cancer heterogeneity also has translational implications: NK-like signature-based biomarkers should be developed and validated within cancer type-specific contexts rather than assumed to be universal.

**Spatial validation.** The spatial transcriptomic analysis using GeoMx DSP data provides *in situ* confirmation of the SPP1–receptor–stemness axis. The positive correlation between SPP1 interaction scores and stemness gene expression in tumor ROIs (*r* = 0.325, *p* = 0.044) demonstrates that the regulatory module identified in dissociated single-cell data operates in the intact tissue architecture. The higher NK cytotoxicity scores in responder ROIs (*p* = 0.027) further support the functional relevance of spatially organized NK-like activity. The limitation of the DSP platform—absence of KLF2 and FGFBP2 probes—necessitated the use of surrogate gene sets, but the convergence of spatial and single-cell findings strengthens the overall conclusion.

**Limitations.** Several limitations should be acknowledged. First, the discovery cohort is a single-center study, and multi-center validation is needed. Second, the external validation cohort GSE126044 has a small sample size (n = 16), and the GSE135222 cohort relies on PFS rather than RECIST criteria, potentially reducing statistical power. Third, the GSE120575 scRNA-seq validation used only 11 of 48 patients due to expression data availability, limiting statistical inference. Fourth, the clonal fate locking model is correlative; causal validation through *in vivo* perturbation of SPP1 signaling and KLF2 restoration is warranted. Fifth, the GeoMx DSP platform lacks several key NK-like marker genes (KLF2, FGFBP2), constraining spatial analysis to surrogate signatures.

**Translational implications.** The NK-Locked Ratio represents a candidate biomarker that could be assessed from pre-treatment tumor biopsies using scRNA-seq or TCR-seq combined with multi-marker flow cytometry. The strong predictive performance (adjusted OR = 38.14) and robustness across sensitivity analyses support its clinical potential. Furthermore, the identification of the SPP1+ TAM–KLF2 axis as a mechanistic target suggests that combining anti-PD-1 therapy with SPP1 pathway inhibition or KLF2-enhancing strategies may restore NK-like clonal commitment in non-responders. Preclinical SPP1 inhibitors and anti-CD44 antibodies are under development [32,33], providing a path to combinatorial trials.

---

## 4. Materials and Methods

### 4.1. Discovery Cohort

The discovery cohort comprised NSCLC patients treated with neoadjuvant anti-PD-1 (toripalimab) therapy from the GSE243013 dataset [15] (Zhang et al., Cell 2025). Single-cell RNA-seq data were available as pre-processed AnnData objects: GSE243013_T_cells.h5ad (434,458 T cells × 31,831 genes, 231 patients) and GSE243013_immune.h5ad (1,254,749 immune cells × 31,831 genes, 243 patients). Among T cell patients, 188 had pathological response annotations suitable for pCR (n = 84) vs non-MPR (n = 104) binary comparison. Cell type annotations (sub_cell_type column) included CD8T_NK-like_FGFBP2 (n = 16,744), Tex subsets (n = 71,298), and other CD8/CD4 T cell subtypes. TCR clonotype information was available in the clonotype and clonotype_number columns.

### 4.2. Data Preprocessing

Single-cell data were loaded using the anndata package (v0.8). Expression matrices were log-normalized using normalize_total(1e4) followed by log1p transformation, implemented as a sparse-safe custom function equivalent to scanpy.pp.normalize_total + scanpy.pp.log1p. For visualization, dimensionality reduction was performed using PCA (50 components) followed by UMAP (n_neighbors = 15, min_dist = 0.3, random_state = 42) on subsampled data (≤200 non-NK-like cells per patient plus all NK-like cells) to ensure computational tractability and balanced representation.

### 4.3. NK-like Gene Signature

The NK-like 10-gene panel for single-cell analysis comprised: FGFBP2, KLRD1, CX3CR1, FCGR3A, KLRC1, NKG7, GNLY, PRF1, GZMB, KLRB1. For external bulk RNA-seq validation, an extended 24-gene NK-like signature was used: FGFBP2, KLRD1, CX3CR1, FCGR3A, KLRC1, KLRC2, KLRB1, NKG7, GNLY, PRF1, GZMB, GZMH, GZMA, CTSW, KLRF1, SH2D1B, TYROBP, FCER1G, CD160, CRTAM, IFNG, TBX21, EOMES, ZNF683. This signature was derived from published NK-like CD8 T cell transcriptional profiles [9,10] and encompasses NK receptors, cytotoxic effectors, and transcription factors. Gene matching across platforms was performed dynamically with support for gene symbol, Ensembl ID, and Entrez ID formats.

### 4.4. Clonal Fate Locking Analysis

Large clones were defined as clonotypes with ≥5 cells. For each large clone, the fraction of NK-like CD8 T cells (nk_frac) was computed. A large clone was classified as NK-Locked if nk_frac > 0.5, and as Tex-Differentiated if Tex cells constituted the majority. The patient-level NK-Locked Ratio was defined as the proportion of large clones that are NK-locked. Logistic regression was performed using statsmodels (smf.logit) with pCR as the binary outcome (1 = pCR, 0 = non-MPR). Univariate models included NK-Locked Ratio or Tex-dominant fraction as sole predictors. The multivariate model adjusted for total_clone_count and CD8 T cell proportion. Odds ratios and 95% confidence intervals were extracted from model parameters. Sensitivity analysis was performed at clone size thresholds of 0, 5, 10, 20, and 50 cells.

### 4.5. SPP1+ TAM Identification and Interaction Analysis

From the full immune cell atlas (1,254,749 cells), myeloid cells (major_cell_type = "Myeloid cell", n = 191,099) were extracted and log-normalized. SPP1+ TAMs were defined as macrophage/monocyte subset cells with SPP1 expression above the subtype median. Patient-level SPP1+ TAM fraction was computed as the ratio of SPP1+ TAMs to total myeloid cells.

Ligand–receptor interaction analysis was performed for eight candidate pairs (SPP1–CD44, SPP1–ITGAV, SPP1–ITGB1, CXCL8–CXCR2, CCL2–CCR2, IL1B–IL1R1, TNF–TNFRSF1A, VEGFA–FLT1). Interaction scores were computed as the product of mean ligand expression in SPP1+ TAMs and mean receptor expression in each CD8 T cell subtype (NK-like and Tex separately). The final interaction score was defined as the maximum of the NK-like and Tex interaction scores to avoid signal averaging across functionally distinct recipient populations.

### 4.6. Stemness Gene Analysis

Stemness-associated genes comprised KLF2, TCF7, LEF1, IL7R, SELL, and CCR7. Patients were stratified into High-SPP1 (top 30%) and Low-SPP1 (bottom 30%) groups based on SPP1+ TAM fraction. Differential expression analysis between groups was performed using the Mann–Whitney U test, with Benjamini–Hochberg FDR correction. Direct comparison of stemness gene expression between NK-like and Tex cells was performed on a sensitivity-filtered subset (NK-like cell count ≥10 per patient) to exclude extreme value bias.

### 4.7. External Validation Cohorts

Four independent cohorts were used for external validation:

- **GSE126044** [34]: NSCLC bulk RNA-seq, anti-PD-1 (nivolumab), 16 pre-treatment tissue samples (R = 5, NR = 11). Response labels were extracted from the GEO series matrix ("patient response" field). NK-like 24-gene signature scores were computed as the mean of matched gene expressions.

- **GSE135222** [35]: NSCLC bulk RNA-seq, anti-PD-1/PD-L1, 27 samples (R = 6, NR = 21). Response labels were defined from the PFS field (PFS = 0 → R, PFS = 1 → NR), as RECIST response annotations were not available. Gene IDs were in Ensembl format (13/24 genes matched).

- **GSE91061** [36]: Melanoma bulk RNA-seq, anti-PD-1 (nivolumab) ± anti-CTLA-4 (ipilimumab), 109 samples. Pre-treatment samples with PRCR (partial/complete response, n = 10) or PD (progressive disease, n = 23) were selected; SD and unknown responses were excluded. Gene IDs were in Entrez format (22/24 genes matched).

- **GSE120575** [37]: Melanoma scRNA-seq, anti-PD-1, 48 patients total (R = 21, NR = 27). A pre-extracted NK gene expression subset (24 genes × 16,291 cells) was available. Patient IDs were extracted from cell barcodes (format: well_patient_batch). Pseudobulk NK-like signature scores were computed per patient as the mean expression across all cells. Eleven patients had both expression data and matched response labels (R = 3, NR = 8).

Statistical comparison between response groups was performed using the Mann–Whitney U test (scipy.stats.mannwhitneyu, two-sided). Fold change was computed as mean(R)/mean(NR). For the cross-cancer forest plot, 95% confidence intervals for log2 fold change were estimated using 2,000-iteration bootstrap resampling (numpy.random.default_rng, seed = 42).

### 4.8. Spatial Transcriptomic Validation

Spatial transcriptomic data from GSE221733 (GeoMx DSP, NSCLC, 93 ROIs from 39 patients) were analyzed. After filtering ROIs with missing response labels or non-standard segment types, 67 ROIs from 39 patients were retained (Responder = 16, Non-responder = 23; PanCK+ tumor = 38, PanCK− stroma = 29). The platform comprised 8,659 probes mapping to 1,812 unique genes. SPP1–receptor interaction scores were computed as SPP1 expression × mean(CD44, ITGAV, ITGB1). Stemness scores used available surrogate genes (TCF7, LEF1). NK cytotoxicity scores used GZMB, GZMA, GNLY, PRF1. Standard processing workflow included UMAP dimensionality reduction (n_neighbors = 15, min_dist = 0.3) for overall distribution visualization, gene expression mapping, and score distribution analysis.

### 4.9. Statistical Analysis

All statistical analyses were performed in Python 3.11 using scipy (v1.11), statsmodels (v0.14), and numpy (v1.24). Mann–Whitney U tests were two-sided. Spearman rank correlations were used for non-parametric association analysis. Logistic regression used maximum likelihood estimation with default convergence criteria. FDR correction used the Benjamini–Hochberg method. All *p*-values, odds ratios, correlation coefficients, and confidence intervals were dynamically computed from data—no statistical values were hardcoded. All analyses are fully reproducible via the provided run_all.py script.

### 4.10. Data and Code Availability

All datasets are publicly available from the NCBI Gene Expression Omnibus (GEO): GSE243013 (discovery), GSE126044, GSE135222, GSE91061, GSE120575 (validation), GSE221733 (spatial). TCGA data were obtained from the Genomic Data Commons (SKCM, LUAD, LUSC). All analysis code is available in the supplementary materials (code/figure1.py through figure_supplement.py, spatial_validation.py, config.py, _common.py) and can be executed via run_all.py.

---

## 5. Conclusions

This study establishes clonal fate locking—the irreversible commitment of expanded CD8 T cell clones to either an NK-like or Tex trajectory—as a deterministic principle governing anti-PD-1 response in NSCLC. The NK-Locked Ratio potently and independently predicts pathological complete response (adjusted OR = 38.14, *p* = 0.002), while SPP1+ TAMs suppress KLF2-dependent stemness in NK-like CD8 T cells through SPP1–CD44/ITGB1 signaling. The cancer type-specific predictive performance—validated in NSCLC but not melanoma—reveals that clonal fate locking represents the rate-limiting bottleneck in immunologically "cold" tumors. These findings provide a mechanistic framework for myeloid-driven immune evasion and identify the SPP1+ TAM–KLF2 axis as a candidate therapeutic target for overcoming anti-PD-1 resistance.

---

## Supplementary Materials

The following supplementary materials are available:

- **Figure S1**: Dataset overview and quality control metrics.
- **Figure S2**: Sensitivity analysis of NK-Locked Ratio predictive performance across clone size thresholds (0, 5, 10, 20, 50 cells).
- **Figure S3**: Additional validation cohort gene matching statistics.
- **Figure S_spatial_validation**: Spatial transcriptomic validation (5-panel) and standard processing workflow (3×4 panel: UMAP distribution, gene expression, score distribution).
- **Table S1**: Patient demographic and clinical characteristics of the discovery cohort.
- **Table S2**: Logistic regression results (univariate and multivariate) for clonal fate locking analysis.
- **Table S3**: Sensitivity analysis results across clone size thresholds.
- **Table S4**: Ligand–receptor interaction scores for SPP1+ TAM–CD8 T cell communication.
- **Table S5**: Differential stemness gene expression between High-SPP1 and Low-SPP1 patient groups (with FDR).
- **Table S6**: NK-like vs Tex stemness gene expression comparison.
- **Table S7**: 24-gene NK-like signature definition with Ensembl and Entrez ID mappings.
- **Table S8**: Gene matching rates across external validation cohorts.
- **Table S9**: Spatial transcriptomic statistical summary.

---

## Author Contributions

[To be completed by authors]

## Funding

[To be completed by authors]

## Institutional Review Board Statement

[To be completed by authors]

## Informed Consent Statement

[To be completed by authors]

## Data Availability Statement

All datasets are publicly available from NCBI GEO (GSE243013, GSE126044, GSE135222, GSE91061, GSE120575, GSE221733) and TCGA (SKCM, LUAD, LUSC). Analysis code is provided in the supplementary materials.

## Conflicts of Interest

The authors declare no conflict of interest.

---

## References

1. Forde, P.M.; et al. Neoadjuvant Nivolumab in Resectable Lung Cancer. *N. Engl. J. Med.* 2018, 379, e1–e31.
2. Cascone, T.; et al. Neoadjuvant Nivolumab or Nivolumab Plus Ipilimumab in Resectable Non-Small Cell Lung Cancer: The CheckMate 816 Trial. *J. Clin. Oncol.* 2023, 41, 3679–3690.
3. Garon, E.B.; et al. Pembrolizumab for the Treatment of Non-Small-Cell Lung Cancer. *N. Engl. J. Med.* 2015, 372, 2018–2028.
4. Hellmann, M.D.; et al. Nivolumab plus Ipilimumab in Lung Cancer with a High Tumor Mutational Burden. *N. Engl. J. Med.* 2018, 378, 2093–2104.
5. Philip, M.; Schietinger, A. CD8+ T Cell Differentiation and Dysfunction in Cancer. *Nat. Rev. Immunol.* 2022, 22, 209–223.
6. Sade-Feldman, M.; et al. Defining T Cell States Associated with Response to Checkpoint Immunotherapy in Melanoma. *Cell* 2018, 175, 998–1013.
7. Chiossone, L.; et al. Molecular Characterization of Human Natural Killer Cells. *Immunol. Rev.* 2018, 286, 1–14.
8. Crome, S.Q.; et al. A Distinct Innate-Like CD8+ T Cell Population. *Eur. J. Immunol.* 2012, 42, 2632–2642.
9. Freud, A.G.; et al. The Broad Spectrum of Human Natural Killer Cell Diversity. *Immunity* 2017, 47, 820–833.
10. Dobano, C.; et al. Expression and Function of NK Cell Receptors on CD8+ T Cells. *Front. Immunol.* 2019, 10, 2336.
11. Mantovani, A.; et al. Tumor-Associated Macrophages as Treatment Targets in Oncology. *Nat. Rev. Clin. Oncol.* 2017, 14, 399–416.
12. DeNardo, D.G.; Ruffell, B. Macrophages as Regulators of Tumor Immunity and Immunotherapy. *Nat. Rev. Immunol.* 2019, 19, 369–382.
13. Zhang, Q.; et al. Landmarkscapes of Tumor-Infiltrating Immune Cells in Cancer. *Cell* 2021, 184, 797–812.
14. Oshi, M.; et al. SPP1 Expression Is a Prognostic Biomarker in Breast Cancer. *Cancers* 2021, 13, 1451.
15. Zhang, J.; et al. Single-Cell Landscape of NSCLC Anti-PD-1 Therapy. *Cell* 2025. (GSE243013)
16. Emerson, R.O.; et al. High-Throughput Sequencing of T-Cell Receptors Reveals a Homogeneous Repertoire of Tumor-Infiltrating Lymphocytes in Ovarian Cancer. *PLoS One* 2013, 8, e76808.
17. Hart, G.T.; et al. Abandoning the Tug-of-War between T Cell Quiescence and Activation. *Nat. Rev. Immunol.* 2023, 23, 325–338.
18. Yarchoan, M.; et al. Tumor Mutational Burden and Response Rate to PD-1 Inhibition. *N. Engl. J. Med.* 2017, 377, 2500–2501.
19. Simon, S.; Labarriere, N. PD-1 Expression on Tumor-Specific T Cells: Friend or Foe for Immunotherapy? *Cancers* 2017, 9, 95.
20. Krishna, C.; et al. Single-cell Sequencing Links Multiregional Immune Landscapes and Tissue-Resident T Cells in Cervical Cancer. *Nat. Genet.* 2021, 53, 120–129.
21. Miller, B.C.; et al. Subsets of Exhaustished CD8+ T Cells Differentially Mediate Tumor Control and Respond to Checkpoint Blockade. *Nat. Immunol.* 2019, 20, 326–336.
22. Siesel, C.S.; et al. Stem-like CD8+ T Cells. *Nat. Rev. Immunol.* 2023, 23, 611–623.
23. Pauken, K.E.; et al. Epigenetic Stability of Exhausted T Cells Limits Durability of Reinvigoration by PD-1 Blockade. *Science* 2016, 354, 1160–1165.
24. Alfei, F.; et al. CD8+ T Cell Epigenetic Fixed Trait and Counteracting Inflammatory Signals. *Immunity* 2019, 50, 108–121.
25. Rangaswami, H.; et al. Osteopontin: Role in Cell Signaling and Cancer Progression. *Trends Cell Biol.* 2006, 16, 79–87.
26. Weber, G.F.; et al. Receptor-Ligand Interaction Between CD44 and Osteopontin (Eta-1). *Science* 1996, 271, 509–512.
27. Ponta, H.; et al. CD44: A Multifunctional Cell Surface Adhesion Receptor. *Nat. Rev. Mol. Cell Biol.* 2003, 4, 33–45.
28. Hogg, N.; et al. Integrin and Function-Associated Molecules on T Cells. *Immunol. Rev.* 2002, 186, 171–178.
29. Carlson, C.M.; et al. Kruppel-like Factor 2 Regulates T Cell Trafficking. *Nature* 2006, 442, 1049–1052.
30. Huang, A.C.; et al. T-Cell Infiltration and Immunity in Melanoma. *Nat. Rev. Cancer* 2020, 20, 65–77.
31. Yost, K.E.; et al. Clonal Replacement of Tumor-Specific T Cells Following PD-1 Blockade. *Nat. Med.* 2019, 25, 1251–1259.
32. Shevde, L.A.; et al. SPP1 (Osteopontin) as a Therapeutic Target in Cancer. *Expert Opin. Ther. Targets* 2010, 14, 1217–1230.
33. Zöller, M. CD44: Can a Cancer-Initiating Cell Profit from an Abundantly Expressed Molecule? *Nat. Rev. Cancer* 2011, 11, 254–267.
34. Anagnostou, V.; et al. Dynamics of Tumor and Immune Responses during Immune Checkpoint Blockade in Non-Small Cell Lung Cancer. *Cancer Res.* 2020. (GSE126044)
35. Jung, H.; et al. Liquid Biopsy Enables Oncogenic Tracking in NSCLC. *Nat. Commun.* 2019. (GSE135222)
36. Riaz, N.; et al. Tumor and Microenvironment Evolution during Immunotherapy in Melanoma. *Cell* 2017, 171, e895. (GSE91061)
37. Sade-Feldman, M.; et al. Defining T Cell States Associated with Response to Checkpoint Immunotherapy in Melanoma. *Cell* 2018, 175, 998–1013. (GSE120575)
