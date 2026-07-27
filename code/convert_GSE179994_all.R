# convert_GSE179994_all.R
# 将 GSE179994_all.Tcell.rawCounts.rds.gz 转换为 Python 可读格式
# 只提取 NK-like 24 基因签名的表达数据，避免导出全量矩阵
suppressPackageStartupMessages({
  library(Matrix)
})

cat("[R] Loading GSE179994 all T cell raw counts...\n")
rds_path <- "adata/GSE179994_all.Tcell.rawCounts.rds.gz"
counts <- readRDS(gzfile(rds_path))

cat("[R] Object class:", class(counts), "\n")
cat("[R] Object type:", typeof(counts), "\n")

if (inherits(counts, "dgCMatrix") || inherits(counts, "Matrix")) {
  cat("[R] Sparse matrix confirmed\n")
} else if (is.matrix(counts)) {
  cat("[R] Dense matrix\n")
} else if (inherits(counts, "list")) {
  cat("[R] List, names:\n")
  print(names(counts))
  for (nm in names(counts)) {
    cat(sprintf("  %s: %s\n", nm, class(counts[[nm]])))
  }
}

# 检查维度
cat("[R] Dimensions:", nrow(counts), "genes x", ncol(counts), "cells\n")

# 检查行名（基因名）和列名（细胞名）
cat("[R] First 10 gene names:\n")
print(head(rownames(counts), 10))
cat("[R] First 10 cell names:\n")
print(head(colnames(counts), 10))

# 分析细胞名格式，提取患者ID
cell_names <- colnames(counts)
cat("[R] Sample cell name patterns:\n")
for (i in 1:min(20, length(cell_names))) {
  cat(sprintf("  %s\n", cell_names[i]))
}

# NK-like 24 基因签名
sig_genes <- c("FGFBP2","KLRD1","CX3CR1","FCGR3A","KLRC1","KLRC2","KLRB1","NKG7",
               "GNLY","PRF1","GZMB","GZMH","GZMA","CTSW","KLRF1","SH2D1B","TYROBP",
               "FCER1G","CD160","CRTAM","IFNG","TBX21","EOMES","ZNF683")

sig_found <- intersect(sig_genes, rownames(counts))
cat(sprintf("\n[R] Signature genes found: %d/%d\n", length(sig_found), length(sig_genes)))
cat("[R] Missing genes:", setdiff(sig_genes, sig_found), "\n")

if (length(sig_found) > 0) {
  # 提取签名基因表达矩阵
  sig_expr <- counts[sig_found, , drop = FALSE]
  cat(sprintf("[R] Signature expression matrix: %d genes x %d cells\n", nrow(sig_expr), ncol(sig_expr)))

  # 计算每列的签名得分（log1p 均值）
  cat("[R] Calculating per-cell NK-like signature scores...\n")
  # 对稀疏矩阵使用列均值运算，然后 log1p 均值
  # 先计算每列均值，再 log1p - 但这与 log1p 后均值不同
  # 正确做法：colMeans(log1p(x)) 对稀疏矩阵优化
  log_sig_expr <- log1p(sig_expr)
  sig_scores <- colMeans(log_sig_expr)
  names(sig_scores) <- colnames(sig_expr)
  cat("[R] Signature scores calculated\n")

  # 保存签名得分
  cat("[R] Writing signature scores...\n")
  write.csv(data.frame(cell = names(sig_scores), score = sig_scores),
            file = "adata/GSE179994/GSE179994_all_signature_scores.csv",
            row.names = FALSE)
  cat("[R] Signature scores saved\n")

  # 从细胞名提取患者ID和时间点
  cat("[R] Parsing cell names for patient info...\n")
  patient_ids <- character(length(cell_names))
  timepoints <- character(length(cell_names))
  sample_types <- character(length(cell_names))

  for (i in seq_along(cell_names)) {
    cn <- cell_names[i]
    # 格式如: P19.blood.post.AAACCTGAGATGTTAG-1 或 P1.ut.AAACCTGAGATGTTAG-1
    parts <- strsplit(cn, "\\.")[[1]]
    if (length(parts) >= 2) {
      patient_ids[i] <- parts[1]
      # 判断时间点
      if (any(grepl("pre|ut|untreated", parts, ignore.case = TRUE))) {
        timepoints[i] <- "pre"
      } else if (any(grepl("post|on", parts, ignore.case = TRUE))) {
        timepoints[i] <- "post"
      } else {
        timepoints[i] <- "unknown"
      }
      # 判断样本类型
      if (any(grepl("blood|pbmc", parts, ignore.case = TRUE))) {
        sample_types[i] <- "blood"
      } else if (any(grepl("tumor|tumour", parts, ignore.case = TRUE))) {
        sample_types[i] <- "tumor"
      } else {
        sample_types[i] <- "unknown"
      }
    } else {
      patient_ids[i] <- "unknown"
      timepoints[i] <- "unknown"
      sample_types[i] <- "unknown"
    }
  }

  cat("[R] Unique patients:", length(unique(patient_ids)), "\n")
  cat("[R] Patient list:", sort(unique(patient_ids)), "\n")
  cat("[R] Timepoint distribution:\n")
  print(table(timepoints))
  cat("[R] Sample type distribution:\n")
  print(table(sample_types))

  # 保存细胞元信息
  cell_meta <- data.frame(
    cell = cell_names,
    patient = patient_ids,
    timepoint = timepoints,
    sample_type = sample_types,
    nklike_score = sig_scores,
    stringsAsFactors = FALSE
  )
  write.csv(cell_meta, file = "adata/GSE179994/GSE179994_all_cell_metadata.csv", row.names = FALSE)
  cat("[R] Cell metadata saved\n")

  # 按患者 × 时间点计算汇总
  cat("[R] Aggregating by patient x timepoint...\n")
  df <- data.frame(
    patient = patient_ids,
    timepoint = timepoints,
    sample_type = sample_types,
    score = sig_scores
  )

  # 患者水平统计
  pat_stats <- aggregate(score ~ patient + timepoint + sample_type, data = df,
                         FUN = function(x) c(mean = mean(x), median = median(x), sd = sd(x), n = length(x)))
  pat_df <- do.call(data.frame, pat_stats)
  colnames(pat_df) <- c("patient", "timepoint", "sample_type", "mean_score", "median_score", "sd_score", "n_cells")

  cat("[R] Patient x timepoint groups:", nrow(pat_df), "\n")
  print(head(pat_df, 20))

  write.csv(pat_df, file = "adata/GSE179994/GSE179994_patient_timepoint_scores.csv", row.names = FALSE)
  cat("[R] Patient-timepoint scores saved\n")

  cat("\n[R] All outputs saved successfully.\n")
} else {
  cat("[R] ERROR: No signature genes found!\n")
  quit(status = 1)
}
