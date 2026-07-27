# convert_GSE179994.R
# 将 GSE179994 的 Seurat RDS 对象转换为 Python 可读格式
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

cat("[R] Loading GSE179994 Seurat object...\n")
rds_path <- "adata/GSE179994/GSM5444629_P019.blood.post.Tcell.Seurat.object.rds.gz"
obj <- readRDS(gzfile(rds_path))

cat("[R] Object class:", class(obj), "\n")
cat("[R] Object type:", typeof(obj), "\n")

# 检查 Seurat 对象结构
if (inherits(obj, "Seurat")) {
  cat("[R] Seurat object loaded\n")
  cat("[R] Assays:", names(obj@assays), "\n")
  cat("[R] Default assay:", DefaultAssay(obj), "\n")
  cat("[R] Number of cells:", ncol(obj), "\n")
  cat("[R] Number of features:", nrow(obj), "\n")
  cat("[R] Reductions:", names(obj@reductions), "\n")

  # 检查 metadata
  cat("\n[R] Metadata columns:\n")
  print(colnames(obj@meta.data))
  cat("\n[R] Metadata head:\n")
  print(head(obj@meta.data, 3))

  # 检查 metadata 中的患者/样本信息
  meta_cols <- colnames(obj@meta.data)
  for (col in meta_cols) {
    if (grepl("patient|sample|orig|group|response|treatment", col, ignore.case = TRUE)) {
      cat(sprintf("\n[R] Column '%s' unique values (%d):\n", col, length(unique(obj@meta.data[[col]]))))
      print(head(unique(obj@meta.data[[col]]), 20))
    }
  }

  # 提取表达矩阵（counts）
  cat("\n[R] Extracting counts matrix...\n")
  counts <- GetAssayData(obj, assay = DefaultAssay(obj), layer = "counts")
  cat("[R] Counts matrix:", nrow(counts), "genes x", ncol(counts), "cells\n")
  cat("[R] Counts class:", class(counts), "\n")
  cat("[R] Counts sparse:", inherits(counts, "dgCMatrix"), "\n")

  # 检查细胞名（barcodes）
  cat("\n[R] First 5 cell names:\n")
  print(head(colnames(counts), 5))

  # 检查基因名
  cat("\n[R] First 5 gene names:\n")
  print(head(rownames(counts), 5))

  # 导出为 CSV（仅签名基因 + metadata）
  # NK-like 24 基因签名
  sig_genes <- c("FGFBP2","KLRD1","CX3CR1","FCGR3A","KLRC1","KLRC2","KLRB1","NKG7",
                 "GNLY","PRF1","GZMB","GZMH","GZMA","CTSW","KLRF1","SH2D1B","TYROBP",
                 "FCER1G","CD160","CRTAM","IFNG","TBX21","EOMES","ZNF683")

  sig_found <- intersect(sig_genes, rownames(counts))
  cat(sprintf("\n[R] Signature genes found: %d/%d\n", length(sig_found), length(sig_genes)))

  if (length(sig_found) > 0) {
    # 导出签名基因表达矩阵
    sig_expr <- counts[sig_found, , drop = FALSE]
    cat("[R] Writing signature gene expression matrix...\n")
    write.csv(as.matrix(sig_expr), "adata/GSE179994/GSE179994_P019_signature_expr.csv")
    cat("[R] Signature expression saved: adata/GSE179994/GSE179994_P019_signature_expr.csv\n")
  }

  # 导出完整 metadata
  cat("[R] Writing metadata...\n")
  write.csv(obj@meta.data, "adata/GSE179994/GSE179994_P019_metadata.csv")
  cat("[R] Metadata saved: adata/GSE179994/GSE179994_P019_metadata.csv\n")

  # 导出完整 counts 矩阵为 Matrix Market 格式（更高效）
  cat("[R] Writing full counts matrix (Matrix Market format)...\n")
  dir.create("adata/GSE179994/GSE179994_P019_mtx", showWarnings = FALSE)
  writeMM(counts, "adata/GSE179994/GSE179994_P019_mtx/matrix.mtx")
  writeLines(rownames(counts), "adata/GSE179994/GSE179994_P019_mtx/genes.tsv")
  writeLines(colnames(counts), "adata/GSE179994/GSE179994_P019_mtx/barcodes.tsv")
  gzip_cmd <- "gzip -f adata/GSE179994/GSE179994_P019_mtx/matrix.mtx"
  system(gzip_cmd)
  cat("[R] Full counts matrix saved: adata/GSE179994/GSE179994_P019_mtx/\n")

  cat("\n[R] All outputs saved successfully.\n")

} else if (is.list(obj)) {
  # 可能是旧版 Seurat 对象或其他格式
  cat("[R] Object is a list, names:\n")
  print(names(obj))

  # 尝试提取表达矩阵
  if ("counts" %in% names(obj)) {
    counts <- obj$counts
    cat("[R] Found 'counts':", nrow(counts), "x", ncol(counts), "\n")
  } else if ("data" %in% names(obj)) {
    counts <- obj$data
    cat("[R] Found 'data':", nrow(counts), "x", ncol(counts), "\n")
  }
} else {
  cat("[R] Unknown object format. Trying as.matrix...\n")
  mat <- as.matrix(obj)
  cat("[R] Matrix:", nrow(mat), "x", ncol(mat), "\n")
  print(head(rownames(mat), 5))
  print(head(colnames(mat), 5))
}
