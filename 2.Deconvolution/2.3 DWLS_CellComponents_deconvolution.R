suppressMessages(library(energy))
suppressMessages(library(dplyr))
suppressMessages(library(SeuratObject))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(remotes))
suppressMessages(library(devtools))
suppressMessages(library(Biobase))
suppressMessages(library(Seurat))
suppressMessages(library(future.apply))
suppressMessages(library(DWLS))
suppressMessages(library(data.table))

# === 参数解析封装修改：接收来自 Bash 脚本的输入输出参数 ===
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("错误: 必须提供输入 RDS 路径和输出目录路径！", call. = FALSE)
}
input_rds <- args[1]   # 对应 Bash 中的 $DEFAULT_REF_RDS (单细胞参考集)
output_dir <- args[2]  # 对应 Bash 中的 $DEFAULT_OUTPUT_DIR (统一输出根目录，承接前文输出)
# ========================================================

# 修改你原本代码中的硬编码读取和写入位置：
sc <- readRDS(input_rds)

##prepare single cell reference
# sc <- readRDS(list1[[i]])
Idents(sc) <- "CellType"
celltype <- unique(sc$CellType)
if ("Unknown cells" %in% celltype == TRUE) {
  sc <- subset(sc, idents = "Unknown cells", invert = TRUE)   
}
counts <- data.frame(sc@assays$RNA@counts)
meta <- data.frame(sc@meta.data)
phenoData <- data.frame(cbind(meta$CellType, meta$sample))
rownames(phenoData) <- rownames(meta)
names(phenoData) <- c("CellType","sample")

# === 路径修改：将 MAST 签名矩阵临时存储路径指定到统一输出目录下 ===
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
Signature <- buildSignatureMatrixMAST(scdata = counts, id = phenoData[,"CellType"], path = output_dir, diff.cutoff = 0.5, pval.cutoff = 0.01)
# ====================================================================

##input bulk data (tpm)
# === 路径修改：动态流转，读取 2.1 模块在 Uniform_Fraction 目录下生成的伪大样 pseudo.txt ===
pseudo_tpm_path <- file.path(output_dir, "Uniform_Fraction", "random_pseudo", "pseudo.txt")
tpm <- fread(pseudo_tpm_path)
# ====================================================================================
tpm <- na.omit(tpm)

##run DWLS
res <- data.frame()
for (j in 1:ncol(tpm)) {
  b = setNames(tpm[,j], rownames(tpm))
  tr <- trimData(Signature, b)
  RES <- data.frame(t(solveDampenedWLS(tr$sig, tr$tpm)))
  res <- rbind(res,RES)
}
rownames(res) <- colnames(tpm)
res[res < 10^-5] <- 0

# === 路径修改：使用 Bash 传入的 output_dir 动态创建结果子目录并保存预测比例结果 ===
dwls_out_dir <- file.path(output_dir, "Results_DWLS", "Results")
dir.create(dwls_out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(res, file.path(dwls_out_dir, "Predict_DWLS.csv"))
# ===============================================================================
