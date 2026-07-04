options(stringsAsFactors = FALSE)
library(edgeR)
library(preprocessCore)
library(RNOmni)
library(data.table)
library(R.utils)
library(SNPRelate)
library(dplyr)
library(coloc)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(DESeq2)
library(clusterProfiler)
library(org.Bt.eg.db)
library(enrichplot)

# ====================================================================
# 【原生参数接管】抓取外部传入的项目根目录变量
# ====================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("错误: 必须提供项目根目录路径 (PROJECT_DIR)！", call. = FALSE)
}
PROJECT_DIR <- args[1]
# ====================================================================

# 路径解耦：使用 PROJECT_DIR 动态外延
setwd(paste0(PROJECT_DIR, "/CattleGTEx/Selection_region/share/home/zju_zhaopj/01-PanCattle-RNA/99-HC/02-Merge"))
all_fst <- data.frame()
all_XPclr <- data.frame()
all_XPEHH <- data.frame()
for (chr in 1:29) {
    fst <- fread(paste0(chr, ".Fst.windowed.txt"))
    fst <- na.omit(fst)
    fst$chrm <- paste0("chr", chr)
    all_fst <- rbind(all_fst, fst)
}
all_fst$decile <- ntile(all_fst$WEIGHTED_FST, 10) 


path <- paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/eQTL/")
tissue <- list.dirs(path, full.names = TRUE, recursive = FALSE)
tissue <- sapply(tissue, function(x) unlist(strsplit(x, "\\/"))[9])
tissue <- data.frame(tissue)
tissue <- tissue$tissue

setwd(paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/MASHR_Celltype/Tissue_sharing/Bulk/Final_output1/output_top_paris_across_all"))
lfsr_threshold <- 0.05
lfsr_result  <- readRDS("lfsr_m.s.RDS")
significant_matrix <- lfsr_result < lfsr_threshold
num_significant_ct <- rowSums(significant_matrix)
specific_eqtl_matrix <- significant_matrix[num_significant_ct < 3, ]
specific_eqtl_matrix <- data.frame(specific_eqtl_matrix)
sharing_eqtl_matrix <- significant_matrix[num_significant_ct > 2, ]
sharing_eqtl_matrix <- data.frame(sharing_eqtl_matrix)

all_OR <- data.frame()
for (k in 1:length(tissue)) {
  finemap <- fread(paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/SusieR/bulk_eQTL/", tissue[[k]], "_finemapping.txt"))
  finemap$chr <- sub("_.*","",finemap$SNP)
  finemap$position <- sub(".*_","",finemap$SNP)
  finemap <- finemap[!duplicated(finemap$SNP),]
  tissue_specific <- specific_eqtl_matrix[,colnames(specific_eqtl_matrix) == tissue[[k]]]
  tissue_specific <- data.frame(tissue_specific)
  rownames(tissue_specific) <- rownames(specific_eqtl_matrix)
  tissue_specific$gene <- sub(".*,", "", rownames(tissue_specific))
  
  tissue_sharing <- sharing_eqtl_matrix[,colnames(sharing_eqtl_matrix) == tissue[[k]]]
  tissue_sharing <- data.frame(tissue_sharing)
  rownames(tissue_sharing) <- rownames(sharing_eqtl_matrix)
  tissue_sharing$gene <- sub(".*,", "", rownames(tissue_sharing))
  
  specific_gene <- tissue_specific$gene[tissue_specific$tissue_specific == "TRUE"]
  sharing_gene <- tissue_sharing$gene[tissue_sharing$tissue_sharing == "TRUE"]
  select_gene <- specific_gene[!specific_gene %in% sharing_gene]
  
  finemap <- finemap[finemap$gene %in% select_gene,]
  
  b = nrow(finemap)
  all_OR1 <- vector()
  for (i in 1:10) {
    select_fst <-all_fst[all_fst$decile == i,]
    chrm <- unique(all_fst$CHROM)
    a = 0
    d = 0
    c = 0
    for (j in chrm) {
      chr_fst <- select_fst[select_fst$CHROM == j,]
      chr_finemap <- finemap[finemap$chr == j,]
      chr_snp <- fread(paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/eQTL/", tissue[[k]], "/", tissue[[k]], "_new_LMM.cis_qtl_pairs.", j, ".txt.gz"))
      chr_snp$position <- sub(".*_","",chr_snp$variant_id)
      chr_snp <- chr_snp[!duplicated(chr_snp$variant_id),]
      dt_snp <- data.table(position = as.numeric(chr_snp$position))
      dt_snp[, `:=`(start = position, end = position)]
      dt_fst <- data.table(BIN_START = as.numeric(chr_fst$BIN_START), 
                           BIN_END = as.numeric(chr_fst$BIN_END))
      setkey(dt_fst, BIN_START, BIN_END)
      result <- foverlaps(dt_snp, dt_fst, by.x = c("start", "end"), type = "within", nomatch = 0L)
      all_snp_in_selection <- length(unique(result$position))
      
      if (nrow(chr_finemap) > 0) {
        dt_snp <- data.table(position = as.numeric(chr_finemap$position))
        dt_snp[, `:=`(start = position, end = position)]
        dt_fst <- data.table(BIN_START = as.numeric(chr_fst$BIN_START), 
                             BIN_END = as.numeric(chr_fst$BIN_END))
        setkey(dt_fst, BIN_START, BIN_END)
        result <- foverlaps(dt_snp, dt_fst, by.x = c("start", "end"), type = "within", nomatch = 0L)
        finemap_snp_in_selection <- length(unique(result$position))
        a = a + finemap_snp_in_selection
      }
      c = c + all_snp_in_selection
      d = d +nrow(chr_snp)
    }
    OR = (a/c)/(b/d)
    all_OR1 <- c(all_OR1, OR)
  }
  all_OR1 <- data.frame(all_OR1)
  all_OR1$tissue_name <- tissue[[k]]
  all_OR <- rbind(all_OR, all_OR1)
}

# 动态输出保护
out_path_2 <- paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/eQTL_selection/bulk")
dir.create(out_path_2, recursive = TRUE, showWarnings = FALSE)
write.csv(all_OR, paste0(out_path_2, "/tissue_specific_OR_fst1-2.csv"))
