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

path <- paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/Coloc/coloc/GWAS/")
trait <- list.dirs(path, full.names = TRUE, recursive = FALSE)
trait <- sapply(trait, function(x) unlist(strsplit(x, "\\/"))[11])
trait <- data.frame(trait)
trait <- trait$trait

all_OR <- data.frame()
for (k in 1:length(trait)) {
  snp <- fread(paste0(PROJECT_DIR, "/CattleGTEx/GWAS/America_gwas/", trait[[k]], ".gwa_All.txt"))
  names(snp)[1] <- "chr"
  finemap <- fread(paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/Coloc/coloc/GWAS/", trait[[k]], "/", trait[[k]], "_coloc.bed"))
  finemap$chr <- sub("_.*","",finemap$variant_id)
  finemap$position <- sub(".*_","",finemap$variant_id)
  finemap <- finemap[!duplicated(finemap$variant_id),]
  
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
      chr_snp <- snp[snp$chr == j,]
      chr_snp <- chr_snp[!duplicated(chr_snp$ID),]
      dt_snp <- data.table(position = as.numeric(chr_snp$POS))
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
  all_OR1$trait_name <- trait[[k]]
  all_OR <- rbind(all_OR, all_OR1)
}

# 动态输出保护
out_path_1 <- paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/eQTL_selection/Global_pattern/GWAS")
dir.create(out_path_1, recursive = TRUE, showWarnings = FALSE)
write.csv(all_OR, paste0(out_path_1, "/gwas_OR_fst.csv"))
