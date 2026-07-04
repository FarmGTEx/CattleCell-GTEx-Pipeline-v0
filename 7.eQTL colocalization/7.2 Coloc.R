options(stringsAsFactors = FALSE)
library(edgeR)
library(preprocessCore)
library(RNOmni)
library(data.table)
library(R.utils)
library(SNPRelate)
library(dplyr)
library(coloc)

# ====================================================================
# 抓取外部传入的项目根目录变量
# ====================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("错误: 必须提供项目根目录路径 (PROJECT_DIR)！", call. = FALSE)
}
PROJECT_DIR <- args[1] # 对应 Bash 中的 $BASE_PROJECT 路径
# ====================================================================

# 路径修改：使用 PROJECT_DIR 动态外延寻址
path <- paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/eQTL/")
tissue <- list.dirs(path, full.names = TRUE, recursive = FALSE)
tissue <- sapply(tissue, function(x) unlist(strsplit(x, "\\/"))[9])
tissue <- data.frame(tissue)
tissue <- tissue$tissue

path<-NULL
for(i in tissue){
  path[[i]]<-paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/eQTL/", i, "/ieQTL/")
}

# 保留特定循环起点 (从第3个组织开始)
for (i in 3:length(tissue)) {
    setwd(path[[i]])
    ct <- list.dirs(full.names = TRUE, recursive = FALSE)
    ct <- sapply(ct, function(x) unlist(strsplit(x, "\\/"))[2])
    ct <- data.frame(ct)
    ct <- ct$ct
    
    # 路径修改：使用 PROJECT_DIR 动态指定上游 7.1 生成的 Bulk 输入
    qtl1 <- fread(paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/Coloc/coloc/Bulk/", tissue[[i]], "_coloc.bed"))
    for (j in 1:length(ct)) {
      # 路径修改：使用 PROJECT_DIR 动态指定上游 7.1 生成的 Cell Components 输入
      comp_bed_path <- paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/Coloc/coloc/all_cell_components/", tissue[[i]], "_", ct[[j]], "_coloc.bed")
      if (file.exists(comp_bed)) {
        qtl2 <- fread(comp_bed)
        gene1 <- unique(qtl1$gene_id)
        gene2 <- unique(qtl2$gene_id)
        gene_list <- intersect(gene1,gene2)
        all_summary <- data.frame()
        if (length(gene_list) > 0){
          for (k in 1:length(gene_list)) {
            gene_qtl1 <- qtl1[qtl1$gene_id == gene_list[[k]], ]
            gene_qtl2 <- qtl2[qtl2$gene_id == gene_list[[k]], ]
            input <- merge(gene_qtl1, gene_qtl2, by="rs_id", all=FALSE, suffixes=c("_bulk_eqtl","_celltype_iqtl"))
            if (nrow(input) > 0){
              input <- input[!duplicated(input$rs_id),]
              input$varbeta1 <- input$slope_se_bulk_eqtl ^ 2
              input$varbeta2 <- input$slope_se_celltype_iqtl ^ 2
              #input$beta <- abs(input$beta)
              #input$slope <- abs(input$slope)
              res <- coloc.abf(dataset1=list(snp = input$rs_id, pvalues=input$pval_nominal_bulk_eqtl, beta=input$slope_bulk_eqtl, varbeta=input$varbeta1, N=input$N_bulk_eqtl, type="quant", MAF=input$maf_bulk_eqtl),
                               dataset2=list(snp = input$rs_id, pvalues=input$pval_nominal_celltype_iqtl, beta=input$slope_celltype_iqtl, varbeta=input$varbeta2, N=input$N_celltype_iqtl, type="quant", MAF=input$maf_celltype_iqtl))
              summary <- data.frame(res$summary)
              summary$gene_id <- gene_list[[k]]
              all_summary <- rbind(all_summary, summary)
            }
          }
        }
        
        # 路径修改：输出目标目录保护并使用 PROJECT_DIR 动态切换写入
        out_dir_res <- paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/Coloc/coloc/Results/Bulk_components/")
        dir.create(out_dir_2, recursive = TRUE, showWarnings = FALSE)
        setwd(out_dir_2)
        write.csv(all_summary, paste0(tissue[[i]], "_", ct[[j]], "_coloc.csv"))
      }
    }   
}
