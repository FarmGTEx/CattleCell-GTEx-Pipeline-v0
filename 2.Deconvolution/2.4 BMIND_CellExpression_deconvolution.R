library(MIND)
library(dplyr)
library(Seurat)
library(pheatmap)
library(parallel)
library(SeuratObject)
library(tibble)

# === 参数解析封装修改：引入外部参数，但保留脚本内部原有的批量循环逻辑 ===
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("错误: 必须提供项目根目录路径 (PROJECT_DIR)！", call. = FALSE)
}
PROJECT_DIR <- args[1] # 对应 Bash 中的 $BASE_PROJECT 路径
# ====================================================================

# === 路径修改：将写死的路径前缀替换为由参数传入的 PROJECT_DIR ===
path <- file.path(PROJECT_DIR, "Deconvolution/Cattle/Bulk/Results_DWLS/")
# ============================================================
tissue <- list.dirs(path, full.names = TRUE, recursive = FALSE)
tissue <- sapply(tissue, function(x) unlist(strsplit(x, "\\/"))[10])
tissue <- data.frame(tissue)
tissue_sc <- tissue$tissue
#bulk
# === 路径修改：将写死的路径前缀替换为由参数传入的 PROJECT_DIR ===
path <- file.path(PROJECT_DIR, "CattleGTEx/OmiGA/eQTL/")
# ============================================================
tissue <- list.dirs(path, full.names = TRUE, recursive = FALSE)
tissue <- sapply(tissue, function(x) unlist(strsplit(x, "\\/"))[9])
tissue <- data.frame(tissue)
tissue_bulk <- tissue$tissue

tissue <- intersect(tissue_sc, tissue_bulk)

#bulk
list1<-NULL
for(i in tissue){
  # === 路径修改：使用 PROJECT_DIR 替换硬编码 ===
  list1[[i]]<-paste0(PROJECT_DIR, "/Global_atlas/Global_cell_atlas/celltype_annotation/", i, "_anno.rds")
}
list2<-NULL
for(i in tissue){
  # === 路径修改：使用 PROJECT_DIR 替换硬编码 ===
  list2[[i]]<-paste0(PROJECT_DIR, "/CattleGTEx/Cattle_bulk_exp/Bulk_", i, ".txt")
}
list3<-NULL
for(i in tissue){
  # === 路径修改：使用 PROJECT_DIR 替换硬编码 ===
  list3[[i]]<-paste0(PROJECT_DIR, "/Deconvolution/Cattle/Bulk/Results_DWLS/", i, "/Results/Predict_", i, "_DWLS.csv")
}
list4<-NULL
for(i in tissue){
  # === 路径修改：使用 PROJECT_DIR 替换硬编码 ===
  list4[[i]]<-paste0(PROJECT_DIR, "/Deconvolution/Cattle/Bulk/Result_bMIND/", i, "2_bMIND.RData")
}

# === 路径修改：使用 PROJECT_DIR 替换硬编码 ===
gene_id <- read.csv(file.path(PROJECT_DIR, "reference/cattle_genes.csv"), row.names = 1)
# ========================================================================
gene_id <- na.omit(gene_id)

for (k in 1:length(tissue)) {
  bulk <- read.csv(list2[[k]], sep = "\t", row.names = 1)
  bulk <- na.omit(bulk)
  bulk <- bulk[! rowSums(bulk) == 0,]
  frac <- read.csv(list3[[k]], row.names = 1)
  colnames(bulk) <- gsub("_","",colnames(bulk))
  rownames(frac) <- gsub("_","",rownames(frac))
  rownames(frac) <- gsub("-","\\.",rownames(frac))
  bulk <- bulk[,colnames(bulk) %in= rownames(frac)]
  frac <- frac[rownames(frac) %in% colnames(bulk),]
  bulk <- log2(bulk + 1)
  ct <- data.frame(colMeans(frac))
  CT <- rownames(ct)[ct$colMeans.frac. > 0.05]
  
  frac <- frac[, colnames(frac) %in% CT]
  frac <- data.frame(frac[,order(colnames(frac))])
  CT1 <- colnames(frac)
  
  sc <- readRDS(list1[[k]])
  Idents(sc) <- "CellType"
  counts <- data.frame(sc@assays$RNA@counts)
  counts <- t(counts)
  counts <- data.frame(counts)
  meta <- data.frame(sc@meta.data)
  meta <- droplevels(meta)
  counts <- counts %>% mutate(cluster = meta$CellType)
  group <- list()
  group <- split(counts,counts$cluster)
  name<- rownames(table(counts$cluster))
  fenleimean <- colMeans(group[[1]][,-ncol(group[[1]])])
  for (j in 2:length(group)) {
    flmean<- colMeans(group[[j]][,-ncol(group[[j]])])
    fenleimean <-rbind(fenleimean,flmean)
  }
  rownames(fenleimean) <- name
  fenleimean = t(fenleimean)
  counts <- t(counts[,-ncol(counts)])
  
  fenleimean1 <- fenleimean[grep("ENSBTA", rownames(fenleimean)), ]
  fenleimean1 <- data.frame(fenleimean1)
  fenleimean2 <- fenleimean[grep("ENSBTA", rownames(fenleimean), invert = TRUE), ]
  fenleimean2 <- data.frame(fenleimean2)
  
  fenleimean2 <- fenleimean2 %>%
  rownames_to_column(var = "gtf_df.gene_name")
  fenleimean2 <- fenleimean2 %>%
      left_join(gene_id, by = "gtf_df.gene_name")

  fenleimean2 <- fenleimean2 %>%
      mutate(RowName = if_else(!is.na(gtf_df.gene_id), gtf_df.gene_id, gtf_df.gene_name)) %>%
      select(-gtf_df.gene_id, -gtf_df.gene_id) %>%
      column_to_rownames(var = "RowName")
  fenleimean2 <- fenleimean2[,-1]
  fenleimean <- rbind(fenleimean1, fenleimean2)
  
  
  
  fenleimean <- fenleimean[,colnames(fenleimean) %in% colnames(frac)]
  fenleimean <- data.frame(fenleimean[,order(colnames(fenleimean))])
  fenleimean <- fenleimean[! rowSums(fenleimean) == 0,]
  fenleimean <- fenleimean[rownames(fenleimean) %in% rownames(bulk),]
  bulk <- bulk[rownames(bulk) %in% rownames(fenleimean),]
  colnames(frac) <- colnames(fenleimean)
  fenleimean <- data.frame(fenleimean[order(rownames(fenleimean)),])
  bulk <- data.frame(bulk[order(rownames(bulk)),])
  
  colnames(bulk) = rownames(frac) = paste0('s', 1:nrow(frac))
  colnames(frac) = colnames(fenleimean) = paste0('c', 1:ncol(frac))
  
  # === 路径修改：确保保存输出时自动创建其父目录以防报错 ===
  dir.create(dirname(list4[[k]]), recursive = TRUE, showWarnings = FALSE)
  # =====================================================
  
  deconv = bMIND(bulk, frac = frac, profile = fenleimean, ncore = 12)
  data_list <- lapply(deconv, function(df) {
    colnames(df) <- CT1
    return(df)
  }) 
  save(data_list, file = list4[[k]])
}
