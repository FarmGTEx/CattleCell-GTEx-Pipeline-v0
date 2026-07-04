suppressMessages(library(dplyr))
suppressMessages(library(SeuratObject))
suppressMessages(library(Seurat))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(slingshot))
suppressMessages(library(uwot))
suppressMessages(library(mclust))
suppressMessages(library(MeDuSA))
library(RColorBrewer)

# === 参数解析封装修改：引入外部参数，但保留脚本内部原有的批量循环逻辑 ===
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("错误: 必须提供项目根目录路径 (PROJECT_DIR)！", call. = FALSE)
}
PROJECT_DIR <- args[1] # 对应 Bash 中的 $BASE_PROJECT 路径
# ====================================================================

## cell trajectory
#load sc data
# === 路径修改：将写死的路径前缀替换为由参数传入的 PROJECT_DIR ===
path <- file.path(PROJECT_DIR, "CattleGTEx/OmiGA/eQTL/")
# ============================================================
tissue <- list.dirs(path, full.names = TRUE, recursive = FALSE)
tissue <- sapply(tissue, function(x) unlist(strsplit(x, "\\/"))[9])
tissue <- data.frame(tissue)
tissue <- tissue$tissue

for (k in 1:length(tissue)) { 
  # === 路径修改：将写死的路径前缀替换为由参数传入的 PROJECT_DIR ===
  sc <- readRDS(paste0(PROJECT_DIR, "/Global_atlas/Global_cell_atlas/celltype_annotation/", tissue[[k]], "_anno.rds"))
  # ============================================================
  Idents(sc) <- "CellType"
  ct <- unique(sc$CellType)
  ct_number <- data.frame(table(Idents(sc)))
  for (i_outer in 1:length(ct)) {
    if (ct_number$Freq[ct_number$Var1 == ct[[i_outer]]] < 100) {
      sc <- subset(sc, idents = ct[[i_outer]], invert = TRUE)
    }
  }
  ct <- unique(sc$CellType)
  if ("Unknown cells" %in= ct == TRUE) {
    sc <- subset(sc, idents = "Unknown cells", invert = TRUE)   
  }
  ct <- unique(sc$CellType)
  
  sc1 <- subset(sc, idents = ct[[1]])
  counts <- as.matrix(sc1@assays$RNA@counts)
  meta <- data.frame(sc1@meta.data)
  sce <- SingleCellExperiment(assays = List(counts = counts))
  #Filter gene with low expression
  geneFilter <- apply(assays(sce)$counts,1,function(x){
    sum(x >= 3) >= 10
  })
  sce <- sce[geneFilter, ]
  #normalization and scale
  FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
  }
  assays(sce)$norm <- FQnorm(assays(sce)$counts)
  #Dimensionality Reduction
  pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
  rd1 <- pca$x[,1:2]
  rd2 <- uwot::umap(t(log1p(assays(sce)$norm)))
  colnames(rd2) <- c('UMAP1', 'UMAP2')
  reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)
  #cluster
  cl1 <- Mclust(rd1)$classification
  colData(sce)$GMM <- cl1
  cl2 <- kmeans(rd1, centers = 4)$cluster
  colData(sce)$kmeans <- cl2
  #trajectory
  sce <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA')
  colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
  plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
  #tiff(filename = '/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/Cell_state/Blood/Classicalmonocytes/trajectory_Classicalmonocytes.tiff',width =1200,height=1000,res=300)
  #plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
  #lines(SlingshotDataSet(sce), lwd=2, col='black')
  #dev.off() 
  #select one-dimension tracjectory
  pseudotime <- data.frame(colData(sce)[1])
  for (i_inner in 4:length(names(colData(sce)))) {
    pseudotime1 <- data.frame(colData(sce)[i_inner])
    pseudotime <- cbind(pseudotime, pseudotime1)
  }
  pseudotime <- data.frame(pseudotime[,-1])
  rownames(pseudotime) <- rownames(pseudotime1)
  na_count <- sapply(pseudotime, function(col) sum(is.na(col)))
  min_na_col <- which.min(na_count)
  pseudotime <- pseudotime[, min_na_col, drop = FALSE]
  pseudotime <- na.omit(pseudotime)
  colnames(pseudotime) <- "pseudotime"
  select_counts <- counts[,colnames(counts) %in% rownames(pseudotime)]
  select_meta <- meta[rownames(meta) %in% rownames(pseudotime),]
  sce <- CreateSeuratObject(counts = select_counts, project = "test")
  sce$cell_trajectory <- pseudotime$pseudotime
  sce$cell_type <- select_meta$CellType
  sce$cell_trajectory <- (sce$cell_trajectory - min(sce$cell_trajectory)) / (max(sce$cell_trajectory) - min(sce$cell_trajectory))
  sce@assays$RNA@counts = sweep(as.matrix(sce@assays$RNA@counts),2,colSums(sce@assays$RNA@counts),'/')*1e+3
  for (m in 2:length(ct)) {
    sc1 <- subset(sc, idents = ct[[m]])
    counts <- as.matrix(sc1@assays$RNA@counts)
    meta <- data.frame(sc1@meta.data)
    sce1 <- SingleCellExperiment(assays = List(counts = counts))
    #Filter gene with low expression
    geneFilter <- apply(assays(sce1)$counts,1,function(x){
      sum(x >= 3) >= 10
    })
    sce1 <- sce1[geneFilter, ]
    #normalization and scale
    FQnorm <- function(counts){
      rk <- apply(counts,2,rank,ties.method='min')
      counts.sort <- apply(counts,2,sort)
      refdist <- apply(counts.sort,1,median)
      norm <- apply(rk,2,function(r){ refdist[r] })
      rownames(norm) <- rownames(counts)
      return(norm)
    }
    assays(sce1)$norm <- FQnorm(assays(sce1)$counts)
    #Dimensionality Reduction
    pca <- prcomp(t(log1p(assays(sce1)$norm)), scale. = FALSE)
    rd1 <- pca$x[,1:2]
    rd2 <- uwot::umap(t(log1p(assays(sce1)$norm)))
    colnames(rd2) <- c('UMAP1', 'UMAP2')
    reducedDims(sce1) <- SimpleList(PCA = rd1, UMAP = rd2)
    #cluster
    cl1 <- Mclust(rd1)$classification
    colData(sce1)$GMM <- cl1
    cl2 <- kmeans(rd1, centers = 4)$cluster
    colData(sce1)$kmeans <- cl2
    #trajectory
    sce1 <- slingshot(sce1, clusterLabels = 'GMM', reducedDim = 'PCA')
    colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
    plotcol <- colors[cut(sce1$slingPseudotime_1, breaks=100)]
    #tiff(filename = '/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/Cell_state/Blood/Classicalmonocytes/trajectory_Classicalmonocytes.tiff',width =1200,height=1000,res=300)
    #plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
    #lines(SlingshotDataSet(sce), lwd=2, col='black')
    #dev.off() 
    #select one-dimension tracjectory
    pseudotime <- data.frame(colData(sce1)[1])
    for (i in 4:length(names(colData(sce1)))) {
      pseudotime1 <- data.frame(colData(sce1)[i])
      pseudotime <- cbind(pseudotime, pseudotime1)
    }
    pseudotime <- data.frame(pseudotime[,-1])
    rownames(pseudotime) <- rownames(pseudotime1)
    na_count <- sapply(pseudotime, function(col) sum(is.na(col)))
    min_na_col <- which.min(na_count)
    pseudotime <- pseudotime[, min_na_col, drop = FALSE]
    pseudotime <- na.omit(pseudotime)
    colnames(pseudotime) <- "pseudotime"
    select_counts <- counts[,colnames(counts) %in% rownames(pseudotime)]
    select_meta <- meta[rownames(meta) %in% rownames(pseudotime),]
    sce1 <- CreateSeuratObject(counts = select_counts, project = "test")
    sce1$cell_trajectory <- pseudotime$pseudotime
    sce1$cell_type <- select_meta$CellType
    sce1$cell_trajectory <- (sce1$cell_trajectory - min(sce1$cell_trajectory)) / (max(sce1$cell_trajectory) - min(sce1$cell_trajectory))
    sce1@assays$RNA@counts = sweep(as.matrix(sce1@assays$RNA@counts),2,colSums(sce1@assays$RNA@counts),'/')*1e+3
    sce <- merge(sce,sce1)
  }
  
  #Begin Medusa
  Idents(sce) <- "cell_type"
  ct_number <- data.frame(table(Idents(sce)))
  ct <- unique(sce$cell_type)
  ct1<- gsub(" ","",ct)
  ct1<- gsub("\\+","",ct1)
  ct1<- gsub("γδ","",ct1)
  ct1<- gsub("+","",ct1)
  ct1<- gsub("-","",ct1)
  ct1<- gsub("'","",ct1)
  ct1<- gsub("\\(","",ct1)
  ct1<- gsub("\\)","",ct1)
  ##medusa
  #sce$cell_trajectory <- (sce$cell_trajectory - min(sce$cell_trajectory)) / (max(sce$cell_trajectory) - min(sce$cell_trajectory))
  #sce@assays$RNA@counts = sweep(as.matrix(sce@assays$RNA@counts),2,colSums(sce@assays$RNA@counts),'/')*1e+3
  # === 路径修改：将写死的路径前缀替换为由参数传入的 PROJECT_DIR ===
  bulk <- read.table(paste0(PROJECT_DIR, "/CattleGTEx/Cattle_bulk_exp/Bulk_", tissue[[k]], ".txt"), sep = "\t", row.names = 1, header = T)
  # ============================================================
  filtered_bulk <- na.omit(bulk)
  bulk1 = sweep(as.matrix(filtered_bulk),2,colSums(filtered_bulk),'/')*1e+3
  # === 路径修改：将写死的路径前缀替换为由参数传入的 PROJECT_DIR ===
  setwd(paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/Cell_state/"))
  # ============================================================
  if (!dir.exists(tissue[[k]])) {
    dir.create(tissue[[k]])
  }
  for (j in 1:length(ct)) {
    # === 路径修改：将写死的路径前缀替换为由参数传入的 PROJECT_DIR ===
    setwd(paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/Cell_state/", tissue[[k]]))
    # ============================================================
    if (!dir.exists(ct1[[j]])) {
      dir.create(ct1[[j]])
    }
    # === 路径修改：将写死的路径前缀替换为由参数传入的 PROJECT_DIR ===
    setwd(paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/Cell_state/", tissue[[k]], "/", ct1[[j]]))
    # ============================================================
    MeDuSA_obj = MeDuSA(bulk1,sce,
                      select.ct = ct[[j]], markerGene = NULL,span = 0.35,method = "wilcox",
                      resolution = 10,smooth = TRUE,fractional = TRUE,ncpu = 12)
    result <- data.frame(MeDuSA_obj@Estimation$cell_state_abundance)
    write.csv(result,paste0(tissue[[k]], "_", ct1[[j]], ".csv"))
  }
}

# === 路径修改：将写死的路径前缀替换为由参数传入的 PROJECT_DIR ===
sc <- readRDS(paste0(PROJECT_DIR, "/Global_atlas/Global_cell_atlas/celltype_annotation/", tissue[[k]], "_anno.rds"))
# ============================================================
Idents(sc) <- "CellType"
ct <- unique(sc$CellType)
ct_number <- data.frame(table(Idents(sc)))
for (i_outer in 1:length(ct)) {
  if (ct_number$Freq[ct_number$Var1 == ct[[i_outer]]] < 100) {
    sc <- subset(sc, idents = ct[[i_outer]], invert = TRUE)
  }
}
ct <- unique(sc$CellType)
if ("Unknown cells" %in= ct == TRUE) {
  sc <- subset(sc, idents = "Unknown cells", invert = TRUE)   
}
ct <- unique(sc$CellType)
for (j in 1:length(ct)) {
  sc1 <- subset(sc, idents = ct[[j]])
  counts <- as.matrix(sc1@assays$RNA@counts)
  meta <- data.frame(sc1@meta.data)
  sce <- SingleCellExperiment(assays = List(counts = counts))
  #Filter gene with low expression
  geneFilter <- apply(assays(sce)$counts,1,function(x){
    sum(x >= 3) >= 10
  })
  sce <- sce[geneFilter, ]
  #normalization and scale
  FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
  }
  assays(sce)$norm <- FQnorm(assays(sce)$counts)
  #Dimensionality Reduction
  pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
  rd1 <- pca$x[,1:2]
  rd2 <- uwot::umap(t(log1p(assays(sce)$norm)))
  colnames(rd2) <- c('UMAP1', 'UMAP2')
  reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)
  #cluster
  cl1 <- Mclust(rd1)$classification
  colData(sce)$GMM <- cl1
  cl2 <- kmeans(rd1, centers = 4)$cluster
  colData(sce)$kmeans <- cl2
  #trajectory
  sce <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA')
  colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
  plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
  # === 路径修改：将写死的路径前缀替换为由参数传入的 PROJECT_DIR ===
  # 自动创建父目录以防报错
  mam_gland_dir <- paste0(PROJECT_DIR, '/CattleGTEx/OmiGA/Cell_state/MammaryGland/')
  dir.create(mam_gland_dir, recursive = TRUE, showWarnings = FALSE)
  tiff(filename = paste0(mam_gland_dir, ct[[j]], '.tiff'),width =1200,height=1000,res=300)
  # ============================================================
  plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
  lines(SlingshotDataSet(sce), lwd=2, col='black')
  dev.off()
}
