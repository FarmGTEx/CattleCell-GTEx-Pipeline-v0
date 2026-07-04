options(stringsAsFactors = FALSE)
suppressMessages(library(edgeR))
suppressMessages(library(preprocessCore))
suppressMessages(library(RNOmni))
suppressMessages(library(data.table))
suppressMessages(library(R.utils))
suppressMessages(library(SNPRelate))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))   # 引入命令行解析包

# ====================================================================
# 1. 命令行参数解析 (CLI 封装)
# ====================================================================
option_list = list(
  make_option(c("-p", "--project_dir"), type="character", default=NULL, help="项目大基座根目录路径")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$project_dir)){
  print_help(opt_parser)
  stop("错误: 必须通过 --project_dir 指定项目路径！", call.=FALSE)
}

# ====================================================================
# 2. 扫描组织并构建路径对齐
# ====================================================================
#extract cell components and corresponding expression matrix
path <- file.path(opt$project_dir, "Deconvolution/Cattle/Bulk/Results_DWLS/")
tissue <- list.dirs(path, full.names = TRUE, recursive = FALSE)
tissue <- sapply(tissue, function(x) unlist(strsplit(x, "\\/"))[10])
tissue <- data.frame(tissue)
tissue_sc <- tissue$tissue

#bulk
path <- file.path(opt$project_dir, "CattleGTEx/OmiGA/eQTL/")
tissue <- list.dirs(path, full.names = TRUE, recursive = FALSE)
tissue <- sapply(tissue, function(x) unlist(strsplit(x, "\\/"))[9])
tissue <- data.frame(tissue)
tissue_bulk <- tissue$tissue
tissue <- intersect(tissue_sc, tissue_bulk)

list1<-NULL
for(i in tissue){
  list1[[i]]<-paste0(opt$project_dir, "/CattleGTEx/OmiGA/eQTL/", i, "/expr_tmm_inv.bed.gz")
}
list2<-NULL
for(i in tissue){
  list2[[i]]<-paste0(opt$project_dir, "/CattleGTEx/OmiGA/eQTL/", i, "/ieQTL/")
}
list3<-NULL
for(i in tissue){
  # 动态流转对接统一输出格式
  list3[[i]]<-paste0(opt$project_dir, "/Deconvolution/Cattle/Bulk/Results_DWLS/", i, "/Results/Predict_", i, "_DWLS.csv")
}

# ====================================================================
# 3. 过滤低质量细胞类型与样本大循环 (原始计算逻辑保持不变)
# ====================================================================
#remove low quality cell type and samples
for (i in 1:length(tissue)) {
  frac <- read.csv(list3[[i]], row.names = 1)
  frac <- na.omit(frac)
  bed <- fread(list1[[i]])
  sample <- colnames(bed)[-c(1:4)]
  frac <- frac[rownames(frac) %in% sample, ]
  zero_ratio <- colSums(frac == 0) / nrow(frac)
  frac1 <- frac[, zero_ratio <= 0.8]
  mean_pro <- colMeans(frac1)
  frac2 <- frac1[, mean_pro >= 0.05]
  ct <- colnames(frac2)
  for (j in 1:length(ct)) {
    # 自动创建 ieQTL 目标输出目录
    dir.create(list2[[i]], recursive = TRUE, showWarnings = FALSE)
    setwd(list2[[i]])
    
    bed <- fread(list1[[i]])
    if (!dir.exists(ct[[j]])) {
        dir.create(ct[[j]])
    }
    setwd(paste0(list2[[i]], ct[[j]]))
    select_ct <- data.frame(frac2[,j])
    rownames(select_ct) <- rownames(frac2)
    names(select_ct)[1] <- "fraction"
    sd <- sd(select_ct$fraction)
    mean <- mean(select_ct$fraction)
    select_ct$select_sample = select_ct$fraction < mean + 3*sd &
          select_ct$fraction > mean - 3*sd
    select_ct <- select_ct[select_ct$select_sample == TRUE,]
    select_frac <- data.frame(select_ct[,-2])
    rownames(select_frac) <- rownames(select_ct)
    names(select_frac)[1] <- "fraction"
    bed1 <- data.frame(bed[,-c(1:4)])
    bed1 <- bed1[,colnames(bed1) %in% rownames(select_frac)]
    bed <- data.frame(bed[,1:4],bed1)
    names(bed)[1] <- "#Chr"
    fwrite(bed, paste0(ct[[j]], "_expr_tmm_inv.bed"))
    write.table(select_frac, paste0(ct[[j]],"_frac.txt"), sep = "\t", col.names = F, row.names = T, quote = F)
    message(paste("[成功] 组织:", tissue[[i]], " -> 细胞类型:", ct[[j]], " ieQTL input 构建完毕。"))
  }
}
