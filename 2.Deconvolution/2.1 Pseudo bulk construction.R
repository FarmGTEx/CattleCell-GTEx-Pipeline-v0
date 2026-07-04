suppressMessages(library(Seurat))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(dplyr))
suppressMessages(library(FamilyRank))

# === 参数解析封装修改：接收来自 Bash 脚本的输入输出参数 ===
args <- commandArgs(trailingOnly = TRUE)
input_rds <- args[1]   # 对应 Bash 中的 $DEFAULT_INPUT_RDS
output_dir <- args[2]  # 对应 Bash 中的 $DEFAULT_OUTPUT_DIR
gtf_path <- args[3]    # 对应 Bash 中的 $GTF_REFERENCE
# ========================================================

# === 路径修改：使用 Bash 传入的 input_rds 替代硬编码路径 ===
sc <- readRDS(input_rds)
# ========================================================

#sc <- readRDS("/faststorage/project/cattle_gtexs/Deconvolution/Cattle/Pseudo/Cerebral cortex/Cerebral cortex_pseudo.rds")
Idents(sc) <- "CellType"
#sc <- subset(sc, idents = "Proliferative cells", invert = TRUE)
#sc <- subset(sc, idents = "Unknown cells", invert = TRUE)
counts <- data.frame(sc@assays$RNA@counts)

### TPM  
# === 路径修改：使用 Bash 传入的 gtf_path 替代硬编码路径 ===
gtf = rtracklayer::import(gtf_path)
# ========================================================
class(gtf)
gtf = as.data.frame(gtf);dim(gtf)
table(gtf$type)
exon = gtf[gtf$type=="exon",
           c("start","end","gene_name")]
gle = lapply(split(exon,exon$gene_name),function(x){
  tmp=apply(x,1,function(y){
    y[1]:y[2]
  })
  length(unique(unlist(tmp)))
})
gle=data.frame(gene_name=names(gle),
               length=as.numeric(gle))
le = gle[match(rownames(counts),gle$gene_name),"length"]
counts$Length <- le
counts <- na.omit(counts)
kb <- counts$Length / 1000
x = ncol(counts)
countdata <- counts[,1:x-1]
rpk <- countdata / kb
rpk = rpk[complete.cases(rpk),]
tpm <- t(t(rpk)/colSums(rpk) * 1000000)

tpm <- t(tpm)
tpm <- data.frame(tpm)
meta <- data.frame(sc@meta.data)
meta <- droplevels(meta)
tpm <- tpm %>% mutate(cluster = meta$CellType)
#sc_tpm$cluster <- paste0("cluster", sc_tpm$cluster)
group <- list()
group <- split(tpm,tpm$cluster)
name<- rownames(table(tpm$cluster))
fenleimean <- colMeans(group[[1]][,-ncol(group[[1]])])
for (j in 2:length(group)) {
  flmean<- colMeans(group[[j]][,-ncol(group[[j]])])
  fenleimean <-rbind(fenleimean,flmean)
}
rownames(fenleimean) <- name
fenleimean = t(fenleimean)


### CPM
library(edgeR)  
cpm <- edgeR::cpm(counts)
cpm <- t(cpm)
cpm <- data.frame(cpm)
meta <- data.frame(sc@meta.data)
meta <- droplevels(meta)
cpm <- cpm %>% mutate(cluster = meta$CellType)
#sc_tpm$cluster <- paste0("cluster", sc_tpm$cluster)
group <- list()
group <- split(cpm,cpm$cluster)
name<- rownames(table(cpm$cluster))
fenleimean <- colMeans(group[[1]][,-ncol(group[[1]])])
for (j in 2:length(group)) {
  flmean<- colMeans(group[[j]][,-ncol(group[[j]])])
  fenleimean <-rbind(fenleimean,flmean)
}
rownames(fenleimean) <- name
fenleimean = t(fenleimean)


### Count
counts <- t(counts)
counts <- data.frame(counts)
meta <- data.frame(sc@meta.data)
meta <- droplevels(meta)
counts <- counts %>% mutate(cluster = meta$CellType)
#sc_tpm$cluster <- paste0("cluster", sc_tpm$cluster)
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



## Binormal distribution
#cell <- data.frame(table(Idents(sc)))
#real <- cell$Freq / sum(cell$Freq)
#real <- sort(real)
combFunc <- function(...) {
    mapply('bind_cols', ..., SIMPLIFY=FALSE)
}
best_params <- NULL
best_similarity <- -Inf
sd1 <- c(0.1,0.5,0.9)
sd2 <- c(0.1,0.5,0.9)

for (i in 1:length(sd1)) {
  for (j in 1:length(sd2)) {
    fraction <- array(data = NA,dim = c(20,1000))
    for (k in 1:1000) {
      simulated_data <- abs(rbinorm(20, 0, 1, sd1[[i]], sd2[[j]], 0.5))
      simulated_data <- simulated_data / sum(simulated_data)
      fraction[,k] <- simulated_data
    } 
    rownames(fraction) <- paste0("cell", c(1:20))
    colnames(fraction) <- paste0("sample", c(1:1000))
    # === 路径修改：使用 Bash 传入的 output_dir 动态拼接输出路径 ===
    dir.create(file.path(output_dir, "Binorm_Fraction"), recursive = TRUE, showWarnings = FALSE)
    write.csv(fraction, file.path(output_dir, "Binorm_Fraction", paste0("Fraction_", sd1[[i]], "_", sd2[[j]], ".csv"))) 
    # ===========================================================
  }
}

#produce pseudo bulk expression
# === 路径修改：使用 Bash 传入的 output_dir 动态切换工作目录 ===
setwd(file.path(output_dir, "Binorm_Fraction"))
# ===========================================================
file <- list.files(pattern = ".csv")
list0 <- tools::file_path_sans_ext(file)
list1<-NULL
for(i in list0){
  # === 路径修改：使用 Bash 传入的 output_dir 动态拼接输出路径 ===
  list1[[i]]<-file.path(output_dir, "Binorm_Fraction", "random_pseudo", paste0(i, "_pseudo.txt"))
  # ===========================================================
}

# === 路径修改：自动创建 random_pseudo 目录以防写入失败 ===
dir.create(file.path(output_dir, "Binorm_Fraction", "random_pseudo"), recursive = TRUE, showWarnings = FALSE)
# =====================================================

for (i in 1:length(file)) {
  pseudo <- array(data = NA,dim = c(nrow(fenleimean),1000))
  fraction <- read.csv(file[[i]], row.names = 1)
  for (j in 1:1000) {
    for (k in 1:nrow(fenleimean)) {
      gene <- sum(fenleimean[k,] * fraction[,j])
      pseudo[k,j] <- gene
    }
  }
  rownames(pseudo) <- rownames(fenleimean)
  colnames(pseudo) <- paste0("sample", c(1:1000))
  write.table(pseudo, list1[[i]], sep = "\t")
}
      

## normal distribution
combFunc <- function(...) {
    mapply('bind_cols', ..., SIMPLIFY=FALSE)
}
best_params <- NULL
best_similarity <- -Inf
sd <- c(0.1,0.5,0.9)

for (i in 1:length(sd)) {
  fraction <- array(data = NA,dim = c(20,1000))
  for (k in 1:1000) {
    simulated_data <- abs(rnorm(20, 0, sd[[i]]))
    simulated_data <- simulated_data / sum(simulated_data)
    fraction[,k] <- simulated_data
  }
  rownames(fraction) <- paste0("cell", c(1:20))
  colnames(fraction) <- paste0("sample", c(1:1000))
  # === 路径修改：使用 Bash 传入的 output_dir 动态拼接输出路径 ===
  dir.create(file.path(output_dir, "Normal_Fraction"), recursive = TRUE, showWarnings = FALSE)
  write.csv(fraction, file.path(output_dir, "Normal_Fraction", paste0("Fraction_sd_", sd[[i]], ".csv"))) 
  # ===========================================================
}

#produce pseudo bulk expression
# === 路径修改：使用 Bash 传入的 output_dir 动态切换工作目录 ===
setwd(file.path(output_dir, "Normal_Fraction"))
# ===========================================================
file <- list.files(pattern = ".csv")
list0 <- tools::file_path_sans_ext(file)
list1<-NULL
for(i in list0){
  # === 路径修改：使用 Bash 传入的 output_dir 动态拼接输出路径 ===
  list1[[i]]<-file.path(output_dir, "Normal_Fraction", "random_pseudo", paste0(i, "_pseudo.txt"))
  # ===========================================================
}

# === 路径修改：自动创建 random_pseudo 目录以防写入失败 ===
dir.create(file.path(output_dir, "Normal_Fraction", "random_pseudo"), recursive = TRUE, showWarnings = FALSE)
# =====================================================

for (i in 1:length(file)) {
  pseudo <- array(data = NA,dim = c(nrow(fenleimean),1000))
  fraction <- read.csv(file[[i]], row.names = 1)
  for (j in 1:1000) {
    for (k in 1:nrow(fenleimean)) {
      gene <- sum(fenleimean[k,] * fraction[,j])
      pseudo[k,j] <- gene
    }
  }
  rownames(pseudo) <- rownames(fenleimean)
  colnames(pseudo) <- paste0("sample", c(1:1000))
  write.table(pseudo, list1[[i]], sep = "\t")
}



## uniform distribution
combFunc <- function(...) {
    mapply('bind_cols', ..., SIMPLIFY=FALSE)
}
best_params <- NULL
best_similarity <- -Inf
fraction <- array(data = NA,dim = c(20,1000))
for (k in 1:1000) {
  simulated_data <- runif(20, 1, 99)
  simulated_data <- simulated_data / sum(simulated_data)
  fraction[,k] <- simulated_data
}
rownames(fraction) <- paste0("cell", c(1:20))
colnames(fraction) <- paste0("sample", c(1:1000))
# === 路径修改：使用 Bash 传入的 output_dir 动态拼接输出路径 ===
dir.create(file.path(output_dir, "Uniform_Fraction"), recursive = TRUE, showWarnings = FALSE)
write.csv(fraction, file.path(output_dir, "Uniform_Fraction", "Farction.csv"))
# ===========================================================

#produce pseudo bulk expression
# === 路径修改：使用 Bash 传入的 output_dir 动态切换工作目录与拼接输出路径 ===
setwd(file.path(output_dir, "Uniform_Fraction"))
pseudo <- array(data = NA,dim = c(nrow(fenleimean),1000))
fraction <- read.csv("Farction.csv", row.names = 1) # 统一修正对应上文写入的 Farction.csv
for (j in 1:1000) {
  for (k in 1:nrow(fenleimean)) {
    gene <- sum(fenleimean[k,] * fraction[,j])
    pseudo[k,j] <- gene
  }
}
rownames(pseudo) <- rownames(fenleimean)
colnames(pseudo) <- paste0("sample", c(1:1000))
dir.create(file.path(output_dir, "Uniform_Fraction", "random_pseudo"), recursive = TRUE, showWarnings = FALSE)
write.table(pseudo, file.path(output_dir, "Uniform_Fraction", "random_pseudo", "pseudo.txt"), sep = "\t")
# =======================================================================
