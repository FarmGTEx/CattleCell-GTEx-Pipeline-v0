library(scPagwas)
library(rtracklayer)
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("错误: 必须提供项目根目录路径 (PROJECT_DIR)！", call. = FALSE)
}
PROJECT_DIR <- args[1]
# ====================================================================

#prepare sc data
str <- paste0(PROJECT_DIR, "/Global_atlas/Global_cell_atlas/celltype_annotation/")
setwd(paste0(PROJECT_DIR, "/Global_atlas/Global_cell_atlas/celltype_annotation"))
list0 <- list.files(pattern = ".rds")
list0 <- sapply(list0, function(x) unlist(strsplit(x, "\\."))[1])
list0 <- data.frame(list0)
tissue <- list0$list0
list1<-NULL
for(i in tissue){
  list1[[i]]<-paste0(str, i, ".rds")
}

#prepare genome annotation data
gtf_df<- rtracklayer::import(paste0(PROJECT_DIR, "/reference/Bos_taurus.ARS-UCD1.2.110.gtf"))
gtf_df <- as.data.frame(gtf_df)
gtf_df <- gtf_df[,c("seqnames","start","end","type","gene_name")]
gtf_df <- gtf_df[gtf_df$type=="gene",]
block_annotation<-gtf_df[,c(1,2,3,5)]
colnames(block_annotation)<-c("chrom", "start","end","label")
block_annotation$chrom <- paste0("chr", block_annotation$chrom)
#prepare ld file
#ld <- read.table("produc.ld.ld", sep = "", header = T)
#lapply(unique(ld$CHR_A), function(i){
#a<-data.table(ld[ld$CHR_A == i,])
#file_namePath <- paste0("/faststorage/project/cattle_gtexs/Downstream_analysis/complextraits/",i,".Rds")
#saveRDS(a, file = file_name)
#})
ld<-lapply(as.character(1:29),function(chrom){
chrom_ld_file_path <- paste(PROJECT_DIR, '/Downstream_analysis/complextraits/ld_chrom', '/', chrom, '.Rds', sep = '')
ld_data <- readRDS(chrom_ld_file_path)[, .(SNP_A, SNP_B, R)]
return(ld_data)
})
ld <- setNames(ld, paste0("chr", c(1:29)))

#prepare pathway data
#merge_dat <- data.frame()
#for (i in 1:length(Genes_by_pathway_kegg)) {
#list1 <- data.frame(Genes_by_pathway_kegg[i])
#list1$group <- names(Genes_by_pathway_kegg[i])
#names(list1)[1] <- "gene"
#merge_dat <- rbind(merge_dat, list1)
#}
#transfer the human genes to cattle genes by using Biomart
pathway <- read.csv(paste0(PROJECT_DIR, "/Downstream_analysis/complextraits/scpagwas/cattle_human_pathway.csv"))
pathway <- pathway[,-1]
kegg <- split(pathway$gene, pathway$group)


for (k in 1:length(tissue)) { 
sc <- readRDS(list1[[k]])
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sc)
sc <- ScaleData(sc, features = all.genes)
#run the cycle
#str1 <- "/faststorage/project/cattle_gtexs/Downstream_analysis/complextraits/MAGMA/"
#group <- list.dirs(str1, full.names = FALSE, recursive = FALSE)
#group <- group[c(5,6)]
#for(i in 1:length(group)) {
#str2 <- paste0(str1, group[[i]], "/")
#list2 := list.dirs(str2, full.names = FALSE, recursive = FALSE)
#list3<-NULL
#for(j in list2){
#list3[[j]]<-paste0(str2, j, "/", j, "_china.mlma")
#}
#for(j in 1:length(list2)) {
#prepare gwasdata
path <- paste0(PROJECT_DIR, "/Downstream_analysis/complextraits/05sperm/")
traits <- list.dirs(path, full.names = TRUE, recursive = FALSE)
traits <- sapply(traits, function(x) unlist(strsplit(x, "\\/"))[9])
traits <- data.frame(traits)
traits <- traits$traits
path<-NULL
for(i in traits){
  path[[i]]<-paste0(PROJECT_DIR, "/Downstream_analysis/complextraits/05sperm/", i, "/")
}


for (i in 1:length(traits)) {
gwas <- read.table(paste0(PROJECT_DIR, "/Downstream_analysis/complextraits/05sperm/", traits[[i]], "/", traits[[i]], "_china.mlma") ,sep = "\t")
gwas <- gwas[,c(1,2,3,8,7,6)]
names(gwas) <- c("chrom","pos","rsid","se","beta","maf")
##separate the pipeline
Pagwas <- list()
Pagwas <- Single_data_input(
Pagwas = Pagwas,
assay = "RNA",
Single_data = sc,
#Pathway_list = Genes_by_pathway_kegg
Pathway_list = kegg
)
sc <- sc[, colnames(Pagwas$data_mat)]
  
Pagwas <- Pathway_pcascore_run(
Pagwas = Pagwas,
#Pathway_list = Genes_by_pathway_kegg
Pathway_list = kegg
)

Pagwas <- GWAS_summary_input(
Pagwas = Pagwas,
gwas_data = gwas,
maf_filter = 0.1
)

Pagwas$snp_gene_df <- SnpToGene(
gwas_data = Pagwas$gwas_data,
block_annotation = block_annotation,
marg = 10000
)

Pagwas <- Pathway_annotation_input(
Pagwas = Pagwas,
block_annotation = block_annotation
)
        
Pagwas <- Link_pathway_blocks_gwas(
Pagwas = Pagwas,
chrom_ld = ld,
singlecell = F,
celltype = T,
backingpath=paste0(PROJECT_DIR, "/Downstream_analysis/complextraits/scpagwas/scPagwastest_output"))

Pagwas$lm_results <- Pagwas_perform_regression(Pathway_ld_gwas_data = Pagwas$Pathway_ld_gwas_data)
Pagwas <- Boot_evaluate(Pagwas, bootstrap_iters = 200, part = 0.5)
Pagwas$Pathway_ld_gwas_data <- NULL
results <- data.frame(Pagwas$bootstrap_results)
write.csv(results, paste0(PROJECT_DIR, "/CCA_NG/Qi/Germ/new_", traits[[i]], ".csv"))
}
}
