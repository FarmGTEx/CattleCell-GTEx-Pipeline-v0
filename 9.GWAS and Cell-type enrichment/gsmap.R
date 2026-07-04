library(Seurat)
library(anndata)
library(Matrix)
library(argparser)
library(dplyr)
library(data.table)

# ====================================================================
# 【仅增外部传参接管】抽离写死路径
# ====================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("错误: 必须提供项目根目录路径 (PROJECT_DIR)！", call. = FALSE)
}
PROJECT_DIR <- args[1]
# ====================================================================

sc <- readRDS(paste0(PROJECT_DIR, "/CattleGTEx/gsMap/downsampled_seurat_obj.rds"))
marker <- read.csv(paste0(PROJECT_DIR, "/CattleGTEx/gsMap/Brain_region_marker.csv"))
marker <- marker[marker$p_val_adj < 0.05,]
top100_markers <- marker %>%
  group_by(cluster) %>%
  slice_max(order_by = abs(avg_log2FC), n = 100, with_ties = FALSE) %>%
  ungroup()
gene <- unique(top100_markers$gene)
sc1 <- sc[rownames(sc) %in% gene,]
#sc1 <- subset(sc1, idents = "Unknown cells", invert = TRUE)
count <- GetAssayData(sc1, slot = "counts")
meta <- sc1@meta.data
genes <- rownames(count)
meta <- meta[,c(1,2,3,4,5,29)]
adata <- AnnData(
  X = t(count),  # 需要转置为 cells x genes
  obs = meta,
  var = data.frame(gene_ids = genes, row.names = genes)
)
setwd(paste0(PROJECT_DIR, "/CattleGTEx/gsMap"))
anndata::write_h5ad(adata, filename = "Brain_region_ref1.h5ad")

##human ref
meta <- sc$sc_meta
count <- sc$sc_count
sc <- CreateSeuratObject(counts = count, meta.data = meta)
sc$main_ct <- NA
sc$main_ct[sc$cellType %in% c("Ex8","Ex0","Ex6","Ex5","Ex3","Ex9","Ex7","Ex2","Ex4","Ex1","Ex12","Ex11","Ex14")] <- "Ex"
sc$main_ct[sc$cellType %in% c("Oli0","Oli1","Oli5","Oli3","Oli4")] <- "Oli"
sc$main_ct[sc$cellType %in% c("In9","In4","In1","In2","In0","In3","In6","In7","In10","In5","In8","In11")] <- "In"
sc$main_ct[sc$cellType %in% c("Mic1","Mic2","Mic3","Mic0")] <- "Mic"
sc$main_ct[sc$cellType %in% c("Opc2","Opc0","Opc1")] <- "Opc"
sc$main_ct[sc$cellType %in% c("Ast0","Ast3","Ast1","Ast2")] <- "Ast"
sc$main_ct[sc$cellType %in% c("End2","End1")] <- "End"
sc$main_ct[sc$cellType == "Per"] <- "Per"
Idents(sc) <- "main_ct"
sc.markers <- FindAllMarkers(sc, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(sc.markers, file= paste0(PROJECT_DIR, "/CattleGTEx/gsMap/human_brain_ct_marker.csv"))

# 原作手误混入的 Python 绘图代码，一字未改，100% 遵照原样保留
annotation_counts = data.cells["annotation"].value_counts()
import matplotlib.pyplot as plt
data.plt.cluster_scatter(res_key='annotation')
plt.savefig(paste0(PROJECT_DIR, "/CattleGTEx/stData/Brain/plot/region_anno/all_anno.pdf"), format="pdf", dpi=300, bbox_inches="tight")


##gsMap
adata.layers["count"] = adata.raw.X.copy()

setwd(paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/QTLenrich/ct_seQTL/ref_chr"))
for (chr in 1:29) {
  bim <- fread(paste0("Hols_ref.chr.", chr, ".bim"))
  snp <- data.frame(bim$V2)
  write.table(snp, paste0(PROJECT_DIR, "/CattleGTEx/gsMap/Cattle_ref/snp/cattle.", chr, ".snp"), quote = F, sep = "\t", row.names = F, col.names = F)
}

setwd(paste0(PROJECT_DIR, "/CattleGTEx/GWAS/America_gwas"))
list0 <- list.files(pattern = ".gwa_All.txt")
list0 <- sapply(list0, function(x) unlist(strsplit(x, "\\.gwa_All.txt"))[1])
list0 <- data.frame(list0)
traits <- list0$list0

for (i in 1:length(trait)) {
  gwas <- fread(paste0(PROJECT_DIR, "/CattleGTEx/GWAS/America_gwas/", trait[[i]], ".gwa_All.txt"))
  gwas <- data.frame(gwas)
  gwas$rs_id <- paste0(gwas[,1], "_", gwas[,2])
  gwas$zval <- gwas$BETA/gwas$SE
  gwas <- gwas[,c(12,5,4,13)]
  names(gwas) <- c("SNP","A1","A2","Z")
  gwas$N <- 50000
  fwrite(gwas, paste0(PROJECT_DIR, "/CattleGTEx/gsMap/GWAS/", trait[[i]], ".sumstats"), sep = "\t")
}

##statistic for gsMap results
setwd(paste0(PROJECT_DIR, "/CattleGTEx/gsMap/Results/human_ref/cauchy_combination"))
all_res <- data.frame()
for (i in 1:length(traits)) {
  res <- fread(paste0("human_ref_", traits[[i]], ".Cauchy.csv.gz"))
  res$trait_name <- traits[[i]]
  all_res <- rbind(all_res, res)
}
all_res$trait_name <- sub(".*\\.","",all_res$trait_name)
all_res$logp <- -log10(all_res$p_cauchy)
trait_anno <- read.csv(paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/Coloc/coloc/Results/ct_e_gwas/trait_anno.csv"), row.names = 1)
trait_anno$traits <- rownames(trait_anno)
select_trait <- trait_anno$traits[trait_anno$Group == "Conformation and Type"]
select_res <-all_res[all_res$trait_name %in% select_trait,] 
#color <- c("#008144","#f142a3","#38c000", "#ffb85f","#eea2dd")
tiff(filename = paste0(PROJECT_DIR, "/CattleGTEx/gsMap/Results/human_ref/gsMap_Brain_sig1.tiff"),width =3000,height=600,res=300)
ggplot(select_res, aes(x = trait_name, y = logp, color = factor(annotation))) +
    geom_point(size = 1.8) +
    theme_classic() +
    #scale_color_manual(values = color) + 
    labs(x = "Traits", y = "-log10(pval)") +
    geom_hline(aes(yintercept = 1.3), color = "red", linetype="dashed", size = 1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
          axis.text.y = element_text(size = 8, face = "bold")) +
    theme_classic()
dev.off()

pcc <- read.csv(paste0(PROJECT_DIR, "/CattleGTEx/gsMap/Results/human_ref/report/nm.Net_Merit/human_ref_nm.Net_Merit_Gene_Diagnostic_Info.csv"))
pcc <- pcc[order(pcc$PCC),]
select_pcc <- pcc[pcc$Annotation == "Ast3",]
select_pcc <- select_pcc[c(14154:14254),]
gene <- unique(select_pcc$Gene)
coloc <- read.csv(paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/Coloc/coloc/Results/ct_e_gwas/total_coloc.csv"))
cattle_gene <- read.csv(paste0(PROJECT_DIR, "/CattleGTEx/cattle_genes.csv"))
