library(data.table)

# ====================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("错误: 必须提供项目根目录路径 (PROJECT_DIR)！", call. = FALSE)
}
PROJECT_DIR <- args[1]

path <- paste0(PROJECT_DIR, "/Downstream_analysis/complextraits/GWAS/") 
trait <- list.dirs(path, full.names = TRUE, recursive = FALSE)
trait <- sapply(trait, function(x) unlist(strsplit(x, "\\/"))[9])
trait <- data.frame(trait)
trait <- trait$trait

for (i in 1:length(trait)) {
  gwas <- fread(paste0(path, trait[[i]], "/", trait[[i]], "_china.mlma"))
  gwas <- data.frame(gwas)
  gwas$rs_id <- paste0(gwas[,2], "_", gwas[,3])
  select_gwas <- gwas[,c(13,6,5,7,8,9,10,12)]
  names(select_gwas) <- c("SNP","A1","A2","freq","b","se","p","N")
  write.table(select_gwas, paste0(path, trait[[i]], "/", trait[[i]], "_cojo.txt"), sep = "\t", row.names = F, quote = FALSE)
}
