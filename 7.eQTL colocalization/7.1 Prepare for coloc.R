options(stringsAsFactors = FALSE)
library(edgeR)
library(preprocessCore)
library(RNOmni)
library(data.table)
library(R.utils)
library(SNPRelate)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Error: the project root directory path must be provided (PROJECT_DIR)！", call. = FALSE)
}
PROJECT_DIR <- args[1] # Corresponds to the $BASE_PROJECT path in the Bash script
# ====================================================================

###Data prepare for coloc (bulk eQTL)
# Path update: dynamically construct the path using PROJECT_DIR
setwd(paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/eQTL"))
path <- paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/eQTL/")
tissue <- list.dirs(path, full.names = TRUE, recursive = FALSE)
tissue <- sapply(tissue, function(x) unlist(strsplit(x, "\\/"))[9])
tissue <- data.frame(tissue)
tissue <- tissue$tissue
path<-NULL
for(i in tissue){
  path[[i]]<-paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/eQTL/", i, "/")
}

for (i in 1:length(tissue)) {
    setwd(path[[i]])
    eqtl <- fread(paste0(tissue[[i]], "_new_LMM.cis_qtl.txt.gz")) #produced in 4.QTL mapping
    eqtl$is_eGene = eqtl$pval_g1 < eqtl$pval_g1_threshold &
              eqtl$qval_g1 < 0.05
    eGenes <- eqtl[eqtl$is_eGene == TRUE,]
    snp_freq <- fread(paste0(tissue[[i]], "_snp_info.frq"))
    snp_freq <- snp_freq[,c(2,5)]
    names(snp_freq) <- c("variant_id","maf")
    exp <- fread("expr_tmm_inv.bed.gz")
    sample_num <- ncol(exp) - 4
    all_nom_qtl <- data.frame()
    for (chr in 1:29) {
      chr_qtl <- fread(paste0(tissue[[i]], "_new_LMM.cis_qtl_pairs.", chr, ".txt.gz"))
      #chr_qtl <- chr_qtl[chr_qtl$pval_g1 < 0.05, ]
      all_nom_qtl <- rbind(all_nom_qtl, chr_qtl)
    }
    select_qtl <- all_nom_qtl[all_nom_qtl$pheno_id %in% eGenes$pheno_id , ]
    select_qtl <- select_qtl[,c(2,2,1,3,7,5,6)]
    names(select_qtl) <- c("rs_id","variant_id","gene_id","tss_distance","pval_nominal","slope","slope_se")
    select_qtl <- select_qtl %>% inner_join(snp_freq, by = "variant_id")
    select_qtl$N <- sample_num
    
    # Path update: protect the output target directory and dynamically switch paths using PROJECT_DIR
    out_dir_1 <- paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/Coloc/coloc/Bulk")
    dir.create(out_dir_1, recursive = TRUE, showWarnings = FALSE)
    setwd(out_dir_1)
    fwrite(select_qtl, paste0(tissue[[i]], "_coloc.bed"))
}

##cell components ieqtl
# Path update: dynamically extend the path using PROJECT_DIR
setwd(paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/eQTL"))
path <- paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/eQTL/")
tissue <- list.dirs(path, full.names = TRUE, recursive = FALSE)
tissue <- sapply(tissue, function(x) unlist(strsplit(x, "\\/"))[9])
tissue <- data.frame(tissue)
tissue <- tissue$tissue
path<-NULL
for(i in tissue){
  path[[i]]<-paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/eQTL/", i, "/")
}


for (i in 1:length(tissue)) {
    setwd(path[[i]])
    eqtl <- fread(paste0(tissue[[i]], "_bulk_LMM.cis_qtl.txt.gz"))
    eqtl$is_eGene = eqtl$pval_g1 < eqtl$pval_g1_threshold &
              eqtl$qval_g1 < 0.05
    eGenes <- eqtl[eqtl$is_eGene == TRUE,]
    snp_freq <- fread(paste0(tissue[[i]], "_snp_info.frq"))
    snp_freq <- snp_freq[,c(2,5)]
    names(snp_freq) <- c("variant_id","maf")
    exp <- fread("expr_tmm_inv.bed.gz")
    sample_num <- ncol(exp) - 4
    setwd(paste0(path[[i]], "/ieQTL"))
    ct <- list.dirs(full.names = TRUE, recursive = FALSE)
    ct <- sapply(ct, function(x) unlist(strsplit(x, "\\/"))[2])
    ct <- data.frame(ct)
    ct <- ct$ct
    for (j in 1:length(ct)) {
      setwd(paste0(path[[i]], "/ieQTL/", ct[[j]]))
      if (file.exists(paste0(ct[[j]], "_LMM.cis_qtl.txt.gz"))) {
        ieqtl <- fread(paste0(ct[[j]], "_LMM.cis_qtl.txt.gz"))
        ieqtl$is_ieGene = ieqtl$pval_g2 < ieqtl$pval_g2_threshold &
                ieqtl$qval_g2 < 0.05
        ieGenes <- ieqtl[ieqtl$is_ieGene == TRUE,]
        ieGenes <- ieGenes[ieGenes$pheno_id %in% eGenes$pheno_id,]
        
        all_nom_qtl <- data.frame()
        for (chr in 1:29) {
          chr_qtl <- fread(paste0(tissue[[i]], "_", ct[[j]], "_LMM.cis_qtl_pairs.", chr, ".txt.gz"))
          #chr_qtl <- chr_qtl[chr_qtl$pval_g2 < 0.05, ]
          all_nom_qtl <- rbind(all_nom_qtl, chr_qtl)
        }
        select_qtl <- all_nom_qtl[all_nom_qtl$pheno_id %in= ieGenes$pheno_id , ]
        select_qtl <- select_qtl[,c(2,2,1,3,10,8,9)]
        names(select_qtl) <- c("rs_id","variant_id","gene_id","tss_distance","pval_nominal","slope","slope_se")
        select_qtl <- select_qtl %>% inner_join(snp_freq, by = "variant_id")
        select_qtl$N <- sample_num
        
        # Path update: protect the output target directory and dynamically switch paths using PROJECT_DIR
        out_dir_2 <- paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/Coloc/coloc/all_cell_components")
        dir.create(out_dir_2, recursive = TRUE, showWarnings = FALSE)
        setwd(out_dir_2)
        fwrite(select_qtl, paste0(tissue[[i]], "_", ct[[j]], "_coloc.bed"))
      }
    } 
}


##cell type eqtl
# Path update: dynamically extend the path using PROJECT_DIR
setwd(paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/cell_specific"))
path <- paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/cell_specific/")
tissue <- list.dirs(path, full.names = TRUE, recursive = FALSE)
tissue <- sapply(tissue, function(x) unlist(strsplit(x, "\\/"))[9])
tissue <- data.frame(tissue)
tissue <- tissue$tissue
path<-NULL
for(i in tissue){
  path[[i]]<-paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/cell_specific/", i, "/")
}


for (i in 1:length(tissue)) {
    setwd(paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/eQTL/", tissue[[i]]))
    snp_freq <- fread(paste0(tissue[[i]], "_snp_info.frq"))
    snp_freq <- snp_freq[,c(2,5)]
    names(snp_freq) <- c("variant_id","maf")
    exp <- fread("expr_tmm_inv.bed.gz")
    sample_num <- ncol(exp) - 4
    setwd(path[[i]])
    ct <- list.dirs(full.names = TRUE, recursive = FALSE)
    ct <- sapply(ct, function(x) unlist(strsplit(x, "\\/"))[2])
    ct <- data.frame(ct)
    ct <- ct$ct
    for (j in 1:length(ct)) {
      setwd(paste0(path[[i]], ct[[j]]))
      eqtl <- fread(paste0(tissue[[i]], "_", ct[[j]], "_new_LMM.cis_qtl.txt.gz"))
      eqtl$chr <- sub("_.*", "", eqtl$variant_id)
      chr_num <- unique(eqtl$chr)
      eqtl$is_eGene = eqtl$pval_g1 < eqtl$pval_g1_threshold &
                eqtl$qval_g1 < 0.05
      eGenes <- eqtl[eqtl$is_eGene == TRUE,]
      
      all_nom_qtl <- data.frame()
      for (chr in chr_num) {
        chr_qtl <- fread(paste0(tissue[[i]], "_", ct[[j]], "_new_LMM.cis_qtl_pairs.", chr, ".txt.gz"))
        #chr_qtl <- chr_qtl[chr_qtl$pval_g1 < 0.05, ]
        all_nom_qtl <- rbind(all_nom_qtl, chr_qtl)
      }
      select_qtl <- all_nom_qtl[all_nom_qtl$pheno_id %in% eGenes$pheno_id,]
      select_qtl <- select_qtl[,c(2,2,1,3,7,5,6)]
      names(select_qtl) <- c("rs_id","variant_id","gene_id","tss_distance","pval_nominal","beta","varbeta")
      select_qtl <- select_qtl %>% inner_join(snp_freq, by = "variant_id")
      select_qtl$N <- sample_num
      
      # Path update: protect the output target directory and dynamically switch paths using PROJECT_DIR
      out_dir_3 <- paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/Coloc/coloc/all_cell_specific")
      dir.create(out_dir_3, recursive = TRUE, showWarnings = FALSE)
      setwd(out_dir_3)
      fwrite(select_qtl, paste0(tissue[[i]], "_", ct[[j]], "_coloc.bed"))
    }
}


###cell state ieqtl
# Path update: dynamically extend the path using PROJECT_DIR
setwd(paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/Cell_state"))
path <- paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/Cell_state/")
tissue <- list.dirs(path, full.names = TRUE, recursive = FALSE)
tissue <- sapply(tissue, function(x) unlist(strsplit(x, "\\/"))[9])
tissue <- data.frame(tissue)
tissue <- tissue$tissue
path<-NULL
for(i in tissue){
  path[[i]]<-paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/Cell_state/", i, "/")
}

for (1:length(tissue)) {
    setwd(paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/eQTL/", tissue[[i]]))
    eqtl <- fread(paste0(tissue[[i]], "_bulk_LMM.cis_qtl.txt.gz"))
    eqtl$is_eGene = eqtl$pval_g1 < eqtl$pval_g1_threshold &
              eqtl$qval_g1 < 0.05
    eGenes <- eqtl[eqtl$is_eGene == TRUE,]
    snp_freq <- fread(paste0(tissue[[i]], "_snp_info.frq"))
    snp_freq <- snp_freq[,c(2,5)]
    names(snp_freq) <- c("variant_id","maf")
    exp <- fread("expr_tmm_inv.bed.gz")
    sample_num <- ncol(exp) - 4

    setwd(path[[i]])
    all_dirs <- list.dirs(path[[i]], recursive = TRUE, full.names = TRUE)
    bin_pattern <- "^bin[1-9]$|^bin10$"
    bin_exists <- any(grepl(bin_pattern, basename(all_dirs)))
    if (bin_exists) {
        ct <- list.dirs(full.names = TRUE, recursive = FALSE)
        ct <- sapply(ct, function(x) unlist(strsplit(x, "\\/"))[2])
        ct <- data.frame(ct)
        ct <- ct$ct
        for (j in 1:length(ct)) {
            setwd(paste0(path[[i]], ct[[j]]))
            all_dirs <- list.dirs(recursive = FALSE, full.names = TRUE)
            bin_exists <- any(grepl(bin_pattern, basename(all_dirs)))
            if (bin_exists) {
                bin <- list.dirs(full.names = TRUE, recursive = FALSE)
                bin <- sapply(bin, function(x) unlist(strsplit(x, "\\/"))[2])
                bin <- data.frame(bin)
                bin <- bin$bin
                #all_ieGenes <- data.frame()
                for (k in 1:length(bin)) {
                  setwd(paste0(path[[i]], ct[[j]], "/", bin[[k]]))
                  if (file.exists(paste0(bin[[k]], "_LMM.cis_qtl.txt.gz"))) {
                    qtl <- fread(paste0(bin[[k]], "_LMM.cis_qtl.txt.gz"))
                    if (ncol(qtl) > 15) {
                      qtl$is_ieGene = qtl$pval_g2 < qtl$pval_g2_threshold &
                            qtl$qval_g2 < 0.005
                      ieGenes <- qtl[qtl$is_ieGene == TRUE,]
                      ieGenes <- ieGenes[ieGenes$pheno_id %in% eGenes$pheno_id,]
                      #ieGenes$celltype <- ct[[j]]
                      #ieGenes$tissue <- tissue[[i]]
                      all_nom_qtl <- data.frame()
                      for (chr in 1:29) {
                        chr_qtl <- fread(paste0(bin[[k]], "_LMM.cis_qtl_pairs.", chr, ".txt.gz"))
                        #chr_qtl <- chr_qtl[chr_qtl$pval_g1 < 0.05, ]
                        all_nom_qtl <- rbind(all_nom_qtl, chr_qtl)
                      }
                      select_qtl <- all_nom_qtl[all_nom_qtl$pheno_id %in% ieGenes$pheno_id , ]
                      select_qtl <- select_qtl[,c(2,2,1,3,10,8,9)]
                      names(select_qtl) <- c("rs_id","variant_id","gene_id","tss_distance","pval_nominal","slope","slope_se")

                      #all_ieGenes <- all_ieGenes[!duplicated(all_ieGenes$gene_id), ]
                      all_ieGenes <- select_qtl %>% inner_join(snp_freq, by = "variant_id")
                      all_ieGenes$N <- sample_num
                      
                      # Path update: protect the output target directory and dynamically switch paths using PROJECT_DIR
                      out_dir_4 <- paste0(PROJECT_DIR, "/CattleGTEx/OmiGA/Coloc/coloc/all_cell_state")
                      dir.create(out_dir_4, recursive = TRUE, showWarnings = FALSE)
                      setwd(out_dir_4)
                      fwrite(all_ieGenes, paste0(tissue[[i]], "_", ct[[j]], "_", bin[[k]], "_coloc.bed"))
                    }
                  }
                }
            }
            
        }
    }
}
