options(stringsAsFactors = FALSE)
suppressMessages(library(edgeR))
suppressMessages(library(preprocessCore))
suppressMessages(library(RNOmni))
suppressMessages(library(data.table))
suppressMessages(library(R.utils))
suppressMessages(library(SNPRelate))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))   # 引入命令行解析包
suppressMessages(library(rtracklayer)) # 引入依赖包

# ====================================================================
# 1. Command-line argument parsing (CLI wrapper)
# ====================================================================
option_list = list(
  make_option(c("-p", "--project_dir"), type="character", default=NULL, help="Project root directory path"),
  make_option(c("-g", "--gtf_file"), type="character", default=NULL, help="Path to the species-specific genome GTF annotation file")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$project_dir) || is.null(opt$gtf_file)){
  print_help(opt_parser)
  stop("Error: both --project_dir and --gtf_file must be provided to specify the input path and reference genome！", call.=FALSE)
}

# ====================================================================
# 2. Core data import and critical bug fixes
# ====================================================================

## Critical fix 1: Load the GTF file at the beginning of the workflow to ensure that the downstream variable is properly initialized
message("Loading the reference genome GTF file...")
gtf_obj <- rtracklayer::import(opt$gtf_file)
gtf <- as.data.frame(gtf_obj)

## Critical fix 2: Instead of scanning the current directory, directly locate the project path
bulk_exp_dir <- file.path(opt$project_dir, "CattleGTEx/Cattle_bulk_exp")
setwd(bulk_exp_dir)

file <- list.files(pattern = ".txt")
list0 <- tools::file_path_sans_ext(file)
list0 <- sapply(list0, function(x) unlist(strsplit(x, "\\_"))[2])
list0 <- data.frame(list0)
tissue <- list0$list0

list1<-NULL
for(i in tissue){
  list1[[i]]<-paste0(opt$project_dir, "/CattleGTEx/Cattle_bulk_exp/Bulk_", i, ".txt") ##download in Cattle-CellGTEx website (https://file.kiz.ac.cn/cattlecellgtex/RNA-seq%20and%20Genotype/CattleCell-GTEx_v0.Bulk_TPM.zip)
}
list2<-NULL
for(i in tissue){
  list2[[i]]<-paste0(opt$project_dir, "/CattleGTEx/Cattle_bulk_TPM/Bulk_", i, ".txt") ##download in Cattle-CellGTEx website (https://file.kiz.ac.cn/cattlecellgtex/RNA-seq%20and%20Genotype/CattleCell-GTEx_v0.Bulk_TPM.zip)
}

# ====================================================================
# 3. # Core computation and data-filtering logic
# ====================================================================
for (i in 1:length(tissue)) {
  bulk <- read.csv(list1[[i]], sep = "\t")
  bulk <- na.omit(bulk)
  TPM <- read.csv(list2[[i]], sep = "\t")
  TPM <- na.omit(TPM)
  if (ncol(bulk) >= 40) {
    Counts = counts = bulk
    samids = colnames(Counts) # sample id
    expr_counts = Counts
    expr = DGEList(counts=expr_counts) # counts
    nsamples = length(samids) # sample number
    ngenes = nrow(expr_counts) 
    y = calcNormFactors(expr, method="TMM")
    TMM = cpm(y,normalized.lib.sizes=T)

    count_threshold = 6
    tpm_threshold = 0.1
    sample_frac_threshold = 0.2
    sample_count_threshold = 10
    expr_tpm = TPM[rownames(expr_counts),samids]
    tpm_th = rowSums(expr_tpm >= tpm_threshold)
    count_th = rowSums(expr_counts >= count_threshold)
    ctrl1 = tpm_th >= (sample_frac_threshold * nsamples)
    ctrl2 = count_th >= (sample_frac_threshold * nsamples)
    mask = ctrl1 & ctrl2
    TMM_pass = TMM[mask,]
    rank_qnorm <- function(x) {
      result <- qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
      return(result)
    }
    TMM_inv = t(apply(TMM_pass, MARGIN = 1, FUN = rank_qnorm))

    region_annot <- gtf
    geneid = region_annot$gene_id
    expr_matrix = TMM_inv[rownames(TMM_inv) %in% geneid,]


    bed_annot = region_annot[region_annot$gene_id %in% rownames(expr_matrix),]
    bed = data.frame(bed_annot,expr_matrix[bed_annot$gene_id,])
    bed = bed[bed[,1] %in% as.character(1:30),]
    bed[,1] = as.numeric(bed[,1])
    bed = bed[order(bed[,1],bed[,2]),]
    colnames(bed)[1] = "#Chr"

    bed <- bed[bed$type == "gene",]
    bed <- bed[bed$gene_biotype == "protein_coding" | bed$gene_biotype == "lncRNA",]
    start = bed$start[bed$strand == "-"]
    end = bed$end[bed$strand == "-"]
    bed$start[bed$strand == "-"] = end
    bed$end[bed$strand == "-"] = start
    bed_unique <- bed[!duplicated(bed$gene_id), ]
    rownames(bed_unique) = bed_unique$gene_id
    bed_unique1 <- bed_unique[,-c(1:25)]
    bed_unique <- bed_unique[,-c(4:9,11:25)]
    colnames(bed_unique1) <- colnames(bulk)
    bed <- data.frame(bed_unique[,c(1:4)],bed_unique1)
    names(bed)[1:4] <- c("#Chr","start","end","gene_id")
    
    # Dynamically switch to and save outputs in the specified tissue directory
    out_dir <- file.path(opt$project_dir, "CattleGTEx/OmiGA/eQTL", tissue[[i]])
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    setwd(out_dir)
    
    fwrite(bed, file = "expr_tmm_inv.bed", sep = "\t")
    system("bgzip expr_tmm_inv.bed")
    message(paste("[成功] 组织", tissue[[i]], "处理完毕，数据已压缩落盘。"))
  }
}
