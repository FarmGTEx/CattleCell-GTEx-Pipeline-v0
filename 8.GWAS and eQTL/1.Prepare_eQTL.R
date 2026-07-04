#!/usr/bin/env Rscript

# --------- Load Libraries ----------
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(optparse))

# --------- Define Command-line Options ----------
option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", default = NULL,
              help = "Input base directory for tissues (e.g., /path/to/OmiGA/eQTL)", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = "Output directory for .coloc.bed files", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_dir) || is.null(opt$output_dir)) {
  print_help(opt_parser)
  stop("Both --input_dir and --output_dir must be specified", call. = FALSE)
}

input_dir <- opt$input_dir
output_dir <- opt$output_dir

# --------- Get Tissue Names ----------
tissue <- list.dirs(input_dir, full.names = FALSE, recursive = FALSE)
path <- file.path(input_dir, tissue)

# --------- Loop Through Tissues ----------
for (i in seq_along(tissue)) {
  message("Processing tissue: ", tissue[i])
  setwd(path[i])

  eqtl <- fread(paste0(path[i], "/", tissue[i], ".cis_qtl.txt.gz"))
  eqtl$is_eGene <- eqtl$pval_g1 < eqtl$pval_g1_threshold & eqtl$qval_g1 < 0.05
  eGenes <- eqtl[eqtl$is_eGene == TRUE, ]

  snp_freq <- fread(paste0(path[i], "/", tissue[i], "_snp_info.frq"))
  snp_freq <- snp_freq[, c(2, 5)]
  names(snp_freq) <- c("variant_id", "maf")

  exp <- fread("expr_tmm_inv.bed.gz")
  sample_num <- ncol(exp) - 4

  all_nom_qtl <- data.frame()

  for (chr in 1:29) {
    chr_qtl <- fread(paste0(path[i], "/", tissue[i], ".cis_qtl_pairs.", chr, ".txt.gz"))
    all_nom_qtl <- rbind(all_nom_qtl, chr_qtl)
  }

  select_qtl <- all_nom_qtl[all_nom_qtl$pheno_id %in% eGenes$pheno_id, ]
  select_qtl <- select_qtl[, c("rs_id", "variant_id", "gene_id", "tss_distance", "pval_nominal", "slope", "slope_se")]
  select_qtl <- select_qtl %>%
    inner_join(snp_freq, by = "variant_id")

  out_file <- file.path(output_dir, paste0(tissue[i], "_coloc.bed"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(select_qtl, out_file, sep = "\t")
}
