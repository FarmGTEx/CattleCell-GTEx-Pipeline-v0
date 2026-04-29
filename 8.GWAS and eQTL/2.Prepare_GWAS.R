#!/usr/bin/env Rscript

# ---------- Load libraries ----------
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(optparse))

# ---------- Command line options ----------
option_list <- list(
  make_option(c("-g", "--gwas_dir"), type = "character", default = NULL,
              help = "Directory for GWAS result files", metavar = "character"),
  make_option(c("-c", "--cojo_dir"), type = "character", default = NULL,
              help = "Directory for cojo files", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = "Output directory", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ---------- Check arguments ----------
if (is.null(opt$gwas_dir) || is.null(opt$cojo_dir) || is.null(opt$output_dir)) {
  print_help(opt_parser)
  stop("All arguments (--gwas_dir, --cojo_dir, --output_dir) must be provided", call. = FALSE)
}

# ---------- Get traits ----------
trait_dirs <- list.dirs(path = opt$cojo_dir, full.names = FALSE, recursive = FALSE)
traits <- sapply(trait_dirs, function(x) unlist(strsplit(x, "\\/"))[1])
paths <- file.path(opt$cojo_dir, traits)

# ---------- Loop over traits ----------
for (i in seq_along(traits)) {
  message("Processing: ", traits[i])
  setwd(paths[i])

  gwas_file <- file.path(opt$gwas_dir, paste0(traits[i], "_china.mlma"))
  snp <- fread(gwas_file)
  snp <- snp[,c(2,4,8,9,10,7)]
  names(snp) <- c("chr", "variant_id", "beta", "se", "pval", "freq")
  snp$variant_id <- sub(":", "_", snp$variant_id)
  snp <- data.frame(snp)

  cojo_file <- file.path(paths[i], paste0(traits[i], "_cojo_res.jma.cojo"))
  if (file.exists(cojo_file)) {
    inde_snps <- fread(cojo_file)
    if (nrow(inde_snps) > 0) {
      inde_snps <- inde_snps[inde_snps$p < 1e-5, ]
      significant_snps <- snp[snp$variant_id %in% inde_snps$SNP, ]
      significant_snps$start <- significant_snps$POS - 1000000
      significant_snps$end <- significant_snps$POS + 1000000
      candidate_snps <- data.frame()

      for (k in 1:nrow(significant_snps)) {
        chr <- significant_snps$chr[k]
        start <- significant_snps$start[k]
        end <- significant_snps$end[k]
        snps_in_region <- snp[snp$chr == chr & snp$POS >= start & snp$POS <= end, ]
        candidate_snps <- rbind(candidate_snps, snps_in_region)
      }

      candidate_snps <- unique(candidate_snps)
      candidate_snps <- candidate_snps[, c(12, 8, 9, 11)]
      names(candidate_snps)[1:2] <- c("rs_id", "variant_id")

      out_file <- file.path(opt$output_dir, paste0(traits[i], "_coloc.bed"))
      fwrite(candidate_snps, out_file, sep = "\t")
    }
  }
}
