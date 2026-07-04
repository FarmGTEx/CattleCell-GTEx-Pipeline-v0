#!/usr/bin/env Rscript

# ---------------- Load libraries ----------------
options(stringsAsFactors = FALSE)
suppressMessages(library(edgeR))
suppressMessages(library(preprocessCore))
suppressMessages(library(RNOmni))
suppressMessages(library(data.table))
suppressMessages(library(R.utils))
suppressMessages(library(dplyr))
suppressMessages(library(coloc))
suppressMessages(library(optparse))

# ---------------- Command line options ----------------
option_list <- list(
  make_option(c("-g", "--gwas_dir"), type = "character", help = "Directory of GWAS coloc input"),
  make_option(c("-b", "--bulk_dir"), type = "character", help = "Directory of bulk coloc input"),
  make_option(c("-r", "--result_dir"), type = "character", help = "Directory to save coloc results")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$gwas_dir) || is.null(opt$bulk_dir) || is.null(opt$result_dir)) {
  print_help(opt_parser)
  stop("All three options (--gwas_dir, --bulk_dir, --result_dir) must be specified.", call. = FALSE)
}

# ---------------- Get traits and tissues ----------------
traits <- list.dirs(opt$gwas_dir, full.names = FALSE, recursive = FALSE)
tissues <- gsub("_coloc.bed", "", list.files(opt$bulk_dir, pattern = "_coloc.bed"))

# ---------------- Main loop ----------------
for (trait in traits) {
  message("Processing trait: ", trait)

  gwas_file <- file.path(opt$gwas_dir, trait, paste0(trait, "_coloc.bed"))

  gwas <- fread(gwas_file)
  gwas <- data.frame(gwas)
  gwas$maf <- ifelse(gwas$freq > 0.5, 1 - gwas$freq, gwas$freq)

  for (tissue in tissues) {
    qtl_file <- file.path(opt$bulk_dir, paste0(tissue, "_coloc.bed"))
    qtl <- fread(qtl_file)

    same_snp <- merge(qtl, gwas, by = "variant_id", suffixes = c("_qtl", "_gwas"))
    gene_list <- unique(same_snp$gene_id)

    if (length(gene_list) > 0) {
      all_summary <- data.frame()

      for (gene in gene_list) {
        gene_qtl <- qtl[qtl$gene_id == gene, ]
        input <- merge(gene_qtl, gwas, by = "variant_id", suffixes = c("_qtl", "_gwas"))

        if (nrow(input) > 0) {
          input$varbeta1 <- input$slope_se ^ 2
          input$varbeta2 <- input$standard_error ^ 2
          input <- input[!duplicated(input$variant_id),]

          res <- coloc.abf(
            dataset1 = list(
              snp = input$variant_id,
              pvalues = input$pval_nominal,
              beta = input$slope,
              varbeta = input$varbeta1,
              N = input$N_qtl,
              type = "quant",
              MAF = input$maf_qtl
            ),
            dataset2 = list(
              snp = input$variant_id,
              pvalues = input$P_val,
              beta = input$effect_ALT_allele,
              varbeta = input$varbeta2,
              N = input$N_gwas,
              type = "quant",
              MAF = input$maf_gwas
            )
          )

          summary <- data.frame(res$summary)
          summary$gene_id <- gene
          all_summary <- rbind(all_summary, summary)
        }
      }

      out_dir <- file.path(opt$result_dir, trait)
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
      fwrite(all_summary, file.path(out_dir, paste0(tissue, "_coloc.csv")))
    }
  }
}
