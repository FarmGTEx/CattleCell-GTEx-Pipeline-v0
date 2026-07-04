#!/bin/bash

Rscript 1.Prepare_eQTL.R \
  --input_dir "$eQTL_path" \
  --output_dir "$output_path1"

Rscript 2.Prepare_GWAS \
  --gwas_dir "$GWAS_path" \
  --cojo_dir "$cojo_path" \
  --output_dir "$GWAS_path"

Rscript 3.run_coloc \
  --gwas_dir "$GWAS_path" \
  --bulk_dir "$output_path1" \
  --result_dir "$output_path2"