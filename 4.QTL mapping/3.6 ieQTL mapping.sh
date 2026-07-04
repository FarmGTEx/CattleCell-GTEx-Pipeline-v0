#!/bin/bash

main_dir="/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/eQTL/"
for tissue_dir in "$main_dir"/*/; do
  bfile_prefix="${tissue_dir}$(basename "$tissue_dir")_maf_bfile"
  ct_dir="${tissue_dir}ieQTL"
  for subdir in "$ct_dir"/*/; do
    pheno="${subdir}$(basename "$subdir")_expr_tmm_inv.bed"
    prefix_file="$(basename "$subdir")_LMM"
    interaction_file="${subdir}$(basename "$subdir")_frac.txt"
    omiga --mode cis_interaction --genotype "$bfile_prefix" --phenotype "$pheno" --prefix "$prefix_file" --interaction "$interaction_file" --output-dir "$subdir" --verbose --debug --dprop-pc-covar 0.001 --rm-collinear-covar 0.95 --permutations 1000 --multiple-testing clipper
  done
done