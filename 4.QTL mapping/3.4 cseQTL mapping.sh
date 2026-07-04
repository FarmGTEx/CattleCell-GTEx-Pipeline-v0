#!/bin/bash

main_dir="/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/cell_specific/Group6/"
geno_dir="/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/eQTL/Group6/"
for tissue_dir in "$main_dir"/*/; do
  bfile_prefix="${geno_dir}$(basename "$tissue_dir")/$(basename "$tissue_dir")_maf_bfile"
  for ct_dir in "$tissue_dir"/*/; do
    pheno="${ct_dir}expr_tmm_inv.bed"
    prefix_file="$(basename "$tissue_dir")_$(basename "$ct_dir")_new_LMM"
    omiga --mode cis --genotype "$bfile_prefix" --phenotype "$pheno" --prefix "$prefix_file" --output-dir "$ct_dir" --verbose --debug --dprop-pc-covar 0.001 --rm-collinear-covar 0.95 --permutations 1000 --multiple-testing clipper
  done
done
