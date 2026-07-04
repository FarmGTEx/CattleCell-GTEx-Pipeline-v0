#!/bin/bash
main_dir="/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/eQTL/"

for subdir in "$main_dir"/*/; do
  bfile_prefix="${subdir}/$(basename "$subdir")_maf_bfile"
  pheno="${subdir}expr_tmm_inv.bed.gz"
  prefix_file="$(basename "$subdir")_LMM"
  omiga --mode cis --genotype "$bfile_prefix" --phenotype "$pheno" --prefix "$prefix_file" --output-dir "$subdir" --verbose --debug --dprop-pc-covar 0.001 --rm-collinear-covar 0.95 --permutations 1000 --multiple-testing clipper
done

