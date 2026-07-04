#!/usr/bin/env Rscript

library(data.table)
"%&%" = function(a,b) paste0(a,b)

# 保留原版基于位置的参数接收逻辑
ARGS <- commandArgs(trailingOnly = TRUE)
dir_output = ARGS[1]
subset_size = as.numeric(ARGS[2])

DF_all_pairs = NULL
# 完全保留原版的 1:18 染色体大循环，不做任何外延
for (CHR in 1:18){
    message(CHR)
    tmp = fread(dir_output %&% "/nominal_pairs_Chr" %&% CHR %&% ".MashR_input.txt.gz",sep="\t",nThread=20)
    colnames(tmp) = gsub("_Chr(.)*.nominal_pairs_zval","",colnames(tmp))
    # 完全保留原版的 rbind 堆叠逻辑
    DF_all_pairs = rbind(DF_all_pairs, tmp)
}

if (nrow(DF_all_pairs) < subset_size){
    print("# of rows: " %&% nrow(DF_all_pairs))
    random_subset = DF_all_pairs
} else {
    set.seed(9823)
    random_subset = DF_all_pairs[floor(seq(1,nrow(DF_all_pairs),length.out=subset_size)),]
}
saveRDS(random_subset, dir_output %&% "/MashR.random_subset.RDS")

cat("Done\n")
