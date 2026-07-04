#####
step5 Merge all filtered chromosome-wide VCF files
#####

cd ${VCF_path}
ls -1v   recal.simplified_filter_*.bcf >  impute.txt

GLIMPSE2_ligate \
          --input impute.txt \
          --output recal.final.bcf \
          --threads 12
