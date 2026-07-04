#####
step4 Merge imputed bins into chromosome-wide VCFs and filter variants
#####

cd ${VCF_path}
ls -1v  ${chr}_*.bcf >  ${chr}.txt

        GLIMPSE2_ligate \
            --input ${chr}.txt \
            --output impute${chr}.bcf \
            --threads 12

 python3 recalINFO.py --invcf impute${chr}.bcf > recal_${chr}.vcf

 bgzip -f recal_${chr}.vcf
 tabix -p vcf recal_${chr}.vcf.gz

 bcftools annotate -x ID,FORMAT/GP,FORMAT/DS,INFO/RAF  -O b -o recal.simplified_${chr}.bcf recal_${chr}.vcf --threads 12
 bcftools view -i 'INFO > 0.6 & AF > 0.05 & AF < 0.95' recal.simplified_${chr}.bcf -O b -o recal.simplified_filter_${chr}.bcf --threads 12
 bcftools index recal.simplified_filter_${chr}.bcf