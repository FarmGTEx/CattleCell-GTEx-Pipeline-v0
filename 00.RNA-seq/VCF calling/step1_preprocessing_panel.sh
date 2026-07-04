#####
step1 Preprocessing the genotyping imputation panel
#####
conda activate nf-farmgtex
cd ${ref_path}  #the path to the reference panel
ls *.vcf.gz | while read file
do
tabix -p vcf ${file}
bcftools +fill-tags -Oz -o ${ref_path}/AN_AC.${file} -- $file -t AN,AC
bcftools index /faststorage/project/farmgtex/gtex/cattle_ref/SNPs/AN_AC.${file}
done