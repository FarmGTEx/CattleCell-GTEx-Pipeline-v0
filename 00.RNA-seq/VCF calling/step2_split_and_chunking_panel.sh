#####
step2 Chunking and splitting the phased reference panel for GLIMPSE2 imputation
#####

cat chrlist | while read chr            #for each chromosome
do
GLIMPSE2_chunk --input AN_AC.${chr}.phased.vcf.gz --region ${chr} --window-mb 5.0 --buffer-mb 1.0 --output chunks.${chr}.txt --sequential --threads 12
while IFS="" read -r LINE || [ -n "$LINE" ];
do
printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
    IRG=$(echo $LINE | cut -d" " -f3)
    ORG=$(echo $LINE | cut -d" " -f4)
    GLIMPSE2_split_reference --reference AN_AC.${chr}.phased.vcf.gz --input-region ${IRG} --output-region ${ORG} --output ./split/split
done < chunks.$chr.txt
done

cd ${ref_path}/split/
for bin in `ls -1v split_*.bin` ; do basename $bin .bin ; done | sed 's/^split\_//g' > bins.txt

