nextflow.enable.dsl=2 

params.workpath = "/faststorage/project/farmgtex/gtest"

params.genomepath = "${params.workpath}/cattle_genome/"
params.enhancerbed = "${params.genomepath}/bosTau9_E6.bed"
params.align = "${params.genomepath}/Bos_taurusARS-UCD.fa.align"
params.gtf = "${params.genomepath}/Bos_taurus.ARS-UCD1.2.108.chr.gtf"
params.fa = "${params.genomepath}/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa"
params.cdna = "${params.genomepath}/Bos_taurus.ARS-UCD1.2.cdna.all.fa"
params.vcf = "${params.genomepath}/bos_taurus.vcf"

params.enhancer = "${params.genomepath}/enhancer.saf"
params.dict = "${params.genomepath}/Bos_taurus.ARS-UCD1.2.dna.toplevel.dict"
params.fai = "${params.genomepath}/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.fai"

params.starindex = "${params.genomepath}/indexstar/"
params.salmonindex = "${params.genomepath}/indexsalmon/"
params.TEindex = "${params.genomepath}/indexTE/"
params.toolpath = "${params.workpath}/tools/"
params.nthreads = 16
params.chr = "${params.genomepath}/chr.info"
params.snpdb = "${params.genomepath}/snpdb"
params.species = "cow"

process fastadic {
memory '8 GB'
cpus 2

publishDir path: "${params.genomepath}", pattern: '.'

script:
"""
  samtools dict -o ${params.dict} ${params.fa}
  samtools faidx ${params.fa}
  gatk IndexFeatureFile -I ${params.vcf} 
  bgzip -c ${params.vcf} > ${params.vcf}.gz
  tabix -p vcf ${params.vcf}.gz 
"""

}

process chr {
publishDir path: "${params.genomepath}", pattern: '.'

script:
"""
cat ${params.gtf}|grep -v "#" |cut -f 1|sort|uniq|tr '\n' ' ' > ${params.chr}

"""

}

process snpdb {

publishDir path: "${params.genomepath}", pattern: '.'

script:
"""
cat ${params.vcf}|awk '{print \$1"\\t"\$2}' > ${params.snpdb}

"""

}

process STARINDEX {
memory '50 GB'
cpus 16

output:
publishDir path: "${params.starindex}", pattern: '.'

script:
"""
    STAR --runThreadN ${params.nthreads} \\
        --runMode genomeGenerate \\
        --genomeDir ${params.starindex} \\
        --genomeFastaFiles $params.fa \\
        --sjdbGTFfile $params.gtf
"""
}

process SALMONINDEX {
memory '50 GB'
cpus 16

output:
publishDir path: "${params.salmonindex}", pattern: '.'

script:
"""
    bash ${params.toolpath}/SalmonTools-master/scripts/generateDecoyTranscriptome.sh \\
    -a $params.gtf \\
    -g $params.fa \\
    -j 16 \\
    -t $params.cdna \\
    -o temp/
    
    salmon index -t temp/gentrome.fa -p 16 -d temp/decoys.txt -i ${params.salmonindex} \\
    
"""
}

process RNAstablity {
memory '10 GB'
cpus 2

publishDir mode: 'copy', path: "${params.genomepath}", pattern: '*.gtf'
publishDir mode: 'copy', path: "${params.workpath}", pattern: '*.txt'

output :
path ('*.gtf')
path ('*.txt')

script:

"""
cat ${params.gtf} |grep -v '#' | awk -v FS='\t' '\$3=="exon" { exonName=\$1":"\$4":"\$5":"\$7; split(\$9, fields, ";"); geneName=fields[1]; transcriptName=fields[3]; printf("%s\\t%s\\t%s\\n",exonName,geneName,transcriptName); }' | sort | uniq | awk -v FS='\t' '{ eCount[\$1]++; tCount[\$3]++; exonHost[\$1]=\$2; if(tCount[\$3]==1) gCount[\$2]++; } END { for(i in eCount) if(eCount[i]==gCount[exonHost[i]]) { split(i,fields,":"); printf("%s\\tensembl\\texon\\t%s\\t%s\\t.\\t%s\\t.\\t%s;\\n",fields[1],fields[2],fields[3],fields[4],exonHost[i]); } }' > exon.gtf

cat ${params.gtf} |grep -v '#' | awk -v FS='\t' '\$3=="exon" { exonName=\$1":"\$4":"\$5":"\$7; split(\$9, fields, ";"); geneName=fields[1]; transcriptName=fields[3]; printf("%s\\t%s\\t%s\\n",exonName,geneName,transcriptName); }' | sort | uniq | awk -v FS='\t' '{ eCount[\$1]++; tCount[\$3]++; exonHost[\$1]=\$2; if(tCount[\$3]==1) gCount[\$2]++; } END { for(i in eCount) { split(i,fields,":"); printf("%s\\tensembl\\texon\\t%s\\t%s\\t.\\t%s\\t.\\t%s;\\n",fields[1],fields[2],fields[3],fields[4],exonHost[i]); } }' | bedtools sort -i stdin | awk -v FS='\t' '{ if( last_exon[\$9]==1 && (last_exon_end[\$9]+1)<(\$4-1) ) printf("%s\\t%s\\tintron\\t%i\\t%i\\t%s\\t%s\\t%s\\t%s\\n",\$1,\$2,last_exon_end[\$9]+1,\$4-1,\$6,\$7,\$8,\$9); last_exon[\$9]=1; last_exon_end[\$9]=\$5; }' > intron.gtf

echo -e "Label\tFile\tReadType\tBatch" > metadata_for_RNAstability.txt

"""
}

process UTR {
memory '10 GB'
cpus 2

publishDir mode: 'copy', path: "${params.genomepath}", pattern: '*utr.bed'

output : 
path ('*utr.bed')

script:
"""
${params.toolpath}gtfToGenePred -genePredExt ${params.gtf} cow.genepred
${params.toolpath}genePredToBed cow.genepred cow.genepred.bed
cut -f1,12 cow.genepred.bed > cow.IDmapping.txt
python ${params.toolpath}/DaPars2-master/src/DaPars_Extract_Anno.py -b cow.genepred.bed -s cow.IDmapping.txt -o cow.3utr.bed
"""
}

process enhancer {
publishDir mode: 'copy', path: "${params.genomepath}", pattern: '*.saf'
output :
path ('*.saf')
shell:
'''
grep -v '#' !{params.gtf} | awk -v OFS="\t" '{print $1,$4,$5,$7}' | sort -k1,1 -k2,2n | uniq > cow.sort.bed
bedtools intersect -a !{params.enhancerbed} -b cow.sort.bed -wa -v  | awk -v OFS="\t" '{print $1,$2,$3}' | uniq | awk '{{print "enhancer_"$1":"$2":"$3"\t"$1"\t"$2"\t"$3"\t."}}' | sed '1i GeneID\tChr\tStart\tEnd\tStrand' > enhancer.saf
'''

}

process REPEATMASKER {
memory '100 GB'
cpus 16
time '3d'

output:
path('*align')

script:
"""
    RepeatMasker -pa ${params.nthreads} -species ${params.species} -xsmall -a -s -no_is -cutoff 255 -frag 20000 -dir ./  -gff ${params.fa}
"""
}


process TE_fasta {
memory '20 GB'
cpus 8

input:
path(align)
output:
path('*fa')

shell:
'''
python3 !{params.toolpath}/TEdetectionEvaluation-main/helper_scripts/align_parser.py -a !{align} -o ./

sed -i "s/chr//" *.bed

bedtools getfasta -fi !{params.fa} -bed *.bed -name | awk 'BEGIN{FS="::"}{if($1~">"){print $1"\\n" }else{print $0}}' > repeats.fa

'''
}

process TE_index{
memory '50 GB'
cpus 16

input:
path TE_fa

output:
publishDir path: "${params.TEindex}", pattern: '.'

script:

"""
    salmon index \
    -t ${TE_fa} \
    -i ${params.TEindex} \
    
"""
}



workflow {
  fastadic()
  chr()
  snpdb()
  RNAstablity()
  UTR()
}