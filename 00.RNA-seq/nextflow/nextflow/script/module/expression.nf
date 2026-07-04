process GENE_string {
publishDir mode: 'copy', path:"${params.outpath}/${sampleID}/gene_expression/", pattern: '*.tsv'
conda '/home/huicongz/miniconda3/envs/nf-farmgtex/'
memory = '16 GB'
cpus = 4
time = '2 h'

input:
tuple val(sampleID),path(bam)
output:
path ('*.tsv')

script:

"""
        stringtie \\
        -p 4 \\
        -e -B -G ${params.gtf} \\
        -o temp.gtf \\
        -A ${sampleID}_stringtie_gene.tsv \\
        $bam
"""
}

process GENE_feat {
conda '/home/huicongz/miniconda3/envs/nf-farmgtex/'
memory = '8 GB'
cpus = 4
time = '1 h'

publishDir mode: 'copy', path:"${params.outpath}/${sampleID}/gene_expression/", pattern: '*.summary'
publishDir mode: 'copy', path:"${params.outpath}/${sampleID}/gene_expression/", pattern: '*.gz'
publishDir mode: 'copy', path:"${params.outpath}/${sampleID}/alignment/", pattern: '*.data'
input:
tuple val(sampleID),path(bam)
output:
path ('*.gz')
path ('*.summary')
path ('*.data')

script:
        
if (file("${params.input}/${sampleID}/${sampleID}*.gz").size() > 1) {
"""
       featureCounts \\
        -T 4 \\
        -t exon -p -g gene_id \\
        -a ${params.gtf} \\
        -o ${sampleID}_featcount_gene.tsv \\
        $bam
	gzip ${sampleID}_featcount_gene.tsv 
  touch PE.data
"""
}else { 
"""
       featureCounts \\
        -T 4 \\
        -t exon -g gene_id \\
        -a ${params.gtf} \\
        -o ${sampleID}_featcount_gene.tsv \\
        $bam
	gzip ${sampleID}_featcount_gene.tsv 
  touch SE.data
"""
       }
}

process EXON{
conda '/home/huicongz/miniconda3/envs/nf-farmgtex/'
memory = '8 GB'
cpus = 4
time = '1 h'


publishDir mode: 'copy', path:"${params.outpath}/${sampleID}/exon_expression/", pattern: '*.summary'
publishDir mode: 'copy', path:"${params.outpath}/${sampleID}/exon_expression/", pattern: '*.gz'

input:
tuple val(sampleID),path(bam)

output:
path ('*.gz')
path ('*.summary')

script:
        
if (file("${params.input}/${sampleID}/${sampleID}*.gz").size() > 1) {
"""
       featureCounts \\
        -T 4 \\
        -t exon -p -g gene_id \\
        -f -a ${params.gtf} \\
        -o ${sampleID}_featcount_exon.tsv \\
        $bam
	gzip ${sampleID}_featcount_exon.tsv 
"""
}else { 
"""
       featureCounts \\
        -T 4 \\
        -t exon -g gene_id \\
        -f -a ${params.gtf} \\
        -o ${sampleID}_featcount_exon.tsv \\
        $bam
	gzip ${sampleID}_featcount_exon.tsv 
"""
       }
}


process RNA{
conda '/home/huicongz/miniconda3/envs/nf-farmgtex/'
memory = '4 GB'
cpus = 2
time = '1d'

publishDir mode: 'copy', path:"${params.outpath}/${sampleID}/rna_stablity_HTseq/", pattern: '*'

input:
tuple val(sampleID),path(bam)
output:
path ('*.txt')

script:
"""
htseq-count -m intersection-strict -s no -f bam -t exon -n 2 $bam ${params.exon} > ${sampleID}_Exon_counts.txt
htseq-count -m union -f bam -t intron -s no $bam -n 2 ${params.intron} > ${sampleID}_Introns_counts.txt
"""

}

process enhancer{
conda '/home/huicongz/miniconda3/envs/nf-farmgtex/'
memory = '8 GB'
cpus = 4
time = '1 h'

publishDir mode: 'copy', path:"${params.outpath}/${sampleID}/enhancer/", pattern: '*'

input:
tuple val(sampleID),path(bam)
output:
path ('*.tsv')
path ('*.summary')
script:
        
if (file("${params.input}/${sampleID}/${sampleID}*.gz").size() > 1) {
"""
       featureCounts \\
        -T 4 \\
        -p -F SAF -a ${params.enhancer} \\
        -o ${sampleID}_enhancer.tsv \\
        $bam
"""
}else { 
"""
       featureCounts \\
        -T 4 \\
        -F SAF -a ${params.enhancer} \\
        -o ${sampleID}_enhancer.tsv \\
        $bam
"""
       }

}

process Bedgraph{
conda '/home/huicongz/miniconda3/envs/nf-farmgtex/'
memory = '16 GB'
cpus = 4
time = '1 h'

publishDir mode: 'copy', path:"${params.outpath}/wig/", pattern: '*.wig'
publishDir mode: 'copy', path:"${params.workpath}/depar2/", pattern: '*mapped_reads.txt'

input:
tuple val(sampleID),path(bam)
path log
path tsv
output:
path ('*.wig')
path ('*.txt')
shell:
'''
        bedtools genomecov -bga -split -trackline -ibam !{bam} > !{sampleID}.bam.wig
        tissue=`grep -w "^!{sampleID}" !{tsv} | cut -f2 | sort -u`
        depth=`grep 'Uniquely mapped reads number' !{log} | cut -f2`
        echo -e "!{sampleID}\t$tissue\t!{params.outpath}/wig/!{sampleID}.bam.wig\t$depth" >  !{sampleID}.mapped_reads.txt
'''
}