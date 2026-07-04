process LNCRNA {
conda '/home/huicongz/miniconda3/envs/nf-farmgtex/'
memory = '16 GB'
cache 'lenient'
cpus = 4
time = '6 h'

publishDir mode: 'copy', path:"${params.outpath}/${sampleID}/lncRNA/", pattern: '*3.gtf'
input:
tuple val(sampleID),path(bam)
output:
path('*3.gtf')

shell:
'''
stringtie \\
        -p !{params.nthreads} \\
        -G !{params.gtf} \\
        -o !{sampleID}.gtf \\
        -A stringtie_gene.tsv \\
        !{bam}
        
        samtools view -h -F 12 -q 60 !{bam} | samtools view -hb -o uniq.bam
        
        samtools index uniq.bam
        
        regtools junctions extract \
        -s XS \
        -o junctions.bed \
        uniq.bam
           
awk -F "\t|," '{print $1"\t"$2+$13+1"\t"$3-$14"\t"$6"\t"$5}' junctions.bed > junctions.5col.bed

awk 'NR==FNR&&$3=="exon"{a[$12]++}NR>FNR&&$3=="exon"&&a[$12]>1{b[$12]++;if (b[$12]==1){c=$5+1;d=$7}else{print $1"\t"c"\t"$4-1"\t"$7"\t"$10"\t"$12;c=$5+1;d=$7}}' !{sampleID}.gtf !{sampleID}.gtf > !{sampleID}.junction

awk 'NR==FNR{a[$1"\t"$2"\t"$3]=$5}NR>FNR{print $0"\t"a[$1"\t"$2"\t"$3]}' junctions.5col.bed !{sampleID}.junction > junction.junction_reads_number

awk '{a[$6]++;if ($7>=3){b[$6]++}}END{for (i in a){if (a[i]==b[i]){print i}}}' junction.junction_reads_number | awk 'NR==FNR{a[$1]++}NR>FNR&&a[$12]{print $0}' - !{sampleID}.gtf > !{sampleID}.transcripts_junction_3.gtf   
        
        
'''

}
