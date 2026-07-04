process RNAEDIT {
conda '/home/huicongz/miniconda3/envs/nf-farmgtex/'
memory '16 GB'
cpus 4

publishDir mode: 'copy', path:"${params.outpath}/${sampleID}/RNAediting/", pattern: '*gz' 
input: 
tuple val(sampleID),path(bam)
path bai

output:
path("**/outTableSig_*")
path("*gz")

script:

if (file("${params.input}/${sample}/${sampleID}*.gz").size() > 1) {
"""
	samtools view -@ ${params.nthreads} -b -h $bam ${params.chrinfo}  > ${sampleID}chr.bam
	samtools sort -@ ${params.nthreads} -o ${sampleID}chr.sorted.bam  ${sampleID}chr.bam
	gatk --java-options "-Xmx8G -XX:ParallelGCThreads=16 -Djava.io.tmpdir=/tmp -XX:-UseGCOverheadLimit" MarkDuplicates --spark-runner LOCAL -I ${sampleID}chr.sorted.bam -O ${sampleID}chrmkdup.bam -M chr.mkdup_metrics.txt --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --SORTING_COLLECTION_SIZE_RATIO 0.0025 --REMOVE_DUPLICATES true
	samtools index -@ ${params.nthreads} ${sampleID}chrmkdup.bam
  python2 ${params.toolpath}/REDItools-master/main/REDItoolDenovo.py -i ${sampleID}chrmkdup.bam -o ./result \\
        -f ${params.fa} -t ${params.nthreads} -c 1 -m 255 -v 1 -q 30 -e -n 0 -u -l -W -O 5 -p -s 2 -g 2 -T 6-0 -E -V 0.05
	cat ./result/denovo*/outTable_* > ${sampleID}_RNAedit.gz
  cat ./result/denovo*/outTableSig_* > ${sampleID}_RNAeditsig.gz
"""
}else { 
"""
	samtools view -@ ${params.nthreads} -b -h $bam ${params.chrinfo}  > ${sampleID}chr.bam
	samtools sort -@ ${params.nthreads} -o ${sampleID}chr.sorted.bam  ${sampleID}chr.bam
	gatk --java-options "-Xmx8G -XX:ParallelGCThreads=16 -Djava.io.tmpdir=/tmp -XX:-UseGCOverheadLimit" MarkDuplicates --spark-runner LOCAL -I ${sampleID}chr.sorted.bam -O ${sampleID}chrmkdup.bam -M chr.mkdup_metrics.txt --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --SORTING_COLLECTION_SIZE_RATIO 0.0025 --REMOVE_DUPLICATES true
	samtools index -@ ${params.nthreads} ${sampleID}chrmkdup.bam
  python2 ${params.toolpath}/REDItools-master/main/REDItoolDenovo.py -i ${sampleID}chrmkdup.bam -o ./result \\
        -f ${params.fa} -t ${params.nthreads} -c 1 -m 255 -v 1 -q 30 -e -n 0 -u -l -W -O 5 -s 2 -g 2 -T 6-0 -E -V 0.05
  cat ./result/denovo*/outTableSig_* > ${sampleID}_RNAeditsig.gz
	cat ./result/denovo*/outTable_* > ${sampleID}_RNAedit.gz
"""
       }

}
