process BQSR{
conda '/home/huicongz/miniconda3/envs/nf-farmgtex/'
memory = '36 GB'
cpus = 8
time = '1 d'

publishDir mode: 'copy', path:"${params.workpath}/bqsr_bamlist/", pattern: '*bqsr.bam'
publishDir mode: 'copy', path:"${params.workpath}/bqsr_bamlist/", pattern: '*bqsr.bai'
publishDir mode: 'copy', path:"${params.workpath}/bqsr_bamlist/md5", pattern: '*.md5'
input:
tuple val(sampleID),path(bam)
path bai

output:
path ('*bqsr.bam') 
path ('*bqsr.bai') 
path ('*.md5') 

script:
  """
        gatk --java-options "-Xmx8G -XX:ParallelGCThreads=8 -Djava.io.tmpdir=/tmp" AddOrReplaceReadGroups -I $bam -O ${sampleID}.addrg.bam -LB ${sampleID} -PL ILLUMINA -PU ${sampleID} -SM ${sampleID} -ID ${sampleID} -FO ${sampleID} --CREATE_INDEX true
	gatk --java-options "-Xmx8G -XX:ParallelGCThreads=4 -Djava.io.tmpdir=/tmp -XX:-UseGCOverheadLimit" MarkDuplicates --spark-runner LOCAL -I ${sampleID}.addrg.bam -O ${sampleID}.mkdup.bam -M mkdup_metrics.txt --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --SORTING_COLLECTION_SIZE_RATIO 0.0025 --REMOVE_DUPLICATES true
	gatk --java-options "-Xmx36G -XX:ParallelGCThreads=8 -Djava.io.tmpdir=/tmp" SplitNCigarReads --spark-runner LOCAL -I ${sampleID}.mkdup.bam -R ${params.fa} -O ${sampleID}.cigar.bam --create-output-bam-index true --max-reads-in-memory 1000
	gatk --java-options "-Xmx36G -XX:ParallelGCThreads=8 -Djava.io.tmpdir=/tmp" BaseRecalibrator --spark-runner LOCAL -I ${sampleID}.cigar.bam --known-sites ${params.vcf} -O ${sampleID}.bqsr.table -R ${params.fa} 
	gatk --java-options "-Xmx36G -XX:ParallelGCThreads=8 -Djava.io.tmpdir=/tmp" ApplyBQSR --spark-runner LOCAL -I ${sampleID}.cigar.bam --bqsr-recal-file ${sampleID}.bqsr.table -O ${sampleID}_bqsr.bam --create-output-bam-index true --create-output-bam-md5 true 

  """
}
