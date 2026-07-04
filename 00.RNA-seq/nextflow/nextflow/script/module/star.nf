process STAR {
conda '/home/huicongz/miniconda3/envs/nf-farmgtex/'
cache 'lenient'
memory '80 GB'
cpus 12
time '2d'

publishDir mode: 'copy', path:"${params.outpath}/${sampleID}/alignment/", pattern: '*Log.final.out'
publishDir mode: 'copy', path:"${params.outpath}/${sampleID}/alignment/", pattern: '*sortedByCoord.out.bam'
input:
tuple val(sampleID),path(reads)
output:
tuple val(sampleID),path('*sortedByCoord.out.bam'), emit: bam
path('*Log.final.out'), emit: log
path('*bai'),emit: bai
path('*Log.out'), emit: log_out
path('*Log.progress.out'), emit: log_progress
path('*.tab'), emit: tab

script:
"""
STAR     --runMode alignReads \\
        --genomeDir  $params.starindex \\
        --sjdbGTFfile $params.gtf \\
        --readFilesIn $reads \\
        --runThreadN 12 \\
        --outSAMunmapped Within \\
        --readFilesCommand zcat \\
        --outSAMtype BAM Unsorted \\
        --outFileNamePrefix ${sampleID} \\
	--outFilterMismatchNmax 999 \\
        --twopassMode Basic \\
	      --chimSegmentMin 10 \\
	      --chimOutJunctionFormat 1 \\
	      --outFilterType BySJout \\
	      --alignSJoverhangMin 8 \\
	      --alignSJDBoverhangMin 1 
  
  samtools sort -@ 5 ${sampleID}Aligned.out.bam \
		-o ${sampleID}Aligned.sortedByCoord.out.bam
	samtools index  *sortedByCoord.out.bam
        
"""

}
