process leafcutter {

memory '50 GB'
cpus 4
time '8h'

publishDir mode: 'copy', path:"${params.outpath}/${sampleID}/splicing2/", pattern: '*.junc'

input:
tuple val(sampleID),path(bam)

output:
tuple val(sampleID),path('*.junc')

script:
"""
        touch ${sampleID}.junc
        sh  ${params.toolpath}/leafcutter-master/scripts/bam2junc.sh $bam ${sampleID}.junc
"""
}

process phaser_snp{
label 'phaser'

memory '80 GB'
cpus 8
time '24h'

publishDir mode: 'copy', path:"${params.outpath}/${sampleID}/phaser2/snp/", pattern: '*haplotypic_counts.txt'
publishDir mode: 'copy', path:"${params.outpath}/${sampleID}/phaser2/snp/", pattern: '*allelic_counts.txt'
publishDir mode: 'copy', path:"${params.outpath}/${sampleID}/phaser2/snp/", pattern: '*haplotypes.txt'
input:
tuple val(sampleID),path(bam)

output:
tuple val(sampleID),path('*haplotypic_counts.txt'), emit:count

script:
  """
  samtools sort -o ${sampleID}_sorted.bam $bam
  samtools index ${sampleID}_sorted.bam
if [ -f "${params.outpath}/${sampleID}/alignment/SE.data" ];then
 python2 ${params.toolpath}/phaser-master/phaser/phaser.py \\
        --vcf ${params.recalvcf} \\
        --bam ${sampleID}_sorted.bam \\
        --paired_end 0 \\
        --mapq 255 \\
        --pass_only 0 \\
        --baseq 10 \\
        --sample ${sampleID}_bqsr \\
        --threads 8 \\
        --o ${sampleID}_phaser \\
        --output_read_ids 0 \\
        --gw_phase_vcf 1 
 else
 python2 ${params.toolpath}/phaser-master/phaser/phaser.py \\
        --vcf ${params.recalvcf} \\
        --bam ${sampleID}_sorted.bam \\
        --paired_end 1 \\
        --mapq 255 \\
        --pass_only 0 \\
        --baseq 10 \\
        --sample ${sampleID}_bqsr \\
        --threads 8 \\
        --o ${sampleID}_phaser \\
        --output_read_ids 0 \\
        --gw_phase_vcf 1 
fi
  """
}

process phaser_gene{
label 'phaser'

memory '40 GB'
cpus 4
time '8h'

publishDir mode: 'copy', path:"${params.outpath}/${sampleID}/phaser2/gene/", pattern: '*_phaser.gene.txt'

input:
tuple val(sampleID),path(count)

output:
tuple val(sampleID),path('*_phaser.gene.txt')

script:
  """

python2 ${params.toolpath}/phaser-master/phaser_gene_ae/phaser_gene_ae.py \\
        --haplotypic_counts $count \\
        --features ${params.genomepath}/Bos_taurus.ARS-UCD1.2.108.phaser.bed \\
        --o ${sampleID}_phaser.gene.txt
  """
}
