process h5_prepare_1{
memory '20 GB'
cpus 2
time '12h'

input:
each chr

output:
path('*vcf'), emit: vcf

script:
"""
bcftools view -r ${chr} -O v -o ${chr}_cattle.vcf $params.recalvcf
"""
}


process h5_prepare_2{
memory '80 GB'
cpus 2
time '24h'

input:
path vcf

output:
path('haplotypes.h5'), emit: haplotype
path('snp_index.h5'), emit: index
path('snp_tab.h5'), emit: tab

script:
"""
 ${params.toolpath}/WASP-master/snp2h5/snp2h5 --chrom ${params.genomepath}/chromInfo.txt \\
    --format vcf --haplotype haplotypes.h5 \\
    --snp_index snp_index.h5 \\
    --snp_tab snp_tab.h5 \\
    ${vcf}
"""
}

process WASP_find {

memory '80 GB'
cpus 4
time '24h'

input:
each sampleID
path(haplotype)
path(index)
path(tab)

output:
tuple val(sampleID),path('*to.remap.bam'), emit: remapbam
tuple val(sampleID),path('*.gz'), emit: fq
tuple val(sampleID),path('*.keep.bam'), emit: keepbam

script:
"""
if [ -f "${params.outpath}/${sampleID}/alignment/SE.data" ];then
    python ${params.toolpath}/WASP-master/mapping/find_intersecting_snps.py  --is_sorted --output_dir ./ \\
    --snp_tab $tab --snp_index $index --haplotype $haplotype \\
    --samples ${sampleID}_bqsr \\
    ${params.outpath}/${sampleID}/alignment/${sampleID}Aligned.sortedByCoord.out.bam
else   
    python ${params.toolpath}/WASP-master/mapping/find_intersecting_snps.py --is_paired_end --is_sorted --output_dir ./ \\
    --snp_tab $tab --snp_index $index --haplotype $haplotype \\
    --samples ${sampleID}_bqsr \\
    ${params.outpath}/${sampleID}/alignment/${sampleID}Aligned.sortedByCoord.out.bam
    rm *remap.single.fq.gz
fi
"""
}

process WASP_remap {

memory '100 GB'
cpus 12
time '12h'

input:
tuple val(sampleID),path(fq)
tuple val(sampleID),path(toremap)
tuple val(sampleID),path(findkeep)

output:
tuple val(sampleID),path('*.sorted.merge.bam'), emit: bam

script:
""" 
    STAR   --runMode alignReads \\
        --genomeDir  $params.starindex \\
        --sjdbGTFfile $params.gtf \\
        --readFilesIn $fq \\
        --runThreadN 12 \\
        --outSAMunmapped Within \\
        --readFilesCommand zcat \\
        --outSAMtype BAM Unsorted \\
        --outFileNamePrefix ${sampleID}_remap \\
        --twopassMode Basic \\
	--outFilterMismatchNmax 3 \\
	      --chimSegmentMin 10 \\
	      --chimOutJunctionFormat 1 \\
	      --outFilterType BySJout \\
	      --alignSJoverhangMin 8 \\
	      --alignSJDBoverhangMin 1 
    samtools sort -@ 5 ${sampleID}_remapAligned.out.bam \
		-o ${sampleID}_remap_sorted.bam
	  samtools index  *remap_sorted.bam
     
    python ${params.toolpath}/WASP-master/mapping/filter_remapped_reads.py $toremap ${sampleID}_remap_sorted.bam ${sampleID}_keep.bam
    samtools merge ${sampleID}.keep.merge.bam $findkeep ${sampleID}_keep.bam
    samtools sort -o ${sampleID}.sorted.merge.bam ${sampleID}.keep.merge.bam
     
"""
}

process WASP_dedup {

memory '50 GB'
cpus 4
time '12h'


input:
tuple val(sampleID),path(bam)

output:
tuple val(sampleID),path('*_WASP.bam'), emit: bam

script:
"""
samtools index $bam
if [ -f "${params.outpath}/${sampleID}/alignment/SE.data" ];then
python ${params.toolpath}/WASP-master/mapping/rmdup.py $bam ${sampleID}_WASP.bam

else
python ${params.toolpath}/WASP-master/mapping/rmdup_pe.py $bam ${sampleID}_WASP.bam
fi
"""
}