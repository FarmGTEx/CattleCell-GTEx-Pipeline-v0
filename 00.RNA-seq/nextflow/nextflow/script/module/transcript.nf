process Salmon{
conda '/home/huicongz/miniconda3/envs/nf-farmgtex/'
memory = '16 GB'
cpus = 10
time = '1 h'

publishDir mode: 'copy', path:"${params.outpath}/${sampleID}/transcript_expression/", pattern: '*matrix'
publishDir mode: 'copy', path:"${params.outpath}/${sampleID}/transcript_expression/THISTLE", pattern: '*thistle*'
  
input:
tuple val(sampleID),path(reads)

output:
path('*matrix')
path('*thistle*')

script:

if (reads[1]) {
"""
salmon quant  \\
        -i ${params.salmonindex}  \\
        -l A \\
        -1 ${reads[0]}  \\
        -2 ${reads[1]}  \\
        -p 10  \\
        --validateMappings \\
        -o ./  \\


        awk '{ temp = \$4; \$4 = \$5; \$5 = temp; print }' quant.sf > ${sampleID}isoform.matrix
        ${params.toolpath}/osca-0.46.1-linux-x86_64/osca-0.46.1 --efile ${sampleID}isoform.matrix  --gene-expression --make-bod --out ${sampleID}_thistle

"""
}else {
"""
salmon quant \\
        -i ${params.salmonindex} \\
        -l A \\
        -r ${reads[0]} \\
        -p ${params.nthreads} \\
        --validateMappings \\
        -o ./  \\

        awk '{ temp = \$4; \$4 = \$5; \$5 = temp; print }' quant.sf > ${sampleID}isoform.matrix
        ${params.toolpath}/osca-0.46.1-linux-x86_64/osca-0.46.1 --efile ${sampleID}isoform.matrix  --gene-expression --make-bod --out ${sampleID}_thistle

"""
        }

}

process SalmonTE{
conda '/home/huicongz/miniconda3/envs/nf-farmgtex/'
memory = '16 GB'
cpus = 10
time = '1 h'

publishDir mode: 'copy', path:"${params.outpath}/${sampleID}/TE_expression", pattern: '*gz'

input:
tuple val(sampleID),path(reads)

output:
path('*.gz')
script:

if (reads[1]) {
"""
	salmon quant  \\
    	-i ${params.TEindex}  \\
    	-l A \\
    	-1 ${reads[0]}  \\
    	-2 ${reads[1]}  \\
    	-p 10 \\
    	--validateMappings \\
    	-o ./  \\
     
	gzip quant.sf -c > ${sampleID}_TE_expression.gz
"""
}else {
"""

salmon quant \\
    	-i ${params.TEindex} \\
    	-l A \\
    	-r ${reads[0]} \\
    	-p ${params.nthreads} \\
    	--validateMappings \\
    	-o ./  \\
     
	gzip quant.sf -c  > ${sampleID}_TE_expression.gz
"""
        }

}

