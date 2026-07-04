process FASTP {
conda '/home/huicongz/miniconda3/envs/nf-farmgtex/'
memory '16 GB'
cpus 8
time '4h'
publishDir mode: 'copy', path:"${params.outpath}/${sampleID}/QC/", pattern: '*json'
input:
each sampleID

output:
tuple val(sampleID),path('*.clean.fastq.gz'), emit: reads
path('*.json'), emit: json
path('*.html'), emit: html
script:

if (file("${params.input}/${sampleID}/${sampleID}*.gz").size() > 1) {
"""
fastp   -i ${params.input}/${sampleID}/*_1${params.suffix} \\
        -o ${sampleID}_1.clean.fastq.gz \\
        -I ${params.input}/${sampleID}/*_2${params.suffix} \\
        -O ${sampleID}_2.clean.fastq.gz \\
	      -f 3 -t 3 -l 36 -c -r \\
	      -W 4 -M 15 --detect_adapter_for_pe --thread 8
"""
}else { 
"""
    
fastp   -i ${params.input}/${sampleID}/*${params.suffix} \
        -o ${sampleID}.clean.fastq.gz \\
	      -f 3 -t 3 -l 36 -r \\
	      -W 4 -M 15 --thread 8 
"""
       }

}

