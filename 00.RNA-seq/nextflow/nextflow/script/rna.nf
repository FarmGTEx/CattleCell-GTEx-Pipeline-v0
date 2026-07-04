params.input = "/home/huicongz/farmgtex/fastqfile/"
params.workpath = "/home/huicongz/farmgtex/gtex"
params.samplelist = "${params.workpath}/samplelist.csv"
params.genomepath = "${params.workpath}/cattle_genome/"
params.enhancerbed = "${params.genomepath}/bosTau9_E6.bed"
params.suffix =".fastq.gz"
params.nthreads = 20
params.sampleIDs = file("$params.samplelist").readLines().collect { it.split('\t')[0] }.toSet()
params.enhancer = "${params.genomepath}/enhancer.saf"
params.gtf = "${params.genomepath}/Bos_taurus.ARS-UCD1.2.108.chr.gtf"
params.fa = "${params.genomepath}/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa"
params.vcf = "${params.genomepath}/bos_taurus.vcf"
params.outpath = "${params.workpath}/result_new/"
params.toolpath = "${params.workpath}/tools/"
params.starindex = "${params.genomepath}/indexstar/"
params.salmonindex = "${params.genomepath}/indexsalmon/"
params.TEindex = "${params.genomepath}/indexTE/"
params.exon = "${params.genomepath}/exon.gtf"
params.intron = "${params.genomepath}/intron.gtf"
params.chr = "${params.genomepath}/chr.info"
params.snpdb = "${params.genomepath}/snpdb"
params.chrinfo = file("${params.chr}").text.trim()

  include { FASTP                } from './module/fastp.nf'
  include { STAR                 } from './module/star.nf'
  include { GENE_string          } from './module/expression.nf'
  include { GENE_feat            } from './module/expression.nf'
  include { EXON                 } from './module/expression.nf'
  include { enhancer             } from './module/expression.nf'
  include { RNA                  } from './module/expression.nf'
  include { Salmon               } from './module/transcript.nf'
  include { SalmonTE             } from './module/transcript.nf'
  include { BQSR                 } from './module/bqsr.nf'
  include { RNAEDIT              } from './module/RNAedit.nf'
  include { Bedgraph             } from './module/expression.nf'
  include { LNCRNA               } from './module/LNCRNA.nf'
  
workflow{
  FASTP(params.sampleIDs)
  STAR(FASTP.out.reads)
  Salmon(FASTP.out.reads)
  Bedgraph(STAR.out.bam,STAR.out.log,params.samplelist)
  GENE_string(STAR.out.bam)
  GENE_feat(STAR.out.bam)
  enhancer(STAR.out.bam)
  EXON(STAR.out.bam)
  RNA(STAR.out.bam)
  BQSR(STAR.out.bam,STAR.out.bai)
  LNCRNA(STAR.out.bam)
}
