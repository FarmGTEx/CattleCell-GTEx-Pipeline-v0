params.input = "/home/huicongz/farmgtex/fastqfile/"
params.workpath = "/home/huicongz/farmgtex/gtex"
params.samplelist = "${params.workpath}/bqsr1.csv"
params.genomepath = "${params.workpath}/cattle_genome/"
params.enhancerbed = "${params.genomepath}/bosTau9_E6.bed"
params.WASP = "${params.genomepath}/h5/"
params.suffix =".fastq.gz"
params.starindex = "${params.genomepath}/indexstar/"
params.recalvcf = "${params.workpath}/VCFfile/chr_ind2/recal_simplified_filtered_all.vcf.gz"
params.sampleIDs = file("$params.samplelist").readLines().collect { it.split('\t')[0] }.toSet()
params.gtf = "${params.genomepath}/Bos_taurus.ARS-UCD1.2.108.chr.gtf"
params.fa = "${params.genomepath}/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa"
params.outpath = "${params.workpath}/result/"
params.toolpath = "${params.workpath}/tools/"
params.chr = "${params.genomepath}/chr.info"
params.chrinfo = file("${params.genomepath}/chromInfo.txt").readLines().collect { it.split('\t')[0] }.toSet()

  include { h5_prepare_1         } from './module/WASP.nf'
  include { h5_prepare_2         } from './module/WASP.nf'
  include { WASP_find            } from './module/WASP.nf'
  include { WASP_remap           } from './module/WASP.nf'
  include { WASP_dedup           } from './module/WASP.nf'
  include { leafcutter           } from './module/mp2.nf'
  include { phaser_snp           } from './module/mp2.nf'
  include { phaser_gene          } from './module/mp2.nf'
workflow{
  h5_prepare_1(params.chrinfo)
  h5_prepare_2(h5_prepare_1.out.vcf.collect())
  WASP_find(params.sampleIDs,h5_prepare_2.out.haplotype,h5_prepare_2.out.index,h5_prepare_2.out.tab)
  WASP_remap(WASP_find.out.fq,WASP_find.out.remapbam,WASP_find.out.keepbam)
  WASP_dedup(WASP_remap.out.bam)
  leafcutter(WASP_dedup.out.bam)
  phaser_snp(WASP_dedup.out.bam)
  phaser_gene(phaser_snp.out.count)    
}
