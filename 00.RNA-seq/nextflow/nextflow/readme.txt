
###############
tools downloading(download.sh)
###############

To run this pipeline, you shall download the tools and some annotation files first. You can simply run the download.sh. Before running this, please make sure that you already specified the path in the download.sh file.

the gtf file we use is Bos_taurus.ARS-UCD1.2.108.chr.gtf(from Ensembl).

We also provide a nf-farmgtex.yml file so you can quickly set up the enviroment.

##############
prepare(prepare.nf)
##############
Then you need to prepare some file such as index for STAR, you can run prepare.nf to do so.And the path shall be specified as the download step. It is unnecessary to make all the folders manually because nextflow will creat them automatically.In this case, you will only need to specify the working directory(workpath).

To do the TE analysis, it will take long to run a repeatmasker, so I provide the file we gonna need named as "Bos_taurusARS-UCD.repeats.fa".

Also, the enhancer annotation file(enhancer.saf) is in file folder.

And before you run prepare.nf, please make sure the fa file and the bed file in file folder are moved to genomepath.


##############
nextflow part1(rna.nf)
##############
Now you can try analysis, you need to prepare a table(samplelist.csv) in your workpath with sample names(with no suffix) in the first column and tissue type in the second column. 

the format of samplelist.csv should be as below
#########
SAMEA4447762	Mammary Gland
SAMEA4447830	Mammary Gland
SAMEA5847541	Liver
SAMEA5847543	Liver
#########

This time you need to specify both the work path and the input path in the rna.nf file, and also the suffix should be specified(.fastq.gz/.fq.gz/....).
For example if your path is /path/to/your/file/${sampleID}/*.fastq.gz, your expected {params.input} should be "/path/to/your/file/", the suffix should be ".fastq.gz". If your data do not have a extra "sampleID" folder, you need to modify the scripts in the module a little, we provide a bash file in module folder to solve this.

Also you need to modify the "clusteroption" in config file called nextflow.config if you need to submmit jobs to server. If you want to change the memory and resouce distributed to every single task, simply edit the scripts in module folder.

Now you can run the step1.sh to start the analysis.
##26.01.2024 We remove the TE analysis and RNA editing in the process so it will be more efficient, but the module is still there.

##############
VCF calling
##############
After finish the first part, we will use the bqsr.bam to call snps and do imputation. You need the reference panel of cattleGTEx to do this.
The code is provided in the vcfcalling folder.
 
##############
nextflow part2(rna2.nf)
##############
After you finish vcf calling, you can use the filtered vcf to remove the mapping bias and finish some other phenotypes using the rna2.nf.
Similar with the previous one,except you shall specify the path to the imputed vcf file in {params.recalvcf}.
This time you can run the nextflow2.sh

#############
result sharing
#############

All folders in your output path should be shared with cattle gtex project except the folder "alignment"， and the folder depar2 should be shared to us for the tissue annotation. 


#############
Farmgtex 2025 huicongz@qgg.au.dk
#############


