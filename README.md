# Analysis pipelines for the pilot phase (V0) of CattleCell-GTEx
### **1. Introduction**
Understanding the genetic and molecular architecture of complex traits and adaptative evolution is crucial for advancing sustainable precision breeding in cattle and other livestock. Yet, how genetic variation affects cellular gene expression remains elusive in cattle. Here, by integrating 8,866 bulk RNA-seq samples and 999,192 single cells of 81 cell types in 22 bovine tissues, we presented a comprehensive atlas of regulatory variants at the cell type resolution in cattle. Compared to standard bulk-tissue expression quantitative trait loci (beQTL), we detected 57,043 novel cell-type eQTL in 18,153 genes, which exhibited a stronger tissue/cell-type specificity. By examining genome-wide associations (GWAS) of 44 complex traits, these cell-type eQTL were colocalized with 505 (24%) additional GWAS loci compared to beQTL. Through integrating this resource with selection signals between dairy and beef cattle as well as among 154 ancient DNA samples, we provided tissue/cell-specific regulatory insights into cattle breeding and domestication. Overall, the current atlas of cell-type-specific regulatory variants will serve as an invaluable resource for cattle genomics, selective breeding, and domestication.

![CattleCellGTEx](https://github.com/FarmGTEx/CattleCell-GTEx-Pipeline-v0/blob/main/CattleCellGTEx.png)

### **2. Analysis pipeline**
### This repository contains analysis pipelines used by the CattleCell-GTEx Consortium, including:

:black_circle:Raw single-cell RNA-seq data processing, quality control and cell-type annotation

:black_circle:Deconvolution for cell components, gene expression and cell states

:black_circle:cis-Heritability estimation in cell-type gene expression

:black_circle:eQTL mapping

:black_circle:Functional enrichment

:black_circle:Tissue and cell type sharing pattern

:black_circle:Colocalation between eQTL

:black_circle:Analysis between eQTL and GWAS

:black_circle:GWAS and eQTL enrichment

:black_circle:Enrichment between selection signal and eQTL


### **3. Citation**

### An atlas of cell type specific regulatory effects in cattle (Preprint)
Houcheng Li1†, Huicong Zhang1†, Pengju Zhao3†, Qi Zhang1,2†, Senlin Zhu4†, Tao Shi5†, Bo Han2, Weijie Zheng1,2, Liu Yang6,7, Victoria Elizabeth Mullin8, Jolijn Erven9, Jicai Jiang10; Li Ma7, Mian Gong1,11, Xiaoning Zhu1,12, Qing Lin1,13, Yang Xi1,14, Di Zhu1,12, Jinyan Teng13, Dailu Guan15, Yali Hou11, Fei Wang16, John F. O’Grady17, David E. Machugh17,18,19, Bingjie Li20, Laurent Frantz21,22, Greger Larson23, Zexi Cai1, Goutam Sahana1, Daniel Bradley8*, Yu Jiang5*, Huizeng Sun4*, Dongxiao Sun2*, Geroge E. Liu6*, Lingzhao Fang1*
