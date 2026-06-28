单细胞 RNA-Seq 分析流水线 (Cellranger + Seurat)

本流水线专为处理单细胞测序数据而设计，集成了 Cellranger 进行比对与计数，并利用 Seurat 进行质量控制、双细胞鉴定、整合及细胞聚类分析。

目录

快速开始

环境准备 (Nextflow 25.04.2)

输入数据集构建 (Samplesheet)

容器运行环境说明

配置文件与资源管理

常见参数说明

1. 快速开始

在准备好环境后(java17)，使用以下命令运行流水线：

./nextflow run main.nf \
    -profile docker \
    --samplesheet /target_dir/samplesheet.csv \
    --ref_fa "target_species_reference_genome/reference_genome.fna" \
    --ref_gtf "target_species_gtf/genomic.gtf" \
    --species_name "species_name" \
    --mito_threshold 10 \
    --pcs_use 30 \
    --outdir ../Single_cell_processing/ \
    -resume


2. 环境准备 (Nextflow 25.04.2)

本流水线推荐使用 Nextflow 25.04.2 版本。请按照以下步骤准备运行环境：

2.1 Java 环境要求

Nextflow 依赖 Java 运行，建议安装 Java 17 或以上版本：

# 检查 Java 版本
java -version

# 若未安装，建议通过 SDKMAN! 安装（推荐）
curl -s "[https://get.sdkman.io](https://get.sdkman.io)" | bash
source "$HOME/.sdkman/bin/sdkman-init.sh"
sdk install java 17.0.x-tem


2.2 安装 Nextflow 25.04.2

直接下载对应版本的可执行文件并设置路径：

# 下载 Nextflow 25.04.2
wget -O nextflow [https://github.com/nextflow-io/nextflow/releases/download/v25.04.2/nextflow-25.04.2-all](https://github.com/nextflow-io/nextflow/releases/download/v25.04.2/nextflow-25.04.2-all)
chmod +x nextflow
sudo mv nextflow /usr/local/bin/

# 验证安装
nextflow -version


2.3 Docker 环境要求

本流水线全程使用 Docker 容器，请确保 Docker 服务已启动且用户拥有权限：

# 确保 Docker 服务运行
sudo systemctl start docker

# 检查权限（确保无需 sudo 即可执行 docker 命令）
docker ps


3. 输入数据集构建 (Samplesheet)

流水线需要一个 CSV 格式的 samplesheet.csv 文件来定位数据。

3.1 文件格式

文件必须包含 sampleid 和 fastq_dir 两列：

sampleid,fastq_dir
sample1,/media/Single_cell_processing/srr_test_project/data/fastqs
sample2,/media/Single_cell_processing/srr_test_project/data/fastqs
sample3,/media/Single_cell_processing/srr_test_project/data/fastqs


3.2 样本名称 (sampleid) 命名规范建议

唯一性: 每个样本必须具有唯一的标识符。

字符限制: 仅使用字母、数字、下划线 _ 或连字符 -。

逻辑命名法: 建议采用 组别_处理方式_重复编号 的格式（例如：Control_WT_Rep1）。

4. 容器运行环境说明

本流水线采用双镜像隔离策略：

4.1 scrna_pipeline_env:latest (默认容器)

职责: 执行 Cellranger 比对步骤。

环境: 封装 Cellranger 官方环境及必要的基础 Linux 系统依赖。

4.2 local_seurat_env:latest (R 分析容器)

职责: 执行 SEURAT_PER_SAMPLE_QC 和 SEURAT_INTEGRATION 步骤。

环境: 封装 R 4.x 环境及所有关键分析包（Seurat v4, Harmony, DoubletFinder 等）。

5. 配置文件与资源管理

5.1 内存与 CPU 资源设置 (conf/base.config)

如需调整单个进程的内存占用，请修改 conf/base.config 中的 withLabel 配额。

5.2 默认路径配置 (nextflow.config)

您可以在 nextflow.config 的 params 区域设置项目范围的默认路径和参数：

params {
    outdir         = "$baseDir/scRNA_results"
    mito_threshold = 10
    pcs_use        = 30
}


6. 常见参数说明

参数

描述

--samplesheet

样本 CSV 文件路径

--ref_fa

参考基因组 FASTA 文件路径

--ref_gtf

参考基因组 GTF 注释文件路径

--species_name

构建索引时的物种名称

--mito_threshold

线粒体基因过滤阈值（默认 10%）

--pcs_use

降维分析使用的 PC 数量（默认 30）

--outdir

结果输出根目录

技术支持

如果运行中遇到图表生成错误（白纸问题），请检查 results/03_Seurat_Integration 目录下是否有 RDS 数据。本流水线已配置自动恢复机制，即使绘图失败也会完好保留核心数据。
