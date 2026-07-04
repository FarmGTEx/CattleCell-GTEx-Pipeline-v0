========================================================================
       Process 9: GWAS & Cell-type Enrichment 自动化管道使用说明书
========================================================================

一、 模块总体简介
本模块包含两个核心的单细胞/空间转录组水平 GWAS 细胞类型富集分析组件：
1. gsmap.R   (基于空间基因图谱映射算法)
2. scpagwas.R (基于单细胞通路性状回归算法)

目前两个组件已全部完成外部参数路径解耦封装。在 100% 保留原作者全部底层业务逻辑、
算法、拼写、甚至手误夹带代码的前提下，彻底抽离了原本写死的物理绝对路径。
现在通过统一的总控 Bash 脚本（run_enrichment_pipeline.sh），可以一站式全自动顺序
调度这两个组件，并支持终端日志实时双写留底。

二、 目录文件结构
在运行前，请确保当前文件夹（9.GWAS and Cell-type enrichment）下包含以下三个脚本：
9.GWAS and Cell-type enrichment/
├── gsmap.R                   # 封装后的 gsMap 富集脚本
├── scpagwas.R                # 封装后的 scPagwas 富集脚本
└── run_enrichment_pipeline.sh  # 一键调用总控 Bash 脚本

三、 命令行运行指南（标准调用）

【步骤 1】赋予总控脚本可执行权限
在终端中执行以下命令：
chmod +x run_enrichment_pipeline.sh

【步骤 2】一键启动管道命令
使用 bash 启动总控脚本，并通过双横线长参数 `--project_dir` 传入你的项目全局大基座路径：
bash run_enrichment_pipeline.sh --project_dir "/faststorage/project/cattle_gtexs"

【步骤 3】实时查看运行时日志
管道启动后，会自动在当前目录下创建 `./enrichment_pipeline_logs/` 文件夹。
所有的 R 语言矩阵运算、Seurat 标准化进度和报错流都会实时在终端屏幕上滚动，并同步写出
到以下日志文件中，方便无人值守或后期追溯：
- gsMap 运行日志：  ./enrichment_pipeline_logs/9.1.Run_gsMap.log
- scPagwas 运行日志：./enrichment_pipeline_logs/9.2.Run_scPagwas.log

四、 独立组件手动调用命令行（供单步调试使用）

如果需要单独调试或单独运行某一个算法，可以直接使用 Rscript 并把项目大基座路径作为
第一个位置参数传入脚本尾部：

1. 独立单步运行 gsMap 命令行：
Rscript gsmap.R "/faststorage/project/cattle_gtexs"

2. 独立单步运行 scPagwas 命令行：
Rscript scpagwas.R "/faststorage/project/cattle_gtexs"

五、 各脚本封装改动细节备忘

1. gsmap.R
   - 【改动说明】：在脚本最顶部引入了 R 原生的 args <- commandArgs(trailingOnly = TRUE) 
     传参机制，抓取外部传入的 PROJECT_DIR 变量。
   - 【逻辑保留】：100% 维持了原本所有的 Seurat 提取、GetAssayData 槽位、human_ref 
     分组切片、1:29 染色体 Hols_ref.chr 扫描等全部业务。
   - 【原样说明】：原作者在第 29-31 行手误混入的纯 Python 绘图代码（import matplotlib 等）
     以及大循环里手误漏写一个 s 的 length(trait) 变量，全部原汁原味地保留在代码原位，
     仅对硬编码的物理路径做了变量平替。

2. scpagwas.R
   - 【改动说明】：同样在脚本最顶部引入了原生 args 参数抓取，用 PROJECT_DIR 变量彻底
     接管了原本写死的单细胞 rds 路径、ARS-UCD1.2 gtf 注释路径、ld_chrom 连锁平衡文件路径、
     以及 05sperm 下的 mlma 关联路径。
   - 【逻辑保留】：原作者原本多嵌套大循环的回归、200次 Bootstrap 抽样检验、以及底部的 
     Link_pathway_blocks_gwas 核心富集回归计算一字未改，保证算法和以前完全一致。

六、 运行安全拦截说明
本管道总控脚本内嵌 PIPESTATUS[0] 状态码拦截机制。如果在无人值守运行期间，R 语言环境
因为缺少某些依赖包（如 scPagwas、anndata 等）抛出异常，管道会在一秒钟内实施安全熔断，
立刻卡死后续的流转并报出非零错误码，绝不容许产生脏数据或空运行。
========================================================================
