========================================================================
         Process 7: eQTL Colocalization 自动化管道使用说明书
========================================================================

一、 模块总体简介
本模块包含 2 个核心的贝叶斯共定位分析组件：
1. 7.1 Prepare for coloc.R (共定位输入表型矩阵提取与清洗)
2. 7.2 Coloc.R             (基于 coloc.abf 算法的贝叶斯共定位核心检验)

本模块的核心任务是顺序串联四个层级的 QTL 信号（常规大样 Bulk eQTL、细胞组分 ieQTL、
细胞特异性 cseQTL、细胞状态 ieQTL），将其与精细映射窗口进行重叠解构，最后喂给 coloc 
算法计算后验概率（PP0-PP4）。
目前已使用 R 原生传参机制完成了全局项目大基座（PROJECT_DIR）路径解耦封装，100% 
保持了原作者原本的嵌套大循环、过滤阈值与统计算法内核。

二、 目录文件结构
在运行前，请确保当前文件夹（7.eQTL colocalization）下包含以下三个核心文件：
7.eQTL colocalization/
├── 7.1 Prepare for coloc.R   # 封装后的共定位数据准备脚本
├── 7.2 Coloc.R               # 封装后的 coloc 核心检验脚本
└── run_coloc_pipeline.sh     # 统一调用总控 Bash 脚本

三、 命令行运行指南（标准一键流转）

【步骤 1】赋予总控脚本可执行权限
在当前目录下执行命令：
chmod +x run_coloc_pipeline.sh

【步骤 2】一键启动管道命令
使用 bash 启动总控脚本，并通过双横线长参数 `--project_dir` 传入你的项目全局大基座路径：
bash run_coloc_pipeline.sh --project_dir "/faststorage/project/cattle_gtexs"

【步骤 3】实时查看运行时日志与排质控
管道启动后，会自动在当前目录下创建 `./coloc_pipeline_logs/` 文件夹。
多层级 QTL 的 1:29 染色体 Nominal pairs 矩阵合并状态、以及 coloc 软件包在拟合
ABF 概率时的迭代日志都会实时在终端滚动，并同步双写输出到以下文件中：
- 数据准备日志：  ./coloc_pipeline_logs/7.1.Prepare_Coloc.log
- 共共定位核心日志：./coloc_pipeline_logs/7.2.Run_Coloc_ABF.log

四、 独立组件手动调用命令行（供单步调试使用）

如果需要单独挑出某一个脚本针对特定的组织或异常块进行调试，可以直接使用 Rscript 并把项目
全局路径作为第一个位置参数传入脚本尾部：

1. 独立运行 coloc 数据准备：
Rscript 7.1\ Prepare\ for\ coloc.R "/faststorage/project/cattle_gtexs"

2. 独立运行 coloc 贝叶斯检验：
Rscript 7.2\ Coloc.R "/faststorage/project/cattle_gtexs"

五、 各脚本封装改动细节与嵌合备忘

1. 7.1 Prepare for coloc.R
   - 【改动说明】：在脚本最顶部引入了 args <- commandArgs(trailingOnly = TRUE) 机制，
     抓取外部传入的 PROJECT_DIR 变量，用以接管原本全部写死的 eQTL、cell_specific、
     Cell_state 绝对路径。在每一块 fwrite 输出前均加入了 dir.create() 动态目录保护。
   - 【原样说明】：原作者在第一块 Bulk QTL 里写死的 "for (i in 7:length(tissue))" 偏置
     （略过前6个组织），以及在第四块 Cell state 里写死的 "for (i in 14:16)" 限制，
     全部原汁原味地保留在代码原位。

2. 7.2 Coloc.R
   - 【改动说明】：同样引入原生 args 参数抓取，用 PROJECT_DIR 动态外延寻址上游 7.1 
     生成的 /Coloc/coloc/Bulk/ 和 /all_cell_components/ 目录下的中间产物 .bed 文件。
   - 【原样说明】：原作者在脚本大循环开头写死的特定循环起点 "for (i in 3:length(tissue))"
     （漏掉前2个组织），以及底部的 merge 杂交规则、coloc.abf() 内部的 quant 定量性状参数
     一字未动，100% 确保产出的 coloc.csv 结果与之前旧版结果无缝一致。

六、 运行安全与阻断说明
本管道总控脚本内嵌高级 PIPESTATUS[0] 错误码拦截。如果因为上游某些特定的 cseQTL 缺失
对应的 *_coloc.bed 物理输入，或者 R 语言环境在 merge 时触发异常，总控会即刻实施安全熔断，
彻底卡死后续流转，并在终端直接报出非零错误码，绝不允许脏数据蔓延。
========================================================================
