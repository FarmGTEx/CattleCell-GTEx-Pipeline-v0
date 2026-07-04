========================================================================
         Process 8: GWAS & eQTL Colocalization 管道使用说明书
========================================================================

一、 模块总体简介
本模块通过串联 GCTA-COJO 条件分析过滤与 COLOC 贝叶斯统计模型，实现复合物状 
GWAS 精细映射核心信号点与大样组织常规 eQTL 映射信号的共定位分析。
全套流程使用项目大基座（PROJECT_DIR）参数完成纯路径解耦封装

二、 目录文件结构
8.GWAS and eQTL/
├── 1.Prepare_eQTL.R              # 提取 eQTL 共定位 BED 文件
├── 2.Prepare_GWAS.R              # 提取 GWAS ±1Mb 候选窗口
├── 3.run_coloc.R                 # 贝叶斯共共定位核心检验
├── GCTA-COJO.R                   # MLMA 转 COJO 格式处理 R 脚本
├── run_gwas_coloc_pipeline.sh    # 统一调用总控 Bash 脚本
└── README.txt                    # 本专项审计说明书

三、 命令行运行指南（标准一键流转）

【步骤 1】赋予总控脚本可执行权限
在当前目录下执行命令：
chmod +x run_gwas_coloc_pipeline.sh

【步骤 2】一键启动管道命令
使用 bash 启动总控脚本，并通过双横线长参数 `--project_dir` 传入你的项目全局大基座路径：
bash run_gwas_coloc_pipeline.sh --project_dir "/faststorage/project/cattle_gtexs"

【步骤 3】实时查看运行时日志与排质控
管道启动后，会自动在当前目录下创建 `./gwas_coloc_logs/` 文件夹。
GCTA 条件信号扫描、QTL 的 1:29 染色体合并状态都会同步双写输出到以下日志中：
- COJO 转换日志：   ./gwas_coloc_logs/8.1.GCTA_Cojo_Format.log
- eQTL 准备日志：   ./gwas_coloc_logs/8.2.Prepare_eQTL.log
- GWAS 准备日志：   ./gwas_coloc_logs/8.3.Prepare_GWAS.log
- 贝叶斯共定位日志：./gwas_coloc_logs/8.4.Run_Coloc.log

四、 独立组件手动单步调用命令行（供单步调试使用）

如果需要单独挑出某个分析模块针对个别表型群体或异常块进行调试，可以直接独立运行各组件：

1. 单步准备 eQTL 共定位块：
Rscript 1.Prepare_eQTL.R --input_dir "/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/eQTL" --output_dir "./output_eQTL"

2. 单步准备 GWAS 候选窗口：
Rscript 2.Prepare_GWAS.R --gwas_dir "/faststorage/project/cattle_gtexs/Downstream_analysis/complextraits/GWAS" --cojo_dir "/faststorage/project/cattle_gtexs/Downstream_analysis/complextraits/GWAS" --output_dir "./output_GWAS"

五、 管道流转与数据闭环嵌合点备忘

总控脚本通过强约束输入输出实现了完全无缝的业务嵌合，在向甲方汇报时可重点强调以下数据链流转：
1. 【格式与群体扫描嵌合】：
   `GCTA-COJO.R` 执行完毕后，在各性状目录下生成的 `*_cojo.txt` 会作为输入无缝喂给外壳 Bash 循环里的 `gcta64 --cojo-file` 条件分析命令，进而扫描并产出独立显著关联信号。
2. 【表型与窗口重叠嵌合】：
   `1.Prepare_eQTL.R` 吐出的各组织 `*_coloc.bed` 文件，与 `2.Prepare_GWAS.R` 提取出的 GWAS 强相关区域候选信息，在业务上是并列的，它们一同作为下游 `3.run_coloc.R` 的并列输入件（即 `--bulk_dir` 与 `--gwas_dir`），通过 `merge(..., by = "variant_id")` 完成了物理位置上的精准咬合。
3. 【下游 Process 10 的基石】：
   本模块最终在结果目录下吐出的各个性状共定位文件，将直接充当后续 Process 10（选择压力分析中 GWAS 富集脚本）的核心输入。

六、 运行安全与阻断说明
本管道总控脚本内嵌高级 PIPESTATUS[0] 错误码拦截。若前序的格式转换、eQTL 染色体纵向合并或 GWAS 区域外延提取中任何一个 Rscript 触发不可逆异常，总控会即刻实施安全熔断，彻底卡死后续流转，绝不允许脏数据或空文件蔓延至下游贝叶斯检验中。
========================================================================
