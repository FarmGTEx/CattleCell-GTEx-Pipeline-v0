========================================================================
       Process 10: Selection Region Analysis 自动化管道使用说明书
========================================================================

一、 模块总体简介
本模块包含 3 个计算选择压力与功能基因组学富集重叠的 R 脚本以及 1 个 vcftools 
全基因组扫描脚本。
目前全套组件已完成基于外部传参的项目大基座（PROJECT_DIR）路径解耦，100% 保留了
原作者原本的多层嵌套大循环逻辑、富集公式算法、以及特定绘图环境配置，彻底去除了
硬编码绑定的固有死路径，并移除了 Slingshot 缓存盘误删数据的隐患。

二、 目录文件结构
进行自动化联动前，请确保同级目录下包含以下 5 个核心文件：
10.Selection region analysis/
├── Enrichment between GWAS and selection regions.R  # GWAS富集组件
├── Enrichment between eQTL and selection regions.R  # eQTL富集组件
├── Strongly selected regions.R                      # 前10%强选择提取组件
├── fst.sh                                           # 修正后的Slurm提交脚本
└── run_selection_pipeline.sh                        # 统一调用总控Bash脚本

三、 命令行运行指南（标准调用）

【步骤 1】赋予总控脚本可执行权限
在当前终端下键入命令：
chmod +x run_selection_pipeline.sh

【步骤 2】一键拉起全套流转
直接通过 bash 启动总控，并在 `--project_dir` 参数后追加你的服务器全局存储基座路径：
bash run_selection_pipeline.sh --project_dir "/faststorage/project/cattle_gtexs"

【步骤 3】运行时终端日志留底
总控启动后会按顺序自动拉起 3 个 R 组件，并在当前工作区内开辟 `./selection_pipeline_logs/` 
目录。多组织 SuSiE 映射文件检索信息、data.table 库高效重叠区间（foverlaps）比对
状态都会同步实时写出到以下日志中：
- GWAS富集日志：  ./selection_pipeline_logs/10.1.GWAS_Selection_Enrich.log
- eQTL富集日志：  ./selection_pipeline_logs/10.2.eQTL_Selection_Enrich.log
- 强选择筛选日志：./selection_pipeline_logs/10.3.Strong_Selected_Enrich.log

四、 独立组件手动调用命令行（单步调试）

如果需要单独挑出某个分析模块针对个别群体进行调试，可以直接使用 Rscript 并把项目
全局大基座路径作为第一个尾巴参数输入：

1. 独立运行 GWAS 富集：
Rscript Enrichment\ between\ GWAS\ and\ selection\ regions.R "/faststorage/project/cattle_gtexs"

2. 独立运行 eQTL 富集：
Rscript Enrichment\ between\ eQTL\ and\ selection\ regions.R "/faststorage/project/cattle_gtexs"

3. 独立运行 10% 强选择区间提取：
Rscript Strongly\ selected\ regions.R "/faststorage/project/cattle_gtexs"

五、 管道运行安全拦截说明
总控内置高级 PIPESTATUS 状态码熔断器。当 R 语言在迭代到某些特定性状、由于上游精细映射
输出（如 finemapping.txt）空缺或格式不符导致 OmiGA 读取终止时，总控会即刻进行安全卡死
并实施运行阻断，绝不允许脏数据向下层泛滥。
========================================================================
