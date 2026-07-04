#!/bin/bash

# ====================================================================
# Selection Region Analysis 自动化流转管道
# ====================================================================

PROJECT_DIR=""
LOG_DIR="./selection_pipeline_logs"

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --project_dir) PROJECT_DIR="$2"; shift ;;
        *) echo -e "\033[31m[错误] 未知参数: $1\033[0m"; exit 1 ;;
    esac
    shift
done

if [ -z "$PROJECT_DIR" ]; then
    echo -e "\033[33m[用法]\033[0m bash run_selection_pipeline.sh --project_dir /path/to/project"
    exit 1
fi

mkdir -p "$LOG_DIR"

run_step() {
    local step_name=$1
    local script_cmd=$2
    local log_file="${LOG_DIR}/${step_name}.log"

    echo -e "\n========================================================================"
    echo "[管道启动] 正在执行: ${step_name}"
    echo "--------------------------- 实时终端输出 ---------------------------"

    eval "${script_cmd}" 2>&1 | tee "${log_file}"
    local exit_code=${PIPESTATUS[0]}

    echo "--------------------------------------------------------------------"
    if [ "$exit_code" -ne 0 ]; then
        echo -e "\033[31m[阻断异常] 组件返回非零状态码 (Exit Code: ${exit_code})\033[0m"
        exit "$exit_code"
    else
        echo -e "\033[32m[运行成功] ${step_name} 执行完毕。\033[0m"
    fi
}

echo ">>>>>>>>>>>>>>>>>>>> 开始运行选择压力与多维信号富集检测管道 <<<<<<<<<<<<<<<<<<<<"
echo "项目基座路径: $PROJECT_DIR"

# Step 1: 运行 GWAS 信号与 Fst 选择压重叠富集分析
run_step "10.1.GWAS_Selection_Enrich" "Rscript Enrichment\\ between\\ GWAS\\ and\\ selection\\ regions.R ${PROJECT_DIR}"

# Step 2: 运行 eQTL 精细映射特异信号与 Fst 选择压重叠富集分析
run_step "10.2.eQTL_Selection_Enrich" "Rscript Enrichment\\ between\\ eQTL\\ and\\ selection\\ regions.R ${PROJECT_DIR}"

# Step 3: 提取前 10% 强选择核心区域重叠矩阵
run_step "10.3.Strong_Selected_Enrich" "Rscript Strongly\\ selected\\ regions.R ${PROJECT_DIR}"

echo -e "\n>>>>>>>>>>>>>>>>>>>> 全套 Selection 管道流转安全结束！ <<<<<<<<<<<<<<<<<<<<"
