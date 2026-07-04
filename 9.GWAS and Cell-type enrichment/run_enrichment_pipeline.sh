#!/bin/bash

# ====================================================================
# GWAS & Cell-type Enrichment 一站式流转管道
# ====================================================================

# 1. 初始变量
PROJECT_DIR=""
LOG_DIR="./enrichment_pipeline_logs"

# 2. 解析双横线命令行参数
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --project_dir) PROJECT_DIR="$2"; shift ;;
        *) echo -e "\033[31m[错误] 未知参数: $1\033[0m"; exit 1 ;;
    esac
    shift
done

if [ -z "$PROJECT_DIR" ]; then
    echo -e "\033[33m[用法]\033[0m bash run_enrichment_pipeline.sh --project_dir /path/to/project"
    exit 1
fi

mkdir -p "$LOG_DIR"

# --------------------------------------------------------------------
# 实时双向流日志打印与拦截器
# --------------------------------------------------------------------
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
        echo -e "\033[31m[错误拦截] 该步骤返回了非零状态码 (Exit Code: ${exit_code})\033[0m"
        exit "$exit_code"
    else
        echo -e "\033[32m[成功完成] ${step_name} 执行完毕。\033[0m"
    fi
}

echo ">>>>>>>>>>>>>>>>>>>> 开始运行 GWAS & 细胞类型富集计算管道 <<<<<<<<<<<<<<<<<<<<"
echo "全局参数基座: $PROJECT_DIR"

# Step 1: 运行 gsMap
run_step "9.1.Run_gsMap" "Rscript gsmap.R ${PROJECT_DIR}"

# Step 2: 运行 scPagwas
run_step "9.2.Run_scPagwas" "Rscript scpagwas.R ${PROJECT_DIR}"

echo -e "\n>>>>>>>>>>>>>>>>>>>> 全套 Enrichment 细胞富集管道运行结束！ <<<<<<<<<<<<<<<<<<<<"
