#!/bin/bash

# ====================================================================
# eQTL Colocalization (coloc) 自动化分析总控脚本
# ====================================================================

PROJECT_DIR=""
LOG_DIR="./coloc_pipeline_logs"

# 解析双横线长参数
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --project_dir) PROJECT_DIR="$2"; shift ;;
        *) echo -e "\033[31m[错误] 未知参数: $1\033[0m"; exit 1 ;;
    esac
    shift
done

if [ -z "$PROJECT_DIR" ]; then
    echo -e "\033[33m[用法]\033[0m bash run_coloc_pipeline.sh --project_dir /path/to/project"
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
        echo -e "\033[31m[致命错误] 该步骤运行崩溃！(Exit Code: ${exit_code})\033[0m"
        exit "$exit_code"
    else
        echo -e "\033[32m[成功] ${step_name} 执行完毕。\033[0m"
    fi
}

echo ">>>>>>>>>>>>>>>>>>>> 开始运行 eQTL Colocalization 管道 <<<<<<<<<<<<<<<<<<<<"

# Step 1: 准备共定位表型 BED 矩阵输入 (将项目根目录作为参数直接悬挂透传)
run_step "7.1.Prepare_Coloc" "Rscript 7.1\\ Prepare\\ for\\ coloc.R ${PROJECT_DIR}"

# Step 2: 启动贝叶斯共定位 ABF 核心检验
run_step "7.2.Run_Coloc_ABF" "Rscript 7.2\\ Coloc.R ${PROJECT_DIR}"

echo -e "\n>>>>>>>>>>>>>>>>>>>> 全套 Colocalization 管道安全跑通！ <<<<<<<<<<<<<<<<<<<<"
