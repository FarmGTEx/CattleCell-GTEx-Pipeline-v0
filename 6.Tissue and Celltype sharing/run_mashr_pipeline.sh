#!/bin/bash

# ====================================================================
# Tissue & Celltype Sharing (MashR) 自动化管道总控脚本
# ====================================================================

# 1. 初始变量
PROJECT_DIR=""
MODE=""
LOG_DIR="./mashr_pipeline_logs"
export SCRIPT_DIR="$(pwd)" # 当前存放 Python 和 R 脚本的目录

# 2. 解析命令行长参数
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --project_dir) PROJECT_DIR="$2"; shift ;;
        --mode)        MODE="$2"; shift ;;
        *) echo -e "\033[31m[错误] 未知参数: $1\033[0m"; exit 1 ;;
    esac
    shift
done

# 3. 校验参数
if [ -z "$PROJECT_DIR" ] || [ -z "$MODE" ]; then
    echo -e "\033[33m[用法]\033[0m bash run_mashr_pipeline.sh --project_dir <路径> --mode <bulk|celltype>"
    exit 1
fi

if [[ "$MODE" != "bulk" && "$MODE" != "celltype" ]]; then
    echo -e "\033[31m[错误] --mode 参数只能是 'bulk' 或 'celltype' \033[0m"
    exit 1
fi

mkdir -p "$LOG_DIR"

# --------------------------------------------------------------------
# 函数：带实时日志与错误码拦截的脚本执行器
# --------------------------------------------------------------------
run_step() {
    local step_name=$1
    local script_cmd=$2
    local log_file="${LOG_DIR}/${step_name}_${MODE}.log"

    echo -e "\n========================================================================"
    echo "[管道启动] 正在执行: ${step_name} (${MODE} 模式)"
    echo "--------------------------- 实时终端输出 ---------------------------"

    eval "${script_cmd}" 2>&1 | tee "${log_file}"
    local exit_code=${PIPESTATUS[0]}

    echo "--------------------------------------------------------------------"
    if [ "$exit_code" -ne 0 ]; then
        echo -e "\033[31m[致命错误] 该步骤运行崩溃！(Exit Code: ${exit_code})\033[0m"
        echo "流程已安全熔断，完整日志请查阅: ${log_file}"
        exit "$exit_code"
    else
        echo -e "\033[32m[成功] ${step_name} 执行完毕。\033[0m"
    fi
}

# --------------------------------------------------------------------
# 管道逻辑与环境变量分发
# --------------------------------------------------------------------
echo ">>>>>>>>>>>>>>>>>>>> 开始运行 MashR Sharing 管道 <<<<<<<<<<<<<<<<<<<<"

if [ "$MODE" == "bulk" ]; then
    # Bulk 模式：发布相应的环境基座路径
    export NOMINAL_DIR="${PROJECT_DIR}/CattleGTEx/OmiGA/eQTL/"
    export RANDOM_OUT="${PROJECT_DIR}/CattleGTEx/OmiGA/MASHR_Celltype/Tissue_sharing/Bulk/Random_output"
    export STRONG_OUT="${PROJECT_DIR}/CattleGTEx/OmiGA/MASHR_Celltype/Tissue_sharing/Bulk/strong_output"
    
    mkdir -p "$RANDOM_OUT" "$STRONG_OUT"

    export OUTPUT_DIR="$RANDOM_OUT"
    export SUBSET_SIZE=1000000
    run_step "1.Bulk_Prepare_Random" "bash Bulk-sharing-Prepare-MashR-Random.sh"

    export OUTPUT_DIR="$STRONG_OUT"
    run_step "2.Bulk_Prepare_Strong" "bash Bulk-sharing-Prepare-MashR-Strong.sh"

    # 设置传入 run_MashR.sh 的变量
    export STRONG_FILE="${STRONG_OUT}/strong_pairs.MashR_input.txt.gz"
    export RANDOM_FILE="${RANDOM_OUT}/nominal_pairs.${SUBSET_SIZE}_subset.MashR_input.txt.gz"
    export MASHR_OUT_DIR="${STRONG_OUT}/output_top_pairs"
    
    run_step "3.Run_MashR_Bulk" "bash run_MashR.sh"

elif [ "$MODE" == "celltype" ]; then
    # Celltype 模式：发布相应的环境基座路径
    export NOMINAL_DIR="${PROJECT_DIR}/CattleGTEx/OmiGA/cell_specific/"
    export RANDOM_OUT="${PROJECT_DIR}/CattleGTEx/OmiGA/MASHR_Celltype/Celltype_sharing/random_output"
    export STRONG_OUT="${PROJECT_DIR}/CattleGTEx/OmiGA/MASHR_Celltype/Celltype_sharing/strong_output"
    
    mkdir -p "$RANDOM_OUT" "$STRONG_OUT"

    export OUTPUT_DIR="$RANDOM_OUT"
    export SUBSET_SIZE=10000  # 修复原脚本中的缺失变量
    run_step "1.Celltype_Prepare_Random" "bash Celltype-sharing-Prepare-MashR-Random.sh"

    export OUTPUT_DIR="$STRONG_OUT"
    run_step "2.Celltype_Prepare_Strong" "bash Celltype-sharing-Prepare-MashR-Strong.sh"

    # 设置传入 run_MashR.sh 的变量
    export STRONG_FILE="${STRONG_OUT}/strong_pairs.MashR_input.txt.gz"
    export RANDOM_FILE="${RANDOM_OUT}/nominal_pairs.${SUBSET_SIZE}_subset.MashR_input.txt.gz"
    export MASHR_OUT_DIR="${STRONG_OUT}/output_top_pairs"
    
    run_step "3.Run_MashR_Celltype" "bash run_MashR.sh"
fi

echo -e "\n>>>>>>>>>>>>>>>>>>>> 全套 MashR 管道运行完毕！ <<<<<<<<<<<<<<<<<<<<"
