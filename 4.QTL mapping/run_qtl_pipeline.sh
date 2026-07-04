#!/bin/bash

# ====================================================================
# QTL Mapping 管道总控脚本 
# ====================================================================

# 1. 初始化变量
PROJECT_DIR=""
GTF_FILE=""
LOG_DIR="./qtl_pipeline_logs"

# 2. 解析命令行长参数
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --project_dir) PROJECT_DIR="$2"; shift ;;
        --gtf_file)    GTF_FILE="$2"; shift ;;
        *) echo -e "\033[31m[错误] 未知参数: $1\033[0m"; exit 1 ;;
    esac
    shift
done

# 3. 参数非空校验
if [ -z "$PROJECT_DIR" ] || [ -z "$GTF_FILE" ]; then
    echo -e "\033[33m[用法示例]\033[0m bash run_qtl_pipeline.sh --project_dir /path/to/project --gtf_file /path/to/gtf"
    echo -e "\033[31m[错误] 缺少必要参数！请同时指定 --project_dir 和 --gtf_file\033[0m"
    exit 1
fi

# 4. 创建日志目录
mkdir -p "$LOG_DIR"

# --------------------------------------------------------------------
# 封装运行与错误捕获函数 (双向输出：屏幕 + 文件)
# --------------------------------------------------------------------
run_step() {
    local step_name=$1
    local script_cmd=$2
    shift 2 # 移除前两个参数，剩余参数全部透传给下游脚本
    
    local log_file="${LOG_DIR}/${step_name}.log"

    echo -e "\n========================================================================"
    echo "[管道启动] 正在执行: ${step_name}"
    echo "--------------------------- 实时终端输出 ---------------------------"

    # 执行命令，2>&1 混合输出，tee 双写，PIPESTATUS[0] 捕获上游真实状态码
    eval "${script_cmd}" "$@" 2>&1 | tee "${log_file}"
    local exit_code=${PIPESTATUS[0]}

    echo "--------------------------------------------------------------------"
    if [ "$exit_code" -ne 0 ]; then
        echo -e "\033[31m[致命错误] 该步骤运行崩溃！ (退出错误码 Exit Code: ${exit_code})\033[0m"
        echo "管道流程已自动触发安全熔断终止。完整日志请查阅: ${log_file}"
        exit "$exit_code"
    else
        echo -e "\033[32m[成功完成] ${step_name} 执行完毕。 (Exit Code: 0)\033[0m"
    fi
}

# --------------------------------------------------------------------
# 管道按序嵌套调用执行区
# --------------------------------------------------------------------
echo ">>>>>>>>>>>>>>>>>>>> 开始运行 QTL Mapping 自动化管道 <<<<<<<<<<<<<<<<<<<<"
echo "启动时间: $(date)"
echo "项目大基座路径: $PROJECT_DIR"

# --- 步骤 1: beQTL 数据准备 (R 脚本带空格转义) ---
run_step "1.Prepare_beQTL" "Rscript 3.1\\ Prepare\\ for\\ beQTL.R" \
  --project_dir "$PROJECT_DIR" \
  --gtf_file "$GTF_FILE"

# --- 步骤 2: beQTL 关联映射 (调用 sh 脚本) ---
# 注意：原作者的 3.2 脚本内部写死了 main_dir，为了让其受外部控制，这里将变量直接通过环境变量 export 出去
export MAIN_DIR="${PROJECT_DIR}/CattleGTEx/OmiGA/eQTL/"
run_step "2.beQTL_Mapping" "bash 3.2\\ beQTL\\ mapping.sh"


# --- 步骤 3: cseQTL 数据准备 (R 脚本带空格转义) ---
run_step "3.Prepare_cseQTL" "Rscript 3.3\\ Prepare\\ for\\ cseQTL.R" \
  --project_dir "$PROJECT_DIR" \
  --gtf_file "$GTF_FILE"

# --- 步骤 4: cseQTL 关联映射 (调用 sh 脚本) ---
export MAIN_DIR_CSE="${PROJECT_DIR}/CattleGTEx/OmiGA/cell_specific/Group6/"
export GENO_DIR_CSE="${PROJECT_DIR}/CattleGTEx/OmiGA/eQTL/Group6/"
run_step "4.cseQTL_Mapping" "bash 3.4\\ cseQTL\\ mapping.sh"


# --- 步骤 5: ieQTL 数据准备 (R 脚本带空格转义) ---
run_step "5.Prepare_ieQTL" "Rscript 3.5\\ Prepare\\ for\\ cell-type\\ ieQTL.R" \
  --project_dir "$PROJECT_DIR"

# --- 步骤 6: ieQTL 细胞互作映射 (调用 sh 脚本) ---
export MAIN_DIR_IE="${PROJECT_DIR}/CattleGTEx/OmiGA/eQTL/"
run_step "6.ieQTL_Mapping" "bash 3.6\\ ieQTL\\ mapping.sh"


echo -e "\n>>>>>>>>>>>>>>>>>>>> 全套 QTL Mapping 管道安全跑通！ <<<<<<<<<<<<<<<<<<<<"
