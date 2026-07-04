#!/bin/bash

# ====================================================================
# 去卷积分析管道总控脚本 
# ====================================================================

# 1. 初始化变量
INPUT_PSEUDO=""
INPUT_REF=""
GTF_FILE=""
OUTPUT_DIR=""
PROJECT_DIR=""
TISSUE_NAME=""
LOG_DIR="./pipeline_logs"

# 2. 解析命令行长参数
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input_pseudo) INPUT_PSEUDO="$2"; shift ;;
        --input_ref)    INPUT_REF="$2"; shift ;;
        --gtf_file)     GTF_FILE="$2"; shift ;;
        --output_dir)   OUTPUT_DIR="$2"; shift ;;
        --project_dir)  PROJECT_DIR="$2"; shift ;;
        --tissue)       TISSUE_NAME="$2"; shift ;;
        *) echo -e "\033[31m[错误] 未知参数: $1\033[0m"; exit 1 ;;
    esac
    shift
done

# 3. 核心参数非空校验
if [ -z "$INPUT_PSEUDO" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$TISSUE_NAME" ]; then
    echo -e "\033[33m[用法示例]\033[0m bash run_deconvolution.sh --input_pseudo <rds> --input_ref <rds> --gtf_file <gtf> --output_dir <dir> --project_dir <dir> --tissue <name>"
    echo -e "\033[31m[错误] 缺少必要参数！请至少提供 --input_pseudo, --output_dir 和 --tissue\033[0m"
    exit 1
fi

# 4. 创建日志与输出目录
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# --------------------------------------------------------------------
# 封装运行与错误捕获函数
# --------------------------------------------------------------------
run_script() {
    local step_name=$1
    local script_cmd=$2 # 接收外面传进来的整条 Rscript 物理硬转义命令
    local log_file="${LOG_DIR}/${step_name}.log"

    echo -e "\n================================================================"
    echo "[启动] 正在执行: ${step_name}"
    echo "------------------------ 实时终端输出 ------------------------"

    # 使用 eval 直接执行外面物理硬转义好的完整命令
    eval "${script_cmd}" 2>&1 | tee "${log_file}"
    local exit_code=${PIPESTATUS[0]}

    echo "--------------------------------------------------------------"
    if [ "$exit_code" -ne 0 ]; then
        echo -e "\033[31m[致命错误] 步骤 ${step_name} 运行崩溃！ (Exit Code: ${exit_code})\033[0m"
        echo "流程已自动熔断终止。完整报错请查阅: ${log_file}"
        exit "$exit_code"
    else
        echo -e "\033[32m[成功] ${step_name} 执行完毕。 (Exit Code: 0)\033[0m"
    fi
}

# --------------------------------------------------------------------
# 管道按序执行区
# --------------------------------------------------------------------
echo ">>>>>>>>>>>>>>>>>>>> 开始运行去卷积分析管道 <<<<<<<<<<<<<<<<<<<<"
echo "当前处理组织: $TISSUE_NAME"

# Step 1: 拟Bulk表达矩阵构建
run_script "Step_1_Pseudo_bulk" "Rscript 2.1\\ Pseudo\\ bulk\\ construction.R --input_rds ${INPUT_PSEUDO} --output_dir ${OUTPUT_DIR} --gtf_path ${GTF_FILE}"

# Step 2: 基准评估与方法筛选
run_script "Step_2_Benchmarking" "Rscript 2.2\\ Deconvolution\\ benchmarking.R --input_rds ${INPUT_REF} --output_dir ${OUTPUT_DIR}"

# Step 3: DWLS 细胞比例反推
run_script "Step_3_DWLS" "Rscript 2.3\\ DWLS_CellComponents_deconvolution.R --input_ref ${INPUT_REF} --input_bulk ${OUTPUT_DIR}/Uniform_Fraction/random_pseudo/pseudo.txt --output_dir ${OUTPUT_DIR} --tissue ${TISSUE_NAME}"

# Step 4: bMIND 细胞特异性表达量反推
run_script "Step_4_bMIND" "Rscript 2.4\\ BMIND_CellExpression_deconvolution.R --project_dir ${PROJECT_DIR}"

# Step 5: MeDuSa 连续细胞状态拟合分析
run_script "Step_5_MeDuSA" "Rscript 2.5\\ MeduSa_CellState_deconvolution.R --project_dir ${PROJECT_DIR}"

echo -e "\n>>>>>>>>>>>>>>>>>>>> 全套去卷积分析管道安全跑通！ <<<<<<<<<<<<<<<<<<<<"
