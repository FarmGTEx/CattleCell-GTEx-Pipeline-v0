#!/bin/bash

# ====================================================================
# Process 8: GWAS & eQTL 联合流转管道总控
# ====================================================================

PROJECT_DIR=""
LOG_DIR="./gwas_coloc_logs"

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --project_dir) PROJECT_DIR="$2"; shift ;;
        *) echo -e "\033[31m[错误] 未知参数: $1\033[0m"; exit 1 ;;
    esac
    shift
done

if [ -z "$PROJECT_DIR" ]; then
    echo -e "\033[33m[用法]\033[0m bash run_gwas_coloc_pipeline.sh --project_dir /path/to/project"
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
        echo -e "\033[31m[阻断拦截] 步骤异常退出 (Exit Code: ${exit_code})\033[0m"
        exit "$exit_code"
    else
        echo -e "\033[32m[成功] ${step_name} 执行完毕。\033[0m"
    fi
}

# ====================================================================
# 基于统一根目录动态推导上游与下游的相互嵌合输入路径
# ====================================================================
EQTL_PATH="${PROJECT_DIR}/CattleGTEx/OmiGA/eQTL"
GWAS_PATH="${PROJECT_DIR}/Downstream_analysis/complextraits/GWAS"
BFILE_FILE="${PROJECT_DIR}/CattleGTEx/panel_Hols"

# 构建中转站输出物理文件夹（确保输出咬合）
COLOC_INPUT_EQTL="${PROJECT_DIR}/CattleGTEx/OmiGA/Coloc/gwas_coloc/Bulk"
COLOC_INPUT_GWAS="${PROJECT_DIR}/CattleGTEx/OmiGA/Coloc/gwas_coloc/GWAS"
COLOC_FINAL_OUT="${PROJECT_DIR}/CattleGTEx/OmiGA/Coloc/gwas_coloc/Results"

echo ">>>>>>>>>>>>>>>>>>>> 开始运行 GWAS & eQTL 共定位全套管道 <<<<<<<<<<<<<<<<<<<<"

# Step 1: 运行 GCTA-COJO 格式转换 (R部分)
run_step "8.1.GCTA_Cojo_Format" "Rscript GCTA-COJO.R ${PROJECT_DIR}"

# Step 2: 运行 GCTA64 群体扫描条件分析 (原代码尾部剥离的 Bash 部分，100%原样保留)
echo -e "\n[执行中] 正在拉起 gcta64 独立信号条件扫描..."
for trait_dir in "$GWAS_PATH"/*/; do
    if [ -d "$trait_dir" ]; then
        gwas_file="${trait_dir}$(basename "$trait_dir")_cojo.txt"
        out_dir="${trait_dir}$(basename "$trait_dir")_cojo_res"
        if [ -f "$gwas_file" ]; then
            gcta64 --bfile "$BFILE_FILE" \
                   --maf 0.05 \
                   --cojo-file "$gwas_file" \
                   --cojo-slct \
                   --cojo-p 1e-5 \
                   --out "$out_dir"
        fi
    fi
done

# Step 3: 准备 eQTL 共定位输入块 (1.Prepare_eQTL.R)
run_step "8.2.Prepare_eQTL" "Rscript 1.Prepare_eQTL.R --input_dir ${EQTL_PATH} --output_dir ${COLOC_INPUT_EQTL}"

# Step 4: 准备 GWAS 窗口外延输入块 (2.Prepare_GWAS.R，使其与上游的 GWAS_PATH 嵌合)
run_step "8.3.Prepare_GWAS" "Rscript 2.Prepare_GWAS.R --gwas_dir ${GWAS_PATH} --cojo_dir ${GWAS_PATH} --output_dir ${COLOC_INPUT_GWAS}"

# Step 5: 启动最终的贝叶斯共定位核心运算 (3.run_coloc.R)
run_step "8.4.Run_Coloc" "Rscript 3.run_coloc.R --gwas_dir ${COLOC_INPUT_GWAS} --bulk_dir ${COLOC_INPUT_EQTL} --result_dir ${COLOC_FINAL_OUT}"

echo -e "\n>>>>>>>>>>>>>>>>>>>> 全套 GWAS-eQTL 嵌合管道运行结束！ <<<<<<<<<<<<<<<<<<<<"
