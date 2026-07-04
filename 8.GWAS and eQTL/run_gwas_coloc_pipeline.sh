#!/bin/bash

# ====================================================================
# Process 8: Master Control Script for the GWAS & eQTL Integration Pipeline
# ====================================================================

PROJECT_DIR=""
LOG_DIR="./gwas_coloc_logs"

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --project_dir) PROJECT_DIR="$2"; shift ;;
        *) echo -e "\033[31m[ERROR] Unknown argument: $1\033[0m"; exit 1 ;;
    esac
    shift
done

if [ -z "$PROJECT_DIR" ]; then
    echo -e "\033[33m[USAGE]\033[0m bash run_gwas_coloc_pipeline.sh --project_dir /path/to/project"
    exit 1
fi

mkdir -p "$LOG_DIR"

run_step() {
    local step_name=$1
    local script_cmd=$2
    local log_file="${LOG_DIR}/${step_name}.log"

    echo -e "\n========================================================================"
    echo "[PIPELINE START] Running: ${step_name}"
    echo "--------------------------- Real-time terminal output ---------------------------"

    eval "${script_cmd}" 2>&1 | tee "${log_file}"
    local exit_code=${PIPESTATUS[0]}

    echo "--------------------------------------------------------------------"
    if [ "$exit_code" -ne 0 ]; then
        echo -e "\033[31m[STOPPED] Step failed with a non-zero exit code (Exit Code: ${exit_code})\033[0m"
        exit "$exit_code"
    else
        echo -e "\033[32m[SUCCESS] ${step_name} completed successfully.\033[0m"
    fi
}

# ====================================================================
# Dynamically derive upstream and downstream input paths from the shared root directory
# ====================================================================
EQTL_PATH="${PROJECT_DIR}/CattleGTEx/OmiGA/eQTL"
GWAS_PATH="${PROJECT_DIR}/Downstream_analysis/complextraits/GWAS"
BFILE_FILE="${PROJECT_DIR}/CattleGTEx/panel_Hols"

# Create intermediate output directories to ensure consistent handoff between steps
COLOC_INPUT_EQTL="${PROJECT_DIR}/CattleGTEx/OmiGA/Coloc/gwas_coloc/Bulk"
COLOC_INPUT_GWAS="${PROJECT_DIR}/CattleGTEx/OmiGA/Coloc/gwas_coloc/GWAS"
COLOC_FINAL_OUT="${PROJECT_DIR}/CattleGTEx/OmiGA/Coloc/gwas_coloc/Results"

echo ">>>>>>>>>>>>>>>>>>>> Starting the full GWAS & eQTL colocalization pipeline <<<<<<<<<<<<<<<<<<<<"

# Step 1: Run GCTA-COJO format conversion using the R script
run_step "8.1.GCTA_Cojo_Format" "Rscript GCTA-COJO.R ${PROJECT_DIR}"

# Step 2: Run GCTA64 conditional analysis for independent signal detection
# This Bash section was extracted from the original script and preserved unchanged
echo -e "\n[RUNNING] Starting gcta64 conditional scan for independent signals..."
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

# Step 3: Prepare eQTL input blocks for colocalization
run_step "8.2.Prepare_eQTL" "Rscript 1.Prepare_eQTL.R --input_dir ${EQTL_PATH} --output_dir ${COLOC_INPUT_EQTL}"

# Step 4: Prepare GWAS window-expanded input blocks and align them with GWAS_PATH
run_step "8.3.Prepare_GWAS" "Rscript 2.Prepare_GWAS.R --gwas_dir ${GWAS_PATH} --cojo_dir ${GWAS_PATH} --output_dir ${COLOC_INPUT_GWAS}"

# Step 5: Run the final Bayesian colocalization analysis
run_step "8.4.Run_Coloc" "Rscript 3.run_coloc.R --gwas_dir ${COLOC_INPUT_GWAS} --bulk_dir ${COLOC_INPUT_EQTL} --result_dir ${COLOC_FINAL_OUT}"

echo -e "\n>>>>>>>>>>>>>>>>>>>> GWAS-eQTL integration pipeline completed successfully! <<<<<<<<<<<<<<<<<<<<"
