#!/bin/bash

# ====================================================================
# Automated Pipeline for Selection Region Analysis
# ====================================================================

PROJECT_DIR=""
LOG_DIR="./selection_pipeline_logs"

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --project_dir) PROJECT_DIR="$2"; shift ;;
        *) echo -e "\033[31m[ERROR] Unknown argument: $1\033[0m"; exit 1 ;;
    esac
    shift
done

if [ -z "$PROJECT_DIR" ]; then
    echo -e "\033[33m[USAGE]\033[0m bash run_selection_pipeline.sh --project_dir /path/to/project"
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
        echo -e "\033[31m[ERROR INTERCEPTED] Component returned a non-zero exit code (Exit Code: ${exit_code})\033[0m"
        exit "$exit_code"
    else
        echo -e "\033[32m[SUCCESS] ${step_name} completed successfully.\033[0m"
    fi
}

echo ">>>>>>>>>>>>>>>>>>>> Starting the selection pressure and multidimensional signal enrichment pipeline <<<<<<<<<<<<<<<<<<<<"
echo "Project root directory: $PROJECT_DIR"

# Step 1: Run enrichment analysis between GWAS signals and Fst-based selection regions
run_step "10.1.GWAS_Selection_Enrich" "Rscript Enrichment\\ between\\ GWAS\\ and\\ selection\\ regions.R ${PROJECT_DIR}"

# Step 2: Run enrichment analysis between fine-mapped eQTL signals and Fst-based selection regions
run_step "10.2.eQTL_Selection_Enrich" "Rscript Enrichment\\ between\\ eQTL\\ and\\ selection\\ regions.R ${PROJECT_DIR}"

# Step 3: Extract the overlap matrix for the top 10% strongly selected regions
run_step "10.3.Strong_Selected_Enrich" "Rscript Strongly\\ selected\\ regions.R ${PROJECT_DIR}"

echo -e "\n>>>>>>>>>>>>>>>>>>>> Selection pipeline completed successfully! <<<<<<<<<<<<<<<<<<<<"
