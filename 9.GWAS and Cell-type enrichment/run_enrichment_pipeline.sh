#!/bin/bash

# ====================================================================
# One-Stop GWAS & Cell-Type Enrichment Pipeline
# ====================================================================

# 1. Initialize variables
PROJECT_DIR=""
LOG_DIR="./enrichment_pipeline_logs"

# 2. Parse command-line long options
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --project_dir) PROJECT_DIR="$2"; shift ;;
        *) echo -e "\033[31m[ERROR] Unknown argument: $1\033[0m"; exit 1 ;;
    esac
    shift
done

if [ -z "$PROJECT_DIR" ]; then
    echo -e "\033[33m[USAGE]\033[0m bash run_enrichment_pipeline.sh --project_dir /path/to/project"
    exit 1
fi

mkdir -p "$LOG_DIR"

# --------------------------------------------------------------------
# Real-time dual-output logging and error interception
# --------------------------------------------------------------------
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
        echo -e "\033[31m[ERROR INTERCEPTED] This step returned a non-zero exit code (Exit Code: ${exit_code})\033[0m"
        exit "$exit_code"
    else
        echo -e "\033[32m[SUCCESS] ${step_name} completed successfully.\033[0m"
    fi
}

echo ">>>>>>>>>>>>>>>>>>>> Starting the GWAS & cell-type enrichment pipeline <<<<<<<<<<<<<<<<<<<<"
echo "Global project root: $PROJECT_DIR"

# Step 1: Run gsMap
run_step "9.1.Run_gsMap" "Rscript gsmap.R ${PROJECT_DIR}"

# Step 2: Run scPagwas
run_step "9.2.Run_scPagwas" "Rscript scpagwas.R ${PROJECT_DIR}"

echo -e "\n>>>>>>>>>>>>>>>>>>>> Enrichment pipeline completed successfully! <<<<<<<<<<<<<<<<<<<<"
