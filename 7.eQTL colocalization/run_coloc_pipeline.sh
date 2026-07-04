#!/bin/bash

# ====================================================================
# Master Control Script for the Automated eQTL Colocalization Pipeline
# ====================================================================

PROJECT_DIR=""
LOG_DIR="./coloc_pipeline_logs"

# Parse command-line long options
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --project_dir) PROJECT_DIR="$2"; shift ;;
        *) echo -e "\033[31m[ERROR] Unknown argument: $1\033[0m"; exit 1 ;;
    esac
    shift
done

if [ -z "$PROJECT_DIR" ]; then
    echo -e "\033[33m[USAGE]\033[0m bash run_coloc_pipeline.sh --project_dir /path/to/project"
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
        echo -e "\033[31m[FATAL ERROR] This step failed! (Exit Code: ${exit_code})\033[0m"
        exit "$exit_code"
    else
        echo -e "\033[32m[SUCCESS] ${step_name} completed successfully.\033[0m"
    fi
}

echo ">>>>>>>>>>>>>>>>>>>> Starting the eQTL Colocalization pipeline <<<<<<<<<<<<<<<<<<<<"

# Step 1: Prepare phenotype BED matrix inputs for colocalization analysis
# Pass the project root directory directly as a positional argument
run_step "7.1.Prepare_Coloc" "Rscript 7.1\\ Prepare\\ for\\ coloc.R ${PROJECT_DIR}"

# Step 2: Run the core Bayesian ABF colocalization analysis
run_step "7.2.Run_Coloc_ABF" "Rscript 7.2\\ Coloc.R ${PROJECT_DIR}"

echo -e "\n>>>>>>>>>>>>>>>>>>>> Colocalization pipeline completed successfully! <<<<<<<<<<<<<<<<<<<<"
