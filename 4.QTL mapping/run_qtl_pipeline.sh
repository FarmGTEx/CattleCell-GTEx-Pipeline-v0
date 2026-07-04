#!/bin/bash

# ====================================================================
# Master Control Script for the QTL Mapping Pipeline
# ====================================================================

# 1. Initialize variables
PROJECT_DIR=""
GTF_FILE=""
LOG_DIR="./qtl_pipeline_logs"

# 2. Parse command-line long options
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --project_dir) PROJECT_DIR="$2"; shift ;;
        --gtf_file)    GTF_FILE="$2"; shift ;;
        *) echo -e "\033[31m[ERROR] Unknown argument: $1\033[0m"; exit 1 ;;
    esac
    shift
done

# 3. Validate required arguments
if [ -z "$PROJECT_DIR" ] || [ -z "$GTF_FILE" ]; then
    echo -e "\033[33m[USAGE EXAMPLE]\033[0m bash run_qtl_pipeline.sh --project_dir /path/to/project --gtf_file /path/to/gtf"
    echo -e "\033[31m[ERROR] Missing required arguments! Please specify both --project_dir and --gtf_file.\033[0m"
    exit 1
fi

# 4. Create log directory
mkdir -p "$LOG_DIR"

# --------------------------------------------------------------------
# Define a function for script execution and error handling
# Dual output: terminal + log file
# --------------------------------------------------------------------
run_step() {
    local step_name=$1
    local script_cmd=$2
    shift 2 # Remove the first two arguments and pass all remaining arguments to the downstream script
    
    local log_file="${LOG_DIR}/${step_name}.log"

    echo -e "\n========================================================================"
    echo "[PIPELINE START] Running: ${step_name}"
    echo "--------------------------- Real-time terminal output ---------------------------"

    # Execute the command, merge stdout and stderr, write output to both terminal and log file,
    # and capture the true exit code of the upstream command using PIPESTATUS[0]
    eval "${script_cmd}" "$@" 2>&1 | tee "${log_file}"
    local exit_code=${PIPESTATUS[0]}

    echo "--------------------------------------------------------------------"
    if [ "$exit_code" -ne 0 ]; then
        echo -e "\033[31m[FATAL ERROR] This step failed! (Exit Code: ${exit_code})\033[0m"
        echo "The pipeline has triggered fail-fast termination. See the full log at: ${log_file}"
        exit "$exit_code"
    else
        echo -e "\033[32m[SUCCESS] ${step_name} completed successfully. (Exit Code: 0)\033[0m"
    fi
}

# --------------------------------------------------------------------
# Sequential pipeline execution
# --------------------------------------------------------------------
echo ">>>>>>>>>>>>>>>>>>>> Starting the automated QTL Mapping pipeline <<<<<<<<<<<<<<<<<<<<"
echo "Start time: $(date)"
echo "Project root directory: $PROJECT_DIR"

# --- Step 1: Prepare beQTL data (R script with escaped spaces) ---
run_step "1.Prepare_beQTL" "Rscript 3.1\\ Prepare\\ for\\ beQTL.R" \
  --project_dir "$PROJECT_DIR" \
  --gtf_file "$GTF_FILE"

# --- Step 2: beQTL association mapping (call shell script) ---
# Note: the original 3.2 script contains a hard-coded main_dir.
# To allow external control, export the variable as an environment variable here.
export MAIN_DIR="${PROJECT_DIR}/CattleGTEx/OmiGA/eQTL/"
run_step "2.beQTL_Mapping" "bash 3.2\\ beQTL\\ mapping.sh"


# --- Step 3: Prepare cseQTL data (R script with escaped spaces) ---
run_step "3.Prepare_cseQTL" "Rscript 3.3\\ Prepare\\ for\\ cseQTL.R" \
  --project_dir "$PROJECT_DIR" \
  --gtf_file "$GTF_FILE"

# --- Step 4: cseQTL association mapping (call shell script) ---
export MAIN_DIR_CSE="${PROJECT_DIR}/CattleGTEx/OmiGA/cell_specific/Group6/"
export GENO_DIR_CSE="${PROJECT_DIR}/CattleGTEx/OmiGA/eQTL/Group6/"
run_step "4.cseQTL_Mapping" "bash 3.4\\ cseQTL\\ mapping.sh"


# --- Step 5: Prepare ieQTL data (R script with escaped spaces) ---
run_step "5.Prepare_ieQTL" "Rscript 3.5\\ Prepare\\ for\\ cell-type\\ ieQTL.R" \
  --project_dir "$PROJECT_DIR"

# --- Step 6: ieQTL cell-interaction mapping (call shell script) ---
export MAIN_DIR_IE="${PROJECT_DIR}/CattleGTEx/OmiGA/eQTL/"
run_step "6.ieQTL_Mapping" "bash 3.6\\ ieQTL\\ mapping.sh"


echo -e "\n>>>>>>>>>>>>>>>>>>>> QTL Mapping pipeline completed successfully! <<<<<<<<<<<<<<<<<<<<"
