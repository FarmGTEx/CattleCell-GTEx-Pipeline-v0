#!/bin/bash

# ====================================================================
# Master Control Script for the Tissue & Cell Type Sharing MashR Pipeline
# ====================================================================

# 1. Initialize variables
PROJECT_DIR=""
MODE=""
LOG_DIR="./mashr_pipeline_logs"
export SCRIPT_DIR="$(pwd)" # Directory containing the Python and R scripts

# 2. Parse command-line long options
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --project_dir) PROJECT_DIR="$2"; shift ;;
        --mode)        MODE="$2"; shift ;;
        *) echo -e "\033[31m[ERROR] Unknown argument: $1\033[0m"; exit 1 ;;
    esac
    shift
done

# 3. Validate arguments
if [ -z "$PROJECT_DIR" ] || [ -z "$MODE" ]; then
    echo -e "\033[33m[USAGE]\033[0m bash run_mashr_pipeline.sh --project_dir <path> --mode <bulk|celltype>"
    exit 1
fi

if [[ "$MODE" != "bulk" && "$MODE" != "celltype" ]]; then
    echo -e "\033[31m[ERROR] The --mode argument must be either 'bulk' or 'celltype'.\033[0m"
    exit 1
fi

mkdir -p "$LOG_DIR"

# --------------------------------------------------------------------
# Function: script executor with real-time logging and exit-code capture
# --------------------------------------------------------------------
run_step() {
    local step_name=$1
    local script_cmd=$2
    local log_file="${LOG_DIR}/${step_name}_${MODE}.log"

    echo -e "\n========================================================================"
    echo "[PIPELINE START] Running: ${step_name} (${MODE} mode)"
    echo "--------------------------- Real-time terminal output ---------------------------"

    eval "${script_cmd}" 2>&1 | tee "${log_file}"
    local exit_code=${PIPESTATUS[0]}

    echo "--------------------------------------------------------------------"
    if [ "$exit_code" -ne 0 ]; then
        echo -e "\033[31m[FATAL ERROR] This step failed! (Exit Code: ${exit_code})\033[0m"
        echo "The pipeline has been safely terminated. See the full log at: ${log_file}"
        exit "$exit_code"
    else
        echo -e "\033[32m[SUCCESS] ${step_name} completed successfully.\033[0m"
    fi
}

# --------------------------------------------------------------------
# Pipeline logic and environment-variable dispatch
# --------------------------------------------------------------------
echo ">>>>>>>>>>>>>>>>>>>> Starting the MashR Sharing pipeline <<<<<<<<<<<<<<<<<<<<"

if [ "$MODE" == "bulk" ]; then
    # Bulk mode: export the corresponding base environment paths
    export NOMINAL_DIR="${PROJECT_DIR}/CattleGTEx/OmiGA/eQTL/"
    export RANDOM_OUT="${PROJECT_DIR}/CattleGTEx/OmiGA/MASHR_Celltype/Tissue_sharing/Bulk/Random_output"
    export STRONG_OUT="${PROJECT_DIR}/CattleGTEx/OmiGA/MASHR_Celltype/Tissue_sharing/Bulk/strong_output"
    
    mkdir -p "$RANDOM_OUT" "$STRONG_OUT"

    export OUTPUT_DIR="$RANDOM_OUT"
    export SUBSET_SIZE=1000000
    run_step "1.Bulk_Prepare_Random" "bash Bulk-sharing-Prepare-MashR-Random.sh"

    export OUTPUT_DIR="$STRONG_OUT"
    run_step "2.Bulk_Prepare_Strong" "bash Bulk-sharing-Prepare-MashR-Strong.sh"

    # Set variables passed to run_MashR.sh
    export STRONG_FILE="${STRONG_OUT}/strong_pairs.MashR_input.txt.gz"
    export RANDOM_FILE="${RANDOM_OUT}/nominal_pairs.${SUBSET_SIZE}_subset.MashR_input.txt.gz"
    export MASHR_OUT_DIR="${STRONG_OUT}/output_top_pairs"
    
    run_step "3.Run_MashR_Bulk" "bash run_MashR.sh"

elif [ "$MODE" == "celltype" ]; then
    # Celltype mode: export the corresponding base environment paths
    export NOMINAL_DIR="${PROJECT_DIR}/CattleGTEx/OmiGA/cell_specific/"
    export RANDOM_OUT="${PROJECT_DIR}/CattleGTEx/OmiGA/MASHR_Celltype/Celltype_sharing/random_output"
    export STRONG_OUT="${PROJECT_DIR}/CattleGTEx/OmiGA/MASHR_Celltype/Celltype_sharing/strong_output"
    
    mkdir -p "$RANDOM_OUT" "$STRONG_OUT"

    export OUTPUT_DIR="$RANDOM_OUT"
    export SUBSET_SIZE=10000  # Fix the missing variable in the original script
    run_step "1.Celltype_Prepare_Random" "bash Celltype-sharing-Prepare-MashR-Random.sh"

    export OUTPUT_DIR="$STRONG_OUT"
    run_step "2.Celltype_Prepare_Strong" "bash Celltype-sharing-Prepare-MashR-Strong.sh"

    # Set variables passed to run_MashR.sh
    export STRONG_FILE="${STRONG_OUT}/strong_pairs.MashR_input.txt.gz"
    export RANDOM_FILE="${RANDOM_OUT}/nominal_pairs.${SUBSET_SIZE}_subset.MashR_input.txt.gz"
    export MASHR_OUT_DIR="${STRONG_OUT}/output_top_pairs"
    
    run_step "3.Run_MashR_Celltype" "bash run_MashR.sh"
fi

echo -e "\n>>>>>>>>>>>>>>>>>>>> MashR pipeline completed successfully! <<<<<<<<<<<<<<<<<<<<"
