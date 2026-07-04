#!/bin/bash

# ====================================================================
# Master Control Script for the Deconvolution Pipeline
# ====================================================================

# 1. Initialize variables
INPUT_PSEUDO=""
INPUT_REF=""
GTF_FILE=""
OUTPUT_DIR=""
PROJECT_DIR=""
TISSUE_NAME=""
LOG_DIR="./pipeline_logs"

# 2. Parse command-line long options
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input_pseudo) INPUT_PSEUDO="$2"; shift ;;
        --input_ref)    INPUT_REF="$2"; shift ;;
        --gtf_file)     GTF_FILE="$2"; shift ;;
        --output_dir)   OUTPUT_DIR="$2"; shift ;;
        --project_dir)  PROJECT_DIR="$2"; shift ;;
        --tissue)       TISSUE_NAME="$2"; shift ;;
        *) echo -e "\033[31m[ERROR] Unknown argument: $1\033[0m"; exit 1 ;;
    esac
    shift
done

# 3. Validate required arguments
if [ -z "$INPUT_PSEUDO" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$TISSUE_NAME" ]; then
    echo -e "\033[33m[USAGE EXAMPLE]\033[0m bash run_deconvolution.sh --input_pseudo <rds> --input_ref <rds> --gtf_file <gtf> --output_dir <dir> --project_dir <dir> --tissue <name>"
    echo -e "\033[31m[ERROR] Missing required arguments! Please provide at least --input_pseudo, --output_dir, and --tissue.\033[0m"
    exit 1
fi

# 4. Create log and output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# --------------------------------------------------------------------
# Define a function for script execution and error handling
# --------------------------------------------------------------------
run_script() {
    local step_name=$1
    local script_cmd=$2 # Receive the complete escaped Rscript command
    local log_file="${LOG_DIR}/${step_name}.log"

    echo -e "\n================================================================"
    echo "[START] Running: ${step_name}"
    echo "---------------------- Real-time terminal output ----------------------"

    # Execute the complete pre-escaped command using eval
    eval "${script_cmd}" 2>&1 | tee "${log_file}"
    local exit_code=${PIPESTATUS[0]}

    echo "------------------------------------------------------------------------"
    if [ "$exit_code" -ne 0 ]; then
        echo -e "\033[31m[FATAL ERROR] Step ${step_name} failed! (Exit Code: ${exit_code})\033[0m"
        echo "The pipeline has been terminated automatically. See the full error log at: ${log_file}"
        exit "$exit_code"
    else
        echo -e "\033[32m[SUCCESS] ${step_name} completed successfully. (Exit Code: 0)\033[0m"
    fi
}

# --------------------------------------------------------------------
# Sequential pipeline execution
# --------------------------------------------------------------------
echo ">>>>>>>>>>>>>>>>>>>> Starting the deconvolution pipeline <<<<<<<<<<<<<<<<<<<<"
echo "Current tissue: $TISSUE_NAME"

# Step 1: Construct pseudo-bulk expression matrices
run_script "Step_1_Pseudo_bulk" "Rscript 2.1\\ Pseudo\\ bulk\\ construction.R --input_rds ${INPUT_PSEUDO} --output_dir ${OUTPUT_DIR} --gtf_path ${GTF_FILE}"

# Step 2: Benchmark and evaluate deconvolution methods
run_script "Step_2_Benchmarking" "Rscript 2.2\\ Deconvolution\\ benchmarking.R --input_rds ${INPUT_REF} --output_dir ${OUTPUT_DIR}"

# Step 3: Estimate cell-type proportions using DWLS
run_script "Step_3_DWLS" "Rscript 2.3\\ DWLS_CellComponents_deconvolution.R --input_ref ${INPUT_REF} --input_bulk ${OUTPUT_DIR}/Uniform_Fraction/random_pseudo/pseudo.txt --output_dir ${OUTPUT_DIR} --tissue ${TISSUE_NAME}"

# Step 4: Estimate cell-type-specific expression using bMIND
run_script "Step_4_bMIND" "Rscript 2.4\\ BMIND_CellExpression_deconvolution.R --project_dir ${PROJECT_DIR}"

# Step 5: Infer continuous cell states using MeDuSA
run_script "Step_5_MeDuSA" "Rscript 2.5\\ MeduSa_CellState_deconvolution.R --project_dir ${PROJECT_DIR}"

echo -e "\n>>>>>>>>>>>>>>>>>>>> Deconvolution pipeline completed successfully! <<<<<<<<<<<<<<<<<<<<"
