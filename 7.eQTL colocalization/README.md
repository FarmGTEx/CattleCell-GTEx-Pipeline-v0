========================================================================
       Process 7: eQTL Colocalization Pipeline User Guide
========================================================================

1. Module Overview

This module consists of two core components for Bayesian colocalization analysis:

1. 7.1 Prepare for coloc.R   (extraction and preprocessing of colocalization input matrices)
2. 7.2 Coloc.R               (core Bayesian colocalization analysis using coloc.abf)

The main purpose of this module is to sequentially integrate four levels of QTL signals, including conventional bulk eQTL, cell-component ieQTL, cell-type-specific cseQTL, and cell-state ieQTL. These signals are overlapped with fine-mapped genomic windows and subsequently analyzed using the coloc algorithm to estimate posterior probabilities (PP0–PP4).

The global project root path (PROJECT_DIR) has been fully decoupled from hard-coded paths using native R command-line argument parsing.


2. Directory and File Structure

Before running the pipeline, ensure that the current directory (`7.eQTL colocalization`) contains the following three core files:

7.eQTL colocalization/
├── 7.1 Prepare for coloc.R   # Wrapped data preparation script for colocalization analysis
├── 7.2 Coloc.R               # Wrapped script for the core coloc analysis
└── run_coloc_pipeline.sh     # Master Bash script for unified pipeline execution


3. Command-Line Usage Guide

[Step 1] Grant execution permission to the master script

Run the following command in the current directory:

chmod +x run_coloc_pipeline.sh


[Step 2] Run the complete pipeline

Launch the master Bash script and specify the global project root directory using the `--project_dir` long option:

bash run_coloc_pipeline.sh --project_dir "/faststorage/project/cattle_gtexs"


[Step 3] Monitor runtime logs and quality-control output

Once the pipeline starts, the `./coloc_pipeline_logs/` directory will be automatically created in the current working directory.

The merging status of chromosome 1–29 nominal-pair matrices across multiple QTL levels, together with runtime output from the coloc package during ABF-based posterior probability estimation, will be displayed in real time in the terminal and simultaneously written to the following log files:

- Data preparation log:
  ./coloc_pipeline_logs/7.1.Prepare_Coloc.log

- Core colocalization analysis log:
  ./coloc_pipeline_logs/7.2.Run_Coloc_ABF.log


4. Manual Execution of Individual Components

For step-by-step debugging of a specific tissue or problematic data block, each script can be executed independently using Rscript, with the global project root directory passed as the first positional argument.

1. Run the colocalization data preparation step independently:

Rscript 7.1\ Prepare\ for\ coloc.R "/faststorage/project/cattle_gtexs"

2. Run the Bayesian colocalization analysis independently:

Rscript 7.2\ Coloc.R "/faststorage/project/cattle_gtexs"
