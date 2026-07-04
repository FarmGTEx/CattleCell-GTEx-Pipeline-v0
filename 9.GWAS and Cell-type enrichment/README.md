========================================================================
     Process 9: GWAS & Cell-Type Enrichment Pipeline User Guide
========================================================================

1. Module Overview

This module consists of two core components for GWAS-based cell-type enrichment analysis using single-cell and spatial transcriptomic data:

1. gsmap.R      (spatial gene mapping-based enrichment analysis)
2. scpagwas.R   (single-cell pathway-based trait association analysis)

Both components have been refactored to support externally specified project paths without relying on hard-coded locations.

A unified master Bash script (run_enrichment_pipeline.sh) now automatically executes the two components sequentially and supports real-time terminal output with simultaneous log recording.


2. Directory and File Structure

Before running the pipeline, ensure that the current directory (`9.GWAS and Cell-type enrichment`) contains the following three scripts:

9.GWAS and Cell-type enrichment/
├── gsmap.R                       # Wrapped gsMap enrichment script
├── scpagwas.R                    # Wrapped scPagwas enrichment script
└── run_enrichment_pipeline.sh    # Master Bash script for unified pipeline execution


3. Command-Line Usage Guide

[Step 1] Grant execution permission to the master script

Run the following command in the terminal:

chmod +x run_enrichment_pipeline.sh


[Step 2] Run the complete pipeline

Launch the master Bash script and specify the global project root directory using the `--project_dir` long option:

bash run_enrichment_pipeline.sh --project_dir "/faststorage/project/cattle_gtexs"


[Step 3] Monitor runtime logs

Once the pipeline starts, the `./enrichment_pipeline_logs/` directory will be automatically created in the current working directory.

All R-based matrix operations, Seurat normalization progress, and error messages will be displayed in real time in the terminal and simultaneously written to the following log files for unattended execution and subsequent troubleshooting:

- gsMap runtime log:
  ./enrichment_pipeline_logs/9.1.Run_gsMap.log

- scPagwas runtime log:
  ./enrichment_pipeline_logs/9.2.Run_scPagwas.log


4. Manual Execution of Individual Components

For independent debugging or execution of a specific method, each component can be run directly using Rscript, with the project root directory passed as the first positional argument.

1. Run gsMap independently:

Rscript gsmap.R "/faststorage/project/cattle_gtexs"

2. Run scPagwas independently:

Rscript scpagwas.R "/faststorage/project/cattle_gtexs"
