========================================================================
       Process 10: Selection Region Analysis Pipeline User Guide
========================================================================

1. Module Overview

This module consists of three R scripts for analyzing enrichment overlap between selection signals and functional genomic features, together with one vcftools-based genome-wide scanning script.

All components have been refactored to support external project-root path specification through PROJECT_DIR, removing hard-coded paths.


2. Directory and File Structure

Before automated execution, ensure that the following five core files are placed in the same directory:

10.Selection region analysis/
├── Enrichment between GWAS and selection regions.R  # GWAS enrichment component
├── Enrichment between eQTL and selection regions.R  # eQTL enrichment component
├── Strongly selected regions.R                      # Top 10% strong-selection extraction component
├── fst.sh                                           # Revised Slurm submission script
└── run_selection_pipeline.sh                        # Master Bash script for unified execution


3. Command-Line Usage Guide

[Step 1] Grant execution permission to the master script

Run the following command in the current terminal:

chmod +x run_selection_pipeline.sh


[Step 2] Run the complete pipeline

Launch the master Bash script and specify the global project root directory using the `--project_dir` option:

bash run_selection_pipeline.sh --project_dir "/faststorage/project/cattle_gtexs"


[Step 3] Monitor runtime logs

After the master script starts, it will sequentially execute the three R components and automatically create the `./selection_pipeline_logs/` directory in the current working directory.

Runtime information, including multi-tissue SuSiE mapping file retrieval and efficient interval overlap using data.table::foverlaps, will be displayed in real time and simultaneously written to the following log files:

- GWAS enrichment log:
  ./selection_pipeline_logs/10.1.GWAS_Selection_Enrich.log

- eQTL enrichment log:
  ./selection_pipeline_logs/10.2.eQTL_Selection_Enrich.log

- Strong-selection filtering log:
  ./selection_pipeline_logs/10.3.Strong_Selected_Enrich.log


4. Manual Execution of Individual Components

For step-by-step debugging of a specific analysis module or population, each component can be executed independently using Rscript, with the global project root directory passed as the first positional argument.

1. Run GWAS enrichment independently:

Rscript Enrichment\ between\ GWAS\ and\ selection\ regions.R "/faststorage/project/cattle_gtexs"

2. Run eQTL enrichment independently:

Rscript Enrichment\ between\ eQTL\ and\ selection\ regions.R "/faststorage/project/cattle_gtexs"

3. Run extraction of the top 10% strongly selected regions independently:

Rscript Strongly\ selected\ regions.R "/faststorage/project/cattle_gtexs"

========================================================================
