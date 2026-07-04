========================================================================
        Tissue & Cell Type Sharing (MashR) Pipeline User Guide
========================================================================

1. Pipeline Overview

This pipeline is designed for QTL-sharing analysis across tissues (bulk level) and cell types (cell-type level) using the MashR empirical Bayes framework.

The complete workflow is managed by a master Bash script (run_mashr_pipeline.sh). According to the specified mode (--mode), the pipeline automatically schedules four data-preparation scripts and one core MashR execution script.

The pipeline supports command-line long options for specifying the base project environment, real-time terminal output with simultaneous log recording through the tee mechanism, and automatic fail-fast termination based on system exit-code monitoring.


2. File Structure

Before deployment and execution, ensure that all scripts and core tools are placed in the same directory on the server:

├── run_mashr_pipeline.sh                       # Master command-line pipeline script
├── Bulk-sharing-Prepare-MashR-Random.sh        # Bulk-level random-set preparation script
├── Bulk-sharing-Prepare-MashR-Strong.sh        # Bulk-level strong-set preparation script
├── Celltype-sharing-Prepare-MashR-Random.sh    # Cell-type-level random-set preparation script
├── Celltype-sharing-Prepare-MashR-Strong.sh    # Cell-type-level strong-set preparation script
├── run_MashR.sh                                # Core MashR execution script
├── MashR-random_subset.R                       # R script for random subset sampling
├── run_MashR.R                                 # R script for MashR empirical Bayes estimation
├── combine_signif_pairs_tjy.py                 # Python tool for combining significant pairs
├── extract_pairs_tjy.py                        # Python tool for extracting association pairs
└── mashr_prepare_input.py                      # Python tool for matrix formatting and conversion

Once the pipeline starts, the following directory will be automatically created in the current working directory:

└── mashr_pipeline_logs/                        # Real-time runtime logs for each module

------------------------------------------------------------------------
[Module 1] Bulk-sharing-Prepare-MashR-Random.sh & Strong.sh
------------------------------------------------------------------------

- Core function:

  The Random script scans all association pairs from conventional beQTL results, calls the Python tools to merge them, and randomly samples 1,000,000 test pairs to construct the null correlation matrix.

  The Strong script filters permutation-adjusted significant eGene association pairs across tissues and extracts the corresponding Z-score matrix as the true-effect set.

------------------------------------------------------------------------
[Module 2] Celltype-sharing-Prepare-MashR-Random.sh & Strong.sh
------------------------------------------------------------------------

- Core function:

  Corresponding to the bulk-level workflow, these scripts perform nested directory scanning of cell-type-specific cseQTL results and construct both the random subset matrix and the strong-effect matrix for cell-type-level effect sharing. The random subset contains 10,000 test pairs by default.

------------------------------------------------------------------------
[Module 3] MashR-random_subset.R
------------------------------------------------------------------------

- Core function:

  During bulk-level merged sampling, this script sequentially loads the chromosome-specific Z-score matrices, merges them, performs random sampling, and serializes the result as a standard `.RDS` object.

------------------------------------------------------------------------
[Module 4] run_MashR.sh & run_MashR.R
------------------------------------------------------------------------

- Core function:

  These scripts load the prepared Strong and Random matrices using multithreading. The random set is first used to estimate the null correlation structure (Vhat). Data-driven covariance matrices are then estimated using PCA and ED, followed by mixture-proportion fitting. Finally, posterior quantities are calculated and exported as RDS files, including lfsr, pm, and seven other core MashR outputs.

4. Troubleshooting and Fail-Fast Termination

1. Real-time terminal and log output:

   The pipeline uses the `tee` command to enable dual output streams. Matrix merging, Python extraction progress, and mixture-proportion fitting iterations from the mashr R package are displayed in real time in the terminal and are simultaneously written to the corresponding log files under `mashr_pipeline_logs/`.

2. Exit-code interception and fail-fast termination:

   The pipeline uses PIPESTATUS[0] to capture the actual exit code of each step. If Python pandas raises an exception or an R script terminates with `stop()` because of mismatched matrix row or column names, the master script immediately triggers fail-fast termination, stops all downstream steps, and prints the non-zero exit code in red in the terminal. This prevents incomplete or invalid intermediate data from propagating downstream.


5. Example Commands

The master script provides the `--mode` option to switch between the two analysis modes.

[Example 1: Run tissue-level sharing analysis in bulk mode]

bash run_mashr_pipeline.sh \
  --project_dir "/faststorage/project/cattle_gtexs" \
  --mode bulk

[Example 2: Run cell-type-level sharing analysis in celltype mode]

bash run_mashr_pipeline.sh \
  --project_dir "/faststorage/project/cattle_gtexs" \
  --mode celltype

========================================================================
