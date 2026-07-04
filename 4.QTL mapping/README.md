========================================================================
                     QTL Mapping Pipeline User Guide
========================================================================

1. Pipeline Overview

This pipeline is a highly integrated workflow for expression quantification and multidimensional QTL mapping, including beQTL, cseQTL, and ieQTL analyses. The complete workflow is controlled by a master Bash script (run_qtl_pipeline.sh), which sequentially executes three R scripts and three Bash scripts for QTL mapping.

The pipeline supports command-line long options (e.g., --project_dir), real-time terminal output with simultaneous log recording, exit-code monitoring for each module, and automatic fail-fast termination upon fatal errors.


2. File Structure

Before deployment and execution, ensure that the following scripts are placed in the same directory on the server:

├── run_qtl_pipeline.sh                         # Master command-line pipeline script
├── 3.1 Prepare for beQTL.R                     # Step 1: Data preparation for beQTL analysis
├── 3.2 beQTL mapping.sh                        # Step 2: beQTL mapping
├── 3.3 Prepare for cseQTL.R                    # Step 3: Data preparation for cseQTL analysis
├── 3.4 cseQTL mapping.sh                       # Step 4: cseQTL mapping
├── 3.5 Prepare for cell-type ieQTL.R           # Step 5: Data preparation for ieQTL analysis
└── 3.6 ieQTL mapping.sh                        # Step 6: ieQTL mapping

Once the pipeline starts, the following directory will be automatically created in the current working directory:

└── qtl_pipeline_logs/                          # Real-time runtime logs for each module


3. Detailed Description, Modifications, and Critical Fixes for Each Module

------------------------------------------------------------------------
[Step 1] 3.1 Prepare for beQTL.R
------------------------------------------------------------------------

- Core function:

  Reads bulk RNA-seq count and TPM matrices and filters low-abundance genes. It then calculates TMM normalization factors and performs inverse normal transformation. Finally, a standard beQTL phenotype BED matrix is constructed using the species-specific GTF annotation file.

- Critical fixes and modifications:

  1. [Critical fix 1] In the original script, the undefined variable `gtf` was unexpectedly used in `region_annot <- gtf` at line 72, although no GTF file had been loaded or defined beforehand. Running the script independently therefore resulted in an `object 'gtf' not found` error. The updated script now imports the genome annotation using `rtracklayer::import(opt$gtf_file)` and converts it to a data frame, resolving this issue.

  2. [Critical fix 2] The original script used `list.files(pattern = ".txt")` to scan the current working directory without restriction. Unrelated text files in the directory could therefore be incorrectly included in the analysis. The updated script now scans the specified project directory directly.

  3. [CLI wrapper] The script has been refactored to use optparse-based command-line arguments, with all hard-coded absolute paths removed.


------------------------------------------------------------------------
[Step 2] 3.2 beQTL mapping.sh
------------------------------------------------------------------------

- Core function:

  Uses the compressed phenotype BED files generated in Step 3.1 together with genotype data to perform conventional cis-eQTL mapping across tissues using the linear mixed model (LMM) implemented in omiga.

- Modifications and workflow integration:

  The original omiga parameters were fully preserved without modification. The master script dynamically exports the `MAIN_DIR` environment variable, replacing the original hard-coded `main_dir` path prefix and allowing the project path to be specified externally without changing the core mapping logic.


------------------------------------------------------------------------
[Step 3] 3.3 Prepare for cseQTL.R
------------------------------------------------------------------------

- Core function:

  Loads the multidimensional bMIND RData containing deconvolution-derived cell-type-specific expression and extracts expression matrices for individual cell types. It then calculates cumulative exon lengths from the GTF annotation, performs TMM normalization and inverse normal transformation, and constructs cell-type-specific phenotype BED files for cseQTL analysis.

- Critical fixes and modifications:

  1. [Critical fix 1] The original script contained an undefined variable in the loop `for (j in 1:dim(dat1)[2])`, although the loaded multidimensional array was named `dat`. This apparent typographical error caused the loop to fail. The variable has been corrected to `dim(dat)[2]`, with multidimensional array alignment retained.

  2. [CLI wrapper] The `--project_dir` and `--gtf_file` options were introduced. All intermediate directories, including `CT_exp`, and downstream output paths are now dynamically constructed.


------------------------------------------------------------------------
[Step 4] 3.4 cseQTL mapping.sh
------------------------------------------------------------------------

- Core function:

  Performs nested traversal of tissue- and cell-type-specific subdirectories, matches the corresponding genotype data, and uses omiga to perform cis cell-type-specific eQTL mapping with a linear mixed model.

- Modifications and workflow integration:

  The original nested directory traversal structure was fully preserved. The master script exports the `MAIN_DIR_CSE` and `GENO_DIR_CSE` environment variables, which dynamically replace all internal hard-coded absolute paths.


------------------------------------------------------------------------
[Step 5] 3.5 Prepare for cell-type ieQTL.R
------------------------------------------------------------------------

- Core function:

  Reads cell-type proportion matrices predicted by DWLS and removes low-quality samples and cell types with more than 80% zero proportions. Outliers are filtered using the mean ± 3 standard deviations (3 × SD) criterion. The script then generates the phenotype BED files and cell-proportion interaction files required for ieQTL analysis.

- Modifications:

  Path handling was comprehensively refactored. The original hard-coded `/faststorage/project/cattle_gtexs` path was replaced with dynamically controlled project paths, while sample IDs across upstream and downstream data matrices are aligned consistently.


------------------------------------------------------------------------
[Step 6] 3.6 ieQTL mapping.sh
------------------------------------------------------------------------

- Core function:

  Traverses the cell-type interaction subdirectories for each tissue and supplies the corresponding cell-type proportions through the `--interaction` parameter to run the `cis_interaction` mode in omiga.

- Modifications and workflow integration:

  No changes to the core mapping script were required. The `MAIN_DIR_IE` environment variable exported by the master script dynamically controls path resolution and enables direct integration with the interaction matrices generated in Step 3.5.


4. Runtime Environment and Troubleshooting

1. Real-time terminal and log output (tee mechanism):

   The pipeline uses the standard Linux dual-output mechanism:

   eval ... 2>&1 | tee

   During execution, detailed logs from edgeR normalization, TMM calculation, bgzip compression, and omiga association mapping are displayed in real time in the terminal. The same output is simultaneously written to the `qtl_pipeline_logs/` directory for troubleshooting and archival purposes.

2. Exit-code monitoring and fail-fast termination:

   The master script uses the PIPESTATUS mechanism to capture the actual exit code of each R script and omiga process. If any module terminates with an error, the Bash controller immediately stops the entire workflow, preventing empty, incomplete, or invalid data from being passed to downstream steps.

   When the pipeline is terminated, the corresponding non-zero exit code is explicitly displayed in the terminal, facilitating troubleshooting, deployment, and maintenance.


5. Example Command for Running the Complete Pipeline

bash run_qtl_pipeline.sh \
  --project_dir "/faststorage/project/cattle_gtexs" \
  --gtf_file "/faststorage/project/cattle_gtexs/reference/Bos_taurus.ARS-UCD1.2.110.gtf"

========================================================================
