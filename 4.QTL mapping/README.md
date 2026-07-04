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

------------------------------------------------------------------------
[Step 2] 3.2 beQTL mapping.sh
------------------------------------------------------------------------

- Core function:

  Uses the compressed phenotype BED files generated in Step 3.1 together with genotype data to perform conventional cis-eQTL mapping across tissues using the linear mixed model (LMM) implemented in omiga.

------------------------------------------------------------------------
[Step 3] 3.3 Prepare for cseQTL.R
------------------------------------------------------------------------

- Core function:

  Loads the multidimensional bMIND RData containing deconvolution-derived cell-type-specific expression and extracts expression matrices for individual cell types. It then calculates cumulative exon lengths from the GTF annotation, performs TMM normalization and inverse normal transformation, and constructs cell-type-specific phenotype BED files for cseQTL analysis.

------------------------------------------------------------------------
[Step 4] 3.4 cseQTL mapping.sh
------------------------------------------------------------------------

- Core function:

  Performs nested traversal of tissue- and cell-type-specific subdirectories, matches the corresponding genotype data, and uses omiga to perform cis cell-type-specific eQTL mapping with a linear mixed model.

------------------------------------------------------------------------
[Step 5] 3.5 Prepare for cell-type ieQTL.R
------------------------------------------------------------------------

- Core function:

  Reads cell-type proportion matrices predicted by DWLS and removes low-quality samples and cell types with more than 80% zero proportions. Outliers are filtered using the mean ± 3 standard deviations (3 × SD) criterion. The script then generates the phenotype BED files and cell-proportion interaction files required for ieQTL analysis.

------------------------------------------------------------------------
[Step 6] 3.6 ieQTL mapping.sh
------------------------------------------------------------------------

- Core function:

  Traverses the cell-type interaction subdirectories for each tissue and supplies the corresponding cell-type proportions through the `--interaction` parameter to run the `cis_interaction` mode in omiga.


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
