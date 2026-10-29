# CattleCell-GTEx RNA-seq Pipeline: Test Dataset

This directory provides a small paired-end RNA-seq dataset for quickly testing the CattleCell-GTEx RNA-seq Nextflow pipeline.

The test FASTQ files were subsampled from the public RNA-seq dataset **SRR6253330** and are intended only to verify that the pipeline, software environment, reference files, and path settings are configured correctly. The test data should **not** be used for biological interpretation.

## Test directory structure

The test directory should contain:

```text
test/
├── README.md
├── samplelist.csv
└── fastq/
    └── SRR6253330/
        ├── SRR6253330_1.fastq.gz
        └── SRR6253330_2.fastq.gz
```

The two FASTQ files are paired-end reads:

```text
SRR6253330_1.fastq.gz    # Read 1
SRR6253330_2.fastq.gz    # Read 2
```

## 1. Set up the software environment

Before running the test, install the required software and dependencies.

The repository provides:

```text
nf-farmgtex.yml
```

which can be used to create the required environment.

For example:

```bash
conda env create -f nf-farmgtex.yml
conda activate nf-farmgtex
```

Alternatively, `mamba` can be used to create the environment.

## 2. Download reference files and tools

Run:

```bash
bash download.sh
```

before running the RNA-seq pipeline.

Please specify the appropriate local paths in `download.sh` before execution.

The pipeline uses the cattle reference genome **ARS-UCD1.2** and the Ensembl gene annotation:

```text
Bos_taurus.ARS-UCD1.2.108.chr.gtf
```

Other required annotation and reference files are described in the main pipeline README.

## 3. Prepare reference files

After downloading the required files, prepare the reference files and indices using:

```bash
nextflow run prepare.nf
```

Please specify the corresponding genome and working paths before running this step.

The STAR index and other files generated during this step are required before running the RNA-seq test.

If required by the current pipeline configuration, move the FASTA and BED files provided in the `file/` directory to the specified genome directory before running `prepare.nf`.

## 4. Prepare the test sample list

The test dataset contains one RNA-seq sample.

Create:

```text
test/samplelist.csv
```

with the following content:

```text
SRR6253330	Testis
```

Please note that although the file is named `samplelist.csv`, the two columns are **tab-delimited**.

The first column contains the sample ID and the second column contains the tissue name.

Do not include a header unless the main pipeline script has been modified to expect one.

## 5. Configure `rna.nf` for the test dataset

Before running the test, modify the paths in `rna.nf`.

For example, if the repository is located at:

```text
/path/to/CattleCell-GTEx-Pipeline-v0/00.RNA-seq/nextflow/nextflow
```

the test input path should be:

```text
/path/to/CattleCell-GTEx-Pipeline-v0/00.RNA-seq/nextflow/nextflow/test/fastq
```

and the test working directory can be:

```text
/path/to/CattleCell-GTEx-Pipeline-v0/00.RNA-seq/nextflow/nextflow/test
```

Accordingly, set:

```text
params.input = "/path/to/CattleCell-GTEx-Pipeline-v0/00.RNA-seq/nextflow/nextflow/test/fastq/"
params.workpath = "/path/to/CattleCell-GTEx-Pipeline-v0/00.RNA-seq/nextflow/nextflow/test"
```

The FASTQ suffix for this test dataset is:

```text
.fastq.gz
```

Therefore, make sure that the FASTQ suffix parameter in `rna.nf` is set accordingly.

The expected input structure is:

```text
${params.input}/${sampleID}/*.fastq.gz
```

which for this test corresponds to:

```text
test/fastq/SRR6253330/SRR6253330_1.fastq.gz
test/fastq/SRR6253330/SRR6253330_2.fastq.gz
```

## 6. Configure Nextflow execution

If the pipeline is being run on an HPC cluster, modify the corresponding cluster settings in:

```text
nextflow.config
```

according to the local job scheduler and computational resources.

Memory, CPUs, execution time, and other process-specific resources can also be adjusted in the corresponding modules if necessary.

For a local test, make sure that the Nextflow configuration is compatible with local execution and that sufficient memory and disk space are available.

## 7. Run the test

After the software environment, reference files, sample list, and paths have been prepared, run:

```bash
bash step1.sh
```

This executes the first RNA-seq processing workflow implemented in:

```text
rna.nf
```

The workflow performs the same major RNA-seq processing steps used for the full CattleCell-GTEx dataset, but using the small test dataset.

## 8. Test output

Output files will be generated under the directory specified by:

```text
params.workpath
```

For the configuration above, the results will therefore be written under:

```text
test/
```

The purpose of this test is to confirm that the workflow can be executed successfully. Because the FASTQ files contain only a small subset of the original sequencing reads, read counts, gene-expression estimates, genotype calls, and other downstream results should **not be interpreted biologically**.

Some downstream analyses may also produce sparse or empty results because of the limited sequencing depth of the test dataset.

## 9. VCF calling and `rna2.nf`

The minimal test dataset is primarily intended to test the first RNA-seq processing workflow (`rna.nf`).

The subsequent genotype-imputation workflow requires the CattleGTEx genotype reference panel and additional reference files. After genotype calling and imputation, the filtered VCF can be supplied to:

```text
rna2.nf
```

through:

```text
params.recalvcf
```

and the second-stage workflow can then be executed using:

```bash
bash nextflow2.sh
```

Users who only want to verify installation and basic execution of the RNA-seq pipeline do not need to run these additional steps.

## Notes

- This test dataset is provided only for pipeline testing.
- The FASTQ files are paired-end and must remain together.
- Use the same ARS-UCD1.2 reference and annotation files as described for the full pipeline.
- Paths in `download.sh`, `prepare.nf`, `rna.nf`, and `nextflow.config` may need to be adapted to the local computing environment.
- Successful completion of the test indicates that the main software dependencies and workflow configuration are functioning correctly, but does not validate biological results from the subsampled dataset.

For the complete workflow and detailed description of all analysis modules, please refer to the main README of the RNA-seq Nextflow pipeline.
