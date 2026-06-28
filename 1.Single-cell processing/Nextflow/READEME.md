# Single-Cell RNA-Seq Analysis Pipeline (Cell Ranger + Seurat)

This pipeline is designed for single-cell RNA sequencing (scRNA-seq) data analysis. It integrates **Cell Ranger** for read alignment and gene counting, followed by **Seurat** for quality control, doublet detection, data integration, and cell clustering.

---

# Table of Contents

* Quick Start
* Environment Setup (Nextflow 25.04.2)
* Input Dataset Preparation (Samplesheet)
* Container Environment
* Configuration and Resource Management
* Common Parameters

---

# 1. Quick Start

After preparing the environment (Java 17 or later), run the pipeline using:

```bash
./nextflow run main.nf \
    -profile docker \
    --samplesheet /target_dir/samplesheet.csv \
    --ref_fa "target_species_reference_genome/reference_genome.fna" \
    --ref_gtf "target_species_gtf/genomic.gtf" \
    --species_name "species_name" \
    --mito_threshold 10 \
    --pcs_use 30 \
    --outdir ../Single_cell_processing/ \
    -resume
```

---

# 2. Environment Setup (Nextflow 25.04.2)

This pipeline is developed and tested with **Nextflow 25.04.2**. Please prepare the environment as follows.

## 2.1 Java Requirement

Nextflow requires Java. We recommend **Java 17 or later**.

```bash
# Check Java version
java -version

# Install Java via SDKMAN! (recommended)
curl -s "https://get.sdkman.io" | bash
source "$HOME/.sdkman/bin/sdkman-init.sh"
sdk install java 17.0.x-tem
```

---

## 2.2 Install Nextflow 25.04.2

Download the corresponding Nextflow release and add it to your system PATH.

```bash
# Download Nextflow 25.04.2
wget -O nextflow https://github.com/nextflow-io/nextflow/releases/download/v25.04.2/nextflow-25.04.2-all

chmod +x nextflow
sudo mv nextflow /usr/local/bin/

# Verify installation
nextflow -version
```

---

## 2.3 Docker Requirement

This pipeline runs entirely inside Docker containers.

Please ensure that Docker is installed, the Docker service is running, and your user has permission to execute Docker commands.

```bash
# Start Docker service
sudo systemctl start docker

# Verify Docker access
docker ps
```

---

# 3. Input Dataset Preparation (Samplesheet)

The pipeline requires a **CSV** file named `samplesheet.csv` to specify sample information.

## 3.1 File Format

The samplesheet must contain the following two columns:

```csv
sampleid,fastq_dir
sample1,/media/Single_cell_processing/srr_test_project/data/fastqs
sample2,/media/Single_cell_processing/srr_test_project/data/fastqs
sample3,/media/Single_cell_processing/srr_test_project/data/fastqs
```

---

## 3.2 Sample Naming Recommendations

Each sample should follow these recommendations:

* **Unique ID:** Each sample must have a unique identifier.
* **Allowed characters:** Use only letters (`A-Z`, `a-z`), numbers (`0-9`), underscores (`_`), or hyphens (`-`).
* **Recommended naming convention:** `Group_Treatment_Replicate`, for example:

```
Control_WT_Rep1
```

---

# 4. Container Environment

The pipeline uses **two independent Docker images**.

## 4.1 `scrna_pipeline_env:latest` (Default Container)

**Purpose**

Runs the **Cell Ranger** alignment and counting steps.

**Environment**

Includes the official Cell Ranger software together with all required Linux dependencies.

---

## 4.2 `local_seurat_env:latest` (R Analysis Container)

**Purpose**

Runs the **SEURAT_PER_SAMPLE_QC** and **SEURAT_INTEGRATION** processes.

**Environment**

Provides an R 4.x environment with all required analysis packages, including:

* Seurat v4
* Harmony
* DoubletFinder
* and other required R packages

---

# 5. Configuration and Resource Management

## 5.1 CPU and Memory Settings (`conf/base.config`)

To adjust CPU or memory allocated to individual processes, modify the corresponding `withLabel` resource settings in:

```
conf/base.config
```

---

## 5.2 Default Parameters (`nextflow.config`)

Project-wide default parameters can be configured in the `params` section of:

```
nextflow.config
```

Example:

```groovy
params {
    outdir         = "$baseDir/scRNA_results"
    mito_threshold = 10
    pcs_use        = 30
}
```

---

# 6. Common Parameters

| Parameter          | Description                                                                                         |
| ------------------ | --------------------------------------------------------------------------------------------------- |
| `--samplesheet`    | Path to the samplesheet CSV file                                                                    |
| `--ref_fa`         | Path to the reference genome FASTA file                                                             |
| `--ref_gtf`        | Path to the reference genome GTF annotation file                                                    |
| `--species_name`   | Species name used for building the Cell Ranger reference                                            |
| `--mito_threshold` | Mitochondrial gene filtering threshold (default: **10%**)                                           |
| `--pcs_use`        | Number of principal components (PCs) used for downstream dimensionality reduction (default: **30**) |
| `--outdir`         | Root directory for output files                                                                     |

---

# Directory Structure

```
Single_cell_processing/
├── 00_Reference/             # Cell Ranger reference files
├── 01_Cellranger_Counts/     # Cell Ranger outputs
├── 02_Seurat_QC/             # Per-sample Seurat QC results
├── 03_Seurat_Integration/    # Integrated Seurat analysis results
├── reference/                # User-provided genome FASTA and GTF files
├── conf/                     # Nextflow configuration files
├── pipeline_info/            # Nextflow execution reports and logs
├── test/                     # Example test dataset
├── work/                     # Nextflow temporary working directory
├── main.nf                   # Main Nextflow workflow
├── nextflow.config           # Pipeline configuration
└── README.md
```

**Directory descriptions**

* **reference/**
  Store the downloaded reference genome FASTA and GTF annotation files. Their locations are specified through `--ref_fa` and `--ref_gtf`.

* **00_Reference/**
  Stores Cell Ranger reference resources generated or used by the pipeline.

* **01_Cellranger_Counts/**
  Contains Cell Ranger alignment and gene-counting results.

* **02_Seurat_QC/**
  Contains quality control, filtering, doublet detection, and per-sample Seurat objects.

* **03_Seurat_Integration/**
  Contains integrated Seurat objects, clustering results, UMAP visualizations, and downstream analyses.

* **pipeline_info/**
  Stores Nextflow execution reports, logs, timeline, DAG, and trace information.

* **work/**
  Temporary working directory automatically created by Nextflow. It stores intermediate files and caches generated during pipeline execution. This directory can become very large. Once the pipeline has successfully completed and the final results have been saved, the `work/` directory can be safely removed if intermediate files are no longer needed.

---

# Technical Support

If the pipeline finishes successfully but generated figures appear blank (white pages), please check whether the RDS files exist in:

```
03_Seurat_Integration/
```

The pipeline includes an automatic recovery mechanism. Even if figure generation fails, the core analysis results and Seurat objects will still be preserved.

