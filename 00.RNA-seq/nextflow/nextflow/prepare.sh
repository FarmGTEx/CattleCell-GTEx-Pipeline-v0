source "/home/huicongz/miniconda3/etc/profile.d/conda.sh"
conda activate nf-farmgtex
nextflow run ./script/prepare.nf -c nextflow.config -resume -w ./farmgtex/work