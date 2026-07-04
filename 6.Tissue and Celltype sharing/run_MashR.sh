#!/bin/bash
#--------------------------------------------------------------------------#
#              Edit Job specifications                                     #
#--------------------------------------------------------------------------#
#SBATCH -p normal                 # Name of the queue
#SBATCH -N 1                       # Number of nodes(DO NOT CHANGE)
#SBATCH -n 24                       # Number of CPU cores
#SBATCH --mem=1024000               # Memory in MiB(10 GiB = 10 * 1024 MiB)
#SBATCH --account cattle_gtexs     #project name
#SBATCH -J MASHR                # Name of the job
#SBATCH --output=slurm_%A.out   # STDOUT
#SBATCH --error=slurm_%A.err    # STDERR
#SBATCH -t 24:00:00              # Job max time - Format = MM or MM:SS or HH:MM:SS or DD-HH or DD-HH:MM

#=========================================================================#
#              Your job script                                            #
#=========================================================================#
export USE_OPENMP=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

strong_file="${STRONG_FILE}"
random_file="${RANDOM_FILE}"

Rscript ${SCRIPT_DIR}/run_MashR.R ${strong_file} ${random_file} 1 ${MASHR_OUT_DIR}
