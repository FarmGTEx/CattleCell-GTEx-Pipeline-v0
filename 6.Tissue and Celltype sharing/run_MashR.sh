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
# 保留原作者的多线程计算加速优化配置
export USE_OPENMP=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# =========================================================================
# 路径修改：由总控 Bash 脚本根据不同模式 (bulk / celltype) 动态分发的环境变量接管
# =========================================================================
strong_file="${STRONG_FILE}"
random_file="${RANDOM_FILE}"

# 路径修改：动态调用指定脚本目录下的 R 脚本，并输出到参数指定的动态路径中
Rscript ${SCRIPT_DIR}/run_MashR.R ${strong_file} ${random_file} 1 ${MASHR_OUT_DIR}

# =========================================================================
# 安全熔断修改：彻底删除原代码末尾无保护的 rm -rf /scratch/$USER 删库死穴
# =========================================================================
