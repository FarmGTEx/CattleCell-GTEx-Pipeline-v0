#!/bin/bash
#--------------------------------------------------------------------------#
#              Edit Job specifications                                     #
#--------------------------------------------------------------------------#
#SBATCH -p normal                 # Name of the queue
#SBATCH -N 1                       # Number of nodes(DO NOT CHANGE)
#SBATCH -n 24                       # Number of CPU cores
#SBATCH --mem=1024000                # Memory in MiB(10 GiB = 10 * 1024 MiB)
#SBATCH --account cattle_gtexs     #project name
#SBATCH -J MASHR                # Name of the job
#SBATCH --output=slurm_%A.out   # STDOUT
#SBATCH --error=slurm_%A.err    # STDERR
#SBATCH -t 10:00:00              # Job max time - Format = MM or MM:SS or HH:MM:SS or DD-HH or DD-HH:MM

# =========================================================================
# 路径修改：由总控 Bash 脚本分发过来的全局环境变量接管，彻底脱离硬编码
# =========================================================================
dir_nominal="${NOMINAL_DIR}"
dir_output="${OUTPUT_DIR}"

### 1. output all top SNP-gene pairs information from permutation results to a .txt file
# (colnames: pheno_id,variant_id,chr,pos)
# file list of permutation results
rm -f -r ${dir_output}/nominal_combined_files.txt
nominal_combined_files=(`find ${dir_nominal} -name "*_new_LMM.cis_qtl_pairs.*.txt.gz"`)
for l in ${nominal_combined_files[*]}
do
{
    echo ${l} >> ${dir_output}/nominal_combined_files.txt
}
done

# 路径修改：动态调用指定目录下的 python 脚本
python3 ${SCRIPT_DIR}/combine_signif_pairs_tjy.py ${dir_output}/nominal_combined_files.txt nominal_pairs -o ${dir_output}
### output file: nominal_pairs.combined_signifpairs.txt.gz
rm -f ${dir_output}/nominal_combined_files.txt


################################################################################################

tis_names=($(find "${dir_nominal}" -name "*_new_LMM.cis_qtl.txt.gz" | \
             sed 's|.*/||' | \
             sed 's/_new_LMM.cis_qtl.txt.gz//' | \
             sort -u))

### 2. extract all nominal pairs from nominal results for each tissue
for ((tis_i=0; tis_i<${#tis_names[@]}; tis_i++))
do
{
    tissue=${tis_names[tis_i]}
    nominal_files2=(`ls ${dir_nominal}/${tissue}/${tissue}_new_LMM.cis_qtl_pairs.*.txt.gz`)
    rm -f ${dir_output}/${tissue}.nominal_files2.txt
    for l in ${nominal_files2[*]}
    do
    {
        echo ${l} >> ${dir_output}/${tissue}.nominal_files2.txt
    }
    done
    # extract_pairs
    # 路径修改：动态调用指定目录下的 python 脚本
    python3 ${SCRIPT_DIR}/extract_pairs_tjy.py ${dir_output}/${tissue}.nominal_files2.txt ${dir_output}/nominal_pairs.combined_signifpairs.txt.gz ${tissue}_nominal_pairs -o ${dir_output}
    #> output file: *_nominal_pairs.extracted_pairs.txt.gz
}
done
wait


################################################################################################

# 变量修改：如果上层没有分发，则使用默认的 1000000 
if [ -z "$SUBSET_SIZE" ]; then
    subset_size=1000000
else
    subset_size="$SUBSET_SIZE"
fi

### 3. prepare random SNP-gene pairs for MashR
nominal_pairs_files=(`ls ${dir_output}/*_nominal_pairs.extracted_pairs.txt.gz`)
rm -f ${dir_output}/nominal_pairs_files.txt
for l in ${nominal_pairs_files[*]}
do
{
    echo ${l} >> ${dir_output}/nominal_pairs_files.txt
}
done
# MashR format file (z-score)
# 路径修改：动态调用指定目录下的 python 脚本
python3 ${SCRIPT_DIR}/mashr_prepare_input.py ${dir_output}/nominal_pairs_files.txt nominal_pairs.${subset_size}_subset -o ${dir_output} --only_zscore --dropna --subset $subset_size --seed 9823
zcat ${dir_output}/nominal_pairs.${subset_size}_subset.MashR_input.txt.gz | sed -e 's/.nominal_pairs_zval//g' | gzip > ${dir_output}/nominal_pairs.${subset_size}_subset.temp.txt.gz
mv ${dir_output}/nominal_pairs.${subset_size}_subset.temp.txt.gz ${dir_output}/nominal_pairs.${subset_size}_subset.MashR_input.txt.gz

# =========================================================================
# 安全熔断修改：彻底删除原代码末尾无保护的 rm -rf /scratch/$USER 删库死穴
# =========================================================================
