#!/bin/bash
#--------------------------------------------------------------------------#
#              Edit Job specifications                                     #
#--------------------------------------------------------------------------#
#SBATCH -p normal                 # Name of the queue
#SBATCH -N 1                       # Number of nodes(DO NOT CHANGE)
#SBATCH -n 24                       # Number of CPU cores
#SBATCH --mem=204800                # Memory in MiB(10 GiB = 10 * 1024 MiB)
#SBATCH --account cattle_gtexs     #project name
#SBATCH -J MASHR                # Name of the job
#SBATCH --output=slurm_%A.out   # STDOUT
#SBATCH --error=slurm_%A.err    # STDERR
#SBATCH -t 10:00:00              # Job max time - Format = MM or MM:SS or HH:MM:SS or DD-HH or DD-HH:MM
# Create a temporary directory for the job in local storage - DO NOT CHANGE #
TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPDIR
mkdir -p $TMPDIR
#=========================================================================#
#         Your job script                                                 #
#=========================================================================#

###############################################################################################

dir_nominal="/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/eQTL/"
dir_output="/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/MASHR_Celltype/Tissue_sharing/Bulk/strong_output"

### 1. output all top SNP-gene pairs information from permutation results to a .txt file
# (colnames: pheno_id,variant_id,chr,pos)
# file list of permutation results
rm -f ${dir_output}/permutation_files.txt
# permutation results
perm_files=(`find ${dir_nominal} -name "*_LMM.cis_eGene.txt"`)
for l in ${perm_files[*]}
do
{
    echo ${l} >> ${dir_output}/permutation_files.txt
}
done
/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/MASHR_Celltype/combine_signif_pairs_tjy.py ${dir_output}/permutation_files.txt strong_pairs -o ${dir_output}


################################################################################################

### 2. extract top pairs from nominal results for each tissue
nFile=`cat ${dir_output}/permutation_files.txt | wc -l`
# echo $nFile # 0..34
perm_files=(`find ${dir_nominal} -name "*_LMM.cis_eGene.txt"`)
NAMEs=(${perm_files[@]/*\//})
NAMEs=(${NAMEs[@]/_new_LMM.cis_eGene.txt/})
for tis_i in `seq 0 $[nFile-1]`
do
{
    name=${NAMEs[tis_i]}
    nominal_files=(`find ${dir_nominal} -name "${name}_LMM.cis_qtl_pairs.*.txt.gz"`)
    rm -f ${dir_output}/${name}.nominal_files.txt
    for l in ${nominal_files[*]}
    do
    {
        echo ${l} >> ${dir_output}/${name}.nominal_files.txt
    }
    done
    # extract_pairs
    /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/MASHR_Celltype/extract_pairs_tjy.py ${dir_output}/${name}.nominal_files.txt ${dir_output}/strong_pairs.combined_signifpairs.txt.gz ${name} -o ${dir_output}
    ### output file: *.extracted_pairs.txt.gz
    # rm -f ${dir_output}/${name}.nominal_files.txt
} &
done
wait


################################################################################################

### 3. prepare strong SNP-gene pairs for MashR
strong_pairs_files=(`ls ${dir_output}/*.extracted_pairs.txt.gz`)
rm -f ${dir_output}/strong_pairs_files.txt
for l in ${strong_pairs_files[*]}
do
{
    echo ${l} >> ${dir_output}/strong_pairs_files.txt
}
done
# MashR format file (z-score)
/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/MASHR_Celltype/mashr_prepare_input.py ${dir_output}/strong_pairs_files.txt strong_pairs -o ${dir_output} --only_zscore 
zcat ${dir_output}/strong_pairs.MashR_input.txt.gz | sed -e 's/_zval//g' | gzip > ${dir_output}/strong_pairs.temp.txt.gz
rm -f ${dir_output}/strong_pairs.MashR_input.txt.gz
mv ${dir_output}/strong_pairs.temp.txt.gz ${dir_output}/strong_pairs.MashR_input.txt.gz

################################################################################################

#=========================================================================#
#         Cleanup  DO NOT REMOVE OR CHANGE                                #
#=========================================================================#
cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID
