#!/bin/bash

#SBATCH --array=1-136
#SBATCH --time=0-1:00:00
#SBATCH --mem-per-cpu=2GB
FILE_NAME=$1
j=`cat $FILE_NAME | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1`
# j=`cat $FILE_NAME | head -n 1 | tail -n 1`
PREFIX=`basename $j | sed 's/.blast.out//g'`
dir=`dirname $j`
source /fs/cbcb-scratch/nidhi/envs/env_python3/bin/activate 
python3 /fs/cbcb-scratch/nidhi/protist_project/latest/protist_reads/outlier_in_BLAST_hits/score_blast.py -q ${dir}/${PREFIX}.fasta -b ${j} -max 500 -out ${dir}/${PREFIX}.outlier.txt -blast ${dir}/${PREFIX}.subset_blast.out