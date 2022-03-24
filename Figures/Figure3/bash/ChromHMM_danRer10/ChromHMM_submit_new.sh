#!/bin/bash
#SBATCH --mem=10G
#SBATCH -c 8
#SBATCH -o slurm-%J.out
# Please state the number of array iteration at sbatch submission

stages=( "$@" )

# careful, this is hardcoded, depends on the order of stage input!!!!

#num_states=( 6 3 5 3 3 6 6 )
state_num=10
echo $state_num

current_stage=${stages[$SLURM_ARRAY_TASK_ID]}
echo $current_stage

cd $current_stage

set -x
java -mx10G -jar /mnt/biggles/opt/ChromHMM/ChromHMM.jar BinarizeBam /mnt/biggles/opt/ChromHMM/CHROMSIZES/danRer10.txt . cell_mark_filetable.txt binary_bam_out

java -mx10G -jar /mnt/biggles/opt/ChromHMM/ChromHMM.jar LearnModel -p 8 -printposterior -printstatebyline binary_bam_out ChromHMM_out_v4 $state_num danRer10 
set +x
