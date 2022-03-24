#!/bin/bash
#SBATCH --mem=10G
#SBATCH -c 8
#SBATCH -o slurm-%J.out
#SBATCH --array=0-4
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
#java -mx10G -jar /mnt/biggles/opt/ChromHMM/ChromHMM.jar BinarizeBam /mnt/biggles/opt/ChromHMM/CHROMSIZES/danRer10.txt . cell_mark_filetable_wATAC.txt binary_bam_out_wATAC
java -mx10G -jar /mnt/biggles/opt/ChromHMM/ChromHMM.jar Reorder \
     -m labelmappingfile_v4.txt -o stateordering_file_v4.txt \
     ChromHMM_out_v4/model_10.txt ChromHMM_out_v4_reordered

#java -mx10G -jar /mnt/biggles/opt/ChromHMM/ChromHMM.jar LearnModel -p 8 -printposterior -printstatebyline binary_bam_out_wATAC ChromHMM_out_v5 $state_num danRer10 
java -mx4000M -jar /mnt/biggles/opt/ChromHMM/ChromHMM.jar MakeSegmentation \
     ChromHMM_out_v4_reordered/model_10.txt binary_bam_out ChromHMM_out_v4_reordered

java -mx4000M -jar /mnt/biggles/opt/ChromHMM/ChromHMM.jar MakeBrowserFiles \
     -c colormappingfile_v5.txt \
     -m labelmappingfile_v4.txt ChromHMM_out_v4_reordered/${current_stage}_10_segments.bed \
     $current_stage ChromHMM_out_v4_reordered/${current_stage}_v4

set +x
