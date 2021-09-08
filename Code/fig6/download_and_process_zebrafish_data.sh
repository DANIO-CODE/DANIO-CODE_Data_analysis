#! /bin/bash

bam_dir=Data/fig6/bam
bw_dir=Data/fig6/bigwig
mkdir -p $bam_dir $bw_dir

# H3K27ac
k27ac_dome=https://danio-code.zfin.org/files/annotated_files/ChIP-seq/DCD006178BS/ChIP-seq_Skarmeta_Lab_H3K27ac_0002AS.DCD000638SQ.USERpanosfirbas.R1.filt.deduplicated.nodup.bam
wget -nc -O $bam_dir/H3K27ac_danRer10_dome.bam $k27ac_dome
samtools index $bam_dir/H3K27ac_danRer10_dome.bam

# H3K27me3
k27me3_prim5=https://danio-code.zfin.org/files/annotated_files/ChIP-seq/DCD006175BS/ChIP-seq_Skarmeta_Lab_H3K27me3_0002AS.DCD000632SQ.USERpanosfirbas.R1.filt.deduplicated.nodup.bam
wget -nc -O $bam_dir/H3K27me3_danRer10_prim5.bam $k27me3_prim5
samtools index $bam_dir/H3K27me3_danRer10_prim5.bam

# ATAC-seq
atac_prim5=https://danio-code.zfin.org/files/annotated_files/ATAC-seq/DCD007260BS/ATAC-seq_Skarmeta_Lab_0001AS.DCD003075SQ.danRer10.USERdanio-user.R1.merged.filt.bam
wget -nc -O $bam_dir/ATAC_danRer10_prim5.bam $atac_prim5
samtools index $bam_dir/ATAC_danRer10_prim5.bam
