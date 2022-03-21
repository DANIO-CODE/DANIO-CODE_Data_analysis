#! /bin/bash

# This script merges bams from different stages to get 'virtual' whole embryo data
bam_dir=Data/fig6/bam
bw_dir=Data/fig6/bigwig

# H3K27me3 E10.5 heart, midbrain, forebrain, limb, hindbrain, embryonic facial prominence (2 replicates each)
mkdir -p $bam_dir/H3K27me3_E10.5
wget -nc -O $bam_dir/H3K27me3_E10.5/ENCFF156KVI.bam https://www.encodeproject.org/files/ENCFF156KVI/@@download/ENCFF156KVI.bam
wget -nc -O $bam_dir/H3K27me3_E10.5/ENCFF831KQT.bam https://www.encodeproject.org/files/ENCFF831KQT/@@download/ENCFF831KQT.bam
wget -nc -O $bam_dir/H3K27me3_E10.5/ENCFF742ORP.bam https://www.encodeproject.org/files/ENCFF742ORP/@@download/ENCFF742ORP.bam
wget -nc -O $bam_dir/H3K27me3_E10.5/ENCFF842YQE.bam https://www.encodeproject.org/files/ENCFF842YQE/@@download/ENCFF842YQE.bam
wget -nc -O $bam_dir/H3K27me3_E10.5/ENCFF854LJH.bam https://www.encodeproject.org/files/ENCFF854LJH/@@download/ENCFF854LJH.bam
wget -nc -O $bam_dir/H3K27me3_E10.5/ENCFF018UFV.bam https://www.encodeproject.org/files/ENCFF018UFV/@@download/ENCFF018UFV.bam
wget -nc -O $bam_dir/H3K27me3_E10.5/ENCFF867KIO.bam https://www.encodeproject.org/files/ENCFF867KIO/@@download/ENCFF867KIO.bam
wget -nc -O $bam_dir/H3K27me3_E10.5/ENCFF582AST.bam https://www.encodeproject.org/files/ENCFF582AST/@@download/ENCFF582AST.bam
wget -nc -O $bam_dir/H3K27me3_E10.5/ENCFF179RSO.bam https://www.encodeproject.org/files/ENCFF179RSO/@@download/ENCFF179RSO.bam
wget -nc -O $bam_dir/H3K27me3_E10.5/ENCFF553SUK.bam https://www.encodeproject.org/files/ENCFF553SUK/@@download/ENCFF553SUK.bam
wget -nc -O $bam_dir/H3K27me3_E10.5/ENCFF790CKS.bam https://www.encodeproject.org/files/ENCFF790CKS/@@download/ENCFF790CKS.bam
wget -nc -O $bam_dir/H3K27me3_E10.5/ENCFF597TOC.bam https://www.encodeproject.org/files/ENCFF597TOC/@@download/ENCFF597TOC.bam
samtools merge $bam_dir/H3K27me3_mm10_multiTissue_E10.5.merged.bam $bam_dir/H3K27me3_E10.5/*bam
samtools index $bam_dir/H3K27me3_mm10_multiTissue_E10.5.merged.bam

# DNase-seq E10.5 neural tube, embryonic facial prominence, forebrain, midbrain, hindbrain, limb, heart (2 replicates each)
mkdir -p $bam_dir/DNase=seq_E10.5
wget -nc -O $bam_dir/DNase-seq_E10.5/ENCFF316WII.bam https://www.encodeproject.org/files/ENCFF316WII/@@download/ENCFF316WII.bam
wget -nc -O $bam_dir/DNase-seq_E10.5/ENCFF193YPR.bam https://www.encodeproject.org/files/ENCFF193YPR/@@download/ENCFF193YPR.bam
wget -nc -O $bam_dir/DNase-seq_E10.5/ENCFF296CIM.bam https://www.encodeproject.org/files/ENCFF296CIM/@@download/ENCFF296CIM.bam
wget -nc -O $bam_dir/DNase-seq_E10.5/ENCFF554OYU.bam https://www.encodeproject.org/files/ENCFF554OYU/@@download/ENCFF554OYU.bam
wget -nc -O $bam_dir/DNase-seq_E10.5/ENCFF609MGC.bam https://www.encodeproject.org/files/ENCFF609MGC/@@download/ENCFF609MGC.bam
wget -nc -O $bam_dir/DNase-seq_E10.5/ENCFF444PVU.bam https://www.encodeproject.org/files/ENCFF444PVU/@@download/ENCFF444PVU.bam
wget -nc -O $bam_dir/DNase-seq_E10.5/ENCFF608FJD.bam https://www.encodeproject.org/files/ENCFF608FJD/@@download/ENCFF608FJD.bam
wget -nc -O $bam_dir/DNase-seq_E10.5/ENCFF011CAK.bam https://www.encodeproject.org/files/ENCFF011CAK/@@download/ENCFF011CAK.bam
wget -nc -O $bam_dir/DNase-seq_E10.5/ENCFF016JFD.bam https://www.encodeproject.org/files/ENCFF016JFD/@@download/ENCFF016JFD.bam
wget -nc -O $bam_dir/DNase-seq_E10.5/ENCFF106HYD.bam https://www.encodeproject.org/files/ENCFF106HYD/@@download/ENCFF106HYD.bam
wget -nc -O $bam_dir/DNase-seq_E10.5/ENCFF298LCK.bam https://www.encodeproject.org/files/ENCFF298LCK/@@download/ENCFF298LCK.bam
wget -nc -O $bam_dir/DNase-seq_E10.5/ENCFF573QHM.bam https://www.encodeproject.org/files/ENCFF573QHM/@@download/ENCFF573QHM.bam
wget -nc -O $bam_dir/DNase-seq_E10.5/ENCFF905VOZ.bam https://www.encodeproject.org/files/ENCFF905VOZ/@@download/ENCFF905VOZ.bam
wget -nc -O $bam_dir/DNase-seq_E10.5/ENCFF781ZLU.bam https://www.encodeproject.org/files/ENCFF781ZLU/@@download/ENCFF781ZLU.bam
samtools merge $bam_dir/DNase_mm10_multiTissue_E10.5.merged.bam $bam_dir/DNase-seq_E10.5/*bam
samtools index $bam_dir/DNase_mm10_multiTissue_E10.5.merged.bam

# # H3K27ac E10.5 midbrain, limb, embryonic facial prominence, hindbrain, forebrain, heart (2 replicates each)
# mkdir -p $bam_dir/H3K27ac_E10.5
# wget -nc -O $bam_dir/H3K27ac_E10.5/ENCFF067LVV.bam https://www.encodeproject.org/files/ENCFF067LVV/@@download/ENCFF067LVV.bam
# wget -nc -O $bam_dir/H3K27ac_E10.5/ENCFF859YQX.bam https://www.encodeproject.org/files/ENCFF859YQX/@@download/ENCFF859YQX.bam
# wget -nc -O $bam_dir/H3K27ac_E10.5/ENCFF720GQP.bam https://www.encodeproject.org/files/ENCFF720GQP/@@download/ENCFF720GQP.bam
# wget -nc -O $bam_dir/H3K27ac_E10.5/ENCFF024MPT.bam https://www.encodeproject.org/files/ENCFF024MPT/@@download/ENCFF024MPT.bam
# wget -nc -O $bam_dir/H3K27ac_E10.5/ENCFF140IYZ.bam https://www.encodeproject.org/files/ENCFF140IYZ/@@download/ENCFF140IYZ.bam
# wget -nc -O $bam_dir/H3K27ac_E10.5/ENCFF112HWN.bam https://www.encodeproject.org/files/ENCFF112HWN/@@download/ENCFF112HWN.bam
# wget -nc -O $bam_dir/H3K27ac_E10.5/ENCFF135OCE.bam https://www.encodeproject.org/files/ENCFF135OCE/@@download/ENCFF135OCE.bam
# wget -nc -O $bam_dir/H3K27ac_E10.5/ENCFF721KJH.bam https://www.encodeproject.org/files/ENCFF721KJH/@@download/ENCFF721KJH.bam
# wget -nc -O $bam_dir/H3K27ac_E10.5/ENCFF722INF.bam https://www.encodeproject.org/files/ENCFF722INF/@@download/ENCFF722INF.bam
# wget -nc -O $bam_dir/H3K27ac_E10.5/ENCFF328STI.bam https://www.encodeproject.org/files/ENCFF328STI/@@download/ENCFF328STI.bam
# wget -nc -O $bam_dir/H3K27ac_E10.5/ENCFF847FKM.bam https://www.encodeproject.org/files/ENCFF847FKM/@@download/ENCFF847FKM.bam
# wget -nc -O $bam_dir/H3K27ac_E10.5/ENCFF601EJB.bam https://www.encodeproject.org/files/ENCFF601EJB/@@download/ENCFF601EJB.bam
# samtools merge $bam_dir/H3K27ac_mm10_multiTissue_E10.5.merged.bam $bam_dir/H3K27ac_E10.5/*bam
# samtools index $bam_dir/H3K27ac_mm10_multiTissue_E10.5.merged.bam

# # H3K4me1
# mkdir -p $bam_dir/H3K4me1_E10.5
# wget -nc -P H3K4me1_E10.5 https://www.encodeproject.org/files/ENCFF655WRL/@@download/ENCFF655WRL.bam
# wget -nc -P H3K4me1_E10.5 https://www.encodeproject.org/files/ENCFF017VRQ/@@download/ENCFF017VRQ.bam
# wget -nc -P H3K4me1_E10.5 https://www.encodeproject.org/files/ENCFF379TJV/@@download/ENCFF379TJV.bam
# wget -nc -P H3K4me1_E10.5 https://www.encodeproject.org/files/ENCFF314WZL/@@download/ENCFF314WZL.bam
# wget -nc -P H3K4me1_E10.5 https://www.encodeproject.org/files/ENCFF053DHE/@@download/ENCFF053DHE.bam
# wget -nc -P H3K4me1_E10.5 https://www.encodeproject.org/files/ENCFF600FIR/@@download/ENCFF600FIR.bam
# wget -nc -P H3K4me1_E10.5 https://www.encodeproject.org/files/ENCFF154SPX/@@download/ENCFF154SPX.bam
# wget -nc -P H3K4me1_E10.5 https://www.encodeproject.org/files/ENCFF889LKY/@@download/ENCFF889LKY.bam
# wget -nc -P H3K4me1_E10.5 https://www.encodeproject.org/files/ENCFF398EKL/@@download/ENCFF398EKL.bam
# wget -nc -P H3K4me1_E10.5 https://www.encodeproject.org/files/ENCFF975XHS/@@download/ENCFF975XHS.bam
# wget -nc -P H3K4me1_E10.5 https://www.encodeproject.org/files/ENCFF152NAF/@@download/ENCFF152NAF.bam
# wget -nc -P H3K4me1_E10.5 https://www.encodeproject.org/files/ENCFF202ZDF/@@download/ENCFF202ZDF.bam
# wget -nc -P H3K4me1_E10.5 https://www.encodeproject.org/files/ENCFF626CWR/@@download/ENCFF626CWR.bam
# wget -nc -P H3K4me1_E10.5 https://www.encodeproject.org/files/ENCFF348TIN/@@download/ENCFF348TIN.bam
# wget -nc -P H3K4me1_E10.5 https://www.encodeproject.org/files/ENCFF226VNN/@@download/ENCFF226VNN.bam
# wget -nc -P H3K4me1_E10.5 https://www.encodeproject.org/files/ENCFF464SGQ/@@download/ENCFF464SGQ.bam
# wget -nc -P H3K4me1_E10.5 https://www.encodeproject.org/files/ENCFF019KCC/@@download/ENCFF019KCC.bam
# wget -nc -P H3K4me1_E10.5 https://www.encodeproject.org/files/ENCFF109KRV/@@download/ENCFF109KRV.bam
# wget -nc -P H3K4me1_E10.5 https://www.encodeproject.org/files/ENCFF064RLB/@@download/ENCFF064RLB.bam
# wget -nc -P H3K4me1_E10.5 https://www.encodeproject.org/files/ENCFF757MGS/@@download/ENCFF757MGS.bam
# wget -nc -P H3K4me1_E10.5 https://www.encodeproject.org/files/ENCFF416MSB/@@download/ENCFF416MSB.bam
# wget -nc -P H3K4me1_E10.5 https://www.encodeproject.org/files/ENCFF170VVI/@@download/ENCFF170VVI.bam
# wget -nc -P H3K4me1_E10.5 https://www.encodeproject.org/files/ENCFF034VML/@@download/ENCFF034VML.bam
# wget -nc -P H3K4me1_E10.5 https://www.encodeproject.org/files/ENCFF295GNC/@@download/ENCFF295GNC.bam
# samtools merge $bam_dir/H3K4me1_mm10_multiTissue_E10.5.merged.bam $bam_dir/H3K4me1_E10.5/*bam
# samtools index $bam_dir/H3K4me1_mm10_multiTissue_E10.5.merged.bam

# # H3K4me3
# mkdir -p $bam_dir/H3K4me3_E10.5
# wget -nc -P H3K4me3_E10.5 https://www.encodeproject.org/files/ENCFF626VEX/@@download/ENCFF626VEX.bam
# wget -nc -P H3K4me3_E10.5 https://www.encodeproject.org/files/ENCFF960QQM/@@download/ENCFF960QQM.bam
# wget -nc -P H3K4me3_E10.5 https://www.encodeproject.org/files/ENCFF031ZLM/@@download/ENCFF031ZLM.bam
# wget -nc -P H3K4me3_E10.5 https://www.encodeproject.org/files/ENCFF232GAW/@@download/ENCFF232GAW.bam
# wget -nc -P H3K4me3_E10.5 https://www.encodeproject.org/files/ENCFF078WCY/@@download/ENCFF078WCY.bam
# wget -nc -P H3K4me3_E10.5 https://www.encodeproject.org/files/ENCFF766EAE/@@download/ENCFF766EAE.bam
# wget -nc -P H3K4me3_E10.5 https://www.encodeproject.org/files/ENCFF249ZVI/@@download/ENCFF249ZVI.bam
# wget -nc -P H3K4me3_E10.5 https://www.encodeproject.org/files/ENCFF437PGG/@@download/ENCFF437PGG.bam
# wget -nc -P H3K4me3_E10.5 https://www.encodeproject.org/files/ENCFF553IOV/@@download/ENCFF553IOV.bam
# wget -nc -P H3K4me3_E10.5 https://www.encodeproject.org/files/ENCFF512HMV/@@download/ENCFF512HMV.bam
# wget -nc -P H3K4me3_E10.5 https://www.encodeproject.org/files/ENCFF974ABK/@@download/ENCFF974ABK.bam
# wget -nc -P H3K4me3_E10.5 https://www.encodeproject.org/files/ENCFF368UHS/@@download/ENCFF368UHS.bam
# wget -nc -P H3K4me3_E10.5 https://www.encodeproject.org/files/ENCFF293MTH/@@download/ENCFF293MTH.bam
# wget -nc -P H3K4me3_E10.5 https://www.encodeproject.org/files/ENCFF251WXG/@@download/ENCFF251WXG.bam
# wget -nc -P H3K4me3_E10.5 https://www.encodeproject.org/files/ENCFF815WJC/@@download/ENCFF815WJC.bam
# wget -nc -P H3K4me3_E10.5 https://www.encodeproject.org/files/ENCFF690BXJ/@@download/ENCFF690BXJ.bam
# wget -nc -P H3K4me3_E10.5 https://www.encodeproject.org/files/ENCFF748DPU/@@download/ENCFF748DPU.bam
# wget -nc -P H3K4me3_E10.5 https://www.encodeproject.org/files/ENCFF163RNG/@@download/ENCFF163RNG.bam
# wget -nc -P H3K4me3_E10.5 https://www.encodeproject.org/files/ENCFF059WDT/@@download/ENCFF059WDT.bam
# wget -nc -P H3K4me3_E10.5 https://www.encodeproject.org/files/ENCFF384HKX/@@download/ENCFF384HKX.bam
# wget -nc -P H3K4me3_E10.5 https://www.encodeproject.org/files/ENCFF224DTN/@@download/ENCFF224DTN.bam
# wget -nc -P H3K4me3_E10.5 https://www.encodeproject.org/files/ENCFF200DJC/@@download/ENCFF200DJC.bam
# wget -nc -P H3K4me3_E10.5 https://www.encodeproject.org/files/ENCFF838GDM/@@download/ENCFF838GDM.bam
# wget -nc -P H3K4me3_E10.5 https://www.encodeproject.org/files/ENCFF512WTX/@@download/ENCFF512WTX.bam
# samtools merge $bam_dir/H3K4me3_mm10_multiTissue_E10.5.merged.bam $bam_dir/H3K4me3_E10.5/*bam
# samtools index $bam_dir/H3K4me3_mm10_multiTissue_E10.5.merged.bam
