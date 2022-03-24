#!/usr/bin/Rscript
#SBATCH --mem 128G
#SBATCH -J fish_plot

library(magrittr)
library(readr)
library(GenomicRanges)

source("fishPlot.R")

peakFiles <- list.files("data/open_devel_strict_refinepeak/", full.names = T)
names(peakFiles) <- c("dome", "epi75", "hpf12", "longPec", "prim5", "shield")

peakGranges <- lapply(peakFiles, function(x){
  x %>% read_tsv(col_names = c("seqnames", "start", "end", "name", "score")) %>%
    makeGRangesFromDataFrame()
})

library(TxDb.Drerio.UCSC.danRer10.ensGene)

fishPlotInputs <- list(
  dome = list(
    peakFile = "data/ATAC_peaks_timecourse/ATAC_dome_rep1-rep2.IDR0.1.narrowPeak.gz",
    assayFiles =  c(H3K4me3 = "/mnt/orca/damir/DANIO-CODE/DataFreeze_02/ChIP-seq/pipeline_output/DCD006178BS/out/align/pooled_rep/ChIP-seq_Skarmeta_Lab_H3K4me3_0002AS.DCD000647SQ.USERpanosfirbas.R1.filt.deduplicated.nodup_pooled.tagAlign.gz",
                    H3K27ac = "/mnt/storage/Damir/histones_marks/Skarmeta_H3K27ac_Dome_DCD000638SQ/align/rep1/ChIP-seq_Skarmeta_Lab_H3K27ac_0002AS.DCD000638SQ.USERpanosfirbas.R1.filt.deduplicated.nodup.tagAlign.gz",
                    H3K4me1 = "/mnt/orca/damir/DANIO-CODE/DataFreeze_02/ChIP-seq/pipeline_output/DCD006178BSk4me1/align/rep1/ChIP-seq_Skarmeta_Lab_H3K4me1_0002AS.DCD000657SQ.USERpanosfirbas.R1.filt.deduplicated.nodup.tagAlign.gz"),
    NucleoATACFiles = "stage_nucleoatac/Dome.nucleoatac_signal.smooth.bedgraph.gz",
    NucleoATACFilesPdre = "nucleoATAC_openChromatin_accross_development/dome.nucleoatac_signal.smooth.bedgraph.gz",
    ATACCutSitesFile = "/mnt/orca/damir/DanRer_ATAC-seq/pipeline_out/atac_dome/align/pooled_rep/ATAC_dome_rep1.R1.trim.PE2SE.nodup.tn5_pooled.tagAlign.gz",
    peakRanges = peakGranges$dome),

  Epi75 = list(
    peakFile = "/mnt/orca/damir/DanRer_ATAC-seq/pipeline_out/atac_80epi/peak/macs2/idr/true_reps/rep1-rep2/ATAC_80epi_rep1-rep2.IDR0.1.narrowPeak.gz",
    assayFiles = c(H3K4me3 = "/mnt/orca/damir/DANIO-CODE/DataFreeze_02/ChIP-seq/pipeline_output/DCD006177BS/out/align/pooled_rep/ChIP-seq_Skarmeta_Lab_H3K4me3_0002AS.DCD000643SQ.USERpanosfirbas.R1.filt.deduplicated.nodup_pooled.tagAlign.gz",
                   H3K27ac = "/mnt/storage/Damir/temp_histones/H3K27ac/Skarmeta_H3K27ac_75epiboly_1e-2/align/pooled_rep/ChIP-seq_Skarmeta_Lab_H3K27ac_0002AS.DCD000633SQ.USERpanosfirbas.R1.filt.deduplicated.nodup_pooled.tagAlign.gz",
                   H3K4me1 = "/mnt/storage/Damir/histones_marks/Skarmeta_H3K4me1_75-epiboly_DCD000656SQ/align/rep1/ChIP-seq_Skarmeta_Lab_H3K4me1_0002AS.DCD000656SQ.USERpanosfirbas.R1.filt.deduplicated.nodup.tagAlign.gz"),
    NucleoATACFiles = "stage_nucleoatac/Epi75.nucleoatac_signal.smooth.bedgraph.gz",
    NucleoATACFilesPdre = "nucleoATAC_openChromatin_accross_development/Epi75.nucleoatac_signal.smooth.bedgraph.gz",
    ATACCutSitesFile = "/mnt/orca/damir/DanRer_ATAC-seq/pipeline_out/atac_80epi/align/pooled_rep/ATAC_80epi_rep1.R1.trim.PE2SE.nodup.tn5_pooled.tagAlign.gz",
    peakRanges = peakGranges$epi75),

  Hpf12 = list(
    peakFile = "/mnt/orca/damir/DanRer_ATAC-seq/pipeline_out/atac_8som/peak/macs2/idr/true_reps/rep1-rep2/ATAC_8som_rep1-rep2.IDR0.1.narrowPeak.gz",
    assayFiles = c(H3K4me3 = "/mnt/orca/damir/DANIO-CODE/DataFreeze_02/ChIP-seq/pipeline_output/DCD006174BS/out/align/pooled_rep/ChIP-seq_Skarmeta_Lab_H3K4me3_0002AS.DCD000641SQ.USERpanosfirbas.R1.filt.deduplicated.nodup_pooled.tagAlign.gz",
                   H3K27ac = "/mnt/storage/Damir/temp_histones/H3K27ac/Skarmeta_H3K27ac_59somites_1e-2/align/pooled_rep/ChIP-seq_Skarmeta_Lab_H3K27ac_0002AS.DCD000634SQ.USERpanosfirbas.R1.filt.deduplicated.nodup_pooled.tagAlign.gz",
                   H3K4me1 = "/mnt/orca/damir/DANIO-CODE/new_samples/caper_test/chip/8113ae36-9790-47e0-a02b-7091c8d4daff/call-bam2ta/shard-0/execution/Som-K4me1_S9_L001_R1_001.merged.nodup.tagAlign.gz"),
    NucleoATACFiles = "stage_nucleoatac/Som8.nucleoatac_signal.smooth.bedgraph.gz",
    NucleoATACFilesPdre = "nucleoATAC_openChromatin_accross_development/Hpf12.nucleoatac_signal.smooth.bedgraph.gz",
    ATACCutSitesFile = "/mnt/orca/damir/DanRer_ATAC-seq/pipeline_out/atac_8som/align/pooled_rep/ATAC_8som_rep1.R1.PE2SE.nodup.tn5_ATAC_8som_rep2.R1.PE2SE.nodup.tn5.tagAlign.gz",
    peakRanges = peakGranges$hpf12),

  prim5 = list(
    peakFile = "data/Skarmeta_Lab_rep1-rep2.IDR0.1.narrowPeak.gz",
    assayFiles = c(H3K4me3 = "/mnt/orca/damir/DanRer_ChIP-seq_H3K4me3/Mueller_Prim-5/out/align/pooled_rep/ChIP-seq_Mueller_lab_H3K4me3_0006AS.DCD001525SQ.USERdunjanik.R1.filt.deduplicated.nodup_pooled.tagAlign.gz",
                   H3K27ac = "/mnt/orca/damir/DANIO-CODE/DataFreeze_02/ChIP-seq/pipeline_output/DCD006175BSk27ac/align/rep2/ChIP-seq_Skarmeta_Lab_H3K27ac_0002AS.DCD000653SQ.USERpanosfirbas.R1.filt.deduplicated.nodup.tagAlign.gz",
                   H3K4me1 = "data/h3k4me1_tagAlign/ChIP-seq_Skarmeta_Lab_H3K4me1_0002AS.DCD000655SQ.USERpanosfirbas.R1.filt.deduplicated.nodup.tagAlign.gz"),
    NucleoATACFiles = "stage_nucleoatac/Prim5.nucleoatac_signal.smooth.bedgraph.gz",
    NucleoATACFilesPdre = "nucleoATAC_openChromatin_accross_development/prim5.nucleoatac_signal.smooth.bedgraph.gz",
    ATACCutSitesFile = "/mnt/biggley/home/damir/pipelines/atac_dnase_pipelines/Skarmeta_Lab/align/pooled_rep/ATAC-seq_Skarmeta_Lab_0001AS.DCD000394SQ.USERpanosfirbas.R1.trim.PE2SE.nodup.tn5_pooled.tagAlign.gz",
    peakRanges = peakGranges$prim5),

  longPec = list(
    peakFile = "/mnt/orca/damir/DanRer_ATAC-seq/pipeline_out/atac_48h/peak/macs2/idr/true_reps/rep1-rep2/ATAC_48h_rep1-rep2.IDR0.1.narrowPeak.gz",
    assayFiles = c(H3K4me3 = "/mnt/orca/damir/DANIO-CODE/DataFreeze_02/ChIP-seq/pipeline_output/DCD006173BS/out/align/rep1/ChIP-seq_Skarmeta_Lab_H3K4me3_0002AS.DCD000640SQ.USERpanosfirbas.R1.filt.deduplicated.nodup.tagAlign.gz",
                   H3K27ac = "/mnt/orca/damir/DANIO-CODE/DataFreeze_02/ChIP-seq/pipeline_output/DCD006173BSk27ac/align/rep1/ChIP-seq_Skarmeta_Lab_H3K27ac_0002AS.DCD000651SQ.USERpanosfirbas.R1.filt.deduplicated.nodup.tagAlign.gz",
                   H3K4me1 = "/mnt/orca/damir/DANIO-CODE/DataFreeze_02/ChIP-seq/pipeline_output/DCD006173BSk4me1/align/rep1/ChIP-seq_Skarmeta_Lab_H3K4me1_0002AS.DCD000650SQ.USERpanosfirbas.R1.filt.deduplicated.nodup.tagAlign.gz"),
    NucleoATACFiles = "stage_nucleoatac/LongPec.nucleoatac_signal.smooth.bedgraph.gz",
    NucleoATACFilesPdre = "nucleoATAC_openChromatin_accross_development/LongPec.nucleoatac_signal.smooth.bedgraph.gz",
    ATACCutSitesFile = "/mnt/orca/damir/DanRer_ATAC-seq/pipeline_out/atac_48h/align/pooled_rep/ATAC_48hpf_rep1.R1.PE2SE.nodup.tn5_ATAC_48hpf_rep2.R1.PE2SE.nodup.tn5.tagAlign.gz",
    peakRanges = peakGranges$longPec)

)

modelMatrices <- lapply(fishPlotInputs, fishPlot2)

saveRDS(modelMatrices, "model_matrices_v3.RDS")


