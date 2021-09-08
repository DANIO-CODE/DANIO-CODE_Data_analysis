## script for filtering for ATAC-seq peaks within TAD domains
suppressMessages(suppressWarnings(library(rtracklayer)))
args <- commandArgs(trailingOnly=T)
atac_peaks <- import.bed(args[1])
tads <- import.bed(args[2])
seqlevels(tads) <- paste0('chr', seqlevels(tads))
atac_peaks_inTAD <- atac_peaks[overlapsAny(atac_peaks, tads)]
export.bed(atac_peaks_inTAD, args[3])
