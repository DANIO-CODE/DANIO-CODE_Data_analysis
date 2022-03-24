# Does a fish plot needs:
# 1. peakFile (ATAC-seq peaks)
# 2. assayFiles (ChIP seq files with H3K4,e3, H3K27ac, H3K4me1)
# 3. ATACCutSites
# 4. NucleoATACFiles


# Import peaks and signal files

pacman::p_load('rtracklayer', 'magrittr', 'GenomicRanges', 'ggplot2', 'stringr', 'tidyr', 'ChIPseeker', 'TFBSTools', 'JASPAR2018', 'BSgenome.Drerio.UCSC.danRer10', 'BiocParallel', 'biomaRt', 'VennDiagram', 'boot', 'pbapply', 'data.table', 'genomation', 'purrr', 'RColorBrewer', 'dplyr', 'readr', 'reticulate', 'rlang', 'spp')
# I deleted tmaptools package
danRer10 <- BSgenome.Drerio.UCSC.danRer10
library(TxDb.Drerio.UCSC.danRer10.refGene)
#TxDB.danRer10 <- TxDb.Drerio.UCSC.danRer10.refGene

library(TxDb.Drerio.UCSC.danRer10.ensGene)
library(org.Dr.eg.db)
library("uwot")

TxDb.danRer10.ENSEMBL <- TxDb.Drerio.UCSC.danRer10.ensGene



importTrim <- function(x, genome = "danRer10", ...){
  x %>% rtracklayer::import(genome = genome, ...) %>% trim() %>% unique()
}


# atacFiles <- c(DCD000394SQ = "data/ATAC-seq_Skarmeta_Lab_0001AS.DCD000394SQ.USERpanosfirbas.R1.trim.PE2SE.nodup.tn5.pf.narrowPeak.gz",
#                DCD000395SQ = "data/ATAC-seq_Skarmeta_Lab_0001AS.DCD000395SQ.USERpanosfirbas.R1.PE2SE.nodup.tn5.pf.narrowPeak.gz",
#                pooled="data/ATAC-seq_Skarmeta_Lab_0001AS.DCD000394SQ.USERpanosfirbas.R1.trim.PE2SE.nodup.tn5_pooled.pf.narrowPeak.gz",
#                idr = "data/ATAC_peaks_timecourse/ATAC_dome_rep1-rep2.IDR0.1.narrowPeak.gz",
#                ppr = "data/Skarmeta_Lab_ppr.IDR0.1.narrowPeak.gz")
# 
# atacPeaks <- lapply(atacFiles, importTrim)
# 
# lapply(atacPeaks, summary)

danRer10SeqInfo <- seqinfo(danRer10)

importAtacCutSites <- function(file, seqinfo = NULL, resize = 1, shift = FALSE, step = 75){
  temp <- fread(file)
  tempRanges <- GRanges(seqnames = temp$V1,
                        ranges = IRanges(start= temp$V2, end = temp$V3),
                        strand = temp$V6,
                        score = temp$V5,
                        name = temp$V4,
                        seqinfo = seqinfo)
  tempRanges <- resize(tempRanges, width = resize, fix = "start", ignore.strand = FALSE)
  
  if (shift){
    mult <- c("+" = 1, "-" = -1)
    tempRanges <- IRanges::shift(tempRanges, shift = mult[as.vector(strand(tempRanges))] * step)
  }
  tempRanges
}

ScoreMatrixWrap <- function(target, window, width, weight = NA, bin.num = NA, fix = 'start', bin.op= 'mean'){
  
  if (fix == 'start'){
    w <- promoters(window, downstream =  width, upstream = width)
    sa = TRUE
  }
  
  args <- list(target = target, windows = w, strand.aware = sa)
  
  if (is.na(weight)){
    args <- append(args, list(rpm = TRUE))
  }else{
    args <- append(args, list(weight.col = weight))
  }
  
  if (!is.na(bin.num)){
    args <- append(args, list(bin.num = bin.num, bin.op = bin.op))
    sm <- invoke(ScoreMatrixBin, args)
  } else {
    sm <- invoke(ScoreMatrix, args)
  }
  
  sm
  
}

scoreMatrixToMatrix <- function(sm, nrow){
  tempMat <- matrix(nrow = nrow, ncol = ncol(sm))
  tempMat[as.numeric(rownames(sm)),] <- sm
  tempMat
}

# 
# system.time(test <- fread("data/atac_tagAlign_ln/ATAC-seq_Skarmeta_Lab_0001AS.DCD000394SQ.USERpanosfirbas.R1.trim.PE2SE.nodup.tn5.tagAlign.gz", col.names = c("seqnames", "start", "end", "name", "score", "strand")) %>% GRanges() %>% resize(width = 1, fix = "start", ignore.strand = FALSE))
# system.time(ATACCutSites <- importAtacCutSites("data/atac_tagAlign_ln/ATAC-seq_Skarmeta_Lab_0001AS.DCD000394SQ.USERpanosfirbas.R1.trim.PE2SE.nodup.tn5.tagAlign.gz", resize = 1, seqinfo = danRer10SeqInfo))


# here the files are links on our server. Can be replaced with DCC links to the files
# assayFiles <- c(H3K4me3 = "/mnt/orca/damir/DANIO-CODE/DataFreeze_02/ChIP-seq/pipeline_output/DCD006178BS/out/align/rep1/ChIP-seq_Skarmeta_Lab_H3K4me3_0002AS.DCD000647SQ.USERpanosfirbas.R1.filt.deduplicated.nodup.tagAlign.gz",
#                      H3K27ac = "/mnt/storage/Damir/histones_marks/Skarmeta_H3K27ac_Dome_DCD000638SQ/align/rep1/ChIP-seq_Skarmeta_Lab_H3K27ac_0002AS.DCD000638SQ.USERpanosfirbas.R1.filt.deduplicated.nodup.tagAlign.gz",
#                      H3K4me1 = "/mnt/orca/damir/DANIO-CODE/DataFreeze_02/ChIP-seq/pipeline_output/DCD006178BSk4me1/align/rep1/ChIP-seq_Skarmeta_Lab_H3K4me1_0002AS.DCD000657SQ.USERpanosfirbas.R1.filt.deduplicated.nodup.tagAlign.gz")
# NucleoATAC signal files doesn't require any links and changes. It is exactly where it is written
#NucleoATACFiles <- "Dome.nucleoatac_signal.smooth.bedgraph.gz"

# here starts the importirng
fishPlot <- function(x){
  print("Import NucleoATAC")
NucleoATACSignal <- rtracklayer::import(x$NucleoATACFiles)
print(NucleoATACSignal)
print("Import ATAC Cut Sites")

ATACCutSites <- importAtacCutSites(x$ATACCutSitesFile, resize = 1, seqinfo = danRer10SeqInfo)
print(ATACCutSites)

print("Importing ChIPS")
Chips <- lapply(x$assayFiles, importAtacCutSites, shift = TRUE)
print(Chips)
# mult <- c("+" = 1, "-" = -1)
# Chips <- lapply(Chips, function(x) shift(x, shift )
exprs <- str_c("Chips$", names(Chips))
names(exprs) <- names(Chips)
chipQuotes <- lapply(exprs, function(x) parse_expr(x))

#Signals <- c(chipQuotes, ATAC = quote(ATACCutSites), Nucleosome = quote(NucleoATACSignal))
#print(Signals)

#peaks <- atacPeaks$idr # [atacPeaks$pooled$signalValue >= 5]
if (x$peakFile != ""){
  print("Importing peaks")
  peaks <- unique(rtracklayer::import(x$peakFile))
  atacSummit <- IRanges::shift(peaks, shift=peaks$peak) %>% resize(width = 1, fix = 'start')
} else if (class(x$peakRanges) == "GRanges"){
  peaks <- x$peakRanges
  atacSummit <- resize(peaks, width = 1, fix = "center")
}else{
  error("Either a narrow peak file or GRanges file must be provided")
}
# Generally, only the window and bins will change in the future, which can be changed with modify with, so I can put now a list of args

smGlobalArgs <- list(H3K4me3 = list(target = Chips$H3K4me3,
                                    window = atacSummit,
                                    width = 750,
                                    weight = NA,
                                    bin.num = 13,
                                    bin.op = 'sum'),
                     H3K27ac = list(target = Chips$H3K27ac,
                                    window = atacSummit,
                                    width = 750,
                                    weight = NA,
                                    bin.num = 13,
                                    bin.op = 'sum'),
                     H3K4me1 = list(target = Chips$H3K4me1,
                                    window = atacSummit,
                                    width = 750,
                                    weight = NA,
                                    bin.num = 13,
                                    bin.op = 'sum'),
                     ATAC = list(target = ATACCutSites,
                                 window = atacSummit,
                                 width = 750,
                                 weight = NA,
                                 bin.num = 13,
                                 bin.op = 'sum'),
                     Nucleosome = list(target = NucleoATACSignal,
                                       window = atacSummit,
                                       width = 750,
                                       weight = 'score',
                                       bin.num = 13,
                                       bin.op = 'mean'))



binNum <- c( H3K4me3 = 13, H3K27ac = 13, H3K4me1 = 13, ATAC = 13, Nucleosome = 13)
featureNames <- c(unlist(mapply(function(x, y) str_c(rep(x, each = y), 1:y, sep = "_"), names(smGlobalArgs), binNum)))

modelNum <- length(atacSummit)
modelMatrix <- smGlobalArgs %>%
  purrr::map(function(x) invoke(ScoreMatrixWrap, x)) %T>%
  {matrixNames <- lapply(., rownames)} %>%
  lapply(scoreMatrixToMatrix, nrow = modelNum) %>%
  {do.call(cbind, .)}  

colnames(modelMatrix) <- featureNames
fallen <- is.na(modelMatrix)
fallenRows <-  apply(fallen, 1, sum) != 0
modelMatrix1 <- modelMatrix[!fallenRows,]

analysedPeaks <- peaks[!fallenRows]
analysedSummits <- atacSummit[!fallenRows]

scaled <- apply(modelMatrix1, 2, function(x) (x - min(x)) / (max(x) - min(x)))
reducedDims <- umap(scaled)
list(umap = reducedDims,
     peaks = analysedPeaks,
     modelmatrix = modelMatrix1,
     discarded = fallenRows)
}

fishPlot2 <- function(x){
  print("Import NucleoATAC")
  NucleoATACSignal <- rtracklayer::import(x$NucleoATACFiles)
  print(NucleoATACSignal)
  NucleoATACSignalPdre <- rtracklayer::import(x$NucleoATACFilesPdre)
  print("Import ATAC Cut Sites")
  
  ATACCutSites <- importAtacCutSites(x$ATACCutSitesFile, resize = 1, seqinfo = danRer10SeqInfo)
  print(ATACCutSites)
  
  print("Importing ChIPS")
  Chips <- lapply(x$assayFiles, importAtacCutSites, shift = TRUE)
  print(Chips)
  # mult <- c("+" = 1, "-" = -1)
  # Chips <- lapply(Chips, function(x) shift(x, shift )
  exprs <- str_c("Chips$", names(Chips))
  names(exprs) <- names(Chips)
  chipQuotes <- lapply(exprs, function(x) parse_expr(x))
  
  #Signals <- c(chipQuotes, ATAC = quote(ATACCutSites), Nucleosome = quote(NucleoATACSignal))
  #print(Signals)
  
  #peaks <- atacPeaks$idr # [atacPeaks$pooled$signalValue >= 5]
  
  print("Importing peaks")
  peaks <- unique(rev(rtracklayer::import(x$peakFile)))
  atacSummit <- IRanges::shift(peaks, shift=peaks$peak) %>% resize(width = 1, fix = 'start')
  
  peaksRanges <- x$peakRanges
  atacSummitRanges <- resize(peaksRanges, width = 1, fix = "center")
  
  # Generally, only the window and bins will change in the future, which can be changed with modify with, so I can put now a list of args
  
  smGlobalArgs <- list(H3K4me3 = list(target = Chips$H3K4me3,
                                      window = atacSummit,
                                      width = 750,
                                      weight = NA,
                                      bin.num = 13,
                                      bin.op = 'sum'),
                       H3K27ac = list(target = Chips$H3K27ac,
                                      window = atacSummit,
                                      width = 750,
                                      weight = NA,
                                      bin.num = 13,
                                      bin.op = 'sum'),
                       H3K4me1 = list(target = Chips$H3K4me1,
                                      window = atacSummit,
                                      width = 750,
                                      weight = NA,
                                      bin.num = 13,
                                      bin.op = 'sum'),
                       ATAC = list(target = ATACCutSites,
                                   window = atacSummit,
                                   width = 750,
                                   weight = NA,
                                   bin.num = 13,
                                   bin.op = 'sum'),
                       Nucleosome = list(target = NucleoATACSignal,
                                         window = atacSummit,
                                         width = 750,
                                         weight = 'score',
                                         bin.num = 13,
                                         bin.op = 'mean'))
  
  
  
  binNum <- c( H3K4me3 = 13, H3K27ac = 13, H3K4me1 = 13, ATAC = 13, Nucleosome = 13)
  featureNames <- c(unlist(mapply(function(x, y) str_c(rep(x, each = y), 1:y, sep = "_"), names(smGlobalArgs), binNum)))
  
  modelNum <- length(atacSummit)
  modelMatrix <- smGlobalArgs %>%
    purrr::map(function(x) invoke(ScoreMatrixWrap, x)) %T>%
    {matrixNames <- lapply(., rownames)} %>%
    lapply(scoreMatrixToMatrix, nrow = modelNum) %>%
    {do.call(cbind, .)}  
  
  colnames(modelMatrix) <- featureNames
  fallen <- is.na(modelMatrix)
  fallenRows <-  apply(fallen, 1, sum) != 0
  modelMatrix1 <- modelMatrix[!fallenRows,]
  
  analysedPeaks <- peaks[!fallenRows]
  analysedSummits <- atacSummit[!fallenRows]
  
  smGlobalArgs2 <- lapply(smGlobalArgs, function(x){
    x$window <- atacSummitRanges
    x
  })
  
  smGlobalArgs2$Nucleosome$target <- NucleoATACSignalPdre
  
  modelNum <- length(atacSummitRanges)
  modelMatrixPdre <- smGlobalArgs2 %>%
    purrr::map(function(x) invoke(ScoreMatrixWrap, x)) %T>%
    {matrixNames <- lapply(., rownames)} %>%
    lapply(scoreMatrixToMatrix, nrow = modelNum) %>%
    {do.call(cbind, .)}  
  
  colnames(modelMatrixPdre) <- featureNames
  fallenPdre <- is.na(modelMatrixPdre)
  fallenRowsPdre <-  apply(fallenPdre, 1, sum) != 0
  print(length(fallenRowsPdre))
  print(sum(fallenRowsPdre))
  modelMatrix2 <- modelMatrixPdre[!fallenRowsPdre,]
  
  print(dim(modelMatrix2))
  print(length(peaksRanges))
  analysedPdre <- peaksRanges[!fallenRowsPdre]
  analysedSummitsPdre <- atacSummitRanges[!fallenRowsPdre]
  
  #scaled <- apply(modelMatrix1, 2, function(x) (x - min(x)) / (max(x) - min(x)))
  #reducedDims <- umap(scaled)
  list(#umap = reducedDims,
       peaks = analysedPeaks,
       modelMatrixStage = modelMatrix1,
       modelMatrixPdre = modelMatrix2,
       discardedStage = fallenRows,
       peaksPdre = analysedPdre,
       discardedPdre = fallenRowsPdre)
}
