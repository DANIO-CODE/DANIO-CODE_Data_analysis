library(rtracklayer)
library(tidyverse)
library(magrittr)
library(GenomicRanges)
load("failsafe.RData")

atacPeakFiles <- list.files("data/ATAC_peaks_timecourse_danRer11/", full.names = TRUE)
names(atacPeakFiles) <- c("Epi30pec", "c256", "c32", "c64", "Dome", 
                          "Epi75pec", "Hpf12", "LongPec",  "Prim5", "Shield" )

orderedStages <- c("c32", "c64", "c256", "Dome", "Epi30pec", "Shield",
                   "Epi75pec", "Hpf12", "Prim5", "LongPec")

atacPeakFiles <- atacPeakFiles[orderedStages]

importUniqueNarrowPeak <- function(file){
  tmp <- import(file)
  tmp <- tmp[order(tmp$pValue, decreasing = TRUE)]
  unique(tmp)
}

allPeaks <- lapply(atacPeakFiles, importUniqueNarrowPeak)
allPeaks <- as(allPeaks, Class = "GRangesList")

openAcrossDevelopment <- IRanges::reduce(unlist(allPeaks))

# I will include only those open in at least 2 consecutive stages

reduceOverlappingNeighbours <- function(x, y){
  tmp <- findOverlaps(x, y)
  IRanges::reduce(c(x[queryHits(tmp)], y[subjectHits(tmp)]))
}

openAccrossDevelopmenFilteredList <- mapply(reduceOverlappingNeighbours, allPeaks[-length(allPeaks)], allPeaks[-1])
openAccrossDevelopmenFilteredList <- as(openAccrossDevelopmenFilteredList, Class = "GRangesList")
openAccrossDevelopmenFiltered <- IRanges::reduce(unlist(openAccrossDevelopmenFilteredList))
cPADRES <- openAccrossDevelopmenFiltered

# ChromHMM

chromHMMfiles <- list.files(path = "data/chromHMM_danRer11/", 
                            pattern = ".+_dense\\.bed", full.names = T)
chromHMMFiles <- chromHMMfiles[c(3, 2, 1, 5,4)]
names(chromHMMFiles) <- names(umapAllStages)

segs <- lapply(chromHMMFiles,rtracklayer::import) 
chromHMMSegmentsList <- lapply(segs, function(x) split(x, x$name))
chromHMMSegmentsList <- lapply(chromHMMSegmentsList, function(x) x[names(x)[c(2, 10, 9, 7, 8, 6, 5, 4, 3, 1)]])

# fTemp <- function(x, y){
#   x[ modelMatrices$prim5$peaks %over% chromHMMSegmentsList[[y]] ] <- y
#   x
# }

fWrap <-  function(peaks, segments){
  function(x, y){
    x[ peaks %over% segments[[y]] ] <- y
    x
  }
}

stagePadres <- allPeaks[c(4, 7, 8, 9, 10)]
names(stagePadres) <- names(chromHMMSegmentsList)

segsStage <- lapply(names(chromHMMSegmentsList), function(nm) {
  purrr::reduce(names(chromHMMSegmentsList[[nm]]), 
                fWrap(stagePadres[[nm]], chromHMMSegmentsList[[nm]]),
                .init = vector(mode = "character",
                               length = length(stagePadres[[nm]])))
  })
names(segsStage) <- names(chromHMMSegmentsList)

segsStagePdre <- lapply(names(chromHMMSegmentsList), 
                        function(nm) {
                          purrr::reduce(names(
                            chromHMMSegmentsList[[nm]]),
                            fWrap(cPADRES, 
                                  chromHMMSegmentsList[[nm]]), 
                            .init = vector(mode = "character",
                                           length = 
                                             length(cPADRES)))
                        })
names(segsStagePdre) <- names(chromHMMSegmentsList)

for (i in names(segsStage)){
  stagePadres[[i]]$name <- segsStage[[i]]
}

# copes and dopes

reduceOverlappingPeaks <- function(x, y) {
  IRanges::reduce(c(x[x %over% y], y[y %over% x]))
}


quiescent_dr11 <- lapply(stagePadres, function(x){
  x %>% as.data.frame() %>% dplyr::filter(name == "10_Quies") %>% 
    GRanges()
})

copes_dr11 <- purrr::reduce(quiescent_dr11, reduceOverlappingPeaks)

pdreSegmentsMat11 <- do.call(cbind, segsStagePdre)
dopeIndex_dr11 <- which(apply(pdreSegmentsMat11, 1,
                         function(x){
                           all(x == "10_Quies")
                         }))

dopes11 <- cPADRES[dopeIndex_dr11]

# cell type inference

scregSeg <- data.table::fread("scregseg-pi_segmentation.tsv")
scregSegGR <- scregSeg %>% makeGRangesFromDataFrame(
  seqnames.field = "V1", start.field = "V2", 
  end.field = "V3", keep.extra.columns = T)

screSegSplit <- lapply(split(scregSegGR, scregSegGR$V4), GenomicRanges::reduce)

state_annotations <- c("Not Assigned", "Epidermis.periderm", "Muscle I", 
                       "Not Assigned", "Epidermis.olfactory.pronephric", 
                       "Muscle.Tailbud", "Spinal cord", "Not Assigned", 
                       "Epidermis", "Not Assigned", "CNS II", "Not Assigned", 
                       "Not Assigned", "Neural crest", "Endothelium", 
                       "Not Assigned", "Blood", "Midbrain", 
                       "Pharyngeal mesoderm.Tailbud", "Tailbud.Spinal cord", 
                       "Not Assigned", "Not Assigned", "Differentiating neurons", 
                       "Muscle II", "Optic vesicle", "CNS I", "Not Assigned", 
                       "Not Assigned", "Not Assigned", "Not Assigned")
names(state_annotations) <- str_c("state_", 0:29) %>% str_sort()

cellType_activity <- lapply(screSegSplit, function(ct){
  subsetByOverlaps(stagePadres$prim5, ct)
})

names(cellType_activity) <- str_c(state_annotations[names(cellType_activity)], 
                                  "_", names(cellType_activity))

names(cellType_activity) <- str_replace(names(cellType_activity), " ", "_")
#names(cellType_activity) <- str_replace(names(cellType_activity), "/", "|")


yavorEnhancers <- rbind(readxl::read_xlsx("EnhancersFinal.xlsx"),
read_xlsx("EnhancersFinal.xlsx",
sheet = "SepandEnhFinal")) %>%
makeGRangesFromDataFrame(keep.extra.columns = T)




# prepare the tracks

paper_colors <- c("#A6CEE3", "#1F78B4",   "#33A02C",   "#B2DF8A",   "#E31A1C", 
                  "#FB9A99" ,  "#FF7F00",   "#6A3D9A",   "#CAB2D6" ,  "#A1A2A3")
names(paper_colors) <- c("1_TssA1", "2_TssA2", "3_TssFlank1", "4_TssFlank2",
                         "5_EnhA1", "6_EnhFlank", "7_EnhWk1", 
                         "8_Pois", "9_ReprPC", "10_Quies" )

stageExport <- lapply(stagePadres, function(x){
  
  tmp <- x
  mcols(tmp) <- DataFrame(name = tmp$name,
                          score = tmp$score,
                          itemRgb = paper_colors[tmp$name],
                          thick = ranges(tmp))
  tmp
})

lapply(names(stageExport), function(name){
  export(stageExport[[name]], 
         str_c("annotations_danRer11/", name, "_PADREs.bed"))
})


cPadreExport <- stageExport <- lapply(segsStagePdre, function(x){
  
  tmp <- cPADRES
  mcols(tmp) <- DataFrame(name = x,
                          score = 1000,
                          itemRgb = paper_colors[x],
                          thick = ranges(tmp))
  tmp
})

lapply(names(cPadreExport), function(name){
  export(cPadreExport[[name]], 
         str_c("annotations_danRer11/", name, "_cPADREs.bed"))
})

stageDopes <- lapply(names(stagePadres), function(name){
  subsetByOverlaps(dopes11, stagePadres[[name]]) %>%
  export(str_c("annotations_danRer11/", name, "_DOPEs.bed"))
  })

export(copes_dr11, "annotations_danRer11/copes_dr11.bed")
export(dopes11, "annotations_danRer11/dopes_all_dr11.bed")

cell_type_names <- names(cellType_activity)
cell_type_names_filtered <- cell_type_names[!str_detect(cell_type_names,
                                                        "Not_Assigned")]
cellType_colours <- c("202,178,214", "51,160,44", "188,128,189", "178,223,138",
                      "252,205,229", "106,61,154", "166,206,227", "190,186,218",
                      "251,154,153", "251,128,114", "217,217,217", "255,237,111",
                      "177,89,40", "204,235,197", "179,222,105", "253,180,98",
                      "31,120,180"
)

names(cellType_colours) <- cell_type_names_filtered
cellType_rgb <- rgb(t(sapply(cellType_colours, function(x) {
  as.numeric(unlist(str_split(x, ",")))
  })), maxColorValue = 255 
  )
names(cellType_rgb) <- cell_type_names_filtered

lapply(cell_type_names_filtered, function(name){
  tmp <- cellType_activity[[name]]
  mcols(tmp) <- DataFrame(name = name,
                          score = 1000,
                          itemRgb = cellType_rgb[name],
                          thick = ranges(tmp))
  export(tmp,
         str_c("annotations_danRer11/scregSeg_", name, ".bed"))
  
})
