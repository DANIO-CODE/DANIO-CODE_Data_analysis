pacman::p_load("BSgenome.Drerio.UCSC.danRer10",
               "GenomicRanges",
               "rtracklayer",
               "data.table",
               "stringr",
               "DESeq2",
               "RColorBrewer",
               "ggplot2",
               "pheatmap",
               "som",
               "beanplot",
               "tidyr",
               "ggsci",
               "dplyr",
               "devtools",
               "biomaRt",
               "GenomicFeatures",
               "ggparallel",
               "glue",
               "PromoterOntology",
               "genomation",
               "seqPattern",
               "heatmaps",
               "TxDb.Drerio.BioMart.ENSEMBLMARTENSEMBL.GRCz10",
               "TxDb.Drerio.UCSC.danRer10.ensGene",
               "org.Dr.eg.db")

danRer10 <- BSgenome.Drerio.UCSC.danRer10

#TxDb.danRer10.ENSEMBL <- makeTxDbFromBiomart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset='drerio_gene_ensembl', host = 'http://dec2017.archive.ensembl.org')
# TxDb.danRer10.ENSEMBL <- TxDb.Drerio.BioMart.ENSEMBLMARTENSEMBL.GRCz10
# snl <- seqlevels(danRer10)
# names(snl) <- seqlevels(TxDb.danRer10.ENSEMBL)
# renameSeqlevels(TxDb.danRer10.ENSEMBL, snl)

TxDb.danRer10.ENSEMBL <- TxDb.Drerio.UCSC.danRer10.ensGene

atacTagAlignFiles <- list.files(path = "data/atac_tagAlign_ln/", full.names = TRUE)

samples <- c("Epi-30pec", 
             "prim5_rep1", "prim5_rep2", 
             "c128_rep1", "c128_rep2", "c128_rep3", 
             "c256_rep1", "c256_rep2", 
             "c32_rep1", "c32_rep2", "c32_rep3", 
             "c64_rep1", "c64_rep2", 
             "LongPec_rep1", "LongPec_rep2",
             "epi75_rep1", "epi75_rep2",
             "hpf12_rep1", "hpf12_rep2", 
             "dome_rep1", "dome_rep2", 
             "shield_rep1", "shield_rep2")
names(atacTagAlignFiles) <- samples

danRer10SeqInfo <- seqinfo(danRer10)

importAtacCutSites <- function(file, seqinfo = danRer10SeqInfo, resize = 1){
  temp <- fread(file)
  tempRanges <- GRanges(seqnames = temp$V1,
                        ranges = IRanges(start= temp$V2, end = temp$V3),
                        strand = temp$V6,
                        score = temp$V5,
                        name = temp$V4,
                        seqinfo = seqinfo)
  tempRanges <- resize(tempRanges, width = resize, fix = "start", ignore.strand = FALSE)
  tempRanges
}

atacSites <- lapply(atacTagAlignFiles, importAtacCutSites)

atacPeakFiles <- list.files("data/ATAC_peaks_timecourse/", full.names = TRUE)
names(atacPeakFiles) <- c("c128", "c256", "c32", "c64", "LongPec", "Epi75pec",
                          "Hpf12", "Dome", "Shield", "Epi30pec", "Prim5")

orderedStages <- c("c32", "c64", "c128", "c256", "Dome", "Epi30pec", "Shield",
                   "Epi75pec", "Hpf12", "Prim5", "LongPec")

atacPeakFiles <- atacPeakFiles[orderedStages]

sampleStages <-  c("30pec_epi", "Prim_5", "Prim_5", "128_cell", "128_cell", "128_cell", "256_cell", "256_cell","32_cell", "32_cell", "32_cell",
                   "64_cell", "64_cell", "Long_pec", "Long_pec", "75pec_epi", "75pec_epi", "12hpf", "12hpf", "Dome", "Dome", "Shield", "Shield")

samplePeriod <- c("gast", "phar", "phar", "preMZT", "preMZT", "preMZT", "preMZT", "preMZT", "preMZT", "preMZT", "preMZT", 
                  "preMZT", "preMZT", "phar", "phar", "gast", "gast", "phar", "phar", "gast", "gast", "gast", "gast")

names(sampleStages) <- names(atacTagAlignFiles)

importUniqueNarrowPeak <- function(file){
  tmp <- import(file)
  tmp <- tmp[order(tmp$pValue, decreasing = TRUE)]
  unique(tmp)
}

allPeaks <- lapply(atacPeakFiles, importUniqueNarrowPeak)
allPeaks <- as(allPeaks, Class = "GRangesList")

openAcrossDevelopment <- reduce(unlist(allPeaks))

# I will include only those open in at least 2 consecutive stages

reduceOverlappingNeighbours <- function(x, y){
  tmp <- findOverlaps(x, y)
  reduce(c(x[queryHits(tmp)], y[subjectHits(tmp)]))
}

openAccrossDevelopmenFilteredList <- mapply(reduceOverlappingNeighbours, allPeaks[-length(allPeaks)], allPeaks[-1])
openAccrossDevelopmenFilteredList <- as(openAccrossDevelopmenFilteredList, Class = "GRangesList")
openAccrossDevelopmenFiltered <- reduce(unlist(openAccrossDevelopmenFilteredList))

distalAcrossDevelopment <- subsetByOverlaps(openAcrossDevelopment,
                                            promoters(TxDb.danRer10.ENSEMBL,
                                                      downstream = 1000,
                                                      upstream = 500),
                                            invert = TRUE)

distalDevelFilt <- subsetByOverlaps(openAccrossDevelopmenFiltered,
                                    promoters(TxDb.danRer10.ENSEMBL,
                                              downstream = 1000,
                                              upstream = 500),
                                    invert = TRUE)
wrapper <- function(regions){
  function(signal){
    print("doing")
    countOverlaps(regions, signal)
  }
}

wrapperDistal <- wrapper(distalAcrossDevelopment)
wrapperFilter <- wrapper(distalDevelFilt)

accessibilityCounts <- lapply(atacSites, wrapperDistal)
accessibilityCountsFilt <- lapply(atacSites, wrapperFilter)

saveRDS(accessibilityCounts, "accessibility_across_development.RDS")
accessibilityCounts <- readRDS("accessibility_across_development.RDS")
saveRDS(accessibilityCountsFilt, "accessibility_development_filtered.RDS")
accessibilityCountsFilt <- readRDS("accessibility_development_filtered.RDS")


accessibilityCountMat <- matrix(unlist(accessibilityCountsFilt), ncol = 23, byrow = FALSE)
colnames(accessibilityCountMat) <- names(atacTagAlignFiles)

  

perMilFactors <- colSums(accessibilityCountMat) / 1e6
accessibilityCPMMat <- t(t(accessibilityCountMat) / perMilFactors)
colnames(accessibilityCPMMat) <- colnames(accessibilityCountMat)

pca <- prcomp(log(accessibilityCPMMat + .001))
pcaDf <- data.frame(samples = colnames(accessibilityCPMMat),
                    stage = c("30pec_epi", "Prim_5", "Prim_5", "128_cell", "128_cell", "128_cell", "256_cell", "256_cell","32_cell", "32_cell", "32_cell",
                              "64_cell", "64_cell", "Long_pec", "Long_pec", "75pec_epi", "75pec_epi", "12hpf", "12hpf", "Dome", "Dome", "Shield", "Shield"),
                    PC1 = pca$rotation[,"PC1"],
                    PC2 = pca$rotation[,"PC2"])
ggplot(pcaDf, aes(PC1, PC2, color = stage)) + geom_point(size = 3) + scale_color_brewer(palette = "Set3") + theme_bw() + xlab(paste0("PC1: ", format((pca$sdev[1] /sum(pca$sdev)) * 100, digits = 4), "%")) + ylab(paste0("PC2: ", format((pca$sdev[2] /sum(pca$sdev)) * 100, digits = 4), "%"))

sampleDistsCPM <- dist(t(accessibilityCPMMat))
library("RColorBrewer")
library(pheatmap)
sampleDistMatrixCPM <- as.matrix(sampleDistsCPM)
rownames(sampleDistMatrixCPM) <- colnames(accessibilityCPMMat)
colnames(sampleDistMatrixCPM) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrixCPM,
         clustering_distance_rows=sampleDistsCPM,
         clustering_distance_cols=sampleDistsCPM,
         col=colors)

colData <- pcaDf[,1:2]
colData$lab <- c("fm", "jls", "jls", "fm", "fm", "fm", "fm", "fm", "fm", "fm", "fm", "fm", "fm", "jls", "jls", "jls", "jls", "jls", "jls", "jls", "jls", "jls", "jls")
dds <- DESeqDataSetFromMatrix(accessibilityCountMat, colData = colData, design = ~stage)
dds <- DESeq(dds, test = "LRT", reduced = ~1)
vsd <- vst(dds)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$stage, vsd$lab, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "YlOrRd")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col = colors)

pcaData <- plotPCA(vsd, intgroup=c("stage", "lab"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData <- mutate(pcaData, stage = factor(stage, levels =c( 
                                            "32_cell","64_cell", "128_cell", "256_cell", 
                                            "Dome","30pec_epi",  "Shield","75pec_epi",  "12hpf","Prim_5","Long_pec")))
ggplot(pcaData, aes(PC1, PC2, color=stage, shape=lab)) +
  geom_point(size=6) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + scale_color_brewer(palette = "RdBu") +
  theme_bw() + ggtitle("Open chromatin during development PCA")

library(kohonen)

matForSom <- t(base::scale(t(assay(vsd)), center = FALSE))
openSom <- som(matForSom,
               grid = somgrid(4, 5, 'hexagonal'),
               mode = 'pbatch',
               cores = 4)

somDf <- cbind(data.frame(elID = 1:nrow(matForSom)),
               as.data.frame(matForSom))

somDf$som <- openSom$unit.classif
somDf <- gather(somDf, "sample", "value", -c("elID", "som"))
somDf$stage <- sampleStages[somDf$sample]
library(forcats)
somDf <- somDf %>% group_by(stage, elID, som) %>% summarize(avg = mean(value)) %>%
  ungroup() %>%
  mutate(stage = factor(stage, levels =c( 
    "32_cell","64_cell", "128_cell", "256_cell", 
    "Dome","30pec_epi",  "Shield","75pec_epi",  "12hpf","Prim_5","Long_pec")))
head(somDf)
somStats <- somDf %>% group_by(som, stage) %>% summarize(avg = median(avg))

ggplot(somDf, aes(stage, avg, fill = stage)) + geom_violin(alpha = .8) +
  facet_wrap(~ som, scales = "free") + geom_point(data = somStats) + geom_line(data = somStats, aes(group = 1)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_npg() + ggtitle("Distal elements") 

###############################################################################
# Here starts the code to normalize each period separatelly                   #
###############################################################################

gastIdx <- which(samplePeriod == "gast")
pharIdx <- which(samplePeriod == "phar")
preMZTIdx <- which(samplePeriod == "preMZT")

gastAcc <- accessibilityCountMat[,gastIdx]
pharAcc <- accessibilityCountMat[,pharIdx]
preMZTAcc <- accessibilityCountMat[,preMZTIdx]

splitColData <- split(colData, samplePeriod)

ddsGast <- DESeqDataSetFromMatrix(gastAcc, colData = splitColData$gast, design = ~stage)
ddsGast <- DESeq(ddsGast, test = "LRT", reduced = ~1)
vsdGast <- vst(ddsGast)

ddsPhar <- DESeqDataSetFromMatrix(pharAcc, colData = splitColData$phar, design = ~stage)
ddsPhar <- DESeq(ddsPhar, test = "LRT", reduced = ~1)
vsdPhar <- vst(ddsPhar)

ddsPreMZT <- DESeqDataSetFromMatrix(preMZTAcc, colData = splitColData$preMZT, design = ~stage)
ddsPreMZT <- DESeq(ddsPreMZT, test = "LRT", reduced = ~1)
vsdPreMZT <- vst(ddsPreMZT)

vsdMat <- cbind(assay(vsdPreMZT), assay(vsdGast), assay(vsdPhar))

sampleDists <- dist(t(vsdMat))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsdMat)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "YlOrRd")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col = colors)

pcaData <- prcomp(vsdMat)
percentVar <- round(100 * pcaData$sdev^2 / sum(pcaData$sdev^2))
pcaDf <- data.frame(samples = colnames(vsdMat),
                    stage =sampleStages[colnames(vsdMat)],
                    PC1 = pcaData$rotation[,"PC1"],
                    PC2 = pcaData$rotation[,"PC2"])
pcaDf <- mutate(pcaDf, stage = factor(stage, levels =c( 
  "32_cell","64_cell", "128_cell", "256_cell", 
  "Dome","30pec_epi",  "Shield","75pec_epi",  "12hpf","Prim_5","Long_pec")))
ggplot(pcaDf, aes(PC1, PC2, color=stage)) +
  geom_point(size=6) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))  + scale_color_brewer(palette = "RdBu") +
  theme_bw() + ggtitle("Open chromatin during development PCA")

library(kohonen)

matForSom <- t(base::scale(t(vsdMat), center = FALSE))
openSom <- som(matForSom,
               grid = somgrid(4, 5, 'hexagonal'),
               mode = 'pbatch',
               cores = 4)

somDf <- cbind(data.frame(elID = 1:nrow(matForSom)),
               as.data.frame(matForSom))

somDf$som <- openSom$unit.classif
somDf <- gather(somDf, "sample", "value", -c("elID", "som"))
somDf$stage <- sampleStages[somDf$sample]
library(forcats)
somDf <- somDf %>% group_by(stage, elID, som) %>% summarize(avg = mean(value)) %>%
  ungroup() %>%
  mutate(stage = factor(stage, levels =c( 
    "32_cell","64_cell", "128_cell", "256_cell", 
    "Dome","30pec_epi",  "Shield","75pec_epi",  "12hpf","Prim_5","Long_pec")))
head(somDf)
somStats <- somDf %>% group_by(som, stage) %>% summarize(avg = median(avg))

ggplot(somDf, aes(stage, avg, fill = stage)) + geom_violin(alpha = .8) +
  facet_wrap(~ som, scales = "free") + geom_point(data = somStats) + geom_line(data = somStats, aes(group = 1)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_npg() + ggtitle("Distal elements") 


######################
# try logFC
######################

expMat <- (colSums(accessibilityCountMat) / sum(width(distalDevelFilt))) * width(distalDevelFilt)
lFC <- log(expMat / (accessibilityCountMat + 1))

sampleDists <- dist(t(lFC))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$stage, vsd$lab, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "YlOrRd")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col = colors)


feFiles <- c(c32 = "/mnt/storage/Damir/ATAC-seq_resuts/32C/signal/macs2/pooled_rep/ATAC32C1_TAAGGCGA_L001_R1_001.trim.PE2SE.nodup.tn5_pooled.pf.fc.signal.bigwig",
             c64 = "/mnt/storage/Damir/ATAC-seq_resuts/64C/signal/macs2/pooled_rep/ATAC64C1_AGGCAGAA_L002_R1_001.trim.PE2SE.nodup.tn5_pooled.pf.fc.signal.bigwig",
             c128 = "/mnt/storage/Damir/ATAC-seq_resuts/128C/signal/macs2/pooled_rep/ATAC128C1_TCCTGAGC_L003_R1_001.trim.PE2SE.nodup.tn5_pooled.pf.fc.signal.bigwig",
             c256 = "/mnt/storage/Damir/ATAC-seq_resuts/256C/signal/macs2/pooled_rep/ATAC256CH_NoIndex_L001_R1_001.trim.PE2SE.nodup.tn5_pooled.pf.fc.signal.bigwig", 
             Dome = "/mnt/orca/damir/DanRer_ATAC-seq/pipeline_out/atac_dome/signal/macs2/pooled_rep/ATAC_dome_rep1.R1.trim.PE2SE.nodup.tn5_pooled.pf.fc.signal.bigwig",
             Epi30 = "~/pipelines/atac_dnase_pipelines/Mueller_Lab/signal/macs2/rep1/ATAC-seq_Mueller_lab_0006A_R1.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig",  
             Shield = "/mnt/orca/damir/DanRer_ATAC-seq/pipeline_out/atac_shield/signal/macs2/pooled_rep/ATAC_shield_rep1.R1.trim.PE2SE.nodup.tn5_pooled.pf.fc.signal.bigwig",
             Epi75 = "/mnt/orca/damir/DanRer_ATAC-seq/pipeline_out/atac_80epi/signal/macs2/pooled_rep/ATAC_80epi_rep1.R1.trim.PE2SE.nodup.tn5_pooled.pf.fc.signal.bigwig",  
             Hpf12 = "/mnt/orca/damir/DanRer_ATAC-seq/pipeline_out/atac_8som/signal/macs2/pooled_rep/ATAC_8som_rep1.R1.PE2SE.nodup.tn5_ATAC_8som_rep2.R1.PE2SE.nodup.tn5.pf.fc.signal.bigwig",
             Prim5 = "~/pipelines/atac_dnase_pipelines/Skarmeta_Lab/signal/macs2/pooled_rep/ATAC-seq_Skarmeta_Lab_0001AS.DCD000394SQ.USERpanosfirbas.R1.trim.PE2SE.nodup.tn5_pooled.pf.fc.signal.bigwig",
             LongPec ="/mnt/orca/damir/DanRer_ATAC-seq/pipeline_out/atac_48h/signal/macs2/pooled_rep/ATAC_48hpf_rep1.R1.PE2SE.nodup.tn5_ATAC_48hpf_rep2.R1.PE2SE.nodup.tn5.pf.fc.signal.bigwig")

feATACSignals <- lapply(feFiles, function(x) {
  foo <- rtracklayer::import(x)
  seqlevels(foo) <- seqlevels(danRer10)
  seqlengths(foo) <- seqlengths(danRer10)
  seqinfo(foo) <- seqinfo(danRer10)
  coverage(foo, weight = "score")
  })

seqlevels(pdreStric) <- seqlevels(danRer10)
seqlengths(pdreStric) <- seqlengths(danRer10)
seqinfo(pdreStric) <- seqinfo(danRer10)

pdreFC <- lapply(feATACSignals, function(x){
  binnedAverage(pdreStric, x, varname = "FC")
})

pdreFCmat <- sapply(pdreFC, function(x) x$FC)

pcaPdreFC <- prcomp(pdreFCmat)

pcaVarPec <- (pcaPdreFC$sdev)^2 / sum((pcaPdreFC$sdev)^2) * 100

library(stringr)

pcaPdreFC$rotation %>% as.data.frame() %>%
  tibble::rownames_to_column(var = "stage") %>%
  mutate(stage = factor(stage, levels = rev(rownames(pcaPdreFC$rotation)))) %>%
  ggplot(aes(PC1, PC2, colour = stage)) + geom_point(size = 5) + scale_colour_brewer(palette = "RdYlBu") +
  theme_bw() + theme(text = element_text(size = 24)) + xlab(sprintf("PC1  (%.2f %%)", pcaVarPec[1])) + 
  ylab(sprintf("PC1  (%.2f %%)", pcaVarPec[2]))

library(kohonen)

pdreFCmat
# enhancers <- import("enhancers_v1.bed")
# enhancers$name <- str_c("ENH", formatC(1:length(enhancers), width = 5, flag = '0'))
# ra
# overlapAndNormalize <- function(cutsites){
#   counts <- countOverlaps(enhancers, cutsites)
#   sizeFactor <- sum(counts) / 1000000
#   counts / sizeFactor
# }
# 
# overlapAndNormalizeLength <- function(cutsites){
#   counts <- countOverlaps(enhancers, cutsites)
#   counts <- counts / width(enhancers)
#   sizeFactor <- sum(counts) / 1000000
#   counts / sizeFactor
# }
# 
# enhExpression <- sapply(atacSites, overlapAndNormalize)
# enhExpressionMat <- sapply(atacSites, function(x) countOverlaps(enhancers, x))
# 
# colOrder <- c("dome_rep1", "dome_rep2", "shield_rep1", "shield_rep2", "epi75_rep1", "epi75_rep2", "hpf12_rep1", "hpf12_rep2", "prim5_rep1", "prim5_rep2", "LongPec_rep1", "LongPec_rep2")
# enhExpressionMat <- enhExpressionMat[,colOrder]
# row.names(enhExpressionMat) <- enhancers$name
# 
# set.seed(1234)
# plotEnhExpression <- enhExpression[sample(nrow(enhExpression), size = 30000, replace = FALSE),]
# plotEnhExpressionWin <- Winsorize(plotEnhExpression, probs = c(0, .95))
# #clust <- kmeans(plotEnhExpression, centers = 6)
# #heatmap(plotEnhExpressionWin[order(clust$cluster),], Colv = NA, Rowv = NA, scale = 'none')
# heatmap(plotEnhExpressionWin, Colv = NA, scale = 'none')
# #pdf("enhancer_hclust.pdf"); heatmap(plotEnhExpression, Colv = FALSE); dev.off()
# 
# accessEnhExpression <- enhExpression[rowSums(enhExpression < 10) < 6,]
# accessEnhExpressionWin <- Winsorize(accessEnhExpression, probs = c(0, .95))
# heatmap(accessEnhExpression, Colv = NA, scale = 'none')
# 
# # different normalization, peak length taken into consideration
# 
# enhExprLenNorm <- sapply(atacSites, overlapAndNormalizeLength)
# 
# enhExprLenNorm <- enhExprLenNorm[,colOrder]
# set.seed(1234)
# plotEnhExpression <- enhExprLenNorm[sample(nrow(enhExprLenNorm), size = 30000, replace = FALSE),]
# plotEnhExpressionWin <- Winsorize(plotEnhExpression, probs = c(0, .95))
# #clust <- kmeans(plotEnhExpression, centers = 6)
# #heatmap(plotEnhExpressionWin[order(clust$cluster),], Colv = NA, Rowv = NA, scale = 'none')
# heatmap(plotEnhExpressionWin, Colv = NA, scale = 'none')
# 
# p<-ggplot(df_out,aes(x=PC1,y=PC2,label=stage,color=stage ))
# p<-p+geom_point(aes(size = 5)) + geom_text(size=5) + scale_color_brewer(palette = "Paired") + ggtitle("Enhancer accessibility PCA") + xlab("PC1 (45%)") + ylab("PC2 (22%)")
# p
# 
# enhExprMat <- readRDS("enhancerAccessibility.RDS")
# enhAnno <- annotatePeak(enhancers, TxDb = TxDb.danRer10.ENSEMBL)
# plotAnnoBar(enhAnno)
# head(enhExprMat)
# colData <- data.frame(stage = rep(c("dome", "shield", "epi75", "hpf12", "prim5", "LongPec"), each = 2), row.names = colOrder)
# colData
# 
# enhDds <- DESeqDataSetFromMatrix(countData = enhExprMat,
#                                  colData = colData,
#                                  design = ~stage)
# enhDds
# enhDds <- DESeq(enhDds, test = "LRT", reduced = ~1)
# vsd <- vst(enhDds, blind=FALSE)
# sampleDists <- dist(t(assay(vsd)))
# sampleDistMatrix <- as.matrix(sampleDists)
# rownames(sampleDistMatrix) <- vsd$stage
# colnames(sampleDistMatrix) <- NULL
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# pheatmap(sampleDistMatrix,
#          clustering_distance_rows=sampleDists,
#          clustering_distance_cols=sampleDists,
#          col=colors)
# plotPCA(vsd, intgroup="stage")
# 
# ord <- order(results(enhDds)$padj, decreasing = FALSE)[1:40000]
# pheatmap(assay(vsd)[ord, ], cluster_cols = FALSE, show_rownames = FALSE)
# 
# enhSom <- som(assay(vsd), xdim = 4, ydim = 4, neigh = "gaussian", topol = "hex")
# enhSom1 <- som(assay(vsd), xdim = 5,ydim = 5, neigh = "gaussian", topol = "hex")
# 
# .myColorRamp <- function (vec, color.low="green", color.high="red", color.mid=NULL, alpha=1,
#                           value.low=min(vec), value.high=max(vec), value.mid=(value.low+value.high)/2, ...) {
#   vec.01 <- rep(NA, length(vec))
#   vec.01[vec <= value.low] <- 0
#   vec.01[vec >= value.high] <- 1
#   vec.01[vec>value.low & vec<=value.mid] <-
#     0.5*(vec[vec>value.low & vec<=value.mid]-value.low)/(value.mid-value.low)
#   vec.01[vec>value.mid & vec<value.high] <-
#     0.5+ 0.5*(vec[vec>value.mid & vec<value.high]-value.mid)/(value.high-value.mid)
#   cr <- colorRamp(c(color.low, color.mid, color.high),...)
#   return(apply (cr(vec.01)/255, 1, function(x) rgb(x[1], x[2], x[3], alpha)))
# }
# 
# .myColorMatrix <- function(color.mx, nrows, ncols, ...) {
#   top.row <- .myColorRamp(1:ncols,
#                           color.low=color.mx[1,1],
#                           color.mid=NULL,
#                           color.high=color.mx[1,2], ...)
#   bottom.row <- .myColorRamp(1:ncols,
#                              color.low=color.mx[2,1],
#                              color.mid=NULL,
#                              color.high=color.mx[2,2], ...)
#   sapply(1:ncols, function(i) .myColorRamp(1:nrows,
#                                            color.low=top.row[i],
#                                            color.mid=NULL,
#                                            color.high=bottom.row[i], ...
#                                            
#   )
#   )
#   
# }
# 
# .plot.clusters.beanplots <- function(value.matrix, cl, cl.method, dim.som.x, dim.som.y, ylim = c(0,2), las = 0, labels = colnames(value.matrix), titles = 'number', cex.axis = 1, cex.main = 2, cex.lab = 1, cols = c("red", "gold", "green", "blue")) {
#   
#   color.matrix.solid <- .myColorMatrix(matrix(cols, nrow=2), nrows=dim.som.y, ncols=dim.som.x)
#   color.matrix.solid = matrix(as.vector(color.matrix.solid), ncol = dim.som.x, byrow = F)
#   
#   for (j in (dim.som.y:1)-1) {
#     for (i in (1:dim.som.x)-1) {
#       to.plot <- data.frame(value.matrix[cl[,1] == i & cl[,2] ==j,])
#       if(titles == 'class'){
#         if(dim.som.y == 1){
#           title = i + 1
#         }else{
#           title = paste(i, j, sep=",")
#         }
#       }
#       if(titles == 'number'){
#         if(cl.method == "som"){
#           title = paste(i, "_", j, " (", nrow(to.plot), ")", sep = "")
#         }else if(cl.method == "kmeans"){
#           title = paste(i, " (", nrow(to.plot), ")", sep = "")
#         }
#       }
#       if (nrow(to.plot)>2000) to.plot <- to.plot[sample(nrow(to.plot), 2000),]
#       if(j==0){
#         beanplot(to.plot, ylim=ylim, bw=0.1, las=las, log = "", beanlinewd=0, cex.axis=cex.axis, col=c(color.matrix.solid[j+1,i+1],rgb(0,0,0,0),rgb(0,0,0,0),'white'), border=NA, main=title, yaxt = 'n', names = labels, cex.lab = cex.lab, cex.main = cex.main, col.main = color.matrix.solid[j+1,i+1])
#       }else{
#         beanplot(to.plot, ylim=ylim, bw=0.1, las=las, log = "", beanlinewd=0, cex.axis=cex.axis, col=c(color.matrix.solid[j+1,i+1],rgb(0,0,0,0),rgb(0,0,0,0),'white'), border=NA, main=title, yaxt = 'n', show.names = F, cex.lab = cex.lab, cex.main = cex.main, col.main = color.matrix.solid[j+1,i+1])
#         
#       }
#       
#     }
#   }
#   
# }
# 
# .plot.clusters.beanplots(assay(vsd), enhSom$visual, cl.method = "som", dim.som.x = 4, dim.som.y = 4,
#                          titles = 'class')
# 
# scaledMat <- t(scale(t(assay(vsd))))
# maps <- str_c(enhSom$visual[,1], enhSom$visual[,2], sep = "_")
# enhSomDf <- as.data.frame(assay(vsd))
# enhSomDf <- cbind(name = row.names(enhSomDf), class = maps, enhSomDf)
# enhSomDfTidy <- gather(enhSomDf, key = "sample", value = "log2_norm_expression", 3:ncol(enhSomDf))
# enhSomDfTidy$sample <- factor(as.character(enhSomDfTidy$sample), levels = colOrder)
# enhSomDfTidy$stage <- str_split(enhSomDfTidy$sample, pattern = "_", simplify = TRUE)[,1]
# ggplot(enhSomDfTidy, aes(sample, log2_norm_expression, fill = stage)) + geom_violin(alpha = .8) + 
#   facet_wrap(~ class, ) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_fill_brewer(palette = "Paired") + ggtitle("Enhancers dynamics")
# hpf <- c(dome = 4.3, shield = 6, epi75 = 8, hpf12 = 12, prim5 = 24, LongPec = 48)
# enhSomDfTidy$hpf <- hpf[enhSomDfTidy$stage]
# ggplot(enhSomDfTidy, aes(hpf, log2_norm_expression, group = name)) + geom_line(alpha = .75) + 
#   facet_wrap(~ class)
# 
# enhSomDf <- as.data.frame(scaledMat)
# enhSomDf <- cbind(name = row.names(enhSomDf), class = maps, enhSomDf)
# enhSomDfTidy <- gather(enhSomDf, key = "sample", value = "log2_norm_expression", 3:ncol(enhSomDf))
# enhSomDfTidy$sample <- factor(as.character(enhSomDfTidy$sample), levels = colOrder)
# enhSomDfTidy$stage <- str_split(enhSomDfTidy$sample, pattern = "_", simplify = TRUE)[,1]
# ggplot(enhSomDfTidy, aes(sample, log2_norm_expression, fill = stage)) + geom_violin(alpha = .8) + 
#   facet_wrap(~ class, ) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_fill_brewer(palette = "Paired") + ggtitle("Enhancers dynamics")
# hpf <- c(dome = 4.3, shield = 6, epi75 = 8, hpf12 = 12, prim5 = 24, LongPec = 48)
# 
# enhSomDf <- as.data.frame(scaledMat)
# maps1 <- str_c(enhSom1$visual[,1], enhSom1$visual[,2], sep = "_")
# enhSomDf <- cbind(name = row.names(enhSomDf), class = maps1, enhSomDf)
# enhSomDfTidy <- gather(enhSomDf, key = "sample", value = "log2_norm_expression", 3:ncol(enhSomDf))
# enhSomDfTidy$sample <- factor(as.character(enhSomDfTidy$sample), levels = colOrder)
# enhSomDfTidy$stage <- str_split(enhSomDfTidy$sample, pattern = "_", simplify = TRUE)[,1]
# ggplot(enhSomDfTidy, aes(sample, log2_norm_expression, fill = stage)) + geom_violin(alpha = .8) + 
#   facet_wrap(~ class, ) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_fill_npg() + ggtitle("Enhancers dynamics")
# hpf <- c(dome = 4.3, shield = 6, epi75 = 8, hpf12 = 12, prim5 = 24, LongPec = 48)
# 
# sampleStage <- rep(c("dome", "shield", "epi75", "hpf12", "prim5", "LongPec"), each = 2)
# names(sampleStage) <- colnames(assay(vsd))
# #sampleStage
# tp1 <- as.data.frame(assay(vsd))
# tp2 <- cbind(name = row.names(tp1), tp1)
# tp2 <- gather(tp2,key = "sample", value = "score", 2:ncol(tp2) )
# tp2$stage <- sampleStage[tp2$sample]
# #tp2
# tp3 <- tp2 %>% group_by(name, stage) %>% summarize(score = mean(score))
# tp3
# tp4 <- spread(tp3, stage, score)
# tp5 <- as.matrix(tp4[,2:ncol(tp4)])
# rownames(tp5) <- tp4$name
# tp5
# 
# scaledMat <- t(scale(t(tp5)))
# dim(scaledMat)
# head(scaledMat)
# stages <- c("dome", "shield", "epi75", "hpf12", "prim5", "LongPec")
# scaledMat <- scaledMat[,stages]
# scaledSom <- som(scaledMat, xdim = 5, ydim = 5, neigh = "gaussian", topol = "hex")
# 
# enhSomDf <- as.data.frame(scaledMat)
# maps1 <- str_c(scaledSom$visual[,1], scaledSom$visual[,2], sep = "_")
# enhSomDf <- cbind(name = row.names(enhSomDf), class = maps1, enhSomDf)
# enhSomDfTidy <- gather(enhSomDf, key = "stage", value = "log2_norm_expression", 3:ncol(enhSomDf))
# enhSomDfTidy$stage <- factor(as.character(enhSomDfTidy$stage), levels = c("dome", "shield", "epi75", "hpf12", "prim5", "LongPec"))
# #enhSomDfTidy$stage <- str_split(enhSomDfTidy$sample, pattern = "_", simplify = TRUE)[,1]
# 
# enhSomStats <- enhSomDfTidy %>% group_by(class, stage) %>% summarize(log2_norm_expression = median(log2_norm_expression))
# 
# #ggplot(enhSomDfTidy, aes(stage, log2_norm_expression, fill = stage)) + geom_violin(alpha = .8) + 
# #  facet_wrap(~ class, ) +
# #  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
# #  scale_fill_npg() + ggtitle("Enhancers dynamics")
# 
# ggplot(enhSomDfTidy, aes(stage, log2_norm_expression, fill = stage)) + geom_violin(alpha = .8) +
#   facet_wrap(~ class, ) + geom_point(data = enhSomStats) + geom_line(data = enhSomStats, aes(group = 1)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_fill_npg() + ggtitle("Enhancers dynamics")
# 
# transcripts <- fread("transcriptTableWStatusNCAGE.tsv", 
#                      col.names = c("seqnames", "start", "end", "ensembl_gene_id", "score", "strand", "gene_biotype", "external_gene_name"))
# transcripts <- makeGRangesFromDataFrame(transcripts, keep.extra.columns = T)
# transcripts
# promoters <- promoters(transcripts, downstream = 500, upstream = 500)
# promoters
# 
# load("/mnt/storage/dunja/danRer10CAGE/danRer10CAGEset.RData")
# 
# tcPrim5 <- makeGRangesFromDataFrame(tagClusters(danRer10CAGEset, samples = "prim5"), keep.extra.columns = TRUE)
# 
# CAGEPromoters <- get_CAGE_promoters(tcPrim5, promoters, protein_coding_only = FALSE)
# CpGIslands <- import("CpG.prediction.Christopher.UCSC.track.danRer10.bed")
# 
# 
# dominantTSSFlank <- promoters(resize(CAGEPromoters, width = 1, fix = "center"), upstream = 50, downstream = 150)
# CpGOverlaps <- findOverlaps(dominantTSSFlank, CpGIslands)
# length(unique(queryHits(CpGOverlaps)))
# CAGEPromoters$CGI <- FALSE
# CAGEPromoters$CGI[queryHits(CpGOverlaps)] <- TRUE
# CAGEPromoters
# CAGEPromoters$gene_biotype <- str_remove(CAGEPromoters$gene_biotype, pattern = "^(un)?annotated_")
# ggplot(as.data.frame(mcols(CAGEPromoters)), aes(gene_biotype, fill = CGI)) + geom_bar(position = "stack", stat = "identity")
# 
# 
# # ovo ne radi! trebas salloc i srun
# slurmSubmit <- function(command, partition = 'low', jobName = 'rscript', memory = '8G', cpus = 1){
#   line <- glue("salloc -p {partition} -J {jobName} --mem={memory} -c {cpus} srun")
#   print(str_c("Running: ", line, sep = ''))
#   system(line, wait = FALSE)
# }
# danRer10Fasta <- "/mnt/storage/Damir/data/kundaje_pipelines_genomes/wdl/danRer10.fa"
# 
# export(CAGEPromoters, "CAGE_promoters_Prim5.bed")
# 
# runNucleoATAC <- function(bam, bed, fasta, out = "nucleoatac_output", cores = 1){
#   line <- glue("nucleoatac run --bed {bed} --bam {bam} --fasta {fasta} --out {out} --cores {cores}")
#   print(line)
#   line
# }
# 
# pheatmap(scaledMat[sample.int(nrow(scaledMat), size = 40000, replace = FALSE),], cluster_cols = FALSE, show_rownames = FALSE)
# 
# prim5NucleoATACCommand <- runNucleoATAC(bam = "~/pipelines/atac_dnase_pipelines/Skarmeta_Lab/align/rep1/ATAC-seq_Skarmeta_Lab_0001AS.DCD000394SQ.USERpanosfirbas.R1.trim.PE2SE.nodup.bam",
#                                         bed = "CAGE_promoters_Prim5.bed",
#                                         fasta = danRer10Fasta,
#                                         out = "CAGE_Prim5_NucleoATAC",
#                                         cores = 4)
# slurmSubmit(prim5NucleoATACCommand, cpus = 4)
# 
# CAGEPromotersExprSorted <- CAGEPromoters[order(CAGEPromoters$tpm.dominant_ctss, decreasing = TRUE),]
# 
# importTagAlign <- function(file, seqinfo = danRer10SeqInfo){
#   temp <- fread(file)
#   tempRanges <- GRanges(seqnames = temp$V1,
#                         ranges = IRanges(start= temp$V2, end = temp$V3),
#                         strand = temp$V6,
#                         score = temp$V5,
#                         name = temp$V4,
#                         seqinfo = seqinfo)
#   #tempRanges <- resize(tempRanges, width = 1, fix = "start", ignore.strand = FALSE)
#   tempRanges
# }
# 
# prim5H3k4me3TagAlign <- importTagAlign("data/H3k4me3_tagAlign/ChIP-seq_Mueller_lab_H3K4me3_0006AS.DCD001525SQ.USERdunjanik.R1.filt.deduplicated.nodup.tagAlign.gz")
# nucleoATACSignal <- import("CAGE_Prim5_NucleoATAC.nucleoatac_signal.bedgraph.gz")
# sm <- ScoreMatrix(target = nucleoATACSignal,
#                   windows = CAGEPromotersExprSorted,
#                   strand.aware = TRUE,
#                   weight.col = "score")
# blues <- brewer.pal(n = 5,name = "Blues")
# heatMatrix(sm,
#            col = blues,
#            winsorize = c(0, 95),
#            xcoords = c(-500, 500))
# prim5AssayFiles <- c(H3K4me3 = "data/H3k4me3_tagAlign/ChIP-seq_Mueller_lab_H3K4me3_0006AS.DCD001525SQ.USERdunjanik.R1.filt.deduplicated.nodup.tagAlign.gz",
#                      H3k27ac = "data/h3k27ac_tagAlign_ln/ChIP-seq_Skarmeta_Lab_H3K27ac_0002AS.DCD000652SQ.USERpanosfirbas.R1.filt.deduplicated.nodup.tagAlign.gz",
#                      H3K4me1 = "data/h3k4me1_tagAlign/ChIP-seq_Skarmeta_Lab_H3K4me1_0002AS.DCD000655SQ.USERpanosfirbas.R1.filt.deduplicated.nodup.tagAlign.gz")
# 
# prim5NucleoATACFiles <- c(nucleoATACSignal = "CAGE_Prim5_NucleoATAC.nucleoatac_signal.bedgraph.gz",
#                           nucleoATACsmoothed = "CAGE_Prim5_NucleoATAC.nucleoatac_signal.smooth.bedgraph.gz",
#                           occupancy = "CAGE_Prim5_NucleoATAC.occ.bedgraph.gz",
#                           insertion_density = "CAGE_Prim5_NucleoATAC.ins.bedgraph.gz")
# 
# prim5ATACTa
# 
# prim5NucleoATACSignals <- lapply(prim5NucleoATACFiles, import)
# 
# smList <- ScoreMatrixList(targets = prim5NucleoATACSignals,
#                           windows = CAGEPromotersExprSorted,
#                           strand.aware = TRUE, 
#                           weight.col = "score")
# 
# multiHeatMatrix(smList, col = blues, xcoords = c(-500, 500), winsorize = c(5, 95))
# 
# prim5ATACCutSites <- importAtacCutSites("data/atac_tagAlign_ln/ATAC-seq_Skarmeta_Lab_0001AS.DCD000394SQ.USERpanosfirbas.R1.trim.PE2SE.nodup.tn5.tagAlign.gz", resize = 50)
# smCut <- ScoreMatrix(target = prim5ATACCutSites,
#                      windows = CAGEPromoters,
#                      strand.aware = TRUE,
#                      rpm = TRUE)
# heatMatrix(smCut,
#            col = blues,
#            winsorize = c(0, 95),
#            xcoords = c(-500, 500))
# 
# prim5Chips <- lapply(prim5AssayFiles, import, format = "BED")
# prim5Signals <- c(prim5Chips, ATAC = prim5ATACCutSites)
# 
# smList <- ScoreMatrixList(targets = prim5Signals,
#                           windows = CAGEPromoters,
#                           strand.aware = TRUE,
#                           rpm = TRUE)
# 
# smNuc <- ScoreMatrix(target = prim5NucleoATACSignals$nucleoATACsmoothed,
#                      windows = CAGEPromoters,
#                      strand.aware = TRUE,
#                      weight.col = "score")
# 
# smListComplete <- c(smList, Nucleosome = smNuc)
# multiHeatMatrix(smListComplete,
#                 col = blues,
#                 xcoords = c(-500, 500),
#                 winsorize = c(0, 95))
# 
# plotMeta(smListComplete,
#          winsorize = c(0, 95))
# 
# shft <- rep(50, length(CAGEPromoters))
# shft[as.vector(strand(CAGEPromoters) == '-')] <- -50
# 
# CAGEPromotersFirstNucleosome <- shift(promoters(resize(CAGEPromoters, width = 1, fix = "center"), downstream = 150, upstream = 0), shft)
# head(CAGEPromotersFirstNucleosome)
# 
# CAGEPromotersFirstNucleosomeSeq <- getSeq(danRer10, CAGEPromotersFirstNucleosome)
# plotPatternDensityMap(CAGEPromotersFirstNucleosomeSeq,
#                       patterns = c("SS", "WW"))
# plotPatternOccurrenceAverage(CAGEPromotersFirstNucleosomeSeq,
#                              patterns = c("SS", "WW"))
# 
# smListNonCoding <- lapply(smListComplete, function(x) x[CAGEPromoters$gene_biotype == "non_coding", ])
# smListNonCoding <- as(smListNonCoding, Class ="ScoreMatrixList")
# multiHeatMatrix(smListNonCoding,
#                 col = blues,
#                 xcoords = c(-500, 500),
#                 winsorize = c(0, 95))
# 
# cgiGroups <- list(CpGi = which(CAGEPromoters$CGI), nonCpGi = which(!CAGEPromoters$CGI))
# multiHeatMatrix(smListComplete,
#                 col = blues,
#                 xcoords = c(-500, 500),
#                 winsorize = c(0, 95),
#                 group = cgiGroups)
# 
# pdf("Prim5_CAGE_promoters.pdf"); multiHeatMatrix(smListComplete,
#                                                  col = blues,
#                                                  xcoords = c(-500, 500),
#                                                  winsorize = c(0, 95),
#                                                  group = cgiGroups); dev.off()
# 
# pdf("meat_plot.pdf"); lapply(smListComplete, plotMeta, winsorize = c(0, 95), xcoords = c(-500, 500)); dev.off()
