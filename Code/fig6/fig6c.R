packages <- c('rtracklayer','ggplot2','tidyr','RColorBrewer','gridExtra','GenomicFeatures', 'DESeq2', 'pheatmap', 'vsn', 'grid', 'ggpmisc', 'Gviz', 'ehmm', 'RColorBrewer', 'pheatmap', 'reshape2', 'ggnewscale', 'data.table', 'BSgenome.Drerio.UCSC.danRer10', 'BSgenome.Mmusculus.UCSC.mm10', 'TxDb.Drerio.UCSC.danRer10.refGene', 'org.Dr.eg.db')
for (pkg in packages) suppressMessages(suppressWarnings(library(pkg, character.only=T)))
source('Code/fig6/functions.R')

args <- commandArgs(trailingOnly=T)
if (length(args) != 3) stop('Usage: Rscript fig6c.R tad_projections.proj H3K27me3_zf.bw H3K27me3_mm_onZfCoords.bw')

# prepare tracks for genomebrowser view plot
prepare_genomebrowser_data <- function(gr, bw, from=NULL, to=NULL) {
  suppressMessages({
    # prepare data for genome browser screenshot plot
    chromosome <- seqlevels(gr)[1]
    if (is.null(from)) from <- min(start(gr)) - 5e3
    if (is.null(to)) to   <- max(end(gr)) + 5e3
    gr <- GRanges(seqnames=chromosome, IRanges(start=from, end=to))

    ## genes
    # removed stacking='dense'
    txdb <- TxDb.Drerio.UCSC.danRer10.refGene
    genes <- GeneRegionTrack(txdb, name='refGene', transcriptAnnotation="symbol", background.title="white", fontcolor.title="black", col.axis="black", rotation.title=0, size=.05, fontsize.group=22)
    symbols <- unlist(mapIds(org.Dr.eg.db, gsub('\\..*', '', gene(genes)), column="SYMBOL", keytype='ENTREZID', multiVals="first"))
    symbol(genes) <- symbols[gene(genes)]

    ## axis
    axis <- GenomeAxisTrack(name=chromosome, showTitle=TRUE, background.title="white", fontcolor.title="black", col.axis="black", rotation.title=0)

    ## data
    region <- GRanges(seqnames = Rle(chromosome), IRanges(from, to))
    col.fill <- 'grey60'
    
    ## bigwig
    if (any(sapply(bw, function(x) length(x)==0))) {
      res <- NULL
     } else {
      bwTracks <- sapply(names(bw), function(x) {
        DataTrack(bw[[x]], name=x, size=.2, type='polygon', baseline=0, fill.mountain=c(col.fill, col.fill),
                  background.title="transparent", col.border.title="transparent", background.panel="transparent", fontcolor.title="black", col.axis="black", rotation.title=0)
      }, simplify=F)
    
      tracks <- c(list(axis=axis, genes=genes), bwTracks)
      res <- list(chromosome=chromosome, from=from, to=to, tracks=tracks)
    }
    })
  return(res)
}

# load TAD projection coordinates for the Irx3/5 TAD and compute trackdata

proj_file <- args[1] # 'data/projection/tad_48h_1kb_bins/tad_0017_1kb.proj'
proj <- read.table(proj_file, sep='\t', skip=2)
bwfiles <- list(zf=args[2], # 'data/bigwig/H3K27me3_danRer10_prim5.cpm.200bp.1kb_bins.bw'
                mm=args[3]) # 'data/bigwig/H3K27me3_mm10_multiTissue_E10.5.cpm.200bp.1kb_bins.on_danRer10_multi.bw'
zf <- ehmm:::aggScore(GRanges(proj[,2]), import.bw(bwfiles$zf, which=GRanges(proj[,2])), 'mean')
mm <- ehmm:::aggScore(GRanges(proj[,2]), import.bw(bwfiles$mm, which=GRanges(proj[,2])), 'mean')
gr <- reduce(zf, min.gapwidth=1000)
gr_irx3 <- GRanges('chr7:35693000-35795000')
trackdata <- prepare_genomebrowser_data(gr, c(zf=zf, mm=mm), from=start(gr), to=end(gr))
trackdata_irx3 <- prepare_genomebrowser_data(gr_irx3, c(zf=zf, mm=mm), from=start(gr_irx3), to=end(gr_irx3))

# plot (without subdomains data as quantiles are computed for 1kb-binned data)
w <- 25
h <- 4
# write plot to file
pdffile <- 'Figures/fig6/fig6c.pdf'
cairo_pdf(pdffile, width=w, height=h, onefile = T)
plotTracks(trackdata$tracks, col.mountain='transparent', chromosome=trackdata$chromosome,
           from=trackdata$from, to=trackdata$to, genome="danRer10", window=1000, fontsize=22,
           cex=.7, cex.axis=.6, cex.title=.7, title.width=1.3)
plotTracks(trackdata_irx3$tracks, col.mountain='transparent', chromosome=trackdata_irx3$chromosome,
           from=trackdata_irx3$from, to=trackdata$to, genome="danRer10", window=1000, fontsize=22,
           cex=.7, cex.axis=.6, cex.title=.7, title.width=1.3)
dev.off()
