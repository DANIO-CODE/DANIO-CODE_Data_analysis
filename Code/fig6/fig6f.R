packages <- c('rtracklayer','ggplot2','tidyr','RColorBrewer','gridExtra','GenomicFeatures', 'DESeq2', 'pheatmap', 'vsn', 'grid', 'ggpmisc', 'Gviz', 'ehmm', 'RColorBrewer', 'pheatmap', 'reshape2', 'ggnewscale', 'data.table', 'BSgenome.Drerio.UCSC.danRer10', 'BSgenome.Mmusculus.UCSC.mm10', 'TxDb.Drerio.UCSC.danRer10.refGene', 'org.Dr.eg.db')
for (pkg in packages) suppressMessages(suppressWarnings(library(pkg, character.only=T)))
source('Code/fig6/functions.R')

## parse args
args <- commandArgs(trailingOnly=T)

# load ensembl annotation
TxDb.Drerio.UCSC.danRer10.ensGene <- suppressMessages(suppressWarnings(makeTxDbFromUCSC(genome = "danRer10", tablename = "ensGene")))

# prepare tracks for genomebrowser view plot
prepare_genomebrowser_data <- function(gr, bed, bw, from=NULL, to=NULL) {
  suppressMessages({
    # prepare data for genome browser screenshot plot
    chromosome <- seqlevels(gr)[1]
    if (is.null(from)) from <- min(start(gr)) - 5e3
    if (is.null(to)) to   <- max(end(gr)) + 5e3
    gr <- GRanges(seqnames=chromosome, IRanges(start=from, end=to))

    ## genes
    txdb <- TxDb.Drerio.UCSC.danRer10.ensGene
    genes <- GeneRegionTrack(txdb, stacking='dense', name='ENSEMBL', transcriptAnnotation="symbol", background.title="white", fontcolor.title="black", col.axis="black", rotation.title=0, size=.05, fontsize.group=22)
    symbols <- unlist(mapIds(org.Dr.eg.db, gsub('\\..*', '', gene(genes)), column="SYMBOL", keytype='ENSEMBL', multiVals="first"))
    symbol(genes) <- symbols[sapply(strsplit(gene(genes), '\\.'), function(x) x[[1]])]

    ## axis
    axis <- GenomeAxisTrack(name=chromosome, showTitle=TRUE, background.title="white", fontcolor.title="black", col.axis="black", rotation.title=0)

    ## data
    region <- GRanges(seqnames = Rle(chromosome), IRanges(from, to))
    col.fill <- 'grey60'
    
    ## bed
    bedTracks <- sapply(names(bed), function(x) {
      AnnotationTrack(sort(bed[[x]]), start=from, end=to, chromosome=chromosome, id=x, width=10, stacking='dense', feature=x, name=x,
                      background.title="transparent", col='transparent', fill=color_palette_publication[1], col.border.title="transparent", background.panel="transparent", fontcolor.title="black", 
                      col.axis="black", col="transparent", rotation.title=0, size=.05)
    })
    
    ## bigwig
    if (any(sapply(bw, function(x) length(x)==0))) {
      res <- NULL
     } else {
      bwTracks <- sapply(names(bw), function(x) {
        DataTrack(bw[[x]], name=x, size=.2, type='polygon', baseline=0, fill.mountain=c(col.fill, col.fill),
                  background.title="transparent", col.border.title="transparent", background.panel="transparent", fontcolor.title="black", col.axis="black", rotation.title=0)
      }, simplify=F)
    
      tracks <- c(list(axis=axis, genes=genes), bwTracks, bedTracks)
      res <- list(chromosome=chromosome, from=from, to=to, tracks=tracks)
    }
    })
  return(res)
}

# load TAD projection coordinates for the Irx3/5 TAD and compute trackdata
proj <- read.table(args[1], sep='\t', skip=2)
direct <- GRanges(unlist(strsplit(do.call('paste', c(as.list(gsub('\\(|\\)| ', '', proj[,7])), sep=',')), ',')))
multi <- GRanges(unlist(strsplit(do.call('paste', c(as.list(gsub('\\(|\\)| ', '', proj[,8])), sep=',')), ',')))
open_chromatin_zf <- import.bed(args[2])
open_chromatin_mm_onZfCoords <- GRanges(read.table(args[3], sep='\t', skip=2)[,4])
gr_irx3 <- GRanges('chr7:35693000-35795000')
trackdata_irx3 <- prepare_genomebrowser_data(gr=gr_irx3, bed=c(direct=direct, multi=multi, open_chromatin_zf=open_chromatin_zf, open_chromatin_mm=open_chromatin_mm_onZfCoords),
                                             bw=NULL, from=start(gr_irx3), to=end(gr_irx3))

# write plot to file
w <- 25
h <- 3
pdffile <- 'Figures/fig6/fig6f.pdf'
cairo_pdf(pdffile, width=w, height=h, onefile = T)
plotTracks(trackdata_irx3$tracks, col.mountain='transparent', chromosome=trackdata_irx3$chromosome,
           from=trackdata_irx3$from, to=trackdata_irx3$to, genome="danRer10", window=1000, fontsize=14,
           cex=.7, cex.axis=.6, cex.title=.7, title.width=1.3)
dev.off()
