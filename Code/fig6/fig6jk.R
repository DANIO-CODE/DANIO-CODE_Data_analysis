packages <- c('rtracklayer','ggplot2','tidyr','RColorBrewer','gridExtra','GenomicFeatures', 'DESeq2', 'pheatmap', 'vsn', 'grid', 'ggpmisc', 'Gviz', 'ehmm', 'RColorBrewer', 'pheatmap', 'reshape2', 'ggnewscale', 'data.table', 'BSgenome.Drerio.UCSC.danRer10', 'BSgenome.Mmusculus.UCSC.mm10', 'TxDb.Drerio.UCSC.danRer10.refGene', 'org.Dr.eg.db')
for (pkg in packages) suppressMessages(suppressWarnings(library(pkg, character.only=T)))
source('Code/fig6/functions.R')

## parse arguments
args <- commandArgs(trailingOnly=T)

# load zf enhancer ensembles and TADs
ee <- import.bed(args[1])
tads <- import.bed(args[2])
tads$name <- sprintf('tad_%03d', seq_along(tads))

# specify paths to data files
files <- list(
  danRer10=args[3],
  mm10=args[4]
)
species <- names(files)
ref <- 'danRer10'
qry <- 'mm10'

# load chromosome sizes
sizes <- list(danRer10=ehmm:::readGenomeSize(args[5]),
              mm10=ehmm:::readGenomeSize(args[6]))

# read whole genome count data from 200bp bigwigs.
whole_genome_count_data <- mclapply(species, function(sp) {
  return(import.bw(files[[sp]]))
}, mc.cores=2)
names(whole_genome_count_data) <- species

# load projections
thresh <- .95 # threshold for choosing the multi projection coordinates instead of direct
gr_list <- mclapply(tads$name, function(tad) {
  proj_file <- sprintf('Data/fig6/projection/tad_ESC_1kb_bins/%s_1kb.proj', tad)
  tad_proj <- read.table(proj_file, sep='\t', skip=2)[,2:6]
  colnames(tad_proj) <- c('coords_ref', 'coords_direct', 'coords_multi', 'score_direct', 'score_multi')
  gr_ref <- GenomicRanges::resize(GRanges(tad_proj$coords_ref), 1000, 'center')
  gr_qry <- GenomicRanges::resize(GRanges(apply(tad_proj, 1, function(x) {
    if (x['score_direct'] > thresh | x['score_direct'] > x['score_multi']) {
      return(x['coords_direct'])
    } else {
      return(x['coords_multi'])
    }
  })), 1000, 'center')
  gr_ref$name <- tad
  gr_qry$name <- tad
  return(list(ref=gr_ref, qry=gr_qry))
}, mc.cores=20)
names(gr_list) <- tads$name
tads_ref <- suppressWarnings(do.call('c', unname(sapply(gr_list, function(x) x$ref))))
tads_qry <- suppressWarnings(do.call('c', unname(sapply(gr_list, function(x) x$qry))))

# aggregate K27ac in 1kb bins and quantile-normalize query to reference
tads_ref <- suppressWarnings(ehmm:::aggScore(tads_ref, whole_genome_count_data$danRer10, 'mean'))
tads_qry <- suppressWarnings(ehmm:::aggScore(tads_qry, whole_genome_count_data$mm10, 'mean'))
tads_qry$score_qn <- suppressWarnings(ehmm:::quantileNormalizeToReference(cm.reference=matrix(tads_ref$score, nrow=1), cm.query=matrix(tads_qry$score, nrow=1))$cm.query.normalized)

# classify bins (zf/mm/both/none) with different thresholds
thresh <- quantile(tads_ref$score, 0.8)
class <- factor(mclapply(seq_along(tads_ref), function(i) {
  if (tads_ref$score[i] >= thresh & tads_qry$score_qn[i] < thresh) return('zebrafish')
  else if (tads_ref$score[i] < thresh & tads_qry$score_qn[i] >= thresh) return('mouse')
  else if (tads_ref$score[i] >= thresh & tads_qry$score_qn[i] >= thresh) return('both')
  else return('none')
}, mc.cores=40), levels=c('zebrafish', 'mouse', 'both', 'none'))

label <- ifelse(overlapsAny(tads_ref, ee), 'ensemble', 'non-ensemble')
data <- rbind(cbind(data.frame(table(class[which(label=='ensemble')])), data.frame(label='ensemble')),
              cbind(data.frame(table(class[which(label=='non-ensemble')])), data.frame(label='non-ensemble')))
colnames(data) <- c('signal', 'n_bins', 'label')
data$label <- factor(data$label, levels=c('non-ensemble', 'ensemble'))

w <- 10
h <- 3
set_plot_dimensions(w,h)
pl <- ggplot(data, aes(label, n_bins, fill=signal)) +
  geom_bar(stat = "identity", position = "fill") +
  coord_flip() +
  labs(y='Fraction of Bins', x='') +
  scale_fill_Publication(cols=c(zebrafish='#87ceeb', mouse='#FF0000', both='#006400', none='#CCCCCC'), name='H3K27ac enrichment') +
  theme_Publication(base_size=20)
ggsave('Figures/fig6/fig6j.pdf', pl, width=w, height=h, device=cairo_pdf)
write.table(data, 'Figures/fig6/fig6j.csv', sep='\t', quote=F, row.names=F)

# allocate quantile normalized mouse score to zebrafish coordinates
tads_zf <- tads_mm <- tads_ref[!duplicated(tads_ref)]
tads_mm$score <- tads_qry$score_qn[!duplicated(tads_ref)] # remove the REF duplicate entries from the QRY values as they positionally match
seqlengths(tads_zf) <- seqlengths(tads_mm) <- seqlengths(BSgenome.Drerio.UCSC.danRer10)[seqlevels(tads_zf)]

# compute quantile logratios and amplitudes
tads_zf$tad <- tads_zf$name
df_tads_zf <- cbind(data.frame(tads_zf)[,c('seqnames','start','end','name','score','strand')],
                    data.frame(thickStart=start(tads_zf)), thickEnd=end(tads_zf))
# prepare quantile data
tads_zf$quantile <- ecdf(tads_zf$score)(tads_zf$score)
tads_mm$quantile <- ecdf(tads_mm$score)(tads_mm$score)
tads_zf$quantile_logratio <- log((tads_mm$quantile+1e-2) / (tads_zf$quantile+1e-2))
# cap quantile log ratios to [-1,1]
tads_zf$quantile_logratio[tads_zf$quantile_logratio > 1] <- 1
tads_zf$quantile_logratio[tads_zf$quantile_logratio < -1] <- -1
tads_zf$quantile_amplitude <- apply(data.frame(zf=tads_zf$quantile, mm=tads_mm$quantile), 1, max)
cols <- colorRampPalette(c('blue','blue','#07731b','red','red'))(100)
tads_zf$itemRgb <- cols[sapply((tads_zf$quantile_logratio + 1) * 50, function(x) max(x,1))] # transform [-1,1] to [1,100] for color association

# load zf genes
TxDb.Drerio.UCSC.danRer10.ensGene <- suppressMessages(suppressWarnings(makeTxDbFromUCSC(genome='danRer10', tablename='ensGene')))
genes <- GeneRegionTrack(TxDb.Drerio.UCSC.danRer10.ensGene, name='ensGene', transcriptAnnotation="symbol", background.title="white", fontcolor.title="black", col.axis="black", rotation.title=0, size=.05, fontsize.group=22)
symbols <- unlist(mapIds(org.Dr.eg.db, gsub('\\..*', '', gene(genes)), column="SYMBOL", keytype='ENSEMBL', multiVals="first"))
symbol(genes) <- symbols[sapply(strsplit(gene(genes), '\\.'), function(x) x[[1]])]

# load direct alignments
aln <- import.bed(args[7])
strand(aln) <- '*'

# prepare tracks for genomebrowser view plot
tad_i <- 'tad_066'
gr_tad <- tads_zf[tads_zf$tad == tad_i]

# prepare data for genome browser screenshot plot
chromosome <- as.character(unique(seqnames(gr_tad)))
from <- min(start(gr_tad)) - 5e3
to   <- max(end(gr_tad)) + 5e3
gr <- GRanges(seqnames=chromosome, IRanges(start=from, end=to))

## axis
axis <- GenomeAxisTrack(name=chromosome, showTitle=TRUE, background.title="white", fontcolor.title="black", col.axis="black", rotation.title=0)

## data
region <- GRanges(seqnames = Rle(chromosome), IRanges(from, to))
col.fill <- 'grey60'

## bed
## names(bedfiles) will be printed
bedTracks <- c(
  ensembl=AnnotationTrack(ee, start=from, end=to, chromosome=chromosome, id='EnhancerEnsembles', width=10, stacking='dense', feature='EnhancerEnsembles', name='EnhancerEnsembles', background.title="transparent", col='transparent',
                          fill='black', col.border.title="transparent", background.panel="transparent", fontcolor.title="black", col.axis="black", col="transparent", rotation.title=0, size=.05),
  signal_overlap=AnnotationTrack(sort(gr_tad), start=from, end=to, chromosome=chromosome, id='H3K27acOverlap', width=10, feature='H3K27acOverlap', name='H3K27acOverlap', background.title="transparent", col='transparent',
                                 fill=gr_tad$itemRgb, alpha.title=1, alpha=gr_tad$quantile_amplitude, col.border.title="transparent", fontcolor.title="black", col.axis="black", col="transparent", rotation.title=0, size=.05),
  direct_alignment=AnnotationTrack(sort(aln), start=from, end=to, chromosome=chromosome, id='DirectAlignment', width=10, stacking='dense', feature='DirectAlignment', name='DirectAlignment', background.title="transparent", col='transparent',
                                 fill='black', col.border.title="transparent", background.panel="transparent", fontcolor.title="black", col.axis="black", col="transparent", rotation.title=0, size=.05)
)

## bigwig
## names(bwfiles) will be printed
bw <- c(zebrafish=tads_zf[tads_zf$tad==tad_i], mouse=tads_mm[tads_zf$tad==tad_i])
bw <- sapply(bw, function(x) {
  mcols(x) <- mcols(x)[,c('name','score')]
  return(x)
}, simplify=F)
bwTracks <- suppressWarnings(sapply(names(bw), function(x) {
  DataTrack(bw[[x]], name=x, size=.2, type='polygon', baseline=0, fill.mountain=c(col.fill, col.fill), ylim=c(0,1.2),
            background.title="transparent", col.border.title="transparent", background.panel="transparent", fontcolor.title="black", col.axis="black", rotation.title=0)
}, simplify=F))

tracks <- c(list(axis=axis, genes=genes), bwTracks, bedTracks)
trackdata <- list(chromosome=chromosome, from=from, to=to, tracks=tracks)

w <- 30
h <- 6
pdffile <- 'Figures/fig6/fig6k.pdf'
cairo_pdf(pdffile, width=w, height=h)
plotTracks(trackdata$tracks, col.mountain='transparent', chromosome=trackdata$chromosome, from=trackdata$from, to=trackdata$to, genome='danRer10', window=1000, fontsize=22, cex=.7, cex.axis=.6, cex.title=.7, title.width=1.3)
dev.off()
