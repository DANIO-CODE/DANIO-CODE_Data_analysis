packages <- c('rtracklayer','ggplot2','tidyr','RColorBrewer','gridExtra','GenomicFeatures', 'DESeq2', 'pheatmap', 'vsn', 'grid', 'ggpmisc', 'Gviz', 'seqPattern', 'scales', 'ehmm', 'RColorBrewer', 'pheatmap', 'reshape2', 'ggnewscale', 'data.table', 'BSgenome.Drerio.UCSC.danRer10', 'BSgenome.Mmusculus.UCSC.mm10', 'TxDb.Drerio.UCSC.danRer10.refGene', 'org.Dr.eg.db')
for (pkg in packages) suppressMessages(suppressWarnings(library(pkg, character.only=T)))
source('Code/fig6/functions.R')

## parse arguments
args <- commandArgs(trailingOnly=T)

# function for calculationg information content of a motif
get_IC <- function(m) {
  if (min(m) == 0) m <- m+1 # pseudocount
  m <- m / colSums(m)
  return(sum(m * log2(m/.25)))
}

# load Piotr's curated motifs from paper F4
lines <- readLines(args[1])
wm_list <- c()
ic <- c()
end_of_wm <- F
for (line in lines) {
  if (startsWith(line, 'NA')) {
    name <- gsub('.p3','',gsub('.p2','',strsplit(line, '\t|\\s{1,}')[[1]][2]))
    m <- c()
  } else if (!is.na(suppressWarnings(as.integer(substring(line,1,1))))) {
    end_of_wm <- T
    m <- c(m, as.numeric(strsplit(line, '\t|\\s{1,}')[[1]][2:5])) # delimiter is either tab or a multiple of one space
  } else if (end_of_wm) {
    m <- t(matrix(m, ncol=4, byrow=T)) + 1 # pseudocount
    rownames(m) <- c('A','C','G','T')
    m <- m / colSums(m) # convert counts to weights
    ic[[name]] <- get_IC(m)
    wm_list[[name]] <- m
    end_of_wm <- F
  }
}

# load all ATAC and DNase peaks as well as the ATAC projections
dnase_peaks <- import.bed(args[2])
proj <- read.table(args[3], sep='\t', skip=2)
atac_peaks <- GRanges(proj$V2)
reference <- 'danRer10'
target <- 'mm10'
peaks_on_mm <- list(danRer10=GRanges(proj[,4], score=proj[,6]), # zf ATAC on mm coordinates
                    mm10=dnase_peaks)

# assign TAD ids to peaks
tads <- GRanges(read.table(args[4], sep='\t', header=F, col.names=c('seqnames','start','end')))
tads <- renameSeqlevels(tads, mapSeqlevels(seqlevels(tads), 'UCSC'))
ov_tad_zf <- findOverlaps(peaks_on_mm$danRer10, tads)
peaks_on_mm$danRer10$tad <- NA
peaks_on_mm$danRer10$tad[queryHits(ov_tad_zf)] <- subjectHits(ov_tad_zf)
ov_tad_mm <- findOverlaps(peaks_on_mm$mm10, tads)
peaks_on_mm$mm10$tad <- NA
peaks_on_mm$mm10$tad[queryHits(ov_tad_mm)] <- subjectHits(ov_tad_mm)

# classify zf peaks
aln_zf_mm <- import.bed(args[5])
thresh <- .99
peaks_on_mm$danRer10$name <- 'NC'
peaks_on_mm$danRer10$name[proj[,6] >= thresh] <- 'IC'
suppressWarnings(peaks_on_mm$danRer10$name[overlapsAny(GRanges(proj[,2]), aln_zf_mm)] <- 'DC')

# identify overlaps (max. gap of 200 bp)
ov <- suppressWarnings(findOverlaps(peaks_on_mm$danRer10, peaks_on_mm$mm10, maxgap=200))
# in the following, zf contains again the zf coordinates, not the projected mm coordinates (used for retrieving the sequence)
peaks_ov <- list(danRer10=GRanges(proj[queryHits(ov),2], name=peaks_on_mm$danRer10$name[queryHits(ov)]),
                 mm10=dnase_peaks[subjectHits(ov)], name=peaks_on_mm$danRer10$name[queryHits(ov)], score=peaks_on_mm$danRer10$score[queryHits(ov)])
peaks_ov$danRer10$tad <- peaks_on_mm$danRer10$tad[queryHits(ov)]
peaks_ov$mm10$tad <- peaks_on_mm$mm10$tad[queryHits(ov)]

# create random sample of zebrafish peaks in cis (within the TAD) and trans (in another TAD) to the matching peak
peaks_dr <- GRanges(proj[,2])
peaks_dr$tad <- peaks_on_mm$danRer10$tad
set.seed(111)
peaks_dr_nc_cis <- mclapply(seq_along(peaks_ov$mm10), function(i) {
  gr <- peaks_dr[which(peaks_dr$tad==peaks_ov$mm10$tad[i])]
  if(length(gr) > 0) {
    return(resize(sample(gr, 1), 200, 'center'))
  } else return(NA) # in case no other peak is found within TAD
}, mc.cores=30)
idx <- which(!is.na(peaks_dr_nc_cis))
peaks_dr_nc_cis <- do.call('c', peaks_dr_nc_cis[idx])
peaks_dr_nc_cis$id <- idx
peaks_dr_nc_trans <- mclapply(seq_along(peaks_ov$mm10), function(i) {
  gr <- peaks_dr[which(peaks_dr$tad!=peaks_ov$mm10$tad[i])]
  if(length(gr) > 0) {
    return(sample(gr, 1))
  } else return(NA)
}, mc.cores=50)
idx <- which(!is.na(peaks_dr_nc_trans))
peaks_dr_nc_trans <- resize(do.call('c', peaks_dr_nc_trans[idx]), 200, 'center')
peaks_dr_nc_trans$id <- idx

# retrieve the sequences for the overlapping regions
# for that, use the original locations in zebrafish and mouse, respectively
seqs <- suppressWarnings(list(
  danRer10=getSeq(BSgenome.Drerio.UCSC.danRer10, resize(peaks_ov$danRer10, 200, 'center')),
  mm10=getSeq(BSgenome.Mmusculus.UCSC.mm10, resize(peaks_ov$mm10, 200, 'center')),
  danRer10_nc_cis=getSeq(BSgenome.Drerio.UCSC.danRer10, peaks_dr_nc_cis), # negative control: randomly picked dr peaks within same TAD whose projections don't overlap the mm peaks
  danRer10_nc_trans=getSeq(BSgenome.Drerio.UCSC.danRer10, peaks_dr_nc_trans) # negative control: randomly picked dr peaks *not* within same TAD whose projections don't overlap the mm peaks
))

# compute motif scores in ATAC / DNase peaks
# ~ 10 mins
scores <- lapply(c(reference,target,'danRer10_nc_cis', 'danRer10_nc_trans'), function(sp) {
  res <- do.call('cbind', mclapply(wm_list, function(m) {
    scores <- suppressWarnings(motifScanScores(regionsSeq=seqs[[sp]], motifPWM=m, asPercentage=T))
    return(rowMax(scores))
  }, mc.cores=50))
  names(res) <- names(wm_list)
  return(res)
})
names(scores) <- c(reference,target,'danRer10_nc_cis', 'danRer10_nc_trans')

# compute motif hits for motifs with a given minimum IC
thresh <- 85
minIC <- 0
scores_x <- sapply(names(scores), function(x) scores[[x]][,ic>=minIC], simplify=F)
wm_list_x <- wm_list[ic>minIC]
hits_zf <- apply(scores_x$danRer10, 1, function(row) names(wm_list_x)[which(row > thresh)])
hits_mm <- apply(scores_x$mm10, 1, function(row) names(wm_list_x)[which(row > thresh)])
hits_zf_nc_cis <- apply(scores_x$danRer10_nc_cis, 1, function(row) names(wm_list_x)[which(row > thresh)])
hits_zf_nc_trans <- apply(scores_x$danRer10_nc_trans, 1, function(row) names(wm_list_x)[which(row > thresh)])

# compute motif links between zf and mm peaks. peaks with no cis control are removed (only 75 out of 8,086)
motif_links <- sapply(seq_along(hits_mm), function(i) length(intersect(hits_mm[[i]], hits_zf[[i]])))
motif_links_nc_cis <- sapply(seq_along(hits_zf_nc_cis), function(i) length(intersect(hits_mm[[peaks_dr_nc_cis$id[i]]], hits_zf_nc_cis[[i]])))
motif_links_nc_trans <- sapply(seq_along(hits_zf_nc_trans), function(i) length(intersect(hits_mm[[i]], hits_zf_nc_trans[[i]])))

# prepare data for plotting
motif_links_list <- list(matched=motif_links, cis=motif_links_nc_cis, trans=motif_links_nc_trans)
df <- do.call('rbind', lapply(names(motif_links_list), function(x) {
  do.call('rbind', lapply(c('DC','IC','NC'), function(cl) {
    df <- data.frame(table(motif_links_list[[x]][peaks_ov$danRer$name==cl]))
    colnames(df) <- c('links', 'count')
    df$cumsum <- rev(cumsum(rev(df$count)))
    df$cumsum_normed <- df$cumsum / max(df$cumsum)
    df$class <- cl
    df$set <- x
    return(df)
  }))
}))
df$links <- as.integer(as.character(df$links))

# plot cumulative motif links
w <- 14
h <- 12
cols <- c(matched='orange', cis='#777777', trans='#aaaaaa')

reverselog_trans <- function(base=exp(1)) { # function for reversing AND log-transforming x-axis
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, log_breaks(base = base), domain = c(1e-100, Inf))
}

plot_cumulative_links <- function(title, data) {
  data$class <- factor(data$class, levels=c('DC','IC','NC'))
  data$set <- factor(data$set, levels=c('matched', 'cis', 'trans'))
  data[data$links==0,'links'] <- .3 # add zero as .3 to the plot and change label manually to zero
  pl <- ggplot(data, aes(x=links, y=cumsum, color=set)) +
    geom_line(size=1, alpha=.7) +
    geom_point(shape=20, size=5, alpha=.8) +#, position=position_dodge(width=1)) +
    facet_wrap(~class, ncol=1, scale='free_y') +
    scale_x_continuous(trans=reverselog_trans(10), breaks=c(.3,1,3,10,30), labels=c(0,1,3,10,30)) +
    scale_y_log10() +
    labs(title=title, x='shared motifs', y='cumulative number of elements') +
    scale_colour_Publication(cols=cols) +
    scale_fill_Publication(cols=cols) +
    theme_Publication(base_size=22)
  return(pl)
}
pl <- plot_cumulative_links('', df)
ggsave('Figures/fig6/fig6i.pdf', pl, width=w, height=h, device=cairo_pdf)
