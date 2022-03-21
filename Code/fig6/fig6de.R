packages <- c('rtracklayer','ggplot2','tidyr','RColorBrewer','gridExtra','GenomicFeatures', 'DESeq2', 'pheatmap', 'vsn', 'grid', 'ggpmisc', 'Gviz', 'ehmm', 'RColorBrewer', 'pheatmap', 'reshape2', 'ggnewscale', 'data.table', 'BSgenome.Drerio.UCSC.danRer10', 'BSgenome.Mmusculus.UCSC.mm10', 'TxDb.Drerio.UCSC.danRer10.refGene', 'org.Dr.eg.db')
for (pkg in packages) suppressMessages(suppressWarnings(library(pkg, character.only=T)))
source('Code/fig6/functions.R')

args <- commandArgs(trailingOnly=T)

tads <- import.bed(args[1]) # 'Data/fig6/TADs/tads_top1000_CNEdensity.bed'
countdata_quantiles <- read.table(args[2], header=T, comment.char='') # 'Data/fig6/countdata_quantiles_zf_tads.csv'
countdata_quantiles <- countdata_quantiles[countdata_quantiles$feature=='H3K27me3',]
tads <- tads[tads$name %in% countdata_quantiles$tad] # remove tads for which no countdata_quantiles could be determined

## assign projection scores to countdata_quantiles
proj <- do.call('rbind', mclapply(tads$name, function(tad_id) {
    proj_file <- sprintf('Data/fig6/projection/tad_48h_1kb_bins/%s_1kb.proj', tad_id)
    tad_proj <- read.table(proj_file, sep='\t', skip=2)
    colnames(tad_proj) <- c('id', apply(read.table(proj_file, sep='\t')[1:2,], 2, paste, collapse='_')[-1])
    return(tad_proj)
}, mc.cores=50))

coords <- sapply(apply(countdata_quantiles[,c('ref_chrom','ref_coord')], 1, paste, collapse=':'), function(x) gsub(' ','',x))
countdata_quantiles$projection_score <- 0
scores <- proj$score_dijkstra[proj$coords_ref %in% coords]
countdata_quantiles[which(proj$coords_ref %in% coords), 'projection_score'] <- scores


# prepare domain data
tads_ordered <- tads$name[order(width(tads), decreasing=T)] # ordered by tad width

# countdata_quantiles
data_countdata_quantiles <- do.call('rbind', lapply(seq_along(tads_ordered), function(i) {
  df_i <- countdata_quantiles[countdata_quantiles$tad==tads_ordered[i],]
  mean_coord <- mean(c(min(df_i$ref_coord), max(df_i$ref_coord)))
  df_i <- cbind(df_i, data.frame(start=df_i$ref_coord-499, end=df_i$ref_coord+500, label='domain',
                                 start_centered=(df_i$ref_coord-499-mean_coord)/1000, end_centered=(df_i$ref_coord+500-mean_coord)/1000))
  return(df_i)
}))
data_countdata_quantiles$tad <- factor(data_countdata_quantiles$tad, levels=tads_ordered)
data_countdata_quantiles$quantile_diff <- data_countdata_quantiles$quantile_mm - data_countdata_quantiles$quantile_zf
data_countdata_quantiles$quantile_logratio <- log(data_countdata_quantiles$quantile_mm) - log(data_countdata_quantiles$quantile_zf)
data_countdata_quantiles$quantile_logratio[data_countdata_quantiles$quantile_logratio > 1] <- 1
data_countdata_quantiles$quantile_logratio[data_countdata_quantiles$quantile_logratio < -1] <- -1
data_countdata_quantiles$quantile_amplitude <- apply(data_countdata_quantiles[,c('quantile_zf','quantile_mm')], 1, max)


# define data frame with start / end coordinates of the bins normalized by the TAD lengths
# also add start / end coordinates relative to the TAD center
df <- do.call('rbind', mclapply(unique(data_countdata_quantiles$tad), function(tad) {
  df_i <- data_countdata_quantiles[which(data_countdata_quantiles$tad == tad),]
  df_i$start_norm <- (df_i$start - min(df_i$start)) / (max(df_i$end) - min(df_i$start))
  df_i$end_norm <- (df_i$end - min(df_i$start)) / (max(df_i$end) - min(df_i$start))
  df_i$chrom <- df_i$ref_chrom
  tad_center <- floor(mean(c(min(df_i$start), max(df_i$end))))
  df_i$start_relative <- df_i$start - tad_center
  df_i$end_relative <- df_i$end - tad_center
  return(df_i[,c('tad','chrom','start','end','start_norm','end_norm','start_relative','end_relative','quantile_logratio','quantile_amplitude')])
}, mc.cores=20))


# load genes
g <- suppressMessages(genes(TxDb.Drerio.UCSC.danRer10.refGene))
g$gene_symbols <- mapIds(org.Dr.eg.db, keys=g$gene_id, keytype='ENTREZID', column='SYMBOL') # get gene symbols from ENSEMBL IDs (without point and trailing number)


# prepare genes data
data_genes <- do.call('rbind', mclapply(seq_along(tads_ordered), function(i) {
  gr_i <- GenomicRanges::shift(resize(GRanges(apply(countdata_quantiles[which(countdata_quantiles$tad==tads_ordered[i]),c('ref_chrom','ref_coord')], 1, paste, collapse=':')), 1000, 'center'),1)
  g_i <- sort(g[overlapsAny(g,gr_i)], ignore.strand=T)
  if (length(g_i)==0) {
    df_i <- data.frame()
  } else {
    df_i <- data.frame(start=start(g_i), end=end(g_i), strand=strand(g_i), label='gene', tad=tads_ordered[i])
    df_i[df_i$strand=='-',c('start','end')] <- df_i[df_i$strand=='-',c('end','start')] # invert genes on minus strand for arrow direction in plot
    mean_tad_coord <- start(resize(tads[tads$name==tads_ordered[i]], 1, 'center'))
    df_i$start_relative <- (df_i$start - mean_tad_coord)
    df_i$end_relative <- (df_i$end - mean_tad_coord)
    df_i$gene_symbol <- g_i$gene_symbols
    df_i$gene_symbol[is.na(df_i$gene_symbol)] <- ''
  }
  return(df_i)
}, mc.cores=20))
data_genes$tad <- factor(data_genes$tad, levels=tads_ordered)

# TAD boundaries
data_tad_boundaries <- data.frame(tads)
colnames(data_tad_boundaries)[6] <- 'tad'
center <- apply(data_tad_boundaries[,c('start','end')], 1, mean)
data_tad_boundaries$start_relative <- (data_tad_boundaries$start - center)
data_tad_boundaries$end_relative <- (data_tad_boundaries$end - center)
data_tad_boundaries <- gather(data_tad_boundaries[,c('tad','start_relative','end_relative')], key='coordinate', value='value', -tad)
data_tad_boundaries$label <- 'TAD boundary'

# within TAD sorting
# order the bins within a TAD according to their absolute logratio multiplied by the amplitude for the greenness sorted plot
# make groups (green/blue/red) based on threshold, then sort the groups according to logratio and amplitude
order_within_tads <- function(df, logratio_thresh=.2, amplitude_power=10, amplitude_thresh=.1) {
  df <- do.call('rbind', lapply(tads$name, function(tad) {
    df_i <- df[which(df$tad == tad),]
    x <- df_i$quantile_logratio
    y <- df_i$quantile_amplitude^amplitude_power
       
    # order blue: highest values first
    blue <- which(x < -logratio_thresh & y > amplitude_thresh)
    ordr_blue <- order(y[blue], decreasing=T)
    
    # order green: highest values in the middle
    green <- which(abs(x) <= logratio_thresh & y > amplitude_thresh)
    ordr_green <- order(y[green], decreasing=F)
    if (length(green) >= 2) ordr_green <- c(ordr_green[seq(1,length(ordr_green),2)], rev(ordr_green[seq(2,length(ordr_green),2)]))
    
    # order red: lowest values first
    red <- which(x > logratio_thresh & y > amplitude_thresh)
    ordr_red <- order(y[red], decreasing=F)
    
    # order white: consider setting their amplitude to zero to make them all appear totally white
    # split them in two to add left and right of the color bins
    # give them the length that results in the green bins being centered, i.e. w1 = x/2 - g/2 - b, and w2 = w - w1
    white <- which(is.nan(x) | y < amplitude_thresh)
    len_white1 <- min(max(ceiling(length(x)/2 - length(green)/2 - length(blue)), 0), length(white))
    len_white2 <- length(white) - len_white1

    x_ordered <- unlist(c(rep(0,len_white1), x[blue][ordr_blue], x[green][ordr_green], x[red][ordr_red], rep(0,len_white2)))
    y_ordered <- unlist(c(rep(0,len_white1), y[blue][ordr_blue], y[green][ordr_green], y[red][ordr_red], rep(0,len_white2)))
    
    df_i$quantile_logratio <- x_ordered
    df_i$quantile_amplitude <- y_ordered^(1/amplitude_power)
    return(df_i)
  }))
  return(df)
}

# function for continuous color-sorting of the TADs 
sort_tads_continuously <- function(df, color=c('green','red','blue'), logratio_thresh=.2, amplitude_power=10, amplitude_thresh=.1) {
  color_score <- sapply(unique(df$tad), function(tad_i) {
    df_i <- data_countdata_quantiles[which(data_countdata_quantiles$tad == tad_i),]
    x <- df_i$quantile_logratio
    y <- df_i$quantile_amplitude^amplitude_power
    y[which(y < amplitude_thresh)] <- 0
    if (color=='green') {
        return(sum(abs(x) <= logratio_thresh, na.rm=T) * mean(y[which(abs(x) <= logratio_thresh)], na.rm=T))
    } else if (color=='blue') {
        return(sum(x < -logratio_thresh, na.rm=T) * mean(y[which(x < -logratio_thresh)], na.rm=T))
    } else if (color=='red') {
        return(sum(x > logratio_thresh, na.rm=T) * mean(y[which(x > logratio_thresh)], na.rm=T))
    }
  })
  names(color_score) <- unique(df$tad)
  tads_by_color <- sort(color_score, decreasing=T)
  return(setNames(names(tads_by_color), tads_by_color)) # switch names / values so that names = score, value = TAD
}

# TAD sorting
# classify TADs as green / blue / red according to the number of bins with the logratio in the target range weighted by the mean amplitudes of those bins
# fractions aren't real fractions here but rather weighted scores
sort_tads <- function(df, tad_order=c('green','red','blue','tad','gbr','brg','bgr','rgb','rbg'), logratio_thresh=.2, amplitude_power=10, amplitude_thresh=.1) {
  if (tad_order %in% c('green','blue','red')) { # don't group TADs, sort them continuously by a color
    return(sort_tads_continuously(df, color=tad_order, amplitude_thresh=amplitude_thresh))
  } else {
    fractions <- data.frame(do.call('rbind', lapply(unique(df$tad), function(tad_i) {
      df_i <- data_countdata_quantiles[data_countdata_quantiles$tad == tad_i,]
      x <- df_i$quantile_logratio
      y <- df_i$quantile_amplitude^amplitude_power
      fractions <- c(green = sum(abs(x) <= logratio_thresh, na.rm=T) * mean(y[abs(x) <= logratio_thresh], na.rm=T),
                     blue = sum(x < -logratio_thresh, na.rm=T) * mean(y[x < -logratio_thresh], na.rm=T),
                     red = sum(x > logratio_thresh, na.rm=T) * mean(y[x > logratio_thresh], na.rm=T))
      fractions[which(is.na(fractions))] <- 0
      return(fractions)
    })))
    rownames(fractions) <- unique(df$tad)

    green_tads <- fractions[which(fractions$green > fractions$blue & fractions$green > fractions$red),]
    green_tads <- green_tads[order(green_tads$green, decreasing=T),]
    blue_tads <- fractions[which(fractions$blue > fractions$green & fractions$blue > fractions$red),]
    blue_tads <- blue_tads[order(blue_tads$blue, decreasing=T),]
    red_tads <- fractions[which(fractions$red > fractions$blue & fractions$red > fractions$green),]
    red_tads <- red_tads[order(red_tads$red, decreasing=T),]
    if (tad_order=='tad') return(c(rownames(green_tads), rownames(red_tads), rownames(blue_tads)))
    if (tad_order=='gbr') return(c(rownames(green_tads), rownames(blue_tads), rownames(red_tads)))
    if (tad_order=='brg') return(c(rownames(blue_tads), rownames(red_tads), rownames(green_tads)))
    if (tad_order=='bgr') return(c(rownames(blue_tads), rownames(green_tads), rownames(red_tads)))
    if (tad_order=='rgb') return(c(rownames(red_tads), rownames(green_tads), rownames(blue_tads)))
    if (tad_order=='rbg') return(c(rownames(red_tads), rownames(blue_tads), rownames(green_tads)))
  }
}

# funnel plot function
plot_funnel <- function(df, n=Inf, plot_genes_and_tad_boundaries=F,
                        logratio_thresh=.2, amplitude_power=10, amplitude_thresh=.1,
                        tad_order=c('green','blue','red'), within_tad_order=c('original','sort'),
                        coordinates=c('original','normalized')) {
  cols <- colorRampPalette(c(blue='#1982c4', blue='#1982c4', green='#07731b', red='#ff0000', red='#ff0000'))(100)
  
  # order bins within TADs
  if (within_tad_order == 'original') df_countdata_quantiles <- df
  if (within_tad_order == 'sort') df_countdata_quantiles <- order_within_tads(df, logratio_thresh=logratio_thresh, amplitude_power=amplitude_power)
  if (coordinates == 'original') {
    df_countdata_quantiles$start <- df_countdata_quantiles$start_relative
    df_countdata_quantiles$end <- df_countdata_quantiles$end_relative
    xlab <- 'distance to TAD center'
  } else if (coordinates == 'normalized') {
    df_countdata_quantiles$start <- df_countdata_quantiles$start_norm
    df_countdata_quantiles$end <- df_countdata_quantiles$end_norm
    xlab <- 'relative position in TAD'
  }
  
  # order TADs
  tads_ordered <- tad_order # if a selection is passed
  if (is.character(tad_order) && tad_order %in% c('green','blue','red')) { # if a color is passed 
    tads_ordered <- sort_tads(df_countdata_quantiles, tad_order=tad_order, logratio_thresh=logratio_thresh, amplitude_power=amplitude_power, amplitude_thresh=amplitude_thresh)
  }
  df_countdata_quantiles$tad <- factor(df_countdata_quantiles$tad, levels=tads_ordered)

  # limit to top n TADs
  n <- min(n, length(unique(df$tad)))
  df_countdata_quantiles <- df_countdata_quantiles[which(df_countdata_quantiles$tad %in% tads_ordered[1:n]),]

  # include genes and TAD boundaries in original coordinate / sorting plot
  df_genes <- data_genes[0,]
  df_tad_boundaries <- data_tad_boundaries[0,]
  size_gene_names <- 5 * n / length(unique(df$tad))
  if (plot_genes_and_tad_boundaries==T && within_tad_order == 'original' && coordinates == 'original') {
    df_genes <- data_genes[which(data_genes$tad %in% tads_ordered[1:n]),]
    df_tad_boundaries <- data_tad_boundaries[which(data_tad_boundaries$tad %in% tads_ordered[1:n]),]
    df_genes$tad <- factor(df_genes$tad, levels=tads_ordered)
    df_tad_boundaries$tad <- factor(df_tad_boundaries$tad, levels=tads_ordered)
  }
  
  pl <- ggplot() +
    geom_rect(data=df_countdata_quantiles, mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1, fill=quantile_logratio, alpha=quantile_amplitude^amplitude_power)) +
    facet_wrap(~tad, ncol=1, strip.position='left') +
    labs(x=xlab, y='TADs') +
    scale_y_reverse(limits=c(1,0), expand=c(0,0), breaks=NULL) +
    scale_x_continuous(expand=c(0,0)) +
    scale_alpha_identity(limits=c(0,1)) + # make the alpha actually start at zero instead of 0.1 (goddamnit)
    scale_fill_gradientn(name='H3K27me3 ratio', colours=cols, limits=c(-1,1)) +
    scale_colour_Publication(name='', cols=c(gene='black', `TAD boundary`='black')) +
    theme_Publication() +
    theme(axis.ticks.y=element_blank(), axis.line.y=element_blank(), axis.text.y=element_blank(),
          panel.spacing=unit(0,'cm'), legend.position='none', panel.grid.major=element_blank(), strip.background=element_blank(), strip.text=element_blank())
  if (plot_genes_and_tad_boundaries==T && within_tad_order == 'original' && coordinates == 'original') pl <- pl +
    geom_vline(data=df_tad_boundaries, mapping=aes(xintercept=value, color=label)) +
    geom_segment(data=df_genes, mapping=aes(x=start_relative, xend=end_relative, y=.7, yend=.7, color=label),
                 alpha=1, arrow=arrow(length=unit(.2,'npc'))) +
    geom_text(data=df_genes, mapping=aes(x=(start_relative+end_relative)/2, y=.3, label=gene_symbol), size=size_gene_names)
  return(pl)
}

# 6d: plot and save
w <- 30
h <- 4
tad_selection <- c('tad_0017','tad_0069','tad_0521','tad_0221')
pl_selection <- plot_funnel(df[which(df$tad %in% tad_selection),], tad_order=tad_selection, within_tad_order='original',
            coordinates='original', logratio_thresh=.2, amplitude_power=10, plot_genes_and_tad_boundaries=T) +
  theme(strip.text.y=element_text(size=9), strip.text.y.left=element_text(angle=0))
ggsave('Figures/fig6/fig6d.pdf', pl_selection, width=w, height=h, device=cairo_pdf)

# 6e: plot (~ 1 min)
w <- 5
h <- 8
set_plot_dimensions(w,h)
pl <- plot_funnel(df, tad_order='green', within_tad_order='sort', coordinates='normalized', logratio_thresh=.2, amplitude_power=10)


# 6e: save to file ( ~ 15 sec)
cairo_pdf('Figures/fig6/fig6e.pdf', width=w, height=h, onefile=T)
pl
dev.off()
