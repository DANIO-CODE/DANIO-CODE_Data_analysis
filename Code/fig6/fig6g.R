packages <- c('rtracklayer','ggplot2','tidyr','RColorBrewer','gridExtra','GenomicFeatures', 'DESeq2', 'pheatmap', 'vsn', 'grid', 'ggpmisc', 'Gviz', 'ehmm', 'RColorBrewer', 'pheatmap', 'reshape2', 'ggnewscale', 'data.table', 'ggsignif', 'effsize', 'BSgenome.Drerio.UCSC.danRer10', 'BSgenome.Mmusculus.UCSC.mm10', 'TxDb.Drerio.UCSC.danRer10.refGene', 'org.Dr.eg.db')
for (pkg in packages) suppressMessages(suppressWarnings(library(pkg, character.only=T)))
source('Code/fig6/functions.R')

## parse args
args <- commandArgs(trailingOnly=T)

# compute effect size (in addition to p-values which are not very meaningful because of large sample size)
compute_effect_sizes <- function(df, classes=c('DC','IC','NC')) {
  effect_size <- do.call('rbind', lapply(seq(1,length(classes)-1), function(i) {
    c1 <- classes[i]
    do.call('rbind', lapply(seq(i+1,length(classes)), function(j) {
      c2 <- classes[j]
      chn <- cohen.d(df[df$class==c1, 'score'],
                     df[df$class==c2, 'score'])
      return(data.frame(c1=c1, c2=c2, cohen_d=round(chn$estimate,2), magnitude=chn$magnitude))
    }))
  }))
    effect_size$y <- 10^c(.34,.64,.04) # order according to how connective lines of p-values are drawn
  effect_size$class <- effect_size$c1
  effect_size$cohen_d <- paste0('cohen\'s d: ', effect_size$cohen_d)
  return(effect_size)
}

plot_fnc <- function(data, test=c(NULL,'t.test','wilcox.test')) {
  pc <- 1e-2
  pl <- ggplot(data, aes(class, value+pc, fill=class, color=class)) + 
    geom_boxplot(alpha=.4, notch=T, outlier.shape=NA) + scale_y_log10() + 
    facet_wrap(~feature, nrow=1, strip.position='bottom') + 
    scale_fill_Publication(breaks=label_order, cols=cols) + 
    scale_colour_Publication(breaks=label_order, cols=cols) +
    theme_Publication(base_size=base_size) + 
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(), strip.placement='outside', strip.background=element_blank(), strip.text=element_text(size=16))
  if (!is.null(test)) {
    effect_size <- compute_effect_sizes(feature_data)
    feature_plot <- feature_plot +
      geom_signif(color='black', comparison=list(group=c('IC','NC')), map_signif_level=sigFunc, test.args=list(alternative="two.sided"), y_position=.1, test=test) +
      geom_signif(color='black', comparison=list(group=c('IC','DC')), map_signif_level=sigFunc, test.args=list(alternative="two.sided"), y_position=.3, test=test) +
      geom_signif(color='black', comparison=list(group=c('NC','DC')), map_signif_level=sigFunc, test.args=list(alternative="two.sided"), y_position=.5, test=test) +
      geom_text(data=effect_size, mapping=aes(x=c1, y=y, label=cohen_d, hjust=hjust), color='black')
  }
  return(pl)
}

# load zebrafish atac peaks and projections

width <- 200
proj <- read.table(args[1], sep='\t', skip=2)
atac_peaks_proj <- resize(GRanges(proj$V4), width=width, fix='center')


# classify peaks based on projection scores (scores > .99 are representing indirectly conserved regions, if the region overlaps a direct alignment it is directly conserved)
label_order <- c('DC','IC','NC')
scores <- data.frame(direct=proj$V5, multi=proj$V6)
score_threshold <- .99
scores$class <- apply(scores, 1, function(row) ifelse(row['multi'] >= score_threshold, 'IC', 'NC'))
aln <- import.bed(args[2])
scores$class[overlapsAny(GRanges(proj$V2), aln)] <- 'DC'
scores$class <- factor(scores$class, levels=rev(label_order))
atac_peaks_proj$class <- scores$class

# plot number of peaks in different classes
cols <- c(IC='#386cb0',NC='#fdb462',DC='#7fc97f')
w <- 12
h <- 2
set_plot_dimensions(w,h)
pl <- ggplot(scores, aes(1, fill=class)) +
  geom_bar(alpha=.7) +
  coord_flip() +
  geom_text(stat='count', aes(label=..count..), size=5, position=position_stack(vjust=.5)) +
  scale_colour_Publication(breaks=label_order, cols=cols) +
  scale_fill_Publication(breaks=label_order, cols=cols) +
  theme_Publication(base_size=12) +
  theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank(), panel.grid.major.y=element_blank())
ggsave('Figures/fig6/fig6f_bottom.pdf', pl, width=w, height=h, device=cairo_pdf)

# assess global maximum (97.5th percentile) for normalization of dnase-seq data
# these are the maxima of the mean coverage over 200 bp windows centering on the projections from zebrafish to mouse

dnase_mm <- import.bw(args[3], which=atac_peaks_proj)
atac_peaks_proj$score <- ehmm:::aggScore(atac_peaks_proj, dnase_mm, 'mean')$score
global_maximum <- quantile(atac_peaks_proj$score, .975)
atac_peaks_proj$score <- sapply(atac_peaks_proj$score / global_maximum, function(x) min(x,1)) # normalize to 97.5th percentile
atac_peaks_proj$class <- factor(atac_peaks_proj$class, levels=label_order)


# function to transform p-values to -log10(p)
sigFunc <- function(x) {
  if (x > 1e-5) {paste0('-log10(p): ', round(-log10(x),2))}
  else {'-log10(p): > 5'}
}

# load epigenetic data in target species for all enhancers, compare classes.
w <- 4
h <- 6
pc <- 1e-2
cols <- c(DC='#7fc97f', IC='#fdb462', NC='#8d99ae')
effect_size <- compute_effect_sizes(data.frame(atac_peaks_proj))
pl <- ggplot(data.frame(atac_peaks_proj), aes(class, score+pc, fill=class, color=class)) + 
  geom_boxplot(alpha=.4, notch=T, outlier.shape=NA) + scale_y_log10() + 
  geom_signif(color='black', comparison=list(group=c('IC','NC')), map_signif_level=sigFunc, test.args=list(alternative="two.sided"), y_position=.1, test='t.test') +
  geom_signif(color='black', comparison=list(group=c('IC','DC')), map_signif_level=sigFunc, test.args=list(alternative="two.sided"), y_position=.4, test='t.test') +
  geom_signif(color='black', comparison=list(group=c('NC','DC')), map_signif_level=sigFunc, test.args=list(alternative="two.sided"), y_position=.7, test='t.test') +
  geom_text(data=effect_size, mapping=aes(x=c1, y=y, label=cohen_d, hjust=0), color='black') +
  labs(y='DNase-seq') +
  scale_fill_Publication(breaks=label_order, cols=cols) + 
  scale_colour_Publication(breaks=label_order, cols=cols) +
  theme_Publication(base_size=20) + 
  theme(axis.title.x=element_blank(), legend.position='none', strip.placement='outside', strip.background=element_blank(), strip.text=element_text(size=16)) 
ggsave('Figures/fig6/fig6g.pdf', pl, width=w, height=h, device=cairo_pdf)
