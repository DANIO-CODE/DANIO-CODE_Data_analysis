packages <- c('rtracklayer','ggplot2','tidyr','RColorBrewer','gridExtra','GenomicFeatures', 'DESeq2', 'pheatmap', 'vsn', 'grid', 'ggpmisc', 'Gviz', 'ehmm', 'RColorBrewer', 'pheatmap', 'reshape2', 'ggnewscale', 'data.table', 'BSgenome.Drerio.UCSC.danRer10', 'BSgenome.Mmusculus.UCSC.mm10', 'TxDb.Drerio.UCSC.danRer10.refGene', 'org.Dr.eg.db')
for (pkg in packages) suppressMessages(suppressWarnings(library(pkg, character.only=T)))
source('Code/fig6/functions.R')

filenames <- commandArgs(trailingOnly=T) # sapply(1:1000, function(i) sprintf('data/projection/tad_48h_1kb_bins/tad_%04d_1kb.proj', i))
print(filenames)
df <- do.call('rbind', lapply(filenames, function(fn) {
    if (length(count.fields(fn)) <= 2) {
        return(data.frame())
    } else {
        headers <- scan(fn, nlines=2, what=character())
        proj <- read.table(fn, sep='\t', skip=2)
        colnames(proj) <- c('id', paste(headers[1:(length(headers)/2)], headers[(length(headers)/2+1):length(headers)], sep='_'))
        return(rbind(data.frame(distance=proj$ref_dist_closest_anchor_direct, method='Direct', species='Zebrafish'), data.frame(distance=proj$ref_dist_closest_anchor_dijkstra, method='Multi', species='Zebrafish'),
                     data.frame(distance=proj$qry_dist_closest_anchor_direct, method='Direct', species='Mouse'), data.frame(distance=proj$qry_dist_closest_anchor_dijkstra, method='Multi', species='Mouse')))
    }
}))

# boxplot anchor distances
w <- 10
h <- 5
set_plot_dimensions(w,h)
df$distance <- as.numeric(df$distance)
pl <- ggplot(df, aes(method, log10(distance+1), fill=species)) +
    geom_boxplot(notch=T, size=.6) +
    scale_y_continuous(breaks=c(2,4,6), labels=c('100 bp', '10 Kb', '1 Mb')) +
    labs(y='distance to closest anchor') +
    coord_flip() +
    scale_fill_Publication(cols=rep(c('#1982c4','#ff0000'),2)) +
    theme_Publication(base_size=22) +
    theme(panel.spacing=unit(3.5,'lines'), legend.position=c(.77,.8))
ggsave('Figures/fig6/fig6b.pdf', pl, width=w, height=h, device=cairo_pdf)
