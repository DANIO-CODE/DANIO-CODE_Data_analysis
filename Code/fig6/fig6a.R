packages <- c('rtracklayer','ggplot2','tidyr','RColorBrewer','gridExtra','GenomicFeatures', 'DESeq2', 'pheatmap', 'vsn', 'grid', 'ggpmisc', 'Gviz', 'ehmm', 'RColorBrewer', 'pheatmap', 'reshape2', 'ggnewscale', 'data.table', 'BSgenome.Drerio.UCSC.danRer10', 'BSgenome.Mmusculus.UCSC.mm10', 'TxDb.Drerio.UCSC.danRer10.refGene', 'org.Dr.eg.db')
for (pkg in packages) suppressMessages(suppressWarnings(library(pkg, character.only=T)))
source('Code/fig6/functions.R')

args <- commandArgs(trailingOnly=T)
if (length(args) != 1) stop('Usage: Rscript fig6a.R path_to_tad_projection_file.proj')

# function for centering data on median coordinate and inverting target start/end coords if on negative strand
prepare_data_for_plotting <- function(reference, target) {
  from_reference <- min(c(start(reference), end(reference)))
  to_reference <- max(c(start(reference), end(reference)))
  center_reference <- round(mean(c(from_reference,to_reference)) / 1e6, 2) * 1e6
  from_target <- min(c(start(target), end(target)))
  to_target <- max(c(start(target), end(target)))
  center_target <- round(mean(c(from_target,to_target)) / 1e6, 2) * 1e6
  strand_target <- ifelse(start(target)[1] < rev(end(target))[1], '+', '-')
  
  df <- cbind(data.frame(chrom_reference=seqnames(reference)[1], coord_reference=start(reference)),
              data.frame(chrom_target=seqnames(target)[1], coord_target=start(target)))
  df$coord_reference_centered <- df$coord_reference - center_reference
  df$coord_target_centered <- df$coord_target - center_target
  if (strand_target == '-') df$coord_target_centered <- -df$coord_target_centered
  res <- list(df=df, from_reference=from_reference, to_reference=to_reference, center_reference=center_reference,
              from_target=from_target, to_target=to_target, center_target=center_target, strand_target=strand_target)
  return(res)
}

## function for converting a minus or plus character to an signed integer
sign_to_int <- function(x=c('-','+')) return(c('-'=-1, '+'=1)[x])

# plot function
plot_fnc <- function(data, plot_genes=F) {
  pl <- ggplot() +
  geom_segment(data=data$df, mapping=aes(x=coord_reference_centered, xend=coord_target_centered, y=0, yend=1), alpha=.8, size=1.2) +
  scale_x_continuous(name = paste(reference, data$df$chrom_reference[1], sep=': '),
                     labels = function(x) (x + data$center_reference) / 1e6,
                     sec.axis = dup_axis(name = paste(target, data$df$chrom_target[1], sep=': '),
                                         labels = function(x) {
                                           y <- (x + data$center_target) / 1e6
                                           if (data$strand_target == '-') y <- rev(y)
                                           return(y)
                                         })) +
  theme_Publication() +
  theme(legend.position='none', axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), axis.line.y=element_blank(),
        panel.grid.major=element_blank())
  if (plot_genes) {
    pl <- pl +
      geom_rect(data=genes_reference, mapping=aes(xmin=start_centered,xmax=end_centered,ymin=-.02,ymax=.02), fill='grey', color='black', alpha=.8) +
      geom_rect(data=genes_target, mapping=aes(xmin=start_centered,xmax=end_centered,ymin=.98,ymax=1.02), fill='grey', color='black', alpha=.8)
  }
  return(pl)
}

# load TAD projection coordinates (irx3) and print the scaling of their sizes in mouse and fish
proj_file <- args[1] # 'data/projection/tad_48h_1kb_bins/tad_0017_1kb.proj'
proj <- read.table(proj_file, sep='\t', skip=2)
zf <- GRanges(proj[,2])
mm <- GRanges(proj[,4])

## prepare data for plotting
reference <- 'zebrafish'
target <- 'mouse'
idx <- seq(1, length(zf), 7)
res <- prepare_data_for_plotting(zf[idx], mm[idx])

# genes
genes_reference <- data.frame(start=35770013, end=35773350, name='Irx3a')
genes_reference$start_centered <- genes_reference$start - res$center_reference
genes_reference$end_centered <- genes_reference$end - res$center_reference
genes_target <- data.frame(start=91798525, end=91802067, name='Irx3')
genes_target$start_centered <- (genes_target$start - res$center_target) * sign_to_int(res$strand_target)
genes_target$end_centered <- (genes_target$end - res$center_target) * sign_to_int(res$strand_target)

# plot full GRB (every 7th bin)
w <- 20
h <- 5
set_plot_dimensions(w,h)
pl1 <- plot_fnc(res, plot_genes=T)

# plot zoomed in Irx3a locus (every 10th bin)
irx3_idx <- which(start(zf) >= 35693000 & start(zf) <= 35791000)
zf_irx3 <- zf[irx3_idx]
mm_irx3 <- mm[irx3_idx]
idx_irx3 <- seq(1, length(zf_irx3), 1)
res_irx3 <- prepare_data_for_plotting(zf_irx3[idx_irx3], mm_irx3[idx_irx3])
pl2 <- plot_fnc(res_irx3)

cairo_pdf('Figures/fig6/fig6a.pdf', width=w, height=h, onefile=T)
pl1
pl2
dev.off()
