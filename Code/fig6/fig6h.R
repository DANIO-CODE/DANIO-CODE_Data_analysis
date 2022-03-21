packages <- c('rtracklayer','ggplot2','tidyr','RColorBrewer','gridExtra','GenomicFeatures', 'DESeq2', 'pheatmap', 'vsn', 'grid', 'ggpmisc', 'Gviz', 'ehmm', 'RColorBrewer', 'pheatmap', 'reshape2', 'ggnewscale', 'data.table', 'BSgenome.Drerio.UCSC.danRer10', 'BSgenome.Mmusculus.UCSC.mm10', 'TxDb.Drerio.UCSC.danRer10.refGene', 'org.Dr.eg.db')
for (pkg in packages) suppressMessages(suppressWarnings(library(pkg, character.only=T)))
source('Code/fig6/functions.R')

## parse arguments
args <- commandArgs(trailingOnly=T)

# load mouse segmentation and annotate  states
chromhmm_mm <- import.bed(args[1])
state_names <- c('1_ReprPC', '2_Quies1', '3_EnhPois1', '4_EnhPois2', '5_EnhA', '6_EnhWk', '7_Quies2', '8_TssBiv', '9_TssFlank', '10_TssA')
chromhmm_mm$state <- sapply(chromhmm_mm$name, function(i) state_names[as.integer(i)])
chromhmm_mm$itemRgb <- chromhmm_mm$thick <- chromhmm_mm$score <- chromhmm_mm$name <- NULL

# define chromhmm state colors (adapted from Damir). I chose the two yellows for the EnhPois states.
# paper_colors <- c('#A6CEE3', '#1F78B4', '#33A02C', '#B2DF8A', '#E31A1C', '#FB9A99', '#FF7F00', '#6A3D9A', '#CAB2D6' , '#A1A2A3')
# names(paper_colors) <- c('1_TssA1', '2_TssA2', '3_TssFlank1', '4_TssFlank2', '5_EnhA1', '6_EnhA2', '7_EnhWk1', '8_TssBiv', '9_ReprPC', '10_Quies')
paper_colors <- c('#CAB2D6', '#A1A2A3', '#FFF600', '#FFD900', '#E31A1C', '#FF7F00', '#B6B8BA', '#6A3D9A', '#B2DF8A', '#1F78B4')
names(paper_colors) <- c('1_ReprPC', '2_Quies1', '3_EnhPois1', '4_EnhPois2', '5_EnhA', '6_EnhWk', '7_Quies2', '8_TssBiv', '9_TssFlank', '10_TssA')

# read zebrafish ATAC peaks
f <- args[2]
atac_proj_df <- read.table(f, sep='\t', skip=2)[,1:6]
colnames(atac_proj_df) <- c('id', 'ref_coord', 'direct_coord', 'multi_coord', 'direct_score', 'multi_score')

# classify peaks based on projection scores (scores > .99 are representing indirectly conserved regions, if the region overlaps a direct alignment it is directly conserved)
label_order <- c('DC','IC','NC')
score_threshold <- .99
atac_proj_df$class <- apply(atac_proj_df, 1, function(row) ifelse(as.numeric(row['multi_score']) >= score_threshold, 'IC', 'NC'))
aln <- import.bed(args[3])
atac_proj_df$class[which(overlapsAny(GRanges(atac_proj_df$ref_coord), aln))] <- 'DC'
atac_proj_df$class <- factor(atac_proj_df$class, levels=rev(label_order))

## make list of ATAC peak GRanges by state, then show the mouse state frequencies for the projections of every zf state (faceted by DC/IC/NC)
gr_zf <- GRanges(atac_proj_df$ref_coord)
chromhmm_zf <- import.bed(args[4])
chromhmm_zf$name[chromhmm_zf$name=='8_Pois'] <- '8_TssBiv' # new state name according to figure 5
chromhmm_zf$name[chromhmm_zf$name=='6_EnhFlank'] <- '6_EnhA2' # new state name according to figure 5
ov_zf <- findOverlaps(gr_zf, chromhmm_zf)
atac_proj_df$state_zf <- chromhmm_zf$name[subjectHits(ov_zf)]
gr_mm <- GRanges(atac_proj_df$multi_coord)
ov_mm <- findOverlaps(gr_mm, chromhmm_mm)
atac_proj_df$state_mm <- chromhmm_mm$state[subjectHits(ov_mm)]
df <- do.call('rbind', sapply(c('DC','IC','NC'), function(class) {
  cbind(data.frame(class=class), do.call('rbind', sapply(unique(atac_proj_df$state_zf), function(state_zf) {
    df_i <- data.frame(table(atac_proj_df[atac_proj_df$class==class & atac_proj_df$state_zf==state_zf, 'state_mm']) / sum(atac_proj_df$class==class & atac_proj_df$state_zf==state_zf))
    colnames(df_i) <- c('state_mm', 'freq')
    df_i$state_zf <- state_zf
    return(df_i[,c('state_zf','state_mm','freq')])
  }, simplify=F)))
}, simplify=F))
rownames(df) <- NULL

df$state_zf <- factor(df$state_zf, levels=c('5_EnhA1', '6_EnhA2', '7_EnhWk1', '1_TssA1', '2_TssA2', '3_TssFlank1', '4_TssFlank2', '8_TssBiv', '9_ReprPC', '10_Quies'))
df$state_mm <- factor(df$state_mm, levels=c('5_EnhA', '6_EnhWk', '3_EnhPois1', '4_EnhPois2', '10_TssA', '9_TssFlank', '8_TssBiv', '1_ReprPC', '2_Quies1', '7_Quies2'))
df$class <- factor(df$class, levels=c('NC', 'IC', 'DC'))

w <- 18
h <- 15
cols <- c(rev(brewer.pal(7, 'Greens'))[1:4], rev(brewer.pal(7, 'Blues'))[1:3], rev(brewer.pal(7, 'Greys'))[3:5])
pl <- ggplot() + 
  geom_bar(data=df, aes(freq, class, fill=state_mm), stat='identity', position='stack', alpha=.5) +
  facet_wrap(~state_zf, ncol=1, scale='free_y') +
  scale_x_reverse() +
  scale_fill_Publication(cols=paper_colors) +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 90))
ggsave('Figures/fig6/fig6h.pdf', width=w, height=h, device=cairo_pdf)
