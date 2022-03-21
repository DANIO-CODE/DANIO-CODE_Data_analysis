packages <- c('rtracklayer','ggplot2','tidyr','RColorBrewer','gridExtra','GenomicFeatures', 'DESeq2', 'pheatmap', 'vsn', 'grid', 'ggpmisc', 'Gviz', 'ehmm', 'RColorBrewer', 'pheatmap', 'reshape2', 'ggnewscale', 'data.table', 'BSgenome.Drerio.UCSC.danRer10', 'BSgenome.Mmusculus.UCSC.mm10', 'TxDb.Drerio.UCSC.danRer10.refGene', 'org.Dr.eg.db', 'dplyr')
for (pkg in packages) suppressMessages(suppressWarnings(library(pkg, character.only=T)))
source('Code/fig6/functions.R')

## parse arguments
args <- commandArgs(trailingOnly=T)
nthreads <- args[1]

# load chromosome sizes
sizes <- list(danRer10=ehmm:::readGenomeSize(args[2]),
              mm10=ehmm:::readGenomeSize(args[3]))

# read whole genome count data from 200bp bigwigs.
# why? later, I will compute percentiles on the whole genome count distribution. this distribution is based on the bin size.
# for example, 2kb bins are less likely to be zero, but the many zeros from the 200bp bins will make the average value of a non-zero bin decrease compared to the 200bp bin distribution.
species <- c('danRer10', 'mm10')
features <- c('H3K27ac', 'H3K27me3')
files <- list(
  danRer10=list(
    H3K27me3=args[4],
    H3K27ac=args[5]
  ),
  mm10=list(
    H3K27me3=args[6],
    H3K27ac=args[7]
  )
)
whole_genome_count_data <- mclapply(species, function(sp) {
  sapply(features, function(f) {
    return(import.bw(files[[sp]][[f]]))
  }, simplify=F)
}, mc.cores=min(length(species), nthreads))
names(whole_genome_count_data) <- species

# read projection coordinates
binsize <- 1000
proj_files <- Sys.glob('Data/fig6/projection/tad_48h_1kb_bins/tad*proj')
tads <- gsub('.proj', '', gsub('_1kb', '', basename(proj_files)))
names(proj_files) <- tads
# # some projections in 6 TADs don't have coordinates for `direct` because they are at a breaking point of an inversion and don't have collinear upstream or downstream anchors, i.e. no projection.
# # remove those TADs from the analysis. (alternatively, if really interesting TADs, just work with multi, because in all 6 cases this only affects direct projections.)
gr <- mclapply(tads, function(tad) {
    if (length(count.fields(proj_files[[tad]])) <= 2) {
        return(NA)
    } else {
        tad_proj <- read.table(proj_files[[tad]], sep='\t', skip=2)[,2:4]
        ## remove coordinates with no direct projection coordinates. if they are at the borders, i.e. the remaining reference coordinates are still consecutive  / neighbouring, keep the TAD, otherwise discard. they mean breaking points (e.g. inversions) and are difficult to handle.
        tad_proj <- tad_proj[apply(tad_proj, 1, function(r) !any(r %in% c('')) & !any(is.na(r))),]
        remaining_ref <- sapply(tad_proj$V2, function(x) as.integer(strsplit(as.character(x), ':')[[1]][[2]]))
        if (length(remaining_ref) == 0 | !all(diff(remaining_ref)==binsize)) {
            return(NA)
        } else {
            gr_ref <- resize(GRanges(tad_proj[,1]), width=binsize, fix='center')
            gr_direct <- resize(GRanges(tad_proj[,2]), width=binsize, fix='center')
            gr_multi <- resize(GRanges(tad_proj[,3]), width=binsize, fix='center')
            gr <- list(ref=gr_ref, qry_direct=gr_direct, qry_multi=gr_multi)
            return(gr)
        }
    }
}, mc.cores=min(20,nthreads))
names(gr) <- tads

# remove tads with no direct-projection coordinates within the tad (print how many) and keep those that have them at the borders (17).
valid_tads <- names(gr)[which(!is.na(gr))]
print(sprintf('removing %s/%s (%.1f%%) tads due to missing direct projections within the tad (inversions / breaking points)', length(gr) - length(valid_tads), length(gr),
              (length(gr)-length(valid_tads))/length(gr)*100))
tads <- valid_tads
gr <- gr[valid_tads]
gr_all_tads <- sapply(c('ref','qry_direct','qry_multi'), function(x) suppressWarnings(do.call('c', lapply(tads, function(tad) gr[[tad]][[x]]))), simplify=F)

# aggregate 200bp bigwig data to tad coordinates
# ~ 7 sec per 100 tads
pc <- exp(-5)
counts_all_tads <- suppressWarnings(sapply(features, function(f) {
  res <- list(ref=ehmm:::aggScore(gr_all_tads[['ref']], whole_genome_count_data[['danRer10']][[f]], 'mean'),
              qry_direct=ehmm:::aggScore(gr_all_tads[['qry_direct']], whole_genome_count_data[['mm10']][[f]], 'mean'),
              qry_multi=ehmm:::aggScore(gr_all_tads[['qry_multi']], whole_genome_count_data[['mm10']][[f]], 'mean'))
  for (x in names(res)) res[[x]]$score <- res[[x]]$score + pc # add PSEUDOCOUNT
  return(res)
}, simplify=F))

# quantile normalize query to reference
for (f in names(counts_all_tads)) {
  names(mcols(counts_all_tads[[f]]$qry_direct)) <- names(mcols(counts_all_tads[[f]]$qry_multi)) <- c('score_cpm')
  counts_all_tads[[f]]$qry_direct$score <- ehmm:::quantileNormalizeToReference(cm.reference=matrix(counts_all_tads[[f]]$ref$score, nrow=1), cm.query=matrix(counts_all_tads[[f]]$qry_direct$score_cpm, nrow=1))$cm.query.normalized
  counts_all_tads[[f]]$qry_multi$score <- ehmm:::quantileNormalizeToReference(cm.reference=matrix(counts_all_tads[[f]]$ref$score, nrow=1), cm.query=matrix(counts_all_tads[[f]]$qry_multi$score_cpm, nrow=1))$cm.query.normalized
}

# save normalized counts by tad to counts_list (separate entry for every tad with normalized counts of all features as mcols)
# ~ 18 sec per 100 tads
tad_labels <- do.call('c', lapply(tads, function(tad) rep(tad, length(gr[[tad]]$qry_multi))))
for (x in c('ref','qry_direct','qry_multi')) gr_all_tads[[x]]$tad <- tad_labels # save the tad information in gr_all_tads
counts_split <- sapply(names(counts_all_tads), function(f) {
  sapply(c('ref','qry_direct','qry_multi'), function(x) split(counts_all_tads[[f]][[x]], tad_labels))
}, simplify=F)

counts_list <- sapply(names(counts_all_tads), function(f) {
  res <- mclapply(tads, function(tad) {
    sapply(c('ref','qry_direct','qry_multi'), function(x) {
      gr <- gr[[tad]][[x]]
      gr$score <- counts_split[[f]][[x]][[tad]]$score
      return(gr)
    }, simplify=F)
  }, mc.cores=min(20,nthreads))
  names(res) <- tads
  return(res)
}, simplify=F)

# create counts data frame that can be used for plotting
counts_df <- do.call('rbind', lapply(names(counts_all_tads), function(f) {
  df <- data.frame(sapply(c('ref','qry_direct','qry_multi'), function(x) {
    counts_all_tads[[f]][[x]]$score
  }))
  df$feature <- factor(f)
  df$tad <- factor(tad_labels)
  df$ref_chrom <- factor(as.vector(seqnames(gr_all_tads$ref)))
  df$ref_coord <- ceiling((start(gr_all_tads$ref) + end(gr_all_tads$ref)) / 2)
  return(df)
}))
counts_df <- counts_df[,c('feature','tad','ref_chrom','ref_coord','ref','qry_direct','qry_multi')]

# compute quantiles
# ~ 10 sec
quantiles <- mclapply(features, function(f) {
  sapply(c('ref','qry_direct','qry_multi'), function(x) {
    sapply(seq(0,1,.001), function(y) quantile(counts_all_tads[[f]][[x]]$score, y))
  }, simplify=F)
}, mc.cores=min(2,nthreads))
names(quantiles) <- features

# write a table with the quantile values for zf and mm
# ~ 7 sec
df_counts <- pivot_longer(counts_df, names_to='spcs', values_to='value', -c(feature,tad,ref_chrom,ref_coord))
df_counts$spcs <- factor(df_counts$spcs, levels=c('ref','qry_direct','qry_multi'))
df_quantiles <- df_counts
df_quantiles$value <- mclapply(seq_len(nrow(df_counts)), function(i) {
  as.numeric(gsub('%', '', names(quantiles[[df_counts$feature[i]]][[df_counts$spcs[i]]])[which.min(abs(quantiles[[df_counts$feature[i]]][[df_counts$spcs[i]]] - as.numeric(df_counts$value[i])))]))/100
}, mc.cores=min(20,nthreads))
qntls <- sapply(spread(df_quantiles, spcs, value)[c('ref','qry_multi')], as.numeric)
quantiles_df <- cbind(spread(df_quantiles, spcs, value)[,c('feature','tad','ref_chrom','ref_coord')],
                      data.frame(quantile_zf=qntls[,'ref'], quantile_mm=qntls[,'qry_multi']))
print(args[8])
write.table(quantiles_df, args[8], sep='\t', quote=F, col.names=T, row.names=F)

# save reference and mapped target H3K27me3 values of all TADs to bigwig files for plotting
counts_to_gr <- function(spcs, f) {
    GenomicRanges::shift(data.frame(df_counts[df_counts$feature==f & df_counts$spcs==spcs, c('ref_chrom','ref_coord','ref_coord','value')]) %>%
                         dplyr::rename(seqnames=ref_chrom, start=ref_coord, end=ref_coord.1, score=value) %>% GRanges() %>% resize(width=1000, fix='center'), 1)
}

ref <- 'danRer10'
qry <- 'mm10'
for (f in features) {
    outfile_ref <- gsub('.bw', '1kb_bins.bw', files[[ref]][[f]])
    outfile_qry <- gsub('.bw', '1kb_bins.bw', files[[qry]][[f]])
    gr_reference <- counts_to_gr('ref', f)
    gr_target_direct <- counts_to_gr('qry_direct', f)
    gr_target_multi <- counts_to_gr('qry_multi', f)
    seqlengths(gr_reference) <- seqlengths(gr_target_direct) <- seqlengths(gr_target_multi) <- sizes[[ref]][seqlevels(gr_reference)]
    ## TAD definitions can overlap. select non-duplicated entries (their values of duplicates are the same)
    gr_reference <- gr_reference[!duplicated(gr_reference)]
    gr_target_multi <- gr_target_multi[!duplicated(gr_target_multi)]
    export.bw(gr_reference, outfile_ref)
    export.bw(gr_target_multi, outfile_qry)
}
