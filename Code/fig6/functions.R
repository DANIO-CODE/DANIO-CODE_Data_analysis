library(rtracklayer); library(motifcounter); library(ggplot2); library(RColorBrewer); library(dplyr); library(ehmm); library(depmixS4); library(fANCOVA)
library(TxDb.Hsapiens.UCSC.hg19.knownGene); library(TxDb.Ggallus.UCSC.galGal6.refGene)

# function for computing the percentage overlap of two ranges, i.e. what fraction of each region in gr1 is overlapped by regions of gr2
percent_overlap <- function(gr1, gr2) {
  ov <- suppressWarnings(findOverlaps(gr2, gr1))
  overlaps <- pintersect(gr2[queryHits(ov)], gr1[subjectHits(ov)])
  percentOverlap <- width(overlaps) / width(gr1[subjectHits(ov)])
  gr1$percentOverlap <- 0
  gr1$percentOverlap[subjectHits(ov)] <- percentOverlap
  return(gr1$percentOverlap)
}

# function to load repeats from '/project/wig/tobias/reg_evo/data/repeats'
get_repeats <- function(spcs, repeat_dir='/project/wig/tobias/reg_evo/data/repeats') {
  if (!dir.exists(repeat_dir)) stop('repeat_dir does not exist')
  gr <- GRanges()
  repeats_file <- file.path(repeat_dir, sprintf('repeats_%s.bed', spcs))
  if (file.exists(repeats_file)) {
    gr <- import.bed(repeats_file)
  }
  gr <- tryCatch(keepStandardChromosomes(gr, pruning.mode='coarse'),
                 error=function(cond) return(gr))
  return(gr)
}

# function to retrieve exons from TxDb if available
get_exons <- function(spcs) {
  gr <- tryCatch({
    txdb <- supportedOrganisms()$TxDb[which(grepl(spcs, supportedOrganisms()$TxDb))]
    library(txdb, character.only=T)
    exons(eval(as.name(txdb)))
    }, error=function(cond) return(GRanges()))
  gr <- tryCatch(keepStandardChromosomes(gr, pruning.mode='coarse'),
                 error=function(cond) {
                   message(cond)
                   return(gr)
                   })
  return(gr)
}

# function to clean up seqlevels, especially for spotted gar (lepOcu1), where some net.axt file come with a different style (e.g. AHA12345.1 instead of chrUn_AHA12345)
tidy_seqlevels <- function(s) {
  s <- gsub('MT', 'M', gsub('\\.1', '', gsub('AHA', 'chrUn_AHA', gsub('JH', 'chrUn_JH', s))))
  s[!startsWith(s, 'chr')] <- paste0('chr', s[!startsWith(s, 'chr')])
  return(s)
}

# order coordinates: if start coordinate is bigger than end coordinate, switch them.
# returns a matrix with start coords in first row, end coords in second.
order_coords <- function(s,e) {
  coords <- sapply(seq_along(s), function(j) {
    if (s[j] > e[j]) {
      return(c(e[j], s[j]))
    } else {
      return(c(s[j], e[j]))
    }
  })
  rownames(coords) <- c('start','end')
  return(coords)
}

# prepare dataframes for plotting
prepare_plot_data <- function(proj_list, cne_list, pr) {
  # here, proj_list, cne_list are for one feature and not a list of lists for all features
  data_list <- mclapply(seq_along(proj_list), function(i) {
    # get projected coordinates
    proj <- proj_list[[i]]
    grb_id <- factor(paste0('GRB ', names(proj_list)[i], ': ', grbs_ordered@first$target_gene[i]))
    x <- start(resize(proj@first, width=1, fix='center'))
    proj_df <- rbind(data.frame(species='zebrafish prim-5', x=x, y=proj@first$score),
                data.frame(species='mouse embryo forebrain E10.5', x=x, y=proj@second$score))
    proj_df$grb_id <- grb_id
    # get promoter coordinates
    pr_i <- list(x=pr$x[overlapsAny(pr$x, proj@first)], y=pr$y[overlapsAny(pr$y, proj@second)])
    pr_bins <- list(x=which(overlapsAny(proj@first, pr_i$x)), y=which(overlapsAny(proj@second, pr_i$y)))
    if (length(pr_bins$x) > 0) {
      pr_df <- rbind(data.frame(x=x[pr_bins$x], y=-.1, species='zebrafish prim-5'),
                     data.frame(x=x[pr_bins$y], y=-.2, species='mouse embryo forebrain E10.5'))
      pr_df$grb_id <- grb_id
    } else {
      pr_df <- data.frame()
    }
    # get CNE coordinates
    cne_y <- list(xb=-.3, xy=-.4)
    cne_label <- list(xb='CNEs zf - sg', xy='CNEs zf - mm')
    cne_df <- do.call('rbind', lapply(names(cne_y), function(pair) {
      data.frame(x=start(resize(cne_list[[i]][[pair]]@first, width=1, fix='center')), y=cne_y[[pair]], label=cne_label[[pair]], grb_id=grb_id)
    }))
    return(list(proj_df=proj_df, pr_df=pr_df, cne_df=cne_df))
  }, mc.cores=10)
  proj_data <- do.call('rbind', lapply(data_list, function(x) x$proj_df))
  pr_data <- do.call('rbind', lapply(data_list, function(x) x$pr_df))
  cne_data <- do.call('rbind', lapply(data_list, function(x) x$cne_df))

  # normalize again by maximum value per species
  proj_data[proj_data$species=='zebrafish prim-5', 'y'] <- proj_data[proj_data$species=='zebrafish prim-5', 'y'] / max(proj_data[proj_data$species=='zebrafish prim-5', 'y'])
  proj_data[proj_data$species=='mouse embryo forebrain E10.5', 'y'] <- proj_data[proj_data$species=='mouse embryo forebrain E10.5', 'y'] / max(proj_data[proj_data$species=='mouse embryo forebrain E10.5', 'y'])

  return(list(proj_data=proj_data, pr_data=pr_data, cne_data=cne_data))
}

# compute projections of 1kb-binned GRBs and load HM coverages
# also save all CNEs per GRB
get_coverage_and_cnes <- function(grbs_ordered, cnes, files, feature, species1, species2, bridging_species, window_size=1000, ncores=20) {
  grb_ids <- grbs_ordered@first$id
  res_list <- mclapply(grb_ids, function(grb_id) {
    res <- suppressWarnings(project_grb(grbs_ordered, cnes, grb_id=grb_id, species1=species1, species2=species2, bridging_species=bridging_species, window_size=window_size))
    proj <- res[[1]]$proj$bridge$xy # only one result because every grb_id is unique (there might only be multiple results if `gene_of_interest` is passed)
    proj_range <- list(x=reduce(proj@first, min.gapwidth=1e5),
                       y=reduce(proj@second, min.gapwidth=1e5))
    proj_range <- lapply(proj_range, function(x) resize(x, width=width(x)+window_size, fix='center'))
    h <- list(x=import.bw(files$x[feature], which=proj_range$x),
              y=import.bw(files$y[feature], which=proj_range$y))
    proj <- GRangePairs(resize(proj@first, width=window_size, fix='center'), resize(proj@second, width=window_size, fix='center'))
    proj <- GRangePairs(ehmm:::aggScore(proj@first, h$x, 'mean'), ehmm:::aggScore(proj@second, h$y, 'mean'))
    cne <- res[[1]]$cne
    return(list(proj=proj, cne=cne))
  }, mc.cores=ncores)
  names(res_list) <- grb_ids
  proj_list <- lapply(res_list, function(x) x$proj)
  cne_list <- lapply(res_list, function(x) x$cne)
  return(list(proj_list=proj_list, cne_list=cne_list))
}

# function for plotting the coordinate mapping for the GRB containing a gene of interest
plot_projection <- function(proj, cnes_xy, gene_of_interest, grb_i) {
  data <- rbind(data.frame(x=start(proj$bridge$xy@first), y=start(proj$bridge$xy@second), model='bridging species'),
                data.frame(x=start(proj$direct$xy@first), y=start(proj$direct$xy@second), model='direct'))
  data$chrx <- unique(seqnames(proj$bridge$xy@first))
  data$chry <- unique(seqnames(proj$bridge$xy@second))
  data_cne <- data.frame(x = (start(cnes_xy@first) + end(cnes_xy@first)) / 2,
                         y = (start(cnes_xy@second) + end(cnes_xy@second)) / 2)
  pl <- ggplot() +
    geom_point(data=data, aes(x, y, color=model), shape=20, size=.01) +
    geom_point(data=data_cne, aes(x,y), shape=1) +
    facet_grid(chry~chrx, space='free', switch='both') +
    scale_color_manual(values=c(`bridging species`='red', direct='black')) +
    scale_x_continuous(labels=function(z) round(z/1e6,1),
                       expand=expand_scale(mult=0, add=0)) +
    scale_y_continuous(labels=function(z) round(z/1e6,1),
                       expand=expand_scale(mult=0, add=0)) +
    labs(title=do.call('paste0', list('GRB ', grb_i@first$id, ': ', grb_i@first$target_gene)),
         x='danRer10 coordinate [Mbp]',
         y='mm10 coordinate [Mbp]') +
  #   coord_fixed(ratio = .1) +
    theme_minimal() +
    theme(strip.placement='outside',
          strip.text=element_text(size=10))
  ggsave(sprintf('../EpigeneticConservation/grb_mapping_%s.png', toupper(gene_of_interest)), pl, width=8, height=5)
  return(pl)
}

# function for mapping the GRB of a gene of interest
project_grb <- function(grbs, cnes, gene_of_interest=NA, grb_id=NA, species1, species2, bridging_species, window_size=1000) {
  if ((is.na(gene_of_interest) & is.na(grb_id)) | (!is.na(gene_of_interest) & !is.na(grb_id))) {
    stop('One (not both) of the arguments `gene_of_interest` or `grb_id` must be specified')
  }
  if (is.na(grb_id)) {
    grbs_of_interest <- grbs[grepl(toupper(gene_of_interest), toupper(grbs@first$target_gene))]
  } else {
    grbs_of_interest <- grbs[which(grbs@first$id == grb_id)]
  }
  res <- sapply(seq_along(grbs_of_interest), function(i) {
    grb_i_tiled_x <- resize(tile(grbs_of_interest@first[i], width=window_size)[[1]], width=1, fix='center') # tile() will adjust the window_size so that the tiles are more or less equally sized
    # use this elaborate approach instead of tile() because I want every tile to be exactly the size of window_size (tile() would adjust the window_size so that the tiles are more or less equally sized)
    # this is important as I project the 1bp window centers, and later resize them to window_size. With tile, this could result in partially overlapping windows.
    strt <- seq(start(grbs_of_interest[i]@first), end(grbs_of_interest[i]@first), window_size)
    coord <- strt + window_size/2
    grb_i_tiled_x <- GRanges(unique(seqnames(grbs_of_interest[i]@first)), IRanges(coord,coord))
    res_i <- sapply(c('bridge','direct'), function(mode) project_coordinate(grb_i_tiled_x, grbs_of_interest[i], cnes, mode=mode, species1, species2, bridging_species), simplify=F)
    proj <- sapply(res_i, function(x) x$proj, simplify=F)
    cne <- unlist(sapply(res_i, function(x) sapply(x$cne, function(y) y, simplify=F), simplify=F)) # cnes between the direct and bridging species for the grb of interest
    models <- unlist(sapply(res_i, function(x) sapply(x$models, function(y) y, simplify=F), simplify=F), recursive=F) # cnes between the direct and bridging species for the grb of interest
    names(cne) <- names(models) <- c('xb', 'by', 'xy')
    pl <- plot_projection(proj, cne$xy, gene_of_interest, grbs_of_interest[i])
    return(list(proj=proj, cne=cne, models=models, pl=pl))
  }, simplify=F)
  return(res)
}

# function to fit LOESS models using direct CNEs or a bridging species
project_coordinate <- function(gr, grb_x, cnes, mode=c('bridge','direct'), species1, species2, bridging_species) {
  # gr is a GRanges object containing the genomic locations to be projected (width of regions must be 1)
  # grb_x is a GRangePair object containing the grb of interest in species1 and species2.
  # caution: grb_x@first must refer to species1, grb_x@second to species2.
  if (any(width(gr) != 1)) stop('gr must contain genomic locations of width 1')
  if (length(unique(seqnames(gr))) != 1) stop('Ranges in gr must be on the same chromosome')
  # fit LOESS models
  if (mode=='bridge') {
    cne <- list(xb=get_GRangePairs(grb_x@first, cnes[[species1]][[bridging_species]]@first), # zf - gar
                by=invert_GRangePairs(get_GRangePairs(grb_x@second, cnes[[species2]][[bridging_species]]@first))) # gar - mouse7
    # remove CNE outliers, i.e. CNEs that are unusually far from the median position of CNEs in GRB in only one of both species
    cne$xb <- remove_CNE_outliers(cne$xb, grb_x@first)
    models <- lapply(list(xb=cne$xb,by=cne$by), function(grp) fit_loess(grp))
    y1 <- sapply(start(gr), function(x) predict(models$xb, data.frame(x=x), se=F))
    y2 <- sapply(y1, function(y_i) predict(models$by, data.frame(x=y_i), se=F))
    # coordinates outside the CNE range won't be extrapolated and return NA. remove them from all ranges.
    idx <- which(!is.na(y1) & !is.na(y2))
    y1 <- y1[idx]
    y2 <- y2[idx]
    gr <- gr[idx]
    gr_intermediate <- GRanges(seqnames=unique(seqnames(cne$by@first)), IRanges(y1,y1))
    gr_out <- GRanges(seqnames=unique(seqnames(cne$by@second)), IRanges(y2,y2))
    proj <- list(xb=GRangePairs(gr, gr_intermediate),
                 by=GRangePairs(gr_intermediate, gr_out),
                 xy=GRangePairs(gr, gr_out))
  } else if (mode=='direct') {
    cne <- list(xy=get_GRangePairs(grb_x@first, cnes[[species1]][[species2]]@first)) # zf - mouse
    models <- list(xy=fit_loess(cne$xy))
    y <- sapply(start(gr), function(x) {
      predict(models$xy, data.frame(x=x), se=F)
    })
    idx <- which(!is.na(y))
    y <- y[idx]
    gr <- gr[idx]
    gr_out <- GRanges(seqnames=unique(seqnames(cne$xy@second)), IRanges(y,y))
    proj <- list(xy=GRangePairs(gr, gr_out))
  }
  return(list(proj=proj, cne=cne, models=models))
}

# function to identify and remove outliers by looking at distance to median position in both X and Y. If a CNE is very far in only one genome, then it is identified as an outlier
# threshold: abs(log10((y - med(y)) / (x - med(x)))) = .9 --> 8-fold difference between genomes of distances from CNE to median position of all CNEs.
remove_CNE_outliers <- function(cne_pairs, grb_x_i) {
  cne_pairs_i <- get_GRangePairs(grb_x_i, cne_pairs@first) # this already gets rid of CNEs that don't map to the major chromosome
  if (length(cne_pairs_i) == 0) {
    return(GRangePairs())
  } else {
    center_x <- (start(cne_pairs_i@first) + end(cne_pairs_i@first)) / 2
    center_y <- (start(cne_pairs_i@second) + end(cne_pairs_i@second)) / 2
    dist_x <- abs(center_x - median(center_x)) # distance between CNE and median position of all CNEs
    dist_y <- abs(center_y - median(center_y))
    dist_y <- abs(center_y - median(center_y))
    buffer <- 2e+4
    dist_ratio <- (dist_y+buffer) / (dist_x+buffer)
    idx_within_grb <- which(abs(log10(dist_ratio)) < .9) # I had to set it to .9 instead of 1 because of GRB 581 with two outliers with a dist_ratio of 9.6-fold (14 MB from GRB median position)
    return(cne_pairs_i[idx_within_grb])
  }
}

# function to invert CNE GRanges (ranges to name, name to ranges)
invert_CNE_GRanges <- function(cne) {
  cne_inv <- GRanges(cne$name)
  cne_inv$name <- paste0(seqnames(cne), ':', start(cne), '-', end(cne))
  return(cne_inv)
}

# function to invert GRangePairs (first to second, second to first)
invert_GRangePairs <- function(grp) {
  ordr <- order(grp@second)
  grp_inv <- GRangePairs(grp@second[ordr], grp@first[ordr]) # sort according to the new 'first'
  return(grp_inv)
}

# function to create GRangePairs object for CNEs of a given GRB
# first species of CNEs and GRB has to be the same (e.g. GRB danRer-cIde and CNEs danRer-calMil)
get_GRangePairs <- function(grb_x, cne) {
  cne_x <- sort(cne[overlapsAny(cne, grb_x)])
  if (length(cne_x) == 0) {
    return(GRangePairs())
  }
  grp <- GRangePairs(cne_x, GRanges(cne_x$name))
  tbl <- table(seqnames(grp@second))
  major_chr <- names(tbl)[which.max(tbl)] # filter for CNEs that point to the majority chromosome
  grp <- grp[seqnames(grp@second) == major_chr]
  return(grp)
}

# function to fit LOESS
fit_loess <- function(grp) {
#   data <- data.frame(x = (start(grp@first) + end(grp@first)) / 2,
#                      y = (start(grp@second) + end(grp@second)) / 2)
#   model <- loess(y~x, data=data, span=.3, degree=1, control=loess.control(surface='direct'))
  x <- (start(grp@first) + end(grp@first)) / 2
  y <- (start(grp@second) + end(grp@second)) / 2
  model <- loess.as(x, y, degree=1, criterion='gcv', user.span=NULL, plot=F) # automated parameter selection process (i.e. the value for `span`). GCV = generalized cross validation.
  return(model)
}

# function to compute two-dimensional CNE-counts relative to two genomes
# the function takes binned genomes as variables x and y, and the CNEs in cne_xy
compute_CNE_matrix <- function(x, y, cne_xy, ncores=30) {
  m <- do.call('rbind', mclapply(seq_along(x), function(i) {
    m_i <- rep(NA,length(y))
    ov_x_i <- findOverlaps(x[i], cne_xy)
    ov_y_i <- findOverlaps(y, GRanges(cne_xy[subjectHits(ov_x_i)]$name))
    cne_count_y_i <- table(queryHits(ov_y_i))
    m_i[as.integer(names(cne_count_y_i))] <- cne_count_y_i
    return(m_i)
  }, mc.cores=ncores))
  return(m)
}

# 2-dimensional euclidian between two points p and q
euclidian <- function(px,py,qx,qy){
  return(sqrt((px-qx)^2 + (py-qy)^2))
}

# call GRBs based on 2D CNE data
call_2D_GRB <- function(gr_x, gr_y, cnes) {
  m <- compute_CNE_matrix(gr_x, gr_y, cnes, ncores=20)
  coords <- data.frame(which(!is.na(m), arr.ind=T))
  grb_x <- GRanges()
  grb_y <- GRanges()
  i <- 1
  while (T) {
    neighbours <- coords[sample(seq_len(nrow(coords)),1),]
    while (T) {
      neighbours_new <- rbind(neighbours, coords[which(euclidian(min(neighbours$row), min(neighbours$col), coords$row, coords$col) <= 5 |
                                                       euclidian(max(neighbours$row), max(neighbours$col), coords$row, coords$col) <= 5),])
      neighbours_new <- neighbours_new[order(neighbours_new$row),]
      neighbours_new <- neighbours_new[!duplicated(neighbours_new),]
      if (all(dim(neighbours_new) == dim(neighbours))){
        if (all(neighbours_new == neighbours)){
          break
        }
      }
      neighbours <- neighbours_new
    }
    grb_x <- c(grb_x, GRanges(seqnames=seqnames(gr_x[min(neighbours$row)]), IRanges(start=start(gr_x[min(neighbours$row)]), end=end(gr_x[max(neighbours$row)])), id=i))
    grb_y <- c(grb_y, GRanges(seqnames=seqnames(gr_y[min(neighbours$col)]), IRanges(start=start(gr_y[min(neighbours$col)]), end=end(gr_y[max(neighbours$col)])), id=i))
    coords <- dplyr::setdiff(coords, neighbours) # remove previously used coordinates.
    i <- i + 1
    if (nrow(coords) == 0){
      break
    }
  }
  grb_x$name <- paste0(seqnames(grb_y), ':', start(grb_y), '-', end(grb_y))
  grb_y$name <- paste0(seqnames(grb_x), ':', start(grb_x), '-', end(grb_x))
  return(list(x=sort(grb_x), y=sort(grb_y)))
}

# function to plot CNE densities in 2D
plot_CNEs_2d <- function(m,x,y,genomeX,genomeY,xmin,xmax,ymin,ymax) {
  data <- melt(m[xmin:xmax, ymin:ymax])
  data$chrx <- as.vector(rep(seqnames(x[xmin:xmax]), times=ymax-ymin+1))
  data$chry <- as.vector(rep(seqnames(y[ymin:ymax]), each=xmax-xmin+1))

  options(repr.plot.width=8, repr.plot.height=5)
  pl <- ggplot(data, aes(Var1,Var2, fill=value)) +
    geom_raster() +
    facet_grid(chry~chrx, space="free", scales="free", switch="both") +
    scale_fill_gradient2(low='white', high='red', na.value='black') +
    scale_x_continuous(labels=function(z) genomic_location_labels(z,x,xmin,xmax),
                       expand=expand_scale(mult=0, add=0)) +
    scale_y_reverse(labels=function(z) genomic_location_labels(z,y,ymin,ymax),
                    expand=expand_scale(mult=0, add=0)) +
    labs(x=sprintf('%s coordinate [Mbp]', genomeX),
         y=sprintf('%s coordinate [Mbp]', genomeY)) +
    theme_minimal() +
    theme(strip.placement='outside',
          strip.text=element_text(size=10),
          axis.ticks = element_line(size=.25),
          panel.border=element_blank(),
          panel.grid=element_blank())
  return(pl)
}

# function to return genome coordinates as tick labels (in MB entities)
# this is only complicated because facet_grid splits the plot for multiple chromosomes, and ggplot adds a 'NA' break in between.
genomic_location_labels <- function(breaks,gr,gr_min,gr_max) {
  labels <- sapply(breaks, function(z) {
    if (is.na(z)) {
      return(NA)
    } else {
      return(round((start(gr[gr_min+z])+end(gr[gr_min+z]))/2/1e6, 1))
    }
  })
  return(labels)
}

# tile GRBs in motif-sized tiles and assign motif scores
assign_scores <- function(grb, scores) {
  n <- sum(!is.na(scores))
  motif_length <- width(grb) - n
  start <- seq_len(n) + start(grb) - 1
  grb_tiled <- resize(GRanges(seqnames=seqnames(grb), IRanges(start, start+motif_length)), width=9, fix='start')
  grb_tiled$score <- scores[!is.na(scores)]
  return(grb_tiled)
}

# function to tile GRB into windows of a given width. GRBs are first resized to a multiple of the desired window-size in order for tiling to produce windows of the exactly right size.
get_windows <- function(grb, width) {
  windows <- GRangesList(mclapply(seq_along(grb), function(i) {
    grb_i <- resize(grb[i], width=as.integer(width(grb[i]) / width) * width, fix='center') # resize to a multiple of desired width, so that tiling works exact
    windows_i <- unlist(tile(grb_i, width=width))
    return(windows_i)
  }, mc.cores=50))
  names(windows) <- grb$grb_id
  return(windows)
}

# function to remove CNEs which map to an uncommon chromosome compared to most of the CNEs of that GRB
remove_cne_uncommon_chr <- function(gr_by_grb) {
  grl <- mclapply(gr_by_grb, function(x) {
    chrs <- sapply(strsplit(x$name, ':'), function(x) x[[1]][1])
    most_common_chr <- names(which.max(table(chrs)))
    x <- x[chrs == most_common_chr]
    return(x)
  }, mc.cores=40)
  # remove any GRBs that still contain overlapping CNEs
  grl <- grl[which(unname(mclapply(cnes_by_grb, function(x) length(findOverlaps(x, x)) == length(x), mc.cores=40)) == T)]
  grl <- GRangesList(grl)
  return(grl)
}

# function to read CNEs
read_cnes <- function(name) {
  gr <- import.bed(as.character(specs[name,'filename']))
  gr$id <- 1:length(gr)
  gr <- sortSeqlevels(gr)
  gr <- sort(gr)
  return(gr)
}

# function to put CNEs in a correct order with respect to possible inversions of subsets of CNEs.
# stretches of neighboring CNEs that are on the minor strand must be reversed in their order compared to the reference genome.
# example:
# 1 2 3 4 5 6 --> 1 [3 2] 4 5 6
# + - - + + +
# correction: stretches of neighboring CNES must be SORTED instead of reversed because of the following rare case of independent inversions next to each other:
# 1 2 3 4 5 6 --> 1 [3 2] [5 4] 6
# + - - - - +
order_inverted_CNEs <- function(cne, inv) {
  if (length(unique(cne$nonRef_strand)) == 1) return(cne) # skip the whole operation if all CNEs are on the same strand
  major_strand <- names(which.max(table(cne$nonRef_strand)))
  
  # add labels to distinguish groups of neighbouring cnes on the same strand
  cne$label <- 1
  for (i in 2:length(cne)) {
    if (cne$nonRef_strand[i] == cne$nonRef_strand[i-1]) { cne$label[i] <- cne$label[i-1]
    } else cne$label[i] <- cne$label[i-1] + 1
  }

  cne_by_label <- split(cne, ~label)
  cne_ordered <- do.call('c', unname(lapply(cne_by_label, function(gr) {
    if (gr$nonRef_strand[1] == major_strand) { return(gr)
    } else return(sort(gr, decreasing=inv)) # I used `rev(gr)` first, but in rare cases there are independent inversions next to each other, e.g. regions must be sorted
  })))

  cne_ordered$label <- NULL
  return(cne_ordered)
}

# function to determine interCNEs per GRB simultaneously in both genomes
define_interCNEs <- function(cnes_x_by_grb, cnes_y_by_grb) {
  grb_ids <- names(cnes_x_by_grb)
  res <- mclapply(grb_ids, function(grb_id) {
    cne_x_i <- cnes_x_by_grb[[grb_id]]
    cne_y_i <- cnes_y_by_grb[[grb_id]]

    interCNEs_i <- lapply(list(cne_x_i, cne_y_i), function(cne) {
      grb <- GRanges(seqnames=seqnames(cne[1]), IRanges(start=min(start(cne)), end=max(end(cne))))
      # problem: CNEs might border in one genome, yielding less interCNEs using `setdiff`.
      # solution: resize CNEs, i.e. cut 1 bp on each side in order to force gaps. select interCNEs on a min. width of N > 2 bp and finally resize them, too.
      cne_trim <- resize(cne, width=width(cne)-2, fix='center')
      inv <- all(cne_trim == sort(cne_trim)) # determine if GRB is inverted --> this only works if all outliers that mess with the order have been removed correctly
      inter_cnes <- sort(setdiff(grb, cne_trim), decreasing=!inv) # decreasing order only if GRB is not inverted, otherwise in increasing order.
      anchor1 <- cne$id[follow(inter_cnes, cne_trim)] # anchoring CNEs preceding interCNEs
      anchor2 <- cne$id[precede(inter_cnes, cne_trim)] # anchoring CNEs following interCNEs
      # add inner-GRB interCNE IDs refering to the order of the interCNEs in the GRB.
      # this helps comparing interCNEs between genomes e.g. for their min. length, and solves the problem of different anchors in case of CNE-subset inversions.
      inter_cnes$id <- seq_along(inter_cnes)
      inter_cnes$anchors <- paste(anchor1, anchor2, sep='-')
      inter_cnes$grb_id <- grb_id
      return(inter_cnes)
    })

    # restrict to interCNEs that are min. N bp in BOTH genomes.
    N <- 10
    idx_minwidth <- do.call('intersect', lapply(interCNEs_i, function(x) which(width(x) > N)))
    interCNEs_i_minwidth <- lapply(interCNEs_i, function(x) x[idx_minwidth])
    interCNEs_i_trim <- lapply(interCNEs_i_minwidth, function(x) resize(x, width=width(x)-2, fix='center')) # CNEs were shrinked by 2bp, so now shrink interCNEs, too.
  }, mc.cores=50)
  return(res)
}

# function to exclude repeats, promoters and exons from interCNEs.
# regions that span over exons will be split at the ends of the repeats/promotres/exons that overlap it.
exclude_RePrEx <- function(interCNEs, re_pr_ex, rdsfile) {
  if (file.exists(rdsfile)) {
    res <- readRDS(rdsfile)
  } else {
    res <- GRangesList(mclapply(interCNEs, function(x) {
      res_x <- suppressWarnings(setdiff(x, re_pr_ex, ignore.strand=T))
      mcols(res_x) <- mcols(x)[subjectHits(findOverlaps(res_x,x)),]
      N <- 10
      res_x <- subset(res_x, width > N)
      return(res_x)
    }, mc.cores=min(50,length(interCNEs))))
    saveRDS(res, rdsfile)
  }
  return(res)
}

### Function to allocate CNEs to given GRBs in genome X, define GRBs in genome Y, and refine GRBs in genome X with respect to genome Y.
### The refinement step involves two steps:
### 1. Identify inversions and put them in correct order
### 2. Identify remaining outliers and remove them from both sets.
### 3. Identify remaining incorrectly ordered CNEs (outliers) and remove them from both sets, thus considering them as non-conserved with regard to the current GRB.
### The 2nd step involves calling 'outliers', i.e. CNEs in GRBs in genome X that map to regions in genome Y far away from the actual GRB.
### Outliers are defined the following: I calculate distances from every CNE to the median position of all CNEs in a GRB. 
### CNEs whose distances differ by more than 10 times between genome X and Y are considered as outliers.
### I remove these 'outlier' CNEs in genome X and thus consider them as non-conserved with regard to the same GRB in genome Y.
###
# explanation for the distance buffer:
# because of the outliers, the median of x and y can be shifted relative to each other which might cause problems at regions close to the median.
# Example: in genome X, the median position might be very close to the position of CNE #19, whereas in genome Y it might be CNE #20.
# If CNE #19 and #20 are not very close to each other, this could give us relatively high differences of distance-to-median between X and Y and thus high ratios.
# This can be buffered by a pseudo-distance, which I choose to be 20 KB.
# Reason: in some cases, some CNEs are within tens of BP to the median and the respective CNE in the other genome tens of KB.

refine_GRBs <- function(grb_x, cne_x, cne_y, minor_strand_limit=.05, outlier_limit=.05){
  ov <- findOverlaps(cne_x, grb_x)
  cne_x$grb <- 'not in GRB'
  cne_x$grb[queryHits(ov)] <- subjectHits(ov)
  cne_x_by_grb <- split(cne_x, cne_x$grb) # split into GRBs
  cne_x_by_grb <- cne_x_by_grb[names(cne_x_by_grb) != 'not in GRB'] # discard regions that are not in a GRB
  grb_ids <- names(cne_x_by_grb)
  
  res <- mclapply(grb_ids, function(grb_id) {
    # 1. identify inversions and put them in correct order
    # also take care of GRBs with > 5 CNEs on the minor strand as this is a potential case of multiple true GRBs
    remove_GRB <- F # initiate flag to FALSE
    cne_x_i <- cne_x_by_grb[[grb_id]]
    n_cnes_original <- length(cne_x_i) # for later comparison with idx_within_grb: if the difference is > 5, set the result to NULL as there are too many outliers.
    cne_ids <- cne_x_i$id
    cne_y_i <- cne_y[match(cne_ids, cne_y$id)]
    cne_y_i$grb <- grb_id
    major_chr_x <- names(which.max(table(seqnames(cne_x_i))))
    major_chr_y <- names(which.max(table(seqnames(cne_y_i))))
    idx_same_chr <- seqnames(cne_x_i) == major_chr_x & seqnames(cne_y_i) == major_chr_y
    cne_x_i <- cne_x_i[idx_same_chr]
    cne_y_i <- cne_y_i[idx_same_chr] # only keep CNEs that in both genomes are on the same chromosome as the majority of the CNEs.
    if ((min(table(cne_x_i$nonRef_strand)) / n_cnes_original) > minor_strand_limit) remove_GRB <- T # too many minor-strand CNEs (more than 5% of all CNES), potential case of multiple true GRBs.
    inv <- names(which.max(table(cne_y_i$nonRef_strand))) == '-' # determine if GRB is inverted on genome Y
    cne_y_i <- order_inverted_CNEs(cne_y_i, inv) # identify subset of inverted CNEs and put them in correct order
    
    # 2. remove outliers called by distance ratio --> 13.11.2019: this has been defined as a function without the outlier_limit in order to use for the epigenetic conservation project.
    # note: in case of GRBs consisting of multiple true GRBs, this will somehow lead to non-matching IDs between genomes, which will be set to NULL and later excluded.
    center_x <- (start(cne_x_i) + end(cne_x_i)) / 2 # center position
    center_y <- (start(cne_y_i) + end(cne_y_i)) / 2
    dist_x <- abs(center_x - median(center_x)) # distance between CNE and median position of all CNEs
    dist_y <- abs(center_y - median(center_y))
    buffer <- 2e+4
    dist_ratio <- (dist_y+buffer) / (dist_x+buffer)
    idx_within_grb <- which(abs(log10(dist_ratio)) < 1)
    if ((length(idx_within_grb) / n_cnes_original) < (1-outlier_limit)) remove_GRB <- T # too many outliers (more than 5% of all CNEs), potential case of multiple true GRBs.
    cne_x_i <- cne_x_i[idx_within_grb]
    cne_y_i <- cne_y_i[idx_within_grb]
    
    # 3. identify remaining incorrectly ordered CNEs and remove them from both sets. (this removes remaining 'outliers' that are not so far away)
    # the algorithm takes the indices of the current order and identifies all indices that are not in a sorted order and removes them from the GRanges object.
    idx <- order(cne_y_i, decreasing=inv) # inv takes care of inverted GRBs
    dump <- c() # index of elements that are not in order
    i <- 1
    while (i <= length(idx)) {
      if (idx[i] != i) {
        dump <- c(dump, idx[i])
        # important: if we e.g. remove index 49 that was in position 40, then the new index at position 40 will be 40. good.
        # however, the new index at postion 49 would now be 50, so we need to decrease indices bigger than the removed index by 1.
        idx[idx > idx[i]] <- idx[idx > idx[i]] - 1
        idx <- idx[-i]
      } else i = i + 1
    }
    cne_x_i <- cne_x_i[idx]
    cne_y_i <- cne_y_i[idx]
    if (remove_GRB) { 
      grb_x_i <- grb_y_i <- GRanges()
    } else {
      grb_x_i <- GRanges(seqnames=major_chr_x, IRanges(start=min(start(cne_x_i)), end=max(end(cne_x_i))), grb_id=grb_id)
      grb_y_i <- GRanges(seqnames=major_chr_y, IRanges(start=min(start(cne_y_i)), end=max(end(cne_y_i))), grb_id=grb_id)
      seqlevels(grb_x_i) <- seqlevels(grb_x)
      seqlevels(grb_y_i) <- seqlevels(cne_y)
    }
    res_i <- list(cne_x_i=cne_x_i, cne_y_i=cne_y_i, grb_x_i=grb_x_i, grb_y_i=grb_y_i)
    if (!all(res_i$cne_x_i$id %in% res_i$cne_y_i$id)) remove_GRB <- T # CNE IDs don't match between genomes
    if (remove_GRB) res_i <- NULL 
    return(res_i)
  }, mc.cores=min(length(grb_ids),50))
  names(res) <- grb_ids
  grb_ids <- grb_ids[!is.null(grb_ids)]
  res <- res[which(sapply(res, function(x) length(x$cne_x_i) >= 10))] # remove GRBs with # CNEs < 10
  cne_x_by_grb <- GRangesList(lapply(res, function(x) x$cne_x_i))
  cne_y_by_grb <- GRangesList(lapply(res, function(x) x$cne_y_i))
  grb_x <- sort(do.call('c', unname(sapply(res, function(x) x$grb_x_i))))
  grb_y <- sort(do.call('c', unname(sapply(res, function(x) x$grb_y_i))), by=~as.numeric(grb_id))
  refinedRegions <- list(cne_x_by_grb=cne_x_by_grb, cne_y_by_grb=cne_y_by_grb, grb_x=grb_x, grb_y=grb_y)
  return(refinedRegions)
}
                                   

get_state_colors <- function(label_categories) {
  # only works on already ordered label categories
  if (sum(label_categories=='both_low') <= 5) {
    both_low <- tail(brewer.pal(9, 'Greys')[3:7], sum(label_categories=='both_low'))
  } else both_low <- colorRampPalette(brewer.pal(9, 'Greys')[3:7])(sum(label_categories=='both_low')) # deal with high number of states
  if (sum(label_categories=='both_signal') <= 3) {
    both_signal <- tail(brewer.pal(9, 'Greens')[c(4,6,8)], sum(label_categories=='both_signal'))
  } else both_signal <- colorRampPalette(brewer.pal(9, 'Greens')[5:8])(sum(label_categories=='both_signal'))
  if (sum(label_categories=='diff_ref') <= 3) {
    diff_ref <- tail(brewer.pal(9, 'Blues')[c(4,6,8)], sum(label_categories=='diff_ref'))
  } else diff_ref <- colorRampPalette(brewer.pal(9, 'Blues')[4:8])(sum(label_categories=='diff_ref'))
  if (sum(label_categories=='diff_qry') <= 3) {
    diff_qry <- tail(brewer.pal(9, 'Reds')[c(4,6,8)], sum(label_categories=='diff_qry'))
  } else diff_qry <- colorRampPalette(brewer.pal(9, 'Reds')[5:8])(sum(label_categories=='diff_qry'))
  state_cols <- c(both_low, diff_ref, diff_qry, both_signal)
  names(state_cols) <- seq_along(state_cols)
  return(state_cols)
}

get_plot_colors <- function() return(brewer.pal(9, 'Set1')[2:1])

plot_domains <- function(df_counts, domains, ref, qry, stage, feature, grb, scaling=c('global','local')) {
  scaling <- match.arg(scaling)
  species_dict <- c(danRer10='Zebrafish', mm10='Mouse', xenTro9='Xenopus')
  labels <- c(ref=species_dict[[ref]], qry_direct=paste0(species_dict[qry],' Direct'), qry_dijkstra=paste0(species_dict[qry],' Dijkstra'))
  df <- df_counts[df_counts$qry==qry & df_counts$feature==feature & df_counts$stage==stage & df_counts$grb==grb & df_counts$spcs %in% c('ref','qry_dijkstra'),]
  dom <- domains[domains$qry==qry & domains$stage==stage & domains$feature==feature & domains$grb==grb,]
  if (scaling=='global') {
    scaling_values <- df_counts[df_counts$qry==qry & df_counts$feature==feature,'value']
  } else if(scaling=='local') {
    scaling_values <- df_counts[df_counts$qry==qry & df_counts$feature==feature & df_counts$grb==grb,'value']
  }
  min_value <- min(scaling_values)
  max_value <- max(scaling_values)
  plot_cols <- get_plot_colors()
  state_cols <- c(brewer.pal(9,'Greys')[3], brewer.pal(9,'YlOrBr')[5], brewer.pal(9,'Greens')[7])
  xticks <- seq(round(df$ref_coord[1], -5), round(rev(df$ref_coord)[1], -5), 5e4)
  pl <- ggplot() +
    geom_tile(data=dom, aes(ref_coord,mean(c(min_value,max_value)),fill=domain), height=(max_value-min_value)*1.1, alpha=.4, show.legend=F) +
    scale_x_continuous(breaks=xticks, labels=xticks/1e6, expand=c(0,0)) +
    scale_fill_manual(values=state_cols, name='') +
    labs(x=sprintf('zebrafish coordinate on %s [MB]', df$ref_chrom[1]),
         y=feature,
         title=sprintf('%s Domains', feature)) +
    geom_path(df, mapping=aes(ref_coord,value,color=spcs), size=1) +
    scale_color_manual(values=plot_cols, name='', labels=labels) +
    theme_classic() +
    theme(axis.title.y=element_text(vjust=2),
          axis.title.x=element_text(vjust=-1),
          axis.line=element_blank(),
          text=element_text(size=18),
          legend.position='bottom')
  return(pl)
}

display_colors <- function(palette, n=length(palette)) {
  return({
    plot(.5,.5, xlim=c(0,n), ylim=c(0,1), asp=1, xlab=NA, ylab=NA, type='n', axes=F)
    for (i in 1:n) rect(i-1,0,i,1,density=1000,col=palette[i])
  })
}
# layout for R plots

theme_Publication <- function(base_size=22, base_family="helvetica", plot_title_hjust=0) { # pass plot_title_hjust=0.5 for centered titles
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust=plot_title_hjust),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(), 
               axis.line = element_line(colour="black"),
               axis.ticks = element_line(),
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               legend.box = 'vertical',
               legend.key = element_rect(colour = NA),
               legend.position = "bottom",
               legend.direction = "horizontal",
               legend.key.size= unit(.5, "cm"),
               legend.margin = margin(unit(0, "cm")),
               legend.title = element_text(face="italic"),
               plot.margin=unit(c(5,5,5,5),"mm"),
#                plot.margin=grid::unit(c(0,0,0,0),'mm'),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
          ))
      
}

color_palette_publication <- c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")

scale_fill_Publication <- function(cols=color_palette_publication, ...){
      library(scales)
      discrete_scale("fill","Publication",manual_pal(values = cols), ...)

}


scale_colour_Publication <- function(cols=color_palette_publication, ...){
      library(scales)
      discrete_scale("colour","Publication",manual_pal(values = cols), ...)
}

set_plot_dimensions <- function(width_choice, height_choice) {
     options(repr.plot.width=width_choice, repr.plot.height=height_choice)
}