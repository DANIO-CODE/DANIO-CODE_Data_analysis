library(BSgenome.Drerio.UCSC.danRer10)
library(tools)
library(jsonlite)

ref_stages <- c("egg", "fertilized_egg", "fertilized_egg.1", "1-cell", "2-cell",
  "4-cell", "8-cell", "16-cell", "32-cell", "64-cell", "128-cell", "256-cell",
  "512-cell", "512-cell.1", "1k-cell", "high", "oblong", "sphere", "dome", "dome.1",
  "dome-30per_epiboly", "30per_epiboly", "50per_epiboly", "germ_ring", "shield",
  "75per_epiboly", "90per_epiboly", "bud", "1-4_somites", "5-8_somites", "10-13_somites",
  "14-19_somites", "20-25_somites", "26+_somites", "prim-5", "prim-5.1", "prim-5.2",
  "prim-15", "prim-25", "high-pec", "long-pec", "pec-fin", "protruding_mouth",
  "day_4", "day_5", "day_6", "days_7-13", "days_14-20", "days_21-29", "days_30-44",
  "days_45-89", "90 days-2_years")


shift_files <- function(assay, scalingFactor, shift, file_stages, file_ext) {
  # shifts bigwig files vertically for joy-division-plot track, 
  # `scalingFactor` can be an number or "log"/"log_neg"
  # `shift` can be an array `(step_size,max_value)` or a single number, which makes `1` the `max_value` and the `step_size` is dependent on the number of stages
  jsonfile <- c()
  max_score <- 0
  files <- list.files(path = paste0("non_shifted/", assay), pattern = file_ext,
    full.names = TRUE)
  file_stages <- make.unique(file_stages)
  names(files) <- file_stages
  stages <- file_stages[order(match(file_stages, ref_stages))]
  files <- files[stages]

  
  if (length(shift) == 1) {
      signalIncrease <- seq(0, length(stages), 1)
  } else {
    signalIncrease <- seq(0, shift[1], shift[2])
  }
  print(paste("=============", "starting", assay, "============="))
  for (i in 1:length(files)) {

    stage_name <- stages[i]
    stage_name <- gsub("\\.\\d", "", stage_name)
    print("start")
    print(basename(files[i]))
    print(stage_name)
    temp <- import(files[i])
    print("imported")
    seqlevels(temp) <- seqlevels(BSgenome.Drerio.UCSC.danRer10)
    seqinfo(temp) <- seqinfo(BSgenome.Drerio.UCSC.danRer10)

    # some assays, have gaps between signal, since the browser interpolates these, 
    # I manually set them to zero to avoid weird interpolation lines
    gap_cond <- substring(assay, 1, nchar("RNA")) == "RNA" || substring(assay,
      1, nchar("CAGE")) == "CAGE" || substring(assay,
        1, nchar("new_CAGE")) == "new_CAGE" ||substring(assay,
          1, nchar("nanti_CAGE")) == "nanti_CAGE" || substring(assay, 1, nchar("4C")) == "4C"
    if (gap_cond) {
      temp_gap <- gaps(temp)[which(strand(gaps(temp)) == "*")]
      temp_gap$score <- rep(0, length(temp_gap))
      comb <- c(temp, temp_gap)
      temp <- sort(comb)
    }
    max_score <- max(max(abs(temp$score)), max_score)
    print(log1p(max_score))
    if (length(shift) == 1) {
      if (shift == "log_neg") {
        temp$score <- log1p(abs(temp$score)) + signalIncrease[i] * scalingFactor -
          20/length(files)

      } else {
        temp$score <- log1p(temp$score) + signalIncrease[i] * scalingFactor -
          20/length(files)
      }
    } else {
      temp$score <- temp$score * 0.1 * scalingFactor + signalIncrease[i]
    }
    print("increased")
    lead_index <- sprintf("%02d", i)
    filename <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(files[i]))
    outFileName <- paste0("shifted/", assay, "/", lead_index, "_", stages[i],
      "_", filename, ".shift.bigwig")

    meta <- data.frame(og_file = paste0("non_shifted/", assay, "/", basename(files[i])),
      stage = stage_name, new_file = outFileName, assay = assay)

    jsonfile <- rbind(jsonfile, meta)
    export(temp, outFileName)
    print("done")
  }

  exportJSON <- jsonlite::toJSON(jsonfile, pretty = TRUE, auto_unbox = TRUE)
  write(exportJSON, paste0("shifted/", assay, "/log.json"))
  return()
}

sort_files <- function(assay, file_stages, file_ext) {
  # sorts files based on their stage from earliest stage to latest
  jsonfile <- c()
  max_score <- 0
  files <- list.files(path = paste0("non_shifted/", assay), pattern = file_ext,
    full.names = TRUE)
  file_stages <- make.unique(file_stages)
  names(files) <- file_stages
  stages <- file_stages[order(match(file_stages, ref_stages))]
  files <- files[stages]

  print(paste("=============", "starting", assay, "============="))
  for (i in 1:length(files)) {

    stage_name <- stages[i]
    stage_name <- gsub("\\.\\d", "", stage_name)
    print("start")
    print(basename(files[i]))
    print(stage_name)
    temp <- import(files[i])
    print("imported")
    seqlevels(temp) <- seqlevels(BSgenome.Drerio.UCSC.danRer10)
    seqinfo(temp) <- seqinfo(BSgenome.Drerio.UCSC.danRer10)
    lead_index <- sprintf("%02d", i)
    filename <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(files[i]))
    outFileName <- paste0("shifted/", assay, "/", lead_index, "_", stages[i],
      "_", filename, ".shift.bed")

    meta <- data.frame(og_file = paste0("non_shifted/", assay, "/", basename(files[i])),
      stage = stage_name, new_file = outFileName, assay = assay)

    jsonfile <- rbind(jsonfile, meta)
    export(temp, outFileName)
    print("done")
  }

  exportJSON <- jsonlite::toJSON(jsonfile, pretty = TRUE, auto_unbox = TRUE)
  write(exportJSON, paste0("shifted/", assay, "/log.json"))
  return()
}


shift_files('Methyl',100,c(50,3), c('prim-5','long-pec','128-cell', '1k-cell',
'32-cell', '64-cell', 'egg', 'germ_ring','256-cell', 'sphere' ), '.bigwig')

shift_files('H3K27AC',1,'log', c('75per_epiboly','5-8_somites', 'shield',
'dome', 'long-pec', 'prim-5', '256-cell'), '.bigwig|.bw')

shift_files('H3K27me3',1,c(10,2), c('shield', 'long-pec', 'dome', 'prim-5',
'5-8_somites'), '.bw|.bigwig')

shift_files('RNA_pos',1,'log', c('prim-5',
'75per_epiboly','long-pec','5-8_somites', 'shield', 'dome', '256-cell'),
'.pos.bw')
shift_files('RNA_neg',1,'log_neg', c('prim-5',
'75per_epiboly','long-pec','5-8_somites', 'shield', 'dome', '256-cell'),
'.neg.bw')
shift_files('ATAC',1,'log', c('long-pec',
'75per_epiboly','5-8_somites','dome', 'shield', 'prim-5', '256-cell') ,
'.bigwig')

shift_files('H3K4me1',5,c(10,1), c('256-cell', 'long-pec',
'prim-5','75per_epiboly', 'dome', '5-8_somites'), '.bw|.bigwig')
shift_files('H3K4me3',15,c(300,30), c('prim-5',
'protruding_mouth','oblong','512-cell',
'long-pec','5-8_somites','75per_epiboly', 'dome', 'shield'), '.bw')
shift_files("4C_six2a",8,c(1000,5), c("prim-5",  "long-pec", "75per_epiboly",
  "dome"), ".bw")

shift_files('CAGE_pos',0.25,'log', c( '14-19_somites', '30per_epiboly',
'512-cell', 'dome-30per_epiboly', 'oblong', 'prim-25', 'prim-5', 'shield' ) ,
'tagging_(Prim5|dome_30pc|30pc_epi|Shield|512_cell|Oblong|14_19_somites|Prim25).*\\.bw')
shift_files('CAGE_neg',0.25,'log_neg', c( '14-19_somites', '30per_epiboly',
'512-cell', 'dome-30per_epiboly', 'oblong', 'prim-25', 'prim-5', 'shield' ) ,
'tagging_(Prim5|dome_30pc|30pc_epi|Shield|512_cell|Oblong|14_19_somites|Prim25).*\\.bw')

shift_files('new_CAGE_pos',0.25,'log', c( '14-19_somites', '30per_epiboly',
'512-cell', 'dome-30per_epiboly', 'oblong', 'prim-25', 'prim-5', 'shield' ) ,
'tagging_(Prim5|dome_30pc|30pc_epi|Shield|512_cell|Oblong|14_19_somites|Prim25).*\\.bw')
shift_files('new_CAGE_neg',0.25,'log_neg', c( '14-19_somites', '30per_epiboly',
'512-cell', 'dome-30per_epiboly', 'oblong', 'prim-25', 'prim-5', 'shield' ) ,
'tagging_(Prim5|dome_30pc|30pc_epi|Shield|512_cell|Oblong|14_19_somites|Prim25).*\\.bw')
shift_files('nanti_CAGE_pos',0.25,'log', c('128-cell', '16-cell', '30per_epiboly', '1-4_somites',
'512-cell','fertilized_egg','prim-5' ),
'nanti_(128_cell|16_cell|30pc_epi|4_somites|512_cell|fertilized_egg|LongPec|Prim5).*\\.bw')
shift_files('nanti_CAGE_neg',0.25,'log_neg', c( '128-cell', '16-cell', '30per_epiboly', '1-4_somites',
'512-cell','fertilized_egg','prim-5' ),
'nanti_(128_cell|16_cell|30pc_epi|4_somites|512_cell|fertilized_egg|LongPec|Prim5).*\\.bw')

sort_files('pdre', c( '5-8_somites', '75per_epiboly',
'dome', 'long-pec', 'prim-5') ,
'.bed')
