### For splitting up DAM files into up vs. down accumulated
### Also for putting canopus annotations onto the DAMs

argg <- commandArgs(T)

if (length(argg) != 3){
  stop('ARGS: 1) path to DAMs 2) path to canopus_summary filtered
       3) desired output file')
}

get_up_vs_down <- function(damfiles){
  olist <- list()
  for (fil in damfiles){
    hmm <- unlist(strsplit(fil, '_diff_'))[2]
    cond <- gsub('.tab', '', hmm, fixed = T)
    fildf <- read.table(fil, header = T, sep = '\t', fill = NA, quote = "")
    condmean <- rowMeans(fildf[,which(grepl(cond, colnames(fildf)))])
    ctrlmean <- rowMeans(fildf[,which(grepl('Ctrl', colnames(fildf)))])
    fildf$ctrlmean <- ctrlmean
    fildf$treatmean <- condmean
    fildf$condition <- cond
    fildf <- fildf[, -which(grepl(cond, colnames(fildf)))]
    fildf <- fildf[, -which(grepl('Ctrl', colnames(fildf)))]
    fildf$direction <- ifelse(fildf$nfold_change > 0, 'UP', 'DOWN')
    olist[[cond]] <- fildf
  }
  out <- do.call(rbind.data.frame, olist)
  return(out)
}

damfiles <- list.files(argg[1], pattern = 'FDR', full.names = T)
canopus <- read.table(argg[2], header = T, fill = NA, quote = "", sep = '\t')
canopus_to_merge <- canopus[, c('name', 'class', 'superclass', 'subclass', 'level.5')]

conddf <- get_up_vs_down(damfiles)
odir <- argg[3]
canopus_dap <- merge(conddf, canopus_to_merge, all.x = T, by.x = 'Alignment.ID', by.y = 'name')

## Sorting first on condition then direction
canopus_dap <- canopus_dap[with(canopus_dap, order(condition, direction)),]
write.table(canopus_dap, file = paste0(odir, '_all_updown_DAM_canopus.tsv'),
            sep = '\t', quote = F, row.names = F)

print('Done!')
