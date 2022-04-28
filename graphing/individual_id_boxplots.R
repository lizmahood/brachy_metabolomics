library(ggplot2)
library(reshape2)

get_groups <- function(ct){
  groups <- unlist(lapply(strsplit(ct, '_'), '[', 1))
  return(groups)
}

argg <- commandArgs(T)

if (length(argg) != 4){
  stop('ARGS: 1) Input file of metabolites or transcripts to plot levels of
       2) IDs of interest, separated by comma. 3) metab OR trans 4) 
       what do you want to label the output? Is not the full output name.')
}

features <- read.table(argg[1], sep = '\t', header = T, quote = "", fill = NA)
ids <- unlist(strsplit(argg[2], ',', fixed = T))
if (!(grepl('[A-Za-z]', ids[1]))){
  ids <- as.numeric(ids)
}

if (argg[3] == 'metab'){
  idvar <- 'Alignment.ID'
  if ('Formula' %in% colnames(features)){
    features <- features[, -c(2:32)]
  }
}else{
  idvar <- 'Gene'
  #colnames(features) <- sub('_', '.', colnames(features), fixed = T)
  #colnames(features) <- sub('_', '.', colnames(features), fixed = T)
}

idrows <- which(features[, idvar] %in% ids)
toplot <- features[idrows,]
mtoplot <- melt(toplot, id.vars = idvar)
grps <- get_groups(as.character(mtoplot$variable))
mtoplot$groups <- grps

out_str <- argg[4]
ct <- 1
print(paste0(dirname(argg[1]), '/boxplots_of_', out_str, '.pdf'))
pdf(file = paste0(dirname(argg[1]), '/boxplots_of_', out_str, '.pdf'))
for (id in unique(mtoplot[,1])){
  tmp <- mtoplot[which(mtoplot[,1] == id),]
  print(ggplot(tmp, aes(groups, value, fill = groups)) + geom_boxplot() + theme_bw() +
          theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8)) +
          xlab('') + ylab('') + guides(fill = F) + 
          ggtitle(paste0('Normalized Values for ', id)))
  ct <- ct + 1
}
dev.off()
