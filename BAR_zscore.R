library(tidyverse)

argg <- commandArgs(T)

if (length(argg) != 2){
  stop('ARGS: 1) BAR file 2) folder and prefix for output files')
}

get_groups <- function(df){
  ## df is a dataframe, whose header you want to get groups out of
  return(substr(colnames(df), 1, nchar(colnames(df)) - 2))
}

get_avg_per_group <- function(df, groups){
  ## df is a dataframe
  ## groups is a vector of strings
  tdf <- as.data.frame(t(df))
  tdf$groups <- groups
  
  newdf <- tdf %>%
             group_by(groups) %>%
             summarise_all(mean)
  
  out <- as.data.frame(t(newdf))
  outnames <- as.character(out[1,])
  out <- out[-1,]
  out <- apply(out, 2, as.numeric)
  colnames(out) <- outnames
  
  return(out)
}

inpdf <- read.table(argg[1], sep = '\t', quote = "", header = T, fill = NA)
groups <- get_groups(inpdf[,-1])

avgdf <- get_avg_per_group(inpdf[,-1], groups)

dfmean <- mean(as.numeric(as.matrix(avgdf)))
dfsd <- sd(as.numeric(as.matrix(avgdf)))

zscoredf <- as.data.frame((avgdf - dfmean) / dfsd)

zscoredf$class <- inpdf$class
write.table(zscoredf, file = paste0(argg[2], '_zscore.tsv'), sep = '\t', 
            quote = F, row.names = F)