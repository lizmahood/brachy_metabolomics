## Making heatmap for canopus classes
## Inputs - canopus output and alignment file
library(pheatmap)
library(htmltools)
library(sunburstR)
library(ggplot2)
library(htmltools)
library(matrixStats)
library(RColorBrewer)
library(tidyverse)

argg <- commandArgs(T)

if (length(argg) != 4){
  stop('ARGS: 1) canopus output file 2) Filterd canopus output file 
  3) alignment file 4) prefix of output files')
}


get_groups <- function(df){
  ct <- colnames(df)
  groups <- unlist(lapply(strsplit(ct, '_'), '[', 1))
  return(groups)
}

get_sunburst_data <- function(canopus, ofil){
  tmp <- canopus[,5:9]
  #tmp <- tmp[,c(4,3,2,1)]
  for (col in 1:ncol(tmp)){
    tmp[,col] <- gsub('-', ' ', tmp[,col])
    tmp[,col] <- gsub('[^[:alnum:] ]', '', tmp[,col])
  }
  large <- c()
  for (i in 1:nrow(tmp)){
    row <- as.character(tmp[i,])
    nrow <- c()
    for (j in 1:length(row)){
      if (row[j] != ''){
        nrow <- c(nrow, row[j])
      }
    }
    large <- c(large, paste0(nrow, collapse = '-'))
  }
  table_classes <- table(large)
  numbers <- unname(table_classes)
  classes <- names(table_classes)
  sundata <- data.frame(cbind(classes, numbers))
  save_html(sunburst(sundata), file = paste0(ofil, '_classes_sunburst.html'))
}

make_numbers_heatmap <- function(canopus, align, method){
  align <- align[,-(2:32)]
  merged <- merge(align, canopus, by.x = 'Alignment.ID', by.y = 'name')
  grps <- get_groups(align[,2:ncol(align)])
  grplist <- list()
  for (clas in unique(merged$class)){
    tmp <- merged[which(merged$class == clas),]
    grptmp <- c()
    for (grp in unique(grps)){
      tmp2 <- tmp[,which(grepl(grp, colnames(tmp)))]
      if (method == 'numbers'){
        grptmp <- c(grptmp, length(which(rowMedians(as.matrix(tmp2)) > 0.01)))
      }else if (method == 'abund'){
        grptmp <- c(grptmp, mean(as.matrix(tmp2)))
      }
    }
    names(grptmp) <- unique(grps)
    if (nchar(clas) > 1){
      grplist[[clas]] <- grptmp
    }
  }
  out <- do.call(cbind.data.frame, grplist)
  
  conditions <- unique(grps)
  brokenconds <- strsplit(conditions, '.', fixed = T)
  tissue <- sapply(brokenconds, tail, 1)
  conds <- c()
  for (cond in brokenconds){
    if (length(cond) == 3){
      conds <- c(conds, paste(cond[1], cond[2], sep = '.'))
    }else {
      conds <- c(conds, cond[1])
    }
  }
  sampcol <- data.frame(cbind(conds, tissue))
  row.names(sampcol) <- unique(grps)
  return(list(out, sampcol))
}

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

get_level6 <- function(df){
  df$alltogether <- paste0(df$superclass, df$class, df$subclass, df$level.5, 
                           sep = '_')
  ## only including most specific class if it's a level 6 classification
  df$is_not_level6 <- mapply(grepl, pattern=df$most.specific.class, 
                                     x=df$alltogether)
  df$most.specific.class[df$is_not_level6] <- ''
  df$is_not_level6 <- NULL; df$alltogether <- NULL
  return(df)
}

make_barplot_df <- function(cfilt, cnonfilt){
  levell <- c(); countl <- c(); typel <- c()
  total <- nrow(cnonfilt)
  for (col in c(5:9)){
    tcol <- cfilt[, col]
    levell <- c(levell, colnames(cfilt)[col])
    countl <- c(countl, length(which(tcol != '')) / total)
    typel <- c(typel, 'Filtered')
  }
  ctr <- 1
  for (col in c(8:4)){
    tcol <- cnonfilt[, col]
    levell <- c(levell, colnames(cnonfilt)[col])
    countl <- c(countl, (length(which(tcol != '')) / total))
    typel <- c(typel, 'Nonfiltered')
    ctr <- ctr + 1
  }
  outdf <- data.frame('level' = levell, 'count' = countl, 'type' = typel)
  outdf$level <- factor(outdf$level, levels = c('superclass', 'class', 
                                                'subclass', 'level.5', 
                                                'most.specific.class'))
  return(outdf)
}

canopus_nonfilt <- read.table(argg[1], sep = '\t', header = T, fill = NA, quote = "")
canopusbf <- read.table(argg[2], sep = '\t', header = T, fill = NA, quote = "")
align <- read.table(argg[3], sep = '\t', header = T, fill = NA, quote = "")
ofil <- argg[4]

## Changing all occurrences of "None" to "" in the filtered file
canopus <- canopusbf %>%
  mutate_all(list( ~ str_replace(., "None", "")))

torem <- which(canopus$superclass == '' & canopus$class == '' & canopus$subclass == '' &
                 canopus$level.5 == '' & canopus$most.specific.class == '')
canopus <- canopus[-torem,]

num_heatmap <- make_numbers_heatmap(canopus, align, 'numbers')
data_norm <- apply(num_heatmap[[1]], 2, cal_z_score)
data_norm[is.nan(data_norm)] <- 0

canopus <- get_level6(canopus)
canopus_nonfilt <- get_level6(canopus_nonfilt)
classifications_toplot <- make_barplot_df(canopus, canopus_nonfilt)

## Plotting
pdf(file = paste0(ofil, '_counts_heatmap.pdf'), width = 11)
pheatmap(data_norm, annotation_row = num_heatmap[[2]], show_rownames = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = 'YlGnBu'))) (100))
dev.off()

pdf(file = paste0(ofil, '_abund_heatmap.pdf'), width = 11)
abund_heatmap <- make_numbers_heatmap(canopus, align, 'abund')
pheatmap(abund_heatmap[[1]], annotation_row = abund_heatmap[[2]], show_rownames = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = 'YlGnBu'))) (100))
dev.off()

get_sunburst_data(canopus, ofil)

pdf(file = paste0(ofil, '_classification_barplot.pdf'), width = 11)
ggplot(classifications_toplot, aes(level, count, fill = type)) + 
      geom_bar(stat = 'identity', position = position_dodge()) + theme_bw()
dev.off()

print('Done!')