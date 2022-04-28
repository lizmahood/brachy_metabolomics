library(ggplot2)
library(ggrepel)
library(viridis)
library(tidyverse)

argg <- commandArgs(T)

if(length(argg) != 5){
  stop('ARGS: 1) canopus file 2) peak area file 3) do you want cosine?
        "no" or path to cosine score matrix 4) output directory for plot
       5) what level of canopus hierarchy to use? "class", "subclass", or "level.5"?')
}

## Functions
get_class_corr <- function(canopus, pkarea, classcol){
  ## Need to get unique canopus
  uniq_canopus <- unique(canopus[,classcol])
  ## Now for the correlation of the genes in those canopus
  corrlist <- list()
  for (clas in uniq_canopus){
    metabs <- unique(canopus[grepl(clas, canopus[,classcol], fixed = T), 'name'])
    if (length(metabs) > 1){
      metabs_for_cor <- as.matrix(pkarea[which(pkarea$Alignment.ID %in% metabs), 
                                            33:ncol(pkarea)])
      tmp_corr_vec <- cor(t(metabs_for_cor), method = 'spearman')
      ut <- tmp_corr_vec[upper.tri(tmp_corr_vec)]
      clas_avg_corr <- mean(ut)
      clas_sd_corr <- sd(ut)
      all_tocombine <- c(clas, clas_avg_corr, clas_sd_corr, length(metabs))
      corrlist[[clas]] <- all_tocombine
    }
  }
  out <- do.call(rbind.data.frame, corrlist)
  return(out)
}


get_class_cosine <- function(canopus, cosmat, classcol){
  ## Need to get unique canopus
  uniq_canopus <- unique(canopus[,classcol])
  
  ## Need to get the pairwise cosine scores of the metabs in each class
  coslist <- list()
  for (clas in uniq_canopus){
    metabs <- unique(canopus[grepl(clas, canopus[,classcol], fixed = T), 'name'])
    metabs <- paste0('X', metabs)
    if (length(metabs) > 1){
      select_cosine <- cosmat[metabs, metabs]
      ut <- select_cosine[upper.tri(select_cosine)]
      clas_avg_cos <- mean(ut)
      clas_sd_cos <- sd(ut)
      all_tocombine <- c(clas, clas_avg_cos, clas_sd_cos, length(metabs))
      coslist[[clas]] <- all_tocombine
    }
  }
  out <- do.call(rbind.data.frame, coslist)
  return(out)
}

get_average_correlation <- function(pkarea){
  corr_vec <- c()
  for (i in 1:10000L){
    toget <- sample(nrow(pkarea), 5, replace = F)
    genes_for_cor <- as.matrix(pkarea[toget, 33:(ncol(pkarea))])
    tmp_corr_vec <- cor(t(genes_for_cor))
    ut <- tmp_corr_vec[upper.tri(tmp_corr_vec)]
    corr_vec <- c(corr_vec, mean(ut))
  }
  return(mean(corr_vec))
}

get_average_cosine <- function(cosmat){
  cos_vec <- c()
  for (i in 1:10000L){
    toget <- sample(nrow(cosmat), 5, replace = F)
    cos_tmp <- cosmat[toget, toget]
    ut <- cos_tmp[upper.tri(cos_tmp)]
    cos_vec <- c(cos_vec, mean(ut))
  }
  return(mean(cos_vec))
}


## Using functions
canopus <- read.table(argg[1], sep = '\t', fill = NA, quote = "", header = T)
pkarea <- read.table(argg[2], sep = '\t', header = T, quote = "", fill = NA)
classcol <- which(colnames(canopus) == argg[5])

if(argg[3] != 'no'){
  cosmat <- read.table(argg[3], sep = '\t', header = T, quote = "", fill = NA,
                       row.names = 1)
  rownames(cosmat) <- paste0('X', rownames(cosmat))
  avg_cos <- get_average_cosine(cosmat)
  class_cos <- get_class_cosine(canopus, cosmat, classcol)
  colnames(class_cos) <- c('MetabClass', 'AvgCos', 'SDCos', 'NMetabs')
  class_cos$AvgCos <- as.numeric(class_cos$AvgCos)
  class_cos$SDCos <- as.numeric(class_cos$SDCos)
}


avg_corr <- get_average_correlation(pkarea)
class_corrs <- get_class_corr(canopus, pkarea, classcol)
colnames(class_corrs) <- c('MetabClass', 'AvgCorr', 'SDCorr', 'NMetabs')
class_corrs$AvgCorr <- as.numeric(class_corrs$AvgCorr)
class_corrs$SDCorr <- as.numeric(class_corrs$SDCorr)
class_corrs$NMetabs <- as.numeric(class_corrs$NMetabs)
print(head(class_corrs))

if (argg[3] != 'no'){
  both <- merge(class_corrs, class_cos, by = 'MetabClass')
  badclasses <- c('', 'None')
  both <- both[-which(both$MetabClass %in% badclasses),]
  both$SDCos <- NULL
  print(head(both))
  max_cos <- max(both$AvgCos)
  print(max_cos)

  ##writing out
  write.table(both, file = paste0(argg[4], '_class_correlation_cosine_metabs.tsv'),
              sep = '\t', row.names= F, quote = F)
  
  pdf(paste0(argg[4], '_class_correlation_cosine_bubble_plot.pdf'))
  ggplot(class_corrs, aes(x = AvgCorr, y = log2(NMetabs), size = SDCorr, color = AvgCorr)) + 
    geom_jitter(alpha = 0.7) + 
    geom_vline(xintercept = avg_corr, color = 'turquoise3', linetype = 'dashed', size = 2) + 
    scale_color_viridis(discrete = F, guide = 'none', option = "plasma") + theme_minimal() + 
    scale_size(range = c(.5, 10), name = 'Class\nSD') + 
    geom_text_repel(data= class_corrs %>% filter(AvgCorr >= 0.5), aes(label = MetabClass),
                    size = 3, color = 'saddlebrown')
  
  ggplot(both, aes(x = AvgCorr, y = AvgCos, size = NMetabs.x)) + 
    geom_jitter(alpha = 0.7, height = 0.01, width = 0.01, color = 'coral') + 
    geom_vline(xintercept = avg_corr, color = 'turquoise3', linetype = 'dashed', size = 1) +
    geom_hline(yintercept = avg_cos, color = 'lavender', linetype = 'dashed', size = 1) + 
    theme_minimal() + 
    scale_size(range = c(.5, 10), name = 'Class\nSize') + ylim(0, (max_cos + 0.04)) + 
    geom_text_repel(data = both %>% filter(AvgCorr >= 0.5 & AvgCos >= 0.125), 
                    aes(label = MetabClass), size = 3, color = 'saddlebrown')
  dev.off()
} else {
  
  pdf(paste0(argg[4], '_', argg[5], '_correlation_bubble_plot.pdf'), height = 2.5, width = 2.5)
  print(ggplot(class_corrs, aes(x = AvgCorr, y = log2(NMetabs), size = 0.5, color = AvgCorr)) + 
    geom_jitter(alpha = 0.7) + 
    geom_vline(xintercept = avg_corr, color = 'turquoise3', linetype = 'dashed', size = 1) + 
    scale_color_viridis(discrete = F, guide = 'none', option = "plasma") + theme_minimal() + 
    scale_size_continuous(guide = 'none') + xlab('') + ylab(''))
  dev.off()
}