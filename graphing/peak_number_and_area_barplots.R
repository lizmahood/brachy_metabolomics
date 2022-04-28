## Script for finding total # peaks and total peak area per sample.
## Also plots these

library(reshape2)
library(ggplot2)
library(ggpubr)

argg <- commandArgs(T)

get_groups <- function(vec){
  groups <- unlist(lapply(strsplit(vec, '_'), '[[', 1))
  return(groups)
}

make_df_toplot <- function(df){
  ## switching 0.01 values with 0
  mat <- as.matrix(df)
  mat[mat <= 0.01] <- 0
  df <- as.data.frame(mat)
  
  ## getting numbers of non-zero peaks per sample
  nonzerol <- lapply(df, function(x) nrow(df) - length(which(x == 0)))
  nonzerodf <- do.call(rbind.data.frame, nonzerol)
  nonzerodf$names <- colnames(df)
  nonzerodf$groups <- get_groups(colnames(df))
  mnonzerodf <- melt(nonzerodf, id.vars = c('names', 'groups'))
  mnonzerodf$variable <- NULL
  colnames(mnonzerodf)[3] <- 'num_peaks'
  
  ## total peak area per sample
  sumdf <- as.data.frame(colSums(df))
  mnonzerodf$pkarea_sum <- sumdf[, 1]
  
  ## changing names
  mnonzerodf$names <- unlist(lapply(strsplit(mnonzerodf$names, 'Run'), '[[', 1))
  mnonzerodf$names <- substr(mnonzerodf$names, 1, nchar(mnonzerodf$names)-1)
  return(mnonzerodf)
}

if (length(argg) != 3){
  stop('ARGS: 1) Peak area file (after all normalization and 
       with negatives turned to 0) for negative 2) Same for positive
       3) Output directory and name')
}

neg <- read.table(argg[1], sep = '\t', header = T, quote = "", fill = NA)
pos <- read.table(argg[2], sep = '\t', header = T, quote = "", fill = NA)

neg <- neg[, -c(1:32)]
pos <- pos[, -c(1:32)]

neg_pknum <- make_df_toplot(neg)
pos_pknum <- make_df_toplot(pos)

colpal <- c('black', 'brown2', 'chartreuse3', 'dodgerblue2', 'cyan', 'magenta',
            'goldenrod', 'grey50', 'darkgreen', 'darkblue', 'yellow1', 'violetred4',
            'yellowgreen', 'springgreen', 'pink', 'wheat4', 'lightcoral')

p1 <- ggbarplot(neg_pknum, x = "names", y = "pkarea_sum",
            fill = "groups",          
            color = "white",           
            palette = colpal,            
            sort.val = "none",          
            sort.by.groups = T,     
            ylab = "Total Peak Area",
            xlab = "",
            legend.title = "",
            rotate = TRUE,
            ggtheme = theme_bw()
)

p2 <- ggbarplot(neg_pknum, x = "names", y = "num_peaks",
                fill = "groups",          
                color = "white",          
                palette = colpal,         
                sort.val = "none",        
                sort.by.groups = T,    
                ylab = "Total Peak Number",
                xlab = "",
                legend.title = "",
                rotate = TRUE,
                ggtheme = theme_bw()
)

pdf(paste0(argg[3], 'NEG_peak_number_and_area_per_sample.pdf'), height = 10, width = 8)
ggarrange(p1, p2, labels = c('Total Peak Area', 'Total Peaks'), nrow = 1, ncol = 2, 
          common.legend = TRUE, legend = "right", align = 'hv')
dev.off()

p3 <- ggbarplot(pos_pknum, x = "names", y = "pkarea_sum",
                fill = "groups",          
                color = "white",           
                palette = colpal,            
                sort.val = "none",          
                sort.by.groups = T,     
                ylab = "Total Peak Area",
                xlab = "",
                legend.title = "",
                rotate = TRUE,
                ggtheme = theme_bw()
)

p4 <- ggbarplot(pos_pknum, x = "names", y = "num_peaks",
                fill = "groups",          
                color = "white",          
                palette = colpal,         
                sort.val = "none",        
                sort.by.groups = T,    
                ylab = "Total Peak Number",
                xlab = "",
                legend.title = "",
                rotate = TRUE,
                ggtheme = theme_bw()
)

pdf(paste0(argg[3], 'POS_peak_number_and_area_per_sample.pdf'), height = 10, width = 8)
ggarrange(p3, p4, labels = c('Total Peak Area', 'Total Peaks'), nrow = 1, ncol = 2, 
          common.legend = T, legend = 'right', align = 'hv')
dev.off()
