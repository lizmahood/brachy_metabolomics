## For finding differences among highly correlated classes vs. lowly correlated
## MERGE WITH THE BUBBLE PLOT SCRIPT
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

argg <- commandArgs(T)
if (length(argg) != 5){
  stop('ARGS: 1) class correlation file output by correlation bubble plots
       2) normalized pkarea file 3) filtered canopus file 4) what level of the 
       hierarchy do you want to analyze (class, subclass, or level.5
       5) output directory and prefix')
}

get_groups <- function(df){
  ct <- colnames(df)
  groups <- unlist(lapply(strsplit(ct, '_'), '[', 1))
  return(groups)
}

get_expression_of_class <- function(pkarea, lvl){
  ## get pkarea of all metabs in classnm, get only quant cols, no metadata
  ## transpose that df
  ## get group of each row (sample) as col and also get tissue as col
  ## summarise df -- get conds/tissues where sum of metab expression is > 0
  ## Save these conds/tissues as string, also # (of how many conds/tissues)
  ##   class is expressed in
  lvlcol <- which(colnames(pkarea) == lvl)
  all_classes <- unique(pkarea[,lvlcol])
  outl <- list()
  for (classnm in unique(pkarea[,lvlcol])){
    if (!(is.na(classnm))){
      pkarea_use <- pkarea[which(pkarea[,lvlcol] == classnm), c(33:(ncol(pkarea) - 3))]
      pkarea_use[pkarea_use == 0.01] <- 0
      groups <- get_groups(pkarea_use)
      tissues <- unlist(lapply(strsplit(groups, '.', fixed = T), 
                               function(x) tail(x, n = 1)))
      tpkarea_use <- as.data.frame(t(pkarea_use))
      tpkarea_use$sum <- rowSums(tpkarea_use)
      tpkarea_use$groups <- groups
      tpkarea_use$tissues <- tissues
      summarise_conds <- tpkarea_use %>% group_by(groups) %>% summarise(mean(sum))
      summarise_tiss <- tpkarea_use %>% group_by(tissues) %>% summarise(mean(sum))
      num_tiss <- length(which(summarise_tiss$`mean(sum)` > 0))
      num_conds <- length(which(summarise_conds$`mean(sum)` > 0))
      all_tiss <- paste0(summarise_tiss$tissues[which(summarise_tiss$`mean(sum)` > 0)], 
                         collapse = '_')
      all_conds <- paste0(summarise_conds$groups[which(summarise_conds$`mean(sum)` > 0)], 
                          collapse = '_')
      outl[[classnm]] <- c(num_tiss, all_tiss,num_conds, all_conds)
    }
  }
  out <- do.call(rbind.data.frame, outl)
  colnames(out) <- c('num_tissues', 'all_tissues', 'num_conds', 'all_conds')
  out$classes <- all_classes[-which(is.na(all_classes))]
  out$num_tissues <- as.numeric(out$num_tissues)
  out$num_conds <- as.numeric(out$num_conds)
  return(out)
}

corr_class <- read.table(argg[1], sep = '\t', header = T, quote = "", fill = NA)
pkarea <- read.table(argg[2], sep = '\t', header = T, quote = "", fill = NA)
canopus <- read.table(argg[3], sep = '\t', header = T, quote = "", fill = NA)
hierarchy <- argg[4] # class, subclass, or level.5

## finding highly correlated classes -- hiher than 75th percentile of correlation
thresh <- unname(summary(abs(corr_class$AvgCorr))[5])
corr_class$corr_type <- ifelse(corr_class$AvgCorr >= thresh, 'high', 'low')

## finding tissues/conditions a class is expressed in 
canopus_tomerge <- canopus[, c('name', 'class', 'subclass', 'level.5')]
merged <- merge(pkarea, canopus_tomerge, by.x = 'Alignment.ID', by.y = 'name',
                all.x = T)

## Want to get: #conditions a class is expressed in, #tissues its expressed in
## A string combination of the conditions its in and same for tissues
class_exp <- get_expression_of_class(merged, hierarchy)
toplot <- merge(corr_class, class_exp, by.x = 'MetabClass', by.y = 'classes', all.x = T)
write.table(toplot, file = paste0(argg[5], hierarchy, '_stats_high_vs_low_corr.tsv'),
            sep = '\t', quote = F, row.names = F)

## plotting
mtoplot <- pivot_longer(toplot[,c(3,5:8,10)], cols = c(1:3,5,6))

p <- ggboxplot(mtoplot, x = "corr_type", y = "value",
               color = 'corr_type', pallete = 'jco',
               add = "jitter", short.panel.labs = FALSE, 
               title = paste0(hierarchy, ', cor threshold = ', thresh))

pq <- facet(p, 'name', scales = 'free_y')

pq <- pq + stat_compare_means(label = "p.format",
                             method = "wilcox.test", paired = FALSE,)

pdf(paste0(argg[5], hierarchy, '_high_vs_low_corr.pdf'))
print(pq)
dev.off()

r2a <- round(summary(lm(AvgCos ~ abs(AvgCorr), data= toplot))$r.squared, digits = 3)
corra <- round(cor(toplot$AvgCos, abs(toplot$AvgCorr)), digits = 3)
r2b <- round(summary(lm(NMetabs.x ~ abs(AvgCorr), data= toplot))$r.squared, digits = 3)
corrb <- round(cor(toplot$NMetabs.x, abs(toplot$AvgCorr)), digits = 3)
r2c <- round(summary(lm(NMetabs.x ~ AvgCos, data= toplot))$r.squared, digits = 3)
corrc <- round(cor(toplot$NMetabs.x, toplot$AvgCos), digits = 3)

cosx <- max(toplot$AvgCos)/2

mytheme = theme(axis.title.x = element_text(size = 16),
                axis.text.x = element_text(size = 14),
                axis.title.y = element_text(size = 16),
                axis.text.y = element_text(size = 14))

pdf(paste0(argg[5], hierarchy, '_corr_vs_size_cosine.pdf'))
print(ggplot(toplot, aes(abs(AvgCorr), AvgCos)) + geom_point(alpha = 0.4, size = 2) + 
      geom_smooth(method = 'lm', se = F, color = 'turquoise4') + theme_bw() + 
      annotate('text',x = 0.5, y = cosx, 
               label = paste0(expression(italic(r)^2), ' =', r2a), 
               parse = TRUE, size = 8) + 
        annotate('text', x = 0.5, y = cosx, label = paste0('\ncorr = ', corra), size = 8) + 
        xlab('Average Spearman Correlation') + ylab('Average Cosine Score') + mytheme)

print(ggplot(toplot, aes(abs(AvgCorr), NMetabs.x)) + geom_point(alpha = 0.4, size = 2) + 
  geom_smooth(method = 'lm', se = F, color = 'turquoise4') + theme_bw() + 
    annotate('text',x = 0.5, y = 90, label = paste0(expression(italic(r)^2), ' =', r2b), 
             parse = TRUE, size = 8) + 
    annotate('text', x = 0.5, y = 75, label = paste0('\ncorr = ', corrb), size = 8) + 
    xlab('Average Spearman Correlation') + ylab('Number of Peaks') + mytheme)

print(ggplot(toplot, aes(AvgCos, NMetabs.x)) + geom_point(alpha = 0.4, size = 2) + 
  geom_smooth(method = 'lm', se = F, color = 'turquoise4') + theme_bw() + 
    annotate('text',x = cosx, y = 90, label = paste0(expression(italic(r)^2), ' =', r2c), 
             parse = TRUE, size = 8) + 
    annotate('text', x = cosx, y = 75, 
             label = paste0('\ncorr = ', corrc), size = 8) + 
    xlab('Average Cosine Score') + ylab('Number of Peaks') + mytheme)
dev.off()

