library(ggplot2)
library(ggrepel)
library(ggpubr)

argg <- commandArgs(T)
if (length(argg) != 6){
  stop('ARGS: 1) Path to shannon entropy tab file 2) RDPI file 
       3) RDPI per class file 4) RDPI per class and metabolites 5) output directory
       6) RDPI per class and stress (violin RDPI plots per class')
}

get_groups <- function(df, typ){
  ct <- colnames(df)
  groups <- unlist(lapply(strsplit(ct, '_'), '[', 1))
  return(groups)
}

make_meta <- function(samples){
  ##metabs is a data frame
  conds <- c(); tiss <- c()
  for (samp in samples){
    if (grepl('Leaf', samp)){
      tiss <- c(tiss, 'Leaf')
    }else if (grepl('Root', samp)){
      tiss <- c(tiss, 'Root')
    }else if (grepl('Culm', samp)){
      tiss <- c(tiss, 'Culm')
    }else if (grepl('Spike', samp)){
      tiss <- c(tiss, 'Spike')
    }
    
    name <- strsplit(samp, '_', fixed = T)[[1]][1]
    conds <- c(conds, paste(strsplit(name, '.', fixed = T)[[1]][c(1,2)], collapse = '.'))
  }
  
  tiss <- as.factor(tiss); conds <- as.factor(conds)
  tmetabs <- cbind.data.frame(tiss, conds)
  return(tmetabs)
}

get_highest_classes_per_stress <- function(class_rdpi){
  rdpil <- list()
  pkareal <- list()
  class_rdpi$X <- NULL
  for (stress in unique(class_rdpi$Stress)){
    print(stress)
    tmp <- class_rdpi[class_rdpi$Stress == stress,]
    ## For RDPI list
    rest <- tmp[-which(tmp$MetabClass == 'None'),]
    toappend <- head(rest[order(rest$RDPI, decreasing = T),], 5)
    rdpil[[stress]] <- toappend
    
    ##For pk area list
    highs <- tmp[tmp$ChangePercentile >= 70,]
    if (length(highs) < 1){
      highs <- tmp[order(tmp$ChangePercentile),]
      highs <- highs[1:5,]
    }else if (nrow(highs) > 5){
      highs <- highs[order(highs$ChangePercentile, decreasing = T),]
      highs <- highs[1:5,]
    }
    print(highs)
    lows <- tmp[tmp$ChangePercentile <= 30,]
    if (length(lows) < 1){
      lows <- tmp[order(tmp$ChangePercentile),]
      lows <- head(lows, 5)
    }else if (nrow(lows) > 5){
      lows <- lows[order(lows$ChangePercentile),]
      lows <- head(lows, 5)
    }
    print(lows)
    toappendpk <- rbind(highs, lows)
    pkareal[[stress]] <- toappendpk
  }
  pkareadf <- do.call(rbind.data.frame, pkareal)
  rdpidf <- do.call(rbind.data.frame, rdpil)
  return(list(pkareadf, rdpidf))
}

filter_metabolite_rdpi_file <- function(class_rdpi_avg, class_rdpi_all){
  ## For each stress
  ## Remove classes not in the class_rdpi_avg file 
  ## (which has top 5 most changed classes per stress)
  ## I guess also paste the class peak area average onto this df
  ## And remove classes (even if they're in top 5 most changed) that don't have
  ## >= 5 entries per stress
  for (strs in unique(class_rdpi_avg$Stress)){
    tmpavg <- class_rdpi_avg[class_rdpi_avg$Stress == strs,]
    tmpall <- class_rdpi_all[class_rdpi_all$Stress == strs,]
    allclasses <- as.data.frame(table(tmpall$Class))
    badclasses <- allclasses$Var1[allclasses$Freq < 5]
    goodclasses <- unique(tmpavg$MetabClass)
    if (any(badclasses %in% goodclasses)) {
      goodclasses <- goodclasses[-which(goodclasses %in% badclasses)]
    }
    class_rdpi_all <- class_rdpi_all[-which(class_rdpi_all$Stress == strs & !(class_rdpi_all$Class %in% goodclasses)),]
  }
  tomerge <- class_rdpi_avg[,c(1,2,5)]
  out <- merge(class_rdpi_all, tomerge, all = F, all.x = T, by.x = c('Stress', 'Class'), by.y = c('Stress', 'MetabClass'))
  return(out)
}

shannon <- read.table(argg[1], sep = '\t', stringsAsFactors = F, quote = "", header = T)
rdpi <-  read.table(argg[2], sep = '\t', stringsAsFactors = F, quote = "", header = T)
rdpi_by_class <- read.table(argg[6], sep = '\t', stringsAsFactors = F, quote = "", header = T)
class_rdpi <- read.table(argg[3], sep = '\t', stringsAsFactors = F, quote = "", header = T, fill = NA)
class_rdpi_all <- read.table(argg[4], sep = '\t', stringsAsFactors = F, quote = "", header = T, fill = NA)

## getting only the type of stress
rdpi$test_data <- unlist(lapply(strsplit(rdpi$Names, '_'), function(x) paste(x[c(3,4)], collapse = '_')))
rdpi$test_data <- substr(rdpi$test_data, 0, nchar(rdpi$test_data)-2)

rdpi_by_class$test_data <- unlist(lapply(strsplit(rdpi_by_class$Names, '_'), function(x) paste(x[c(3,4)], collapse = '_')))
rdpi_by_class$test_data <- substr(rdpi_by_class$test_data, 0, nchar(rdpi_by_class$test_data)-2)

smeta <- make_meta(shannon$SampleName)
stoplot <- cbind.data.frame(shannon, smeta)

classl <- get_highest_classes_per_stress(class_rdpi)
pkareadf <- classl[[1]]
rdpidf <- classl[[2]]
rdpidf$RDPI <- rdpidf$RDPI * 100
parsed_rdpi_metab <- filter_metabolite_rdpi_file(pkareadf, class_rdpi_all)

## Beautifying metabolite class names, for berevity in figure
parsed_rdpi_metab$Class <- gsub(' and derivatives', '', parsed_rdpi_metab$Class)

## changing factor levels of rdpi so all leaves are together and all roots 
rdpi$test_data <- factor(rdpi$test_data, levels = c('Hydro.Cop.Leaf', 'Hydro.Heat.Leaf',
                                                       'Hydro.HeatCop.Leaf', 'Sym.Spore.Leaf',
                                                       'Sym.SporeW.Leaf', 'Hydro.Cop.Root',
                                                       'Hydro.Heat.Root', 'Hydro.HeatCop.Root',
                                                       'Sym.Spore.Root', 'Sym.SporeW.Root'))

rdpi_by_class$test_data <- factor(rdpi_by_class$test_data, levels = c('Hydro.Cop.Leaf', 'Hydro.Heat.Leaf',
                                                    'Hydro.HeatCop.Leaf', 'Sym.Spore.Leaf',
                                                    'Sym.SporeW.Leaf', 'Hydro.Cop.Root',
                                                    'Hydro.Heat.Root', 'Hydro.HeatCop.Root',
                                                    'Sym.Spore.Root', 'Sym.SporeW.Root'))

## making colors be green for leaves and brown for roots
rdpi$tiss <- ifelse(grepl('Leaf', rdpi$test_data), 'Leaf', 'Root')
colorss <- c('forestgreen', 'saddlebrown')

## Now for colors for shannon/dj
cvec <- c('salmon', 'mediumslateblue','palegreen', 'goldenrod2',
          'lightblue', 'forestgreen', 'plum3', 'darkorchid4', 'pink',
          'azure4', 'black', 'darkred')

pdf(paste0(argg[5], '_info_theory_plots.pdf'), height = 10, width = 10)
#ggplot(shannon, aes(TotalPeaks, Shannon)) + geom_point(size = 2, color = 'navy') +
#  theme_bw() + xlab('Total Number of Peaks') + ylab('Shannon Entropy') + 
#  geom_text_repel(aes(label = Condition),size = 4.5, color = 'saddlebrown')

ggplot(stoplot, aes(TotalPeaks, Shannon, color = conds, shape = tiss, size = 1)) + geom_point(alpha = 0.8) +
  theme_bw() + guides(size = 'none', label = 'none',shape = guide_legend(override.aes = list(size = 3)),
                      color = guide_legend(override.aes = list(size = 3))) + scale_color_manual(values = cvec)


ggplot(stoplot, aes(di, Shannon, color = conds, shape = tiss, size = 1)) + geom_point(alpha = 0.8) + 
  theme_bw() + guides(size = 'none', label = 'none', shape = guide_legend(override.aes = list(size = 3)),
                                                     color = guide_legend(override.aes = list(size = 3))) + 
  scale_color_manual(values = cvec)
dev.off()

pdf(paste0(argg[5], '_rdpi_all_main.pdf'), height = 3, width = 5)
ggplot(rdpi, aes(test_data, RDPI, fill = tiss)) + geom_violin(alpha = 0.4) + geom_jitter(size = 0.2, width = 0.25) + theme_bw() +
  labs(x = '', y = 'RDPI') + theme(axis.text.x = element_text(angle = 90, size = 10, hjust = 1)) +
  scale_fill_manual(values = colorss)

dev.off() 

## beautifying names of classes
rdpi_by_class$Class <- gsub(' and derivatives', '', rdpi_by_class$Class)
rdpi_by_class$Class <- gsub(' and substituted derivatives', '', rdpi_by_class$Class)
rdpi_by_class$Class <- gsub(' and steroid derivatives', '', rdpi_by_class$Class)

parsed_rdpi_metab$Class <- gsub(' and derivatives', '', parsed_rdpi_metab$Class)
parsed_rdpi_metab$Class <- gsub(' and substituted derivatives', '', parsed_rdpi_metab$Class)
parsed_rdpi_metab$Class <- gsub(' and steroid derivatives', '', parsed_rdpi_metab$Class)

pdf(paste0(argg[5], '_RDPI_per_class_for_paper.pdf'), width = 10, height = 11)
#for (cls in unique(rdpi_by_class$Class)){
#  tmp <- rdpi_by_class[rdpi_by_class$Class == cls,]
#  print(ggplot(tmp, aes(test_data, RDPI)) + geom_violin(fill = 'green', alpha = 0.4) + geom_jitter(size = 2) + theme_bw() + 
#    labs(x = '', y = 'RDPI') + ggtitle(cls) + theme(axis.text.x = element_text(angle = 45, size = 12, hjust = 1)))
#}
ggplot(rdpi_by_class, aes(test_data, RDPI)) + geom_violin(fill = 'green', alpha = 0.4) + 
  geom_jitter(size = 0.4) + facet_wrap(~Class, scales = 'free_y') + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, size = 7, hjust = 1)) + xlab('')
dev.off()

barplotl <- list()
violinplotl <- list()
rdpi_classl <- list()

## Adding flag for outliers in the violin plots. Added 1/23/22
parsed_rdpi_metab$is_out <- ifelse((abs(parsed_rdpi_metab$AvgPkAreaChange - mean(parsed_rdpi_metab$AvgPkAreaChange)/ sd(parsed_rdpi_metab$AvgPkAreaChange))) < 3, 'no', 'yes')

for (strs in unique(pkareadf$Stress)){
  tmpparsed <- parsed_rdpi_metab[parsed_rdpi_metab$Stress == strs,]
  
  
  tmprdpi <- pkareadf[pkareadf$Stress == strs,]
  tmprdpidf <- rdpidf[rdpidf$Stress == strs,]
  
  violinplotl[[strs]] <- ggplot(tmpparsed, aes(reorder(Class, ChangePercentile), AvgPkAreaChange)) + 
    geom_jitter(size = 1.5, aes(colour = AvgPkAreaChange, shape = is_out), width=0.25, height = 0) + 
    geom_violin(aes(fill = ChangePercentile), alpha = 0.35) + theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, size = 9, hjust = 1), legend.position = "none") + 
    scale_fill_viridis_c() + labs(x = '', y = '') + facet_wrap(~ Stress, scales = 'free') + 
    scale_color_viridis_c() + scale_shape_manual(values=c(1, 3))

  barplotl[[strs]] <- ggplot(tmprdpi, aes(reorder(MetabClass, ChangePercentile), Avg_PeakHeight_change, fill = Avg_PeakHeight_change)) + 
    geom_bar(stat = 'identity') + theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, size = 9, hjust = 1), legend.position = "none") + 
    geom_text(aes(label=round(ChangePercentile, 0)), color = 'saddlebrown', 
              vjust = ifelse(tmprdpi$Avg_PeakHeight_change >= 0, 0.3, -0.6), size=5) + 
    scale_fill_viridis_c() + labs(x = '', y = '') + facet_wrap(~ Stress, scales = 'free')
  
  rdpi_classl[[strs]] <- ggplot(tmprdpidf, aes(reorder(MetabClass, RDPI), RDPI, fill = RDPI)) +
    geom_bar(stat = 'identity') + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, size = 9, hjust = 1), legend.position = 'none') + 
    scale_fill_viridis_c() + labs(x = '', y = '') + facet_wrap(~ Stress, scales = 'free')
}

pdf(paste0(argg[5], '_info_theory_per_stress.pdf'), height = 10, width = 12)

ggarrange(barplotl[[1]], barplotl[[2]], barplotl[[3]], barplotl[[4]], barplotl[[5]], barplotl[[6]], barplotl[[7]], 
          barplotl[[8]], barplotl[[9]], barplotl[[10]], ncol = 4, nrow = 3)

ggarrange(violinplotl[[1]], violinplotl[[2]], violinplotl[[3]], violinplotl[[4]], violinplotl[[5]], violinplotl[[6]], 
          violinplotl[[7]], violinplotl[[8]], violinplotl[[9]], violinplotl[[10]], ncol = 4, nrow = 3)

ggarrange(rdpi_classl[[1]], rdpi_classl[[2]], rdpi_classl[[3]], rdpi_classl[[4]], rdpi_classl[[5]], rdpi_classl[[6]], 
          rdpi_classl[[7]], rdpi_classl[[8]], rdpi_classl[[9]], rdpi_classl[[10]], ncol = 4, nrow = 3)

dev.off()

pdf(paste0(argg[5], '_violin_per_stress_for_figure.pdf'), height = 5, width = 10)
ggarrange(violinplotl[[1]], violinplotl[[2]], violinplotl[[3]], violinplotl[[4]], violinplotl[[5]], violinplotl[[6]], 
          violinplotl[[7]], violinplotl[[8]], violinplotl[[9]], violinplotl[[10]], ncol = 5, nrow = 2, align = 'hv')

ggplot(parsed_rdpi_metab, aes(reorder(Class, ChangePercentile), AvgPkAreaChange)) + 
  geom_jitter(size = 1.5, aes(colour = AvgPkAreaChange, shape = is_out), width=0.25, height = 0) + 
  geom_violin(aes(fill = ChangePercentile), alpha = 0.35) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, size = 9, hjust = 1), legend.position = "none") + 
  scale_fill_viridis_c() + labs(x = '', y = '') + facet_wrap(~ Stress, scales = 'free', ncol = 5) + 
  scale_color_viridis_c() + scale_shape_manual(values=c(1, 3))
dev.off()
