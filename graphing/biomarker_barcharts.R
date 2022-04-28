## Read in cond specific metab file and canopus file
## Make stacked barchart of enriched metabs per condition (color by class)
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(dplyr)

argg <- commandArgs(T)
if (length(argg) != 3){
  stop('ARGS: 1) condition-specific metabolite file 2) canopus file
       3) desired output path and prefix')
}

condspec <- read.table(argg[1], sep = '\t', header = T, quote = "", fill = NA)
canopus <- read.table(argg[2], sep = '\t', header = T, quote = "", fill = NA)
canopus_merge <- canopus[, c('name', 'class')]

both <- merge(condspec, canopus_merge, all.x = T, by.x = 'Alignment.ID', by.y = 'name')
both[is.na(both$class), 'class'] <- 'None'

colors <- c(brewer.pal(12, 'Paired'), 'grey28', 'bisque', 'maroon2', 'turquoise4',
            'palegreen', 'navy', 'beige', 'chartreuse', 'orangered',
            'darkgreen', 'darkred', 'lightblue1', 'deeppink4', 'thistle1',
            'lightgrey', 'slateblue1', 'skyblue1', 'wheat4')

both$class <- gsub(' and derivatives', '', both$class)
both$class <- gsub(' and substituted derivatives', '', both$class)
both$class <- gsub(' and steroid derivatives', '', both$class)

pdf(paste0(argg[3], 'class_barplots.pdf'), width = 8, height = 8)
ggplot(both, aes(Condition, fill = class)) + geom_bar(color = 'black') +
  scale_fill_manual(values= colors) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
