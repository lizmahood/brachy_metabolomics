## Script to make barplots of specialization and shannon entropy per condition

library(ggplot2)
library(data.table)
library(stringr)

argg <- commandArgs(T)

## Read in data
## Get groups
## aggregate
## plot (put leaves together and roots together)

argg <- commandArgs(T)

if (length(argg) != 2){
  stop('ARG: 1) shannon and di per sample file 2) output directory and filename
       prefix')
}

shannondi <- read.table(argg[1], sep = '\t', header = T, fill = NA, quote = "")
shannondi$groups <- str_split_fixed(shannondi$SampleName, '_', n = 2)[,1]
shannondi$SampleName <- NULL
shannondi$TotalPeaks <- NULL

ag <- aggregate(. ~ groups, shannondi, function(x) c(mean = mean(x), sd = sd(x)))
ag <- do.call(data.frame, ag)

ag$groups <- factor(ag$groups, levels = c('Hydro.Cop.Leaf', 'Hydro.Ctrl.Leaf', 
                                          'Hydro.Heat.Leaf','Hydro.HeatCop.Leaf', 
                                          'Sym.Ctrl.Leaf', 'Sym.Spore.Leaf',
                                          'Sym.SporeW.Leaf', 'Tis.Leaf', 
                                          'Hydro.Cop.Root', 'Hydro.Ctrl.Root',
                                          'Hydro.Heat.Root', 'Hydro.HeatCop.Root',
                                          'Sym.Ctrl.Root', 'Sym.Spore.Root', 
                                          'Sym.SporeW.Root',
                                          'Tis.Culm', 'Tis.Spike'))

## Getting organs of each sample
ag$Organ <- gsub('.', '_', ag$groups, fixed = T)
broken <- as.data.frame(str_split_fixed(ag$Organ, '_', 3))
ag$Organ <- broken[ ,ncol(broken)]
ag$Organ[c(15:17)] <- broken[15:17, 2]

## Getting min and max of shannon and di
smin <- min(ag$Shannon.mean - ag$Shannon.sd) - 0.01
smax <- max(ag$Shannon.mean + ag$Shannon.sd) + 0.01

dmin <- min(ag$di.mean - ag$di.sd) - 0.01
dmax <- max(ag$di.mean + ag$di.sd) + 0.01

## barplot for shannon
pdf(paste0(argg[2], '_diversity_specialization_barplots.pdf'))
ggplot(ag, aes(as.factor(groups), Shannon.mean, fill = Organ)) + 
  geom_bar(position = position_dodge(), stat = "identity", color = 'black') +
  geom_errorbar(aes(ymin = Shannon.mean - Shannon.sd, 
                    ymax = Shannon.mean + Shannon.sd),
                width=.2,position=position_dodge(.9)) + theme_bw() + 
  coord_cartesian(ylim=c(smin, smax)) + 
  theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1)) + 
  xlab('') + ylab('Diversity')

## and for dj
ggplot(ag, aes(as.factor(groups), di.mean, fill = Organ)) + 
  geom_bar(position = position_dodge(), stat = "identity", color = 'black') +
  geom_errorbar(aes(ymin = di.mean - di.sd, 
                    ymax = di.mean + di.sd),
                width = .2,position = position_dodge(.9)) + theme_bw() + 
  coord_cartesian(ylim=c(dmin, dmax)) + 
  theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1)) + 
  xlab('') + ylab('Specialization')

dev.off()


