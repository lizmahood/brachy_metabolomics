## making fragmentation patterns for paper
library(ggplot2)

argg <- commandArgs(T)

if (length(argg) != 2){
  stop('ARGS: 1) mSDIAL-output frag pattern for one compound
       2) desired output directory and prefix for plot')
}

## pre-processing
lines <- readLines(argg[1])
bad <- c(1:11, length(lines)-1)
peaks <- read.table(text = lines[-c(bad)], sep = ' ', fill = NA)
peaks$V1 <- as.numeric(peaks$V1)


pdf(paste0(argg[2], '_fragmentation_pattern_graph.pdf'), height = 4, width = 6)
ggplot(peaks, aes(V1, V2)) + geom_linerange(aes(ymin = 0, ymax = V2)) + 
  theme_bw() + xlab('m/z') + ylab('Intensity') +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + geom_hline(yintercept = 0)
dev.off()
