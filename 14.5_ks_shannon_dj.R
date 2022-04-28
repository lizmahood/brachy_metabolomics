## Ks tests for seeing if stressed leaves/roots have statistically significantly
## different values of specialization or diversity from controls

argg <- commandArgs(T)

if (length(argg) != 2){
  stop('ARGs: 1) diversity/specialization file 2) desired output file')
}

get_groups <- function(df){
  groups <- unlist(lapply(strsplit(df, '_'), '[', 1))
  return(groups)
}

get_pairwise_ks <- function(df){
  ## df must have groups col, Shannon/diversity col, and dj/specialization col
  shannon_mat <- matrix(nrow = length(unique(df$groups)),
                        ncol = length(unique(df$groups)),
                        dimnames = list(unique(df$groups), unique(df$groups)))
  dj_mat <- matrix(nrow = length(unique(df$groups)),
                   ncol = length(unique(df$groups)),
                   dimnames = list(unique(df$groups), unique(df$groups)))
  stop2 <- length(unique(df$groups))
  stop1 <- stop2 - 1
  for (i in 1:stop1){
    for (j in i:stop2){
      test <- df[which(df$groups == unique(df$groups)[i]),]
      ctrl <- df[which(df$groups == unique(df$groups)[j]),]
      this_ks_shannon <- ks.test(test$Shannon, ctrl$Shannon)
      shannon_mat[i, j] <- this_ks_shannon$p.value
      
      this_ks_dj <- ks.test(test$di, ctrl$di)
      dj_mat[i, j] <- this_ks_dj$p.value
    }
  }
  return(list(shannon_mat, dj_mat))
}


info <- read.table(argg[1], sep = '\t', header = T, quote = "", fill = NA)
odir <- argg[2]
info$groups <- get_groups(info$SampleName)

ks_lists <- get_pairwise_ks(info)

write.table(ks_lists[[1]], file = paste0(odir, '_pairwise_ks_for_shannon.tsv'),
            sep = '\t', quote = F)

write.table(ks_lists[[2]], file = paste0(odir, '_pairwise_ks_for_dj.tsv'),
            sep = '\t', quote = F)
