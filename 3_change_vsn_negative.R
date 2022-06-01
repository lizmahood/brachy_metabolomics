##seeing if all values that were transformed to negative values were 0 originally
##if yes, then changing all negative values to 0.01

library(dplyr)

get_zeroes <- function(fullinp, fullorg){
  ##records alignment ID and sample name of negative values
  ##RETURNS: 2 col df with alignment ID and sample name
  fullinp <- fullinp[,-c(2:32)]
  
  ##getting columns in fullorg that correspond with fullinp
  org <- fullorg[,-c(2:32)]
  colstokeep <- which(colnames(org) %in% colnames(fullinp))
  org <- org[,colstokeep]
  
  ##making sure all columns in same order
  fullinp <- fullinp[,which(colnames(fullinp) %in% colnames(org))]
  fullinp <- fullinp[,order(colnames(fullinp))]
  org <- org[,order(colnames(org))]
  print(colnames(fullinp))
  print(colnames(org))
  
  out <- c()
  for (rw in 1:nrow(fullinp)){
    if (any(fullinp[rw,] < 0)){
      low <- which(fullinp[rw,] < 0)
      id <- fullinp[rw, 1]
      
      values <- as.numeric(org[which(org$Alignment.ID == id), low])
      out <- c(out, values)
    }
  }
 return(out) 
}

change_values <- function(fullinp){
  ##this only takes negative values and changes them to 0
  tochange <- as.matrix(fullinp[,-c(1:32)])
  tochange[tochange < 0] <- 0.01
  out <- data.frame(fullinp[,c(1:32)], tochange)
  return(out)
}

argg <- commandArgs(T)

if (length(argg) != 2){
  stop('ARGS: 1) VSN Normalized peak area file, 2) Filtered peak area file, but without VSN')
}

fullinp <- read.table(argg[1], sep = '\t', header = T, stringsAsFactors = F, quote = "", fill = NA)
fullorg <- read.table(argg[2], sep = '\t', header = T, stringsAsFactors = F, quote = "", fill = NA)

zeros <- get_zeroes(fullinp = fullinp, fullorg = fullorg)

if (length(unique(zeros)) == 1 & unique(zeros) == 0){
  out <- change_values(fullinp = fullinp)
  write.table(out, file = paste0(argg[1], '_no_negative.tab'), sep = '\t', row.names = F, quote = F)
} else {
  print('Some non-zero original values!')
  print(table(zeros))
}

print('Done!')
