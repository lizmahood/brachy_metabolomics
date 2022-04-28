library(matrixStats)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(grid)
library(reshape2)

get_fc <- function(ctrlv, stressv){
  nctm <- rowMeans(ctrlv)
  nctm[nctm <= 0.01] <- 1e-30
  nstm <- rowMeans(stressv)
  ## new 9.30.21
  nstm[nstm <= 0.01] <- 1e-30
  logfc <- log2(nstm/nctm)
  #ifelse(nstm <= 0.01, logfc <- 0, logfc <- log2(nstm/nctm))
  return(logfc)
}

get_groups <- function(df){
  ct <- colnames(df)
  groups <- unlist(lapply(strsplit(ct, '_'), '[', 1))
  return(groups)
}

get_all_fc <- function(combined, df_of_int){
  all_fc <- list()
  for (tissue in combined){
    if (length(tissue) > 0){
      tname <- unique(tissue[grepl('Ctrl.', tissue)])
      ctrlcol <- which(grepl(tname, colnames(df_of_int)))
      print('Control: ')
      print(colnames(df_of_int[ctrlcol][1]))
      treat <- tissue[-which(grepl('Ctrl', tissue))]
      uniq <- unique(treat)
      for (cond in uniq){
        print(cond)
        ttreat <- which(grepl(paste0('^', cond), colnames(df_of_int)))
        if (length(ttreat) > 1 & length(ctrlcol) > 1){
          print(str(df_of_int[,ttreat]))
          fcs <- get_fc(df_of_int[,ctrlcol], df_of_int[,ttreat])
          print(head(fcs))
          all_fc[[cond]] <- fcs
        }
      }
    }
  }
  fcdf <- do.call(cbind.data.frame, all_fc)
  mfcdf <- melt(fcdf)
  return(mfcdf)
}

argg <- commandArgs(T)

if(length(argg) != 3){
  stop('ARGS: 1) File of NORMALIZED peak area 2) Corresponding file of NON-NORM
       peak area 3) All experiments together? yes OR no')
}


npeakarea <- read.table(argg[1], sep = '\t', fill = NA, quote = "", header = T)
pkareaorg <- read.table(argg[2], sep = '\t', fill = NA, quote = "", header = T)
allg <- argg[3]

## Do once for getting all new fold changes
newgroups <- get_groups(npeakarea[,-c(1:32)])
nleaves <- newgroups[which(grepl('Leaf', newgroups))]
nroots <- newgroups[which(grepl('Root', newgroups))]

if (allg == 'yes'){
  nhleaves <- nleaves[which(grepl('Hydro', nleaves))]
  nhroots <- nroots[which(grepl('Hydro', nroots))]
  nsleaves <- nleaves[which(grepl('Sym', nleaves))]
  nsroots <- nroots[which(grepl('Sym', nroots))]
  ncombined <- list(nhleaves, nhroots, nsleaves, nsroots)
} else {
  ncombined <- list(nleaves, nroots)
}

new_fcs <- get_all_fc(ncombined, npeakarea[,-c(1:32)])


## Do again for getting new fold changes
orggroups <- get_groups(pkareaorg[,-c(1:32, ncol(pkareaorg))])
oleaves <- orggroups[which(grepl('Leaf', orggroups))]
oroots <- orggroups[which(grepl('Root', orggroups))]

if (allg == 'yes'){
  ohleaves <- oleaves[which(grepl('Hydro', oleaves))]
  ohroots <- oroots[which(grepl('Hydro', oroots))]
  osleaves <- oleaves[which(grepl('Sym', oleaves))]
  osroots <- oroots[which(grepl('Sym', oroots))]
  ocombined <- list(ohleaves, ohroots, osleaves, osroots)
} else {
  ocombined <- list(oleaves, oroots)
}

org_fcs <- get_all_fc(ocombined, pkareaorg[,-c(1:32, 115)])

## getting indices of metabs with original fc <2 and >-2 (non-significant)
nonsig_org <- which(abs(org_fcs$value) < 2)
nonsig_new <- new_fcs[nonsig_org, 'value']

print(top_cutoff <- quantile(nonsig_new, 0.99))
print(bottom_cutoff <- quantile(nonsig_new, 0.01))
