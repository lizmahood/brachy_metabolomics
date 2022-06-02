library(matrixStats)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(grid)
library(reshape2)
library(stringr)

get_fc <- function(ctrlv, stressv){
  ctrlv[ctrlv <= 0.01] <- 1e-30
  nctm <- rowMeans(ctrlv)
  
  stressv[stressv <= 0.01] <- 1e-30
  nstm <- rowMeans(stressv)
  
  logfc <- log2(nstm/nctm)
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
      if (length(ctrlcol) > 0){
        treat <- tissue[-which(grepl('Ctrl', tissue))]
        if (length(treat) > 0){
          uniq <- unique(treat)
          for (cond in uniq){
            ttreat <- which(grepl(paste0('^', cond), colnames(df_of_int)))
            if (length(ttreat) > 1 & length(ctrlcol) > 1){
              print('Control: ')
              print(str(df_of_int[,ctrlcol]))
              print('Stress: ')
              print(str(df_of_int[,ttreat]))
              fcs <- get_fc(df_of_int[,ctrlcol], df_of_int[,ttreat])
              print(head(fcs))
              all_fc[[cond]] <- fcs
            }
          }
        }
      }
    }
  }
  fcdf <- do.call(cbind.data.frame, all_fc)
  mfcdf <- melt(fcdf)
  return(mfcdf)
}

argg <- commandArgs(T)

if(length(argg) != 2){
  stop('ARGS: 1) File of NORMALIZED peak area, with negatives turned to 0
       2) Corresponding file of NON-NORM peak area')
}


npeakarea <- read.table(argg[1], sep = '\t', fill = NA, quote = "", header = T)
pkareaorg <- read.table(argg[2], sep = '\t', fill = NA, quote = "", header = T)


## Do once for getting all new fold changes
newgroups <- get_groups(npeakarea[,-c(1:32)])
newgroups <- newgroups[which(str_count(newgroups, fixed(".")) == 2)]

tissues <- unlist(lapply(strsplit(newgroups, '.', fixed = T), '[[', 3))
tissues <- unique(str_to_title(tissues))
ncombined <- list()
for (t in tissues){
  tmp <- newgroups[which(grepl(t, newgroups))]
  conds <- unique(lapply(strsplit(tmp, '.', fixed = T), '[[', 1))
  for (cd in conds){
    tmp2 <- tmp[which(grepl(cd, tmp))]
    ncombined[[paste0(t, cd, sep = '.')]] <- tmp2
  }
}


new_fcs <- get_all_fc(ncombined, npeakarea[,-c(1:32)])


## Do again for getting new fold changes
ogroups <- get_groups(pkareaorg[,-c(1:32, ncol(pkareaorg))])
ogroups <- ogroups[which(str_count(ogroups, fixed(".")) == 2)]


otissues <- unlist(lapply(strsplit(ogroups, '.', fixed = T), '[[', 3))
otissues <- unique(str_to_title(otissues))
ocombined <- list()
for (t in otissues){
  otmp <- ogroups[which(grepl(t, ogroups))]
  oconds <- unique(lapply(strsplit(otmp, '.', fixed = T), '[[', 1))
  for (cd in oconds){
    otmp2 <- otmp[which(grepl(cd, otmp))]
    ocombined[[paste0(t, cd, sep = '.')]] <- otmp2
  }
}

org_fcs <- get_all_fc(ocombined, pkareaorg[,-c(1:32, ncol(pkareaorg))])

## getting indices of metabs with original fc <2 and >-2 (non-significant)
nonsig_org <- which(abs(org_fcs$value) < 2)
nonsig_new <- new_fcs[nonsig_org, 'value']

print(top_cutoff <- quantile(nonsig_new, 0.99))
print(bottom_cutoff <- quantile(nonsig_new, 0.01))
