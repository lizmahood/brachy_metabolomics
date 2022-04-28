library(xcms)
library(RColorBrewer)
library(SummarizedExperiment)
library(pander)
library(magrittr)
library(doParallel)

##This script is for getting the TICs of all files in the three experiments

cores<-detectCores()
#Create cluster with desired number of cores, leave one open for the machine         
#core processes
cl <- makeCluster(cores[1]-1)
#Register cluster
registerDoParallel(cl)

get_tic <- function(fils, datt, outp){
  ##fils is character vector of paths to files to read in
  ##datt is a phenodata df
  ##outp is string, out path of plots
  raw <- readMSData(files = fils, pdata = new('NAnnotatedDataFrame', datt), mode = 'onDisk')
  print('Done reading In')
  tics <- chromatogram(raw, aggregationFun = "sum")
  typs <- length(unique(datt$conds))
  group_c <- brewer.pal(typs, 'Set3')
  #group_c <- colorRamps::primary.colors(length(unique(datt$conds)))
  names(group_c) <- unique(datt$conds)
  ##plotting TIC
  pdf(paste0(outp, '/TIC.pdf'))
  plot(tics, col = group_c[datt$conds])
  legend("topright", legend = names(group_c), fill = unname(group_c), ncol = 3)
  dev.off()
  ##plotting zoomed in
  pdf(paste0(outp, '/zoomed_in_TIC.pdf'))
  plot(tics, col = group_c[datt$conds], ylim = c(0, 2.5e+08))
  legend("topright", legend = names(group_c), fill = unname(group_c), ncol = 3)
  dev.off()
  ##plotting deadvol removed
  pdf(paste0(outp, '/no_deadvol_TIC.pdf'))
  plot(tics, col = group_c[datt$conds], xlim = c(90, 960))
  legend('topright', legend = names(group_c), fill = unname(group_c), ncol = 3)
  dev.off()
}

make_sample_nam <- function(nam, pat){
  ##nam is character vector
  ##pat is pattern to find name after
  ##RETURNS, e.g.: Spore-Leaf_4
  out <- c()
  for (i in 1:length(nam)){
    if (pat == 'all'){
      sampn <- paste(strsplit(basename(nam[i]), '_')[[1]][c(1,2)], collapse = '_') 
    }else{
      bignam <- strsplit(nam[i], pat, fixed = T)[[1]][2]
      sampn <- paste(strsplit(bignam, '_', fixed = T)[[1]][c(1,2)], collapse ='_')
    }
    out <- c(out, sampn)
  }
  return(out)
}

make_condition <- function(snams, typ){
  ##snams is a vector of sample names
  ##typ is string -- are we doing sym, copper, or ctrl?
  conds <- c()
  print(typ)
  for (i in 1:length(snams)){
    print(snams[i])
    if (typ == 'Sym-' | (grepl('Sym-', snams[i], fixed =T) & typ == 'all')){
      print('Symbiosis sample!')
      if (grepl('Leaf', snams[i], fixed = T)) tmp_t <- 'Sym_L'
      else if (grepl('Root',snams[i], fixed = T)) tmp_t <- 'Sym_R'
      if (grepl('W',snams[i], fixed = T)) tmp_c <- 'SW'
      else if (grepl('Spore',snams[i], fixed = T)) tmp_c <- 'S'
      if (grepl('Ex', snams[i], fixed = T)) tmp_c <- 'BLK'
      else if (grepl('Ctrl', snams[i], fixed = T)) tmp_c <- 'CTRL'
      out <- paste(tmp_t, tmp_c, sep = '_')
      conds <- c(conds, out)
    }else if (typ == 'Tis-' | (grepl('Tis-', snams[i], fixed = T) & typ == 'all')){
      print('Tissue sample!')
      if (grepl('Culm', snams[i], fixed = T)) tmpt <- 'Tis_Culm'
      else if (grepl('Leaf', snams[i], fixed = T)) tmpt <- 'Tis_Leaf'
      else if (grepl('Spike', snams[i], fixed = T)) tmpt <- 'Tis_Spklt'
      else if (grepl('ExCt', snams[i], fixed = T)) tmpt <- 'Tis_BLK'
      conds <- c(conds, tmpt)
    }else if (typ == 'Hydro' | (grepl('Hydro-', snams[i], fixed = T) & typ == 'all')){
      print('Hydro sample!')
      if (grepl('HeatRoot', snams[i], fixed = T)) tmpt <- 'Hyd_BLK_HR'
      else if (grepl('Blank', snams[i], fixed = T)) tmpt <- 'Hyd_BLK_H'
      else if (grepl('HeatNoCop', snams[i], fixed = T)){
        if (grepl('Leaf', snams[i], fixed = T)) tmpt <- 'Hyd_HCP_L'
        else tmpt <- 'Hyd_HCP_R'
      }else if (grepl('Ctrl', snams[i], fixed = T)){
        if (grepl('Leaf', snams[i], fixed = T)) tmpt <- 'Hyd_CTRL_L'
        else tmpt <- 'Hyd_CTRL_R'
      }else if (grepl('Cop', snams[i], fixed = T)){
        if (grepl('Leaf', snams[i], fixed = T)) tmpt <- 'Hyd_CP_L'
        else tmpt <- 'Hyd_CP_R'
      }else if (grepl('Heat', snams[i], fixed = T)){
        if (grepl('Leaf', snams[i], fixed = T)) tmpt <- 'Hyd_H_L'
        else tmpt <- 'Hyd_H_R'
      }
      conds <- c(conds, tmpt)
    }else if (grepl('ExCtrl', snams[i])){
      print('Blank sample!')
      tmpt <- 'Hyd_BLK_Add'
      conds <- c(conds, tmpt)
    }
  }
  return(conds)
}

make_pd <- function(snams, conds){
  print(conds)
  return(data.frame(snams, conds, stringsAsFactors = F))
}

# 
# argg <- commandArgs(T)
# 
# if (length(argg) != 1){
#   stop('ARG: Is this for All experiments combined? yes OR no')
# }
# 
# if (argg[1] == 'no'){
#   conditions <- c('Copper_Heat', 'Symbiosis', 'Tissue_ctrl')
#   pats <- c('Hydro-', 'Sym-', 'Tis-')
#   modes <- c('FPS', 'NEG', 'POS')
# }else{
#   conditions <- c('All_exps')
#   pats <- c('all')
#   modes <- c('NEG', 'POS')
# }

modes <- c('NEG', 'POS')
#bigpath <- 'E:/MS_Data/BrachyMetabolites/original_files'
bigpath <- 'E:/MS_Data/BrachyMetabolites/cuheat/second_cuheat_metabolomics/original_files'
for (cnds in 1:length(conditions)){
  for (mod in modes){
    #filpth <- paste(bigpath, conditions[cnds], mod, sep = '/')
    filpth <- paste(bigpath, mod, sep = '/')
    fils <- list.files(filpth, full.names = T, pattern = 'mzML')
    print(head(fils))
    #snames <- make_sample_nam(fils, pats[cnds])
    snames <- make_sample_nam(fils, 'all')
    print(head(snames))
    #condts <- make_condition(snames, pats[cnds])
    condts <- make_condition(snames, 'Hydro')
    print(head(condts))
    pddf <- make_pd(snames, condts)
    get_tic(fils, pddf, filpth)
    print(paste0('Done with ', mod))
  }
}

