###Prettier script for analyzing metabolomic data. Finds differentially
###expressed metabolites for treatments. Normalizes metabolites
###first, then finds differentially expressed ones between healthy and
###not. Does FDR correction.

library(matrixStats)
library(vsn)
library(MetaboDiff)
library(dplyr)


                  ###################      
                  #### Functions ####
                  ###################

get_groups <- function(df){
  ct <- colnames(df)
  groups <- unlist(lapply(strsplit(ct, '_'), '[', 1))
  return(groups)
}

checkmissing <- function(asycl, cutoff){
  ##Should not include blank
  bad <- c()
  for (i in 1:ncol(asycl)){
    cl <- asycl[,i]
    cl[cl == 0] <- NA
    print (paste0(colnames(asycl)[i], ' ', (length(which(is.na(cl))) / length(cl))))
    if ((length(which(is.na(cl))) / length(cl)) >= cutoff) bad <- c(bad, i)
  }

  if (length(bad) > 0){
    stop(paste0('ERROR! Samples ', paste0(colnames(asycl)[bad], collapse = ', '), 
                ' had more missing values than you allow!!'))
  }
}

checksn <- function(asycl, snm, snthrsam){
  ##Should not include blank
  snasy <- as.matrix(snm[,asycl])
  bad <- c()
  for (samp in 1:ncol(snasy)){
    print(paste0(colnames(snm)[asycl[samp]], ' ', mean(snasy[,samp])))
    if (mean(snasy[,samp]) < snthrsam) bad <- c(bad, samp)
  }
  
  if (length(bad) > 0){
    stop(paste0('ERROR! Samples ', paste0(colnames(snasy)[bad], collapse = ', '),
                ' had lower S/N than you allow!!'))
  }
}

make_blank_single <- function(asy){
  inblank <- c()
  bcol <- ncol(asy)
  
  for (i in 1:nrow(asy)){
    if (asy[i, bcol] > 50000 & (asy[i, bcol] / min(as.numeric(unname(asy[i,(1:bcol-1)])))) > 0.5){
      inblank <- c(inblank, i)
    } 
  }
  return(inblank)
}

make_blank_mult <- function(pkl, blkc, blkg){
  ##pkl is whole alignment table, not just asy
  ##blkc is numeric vector of blank column indices
  ##blkg is a string -- which blanks go with which columns
  
  ##parsing the blkg string into a list
  blkgroups <- unlist(strsplit(blkg, '[', fixed = T))
  blkgroups <- blkgroups[blkgroups != '']
  inblank <- c()
  
  for (grp in blkgroups){
    ##parsing the string
    grp <- gsub(']', '', grp, fixed = T)
    grplst <- strsplit(grp, '-', fixed = T)[[1]]
    blkcol <- as.numeric(grplst[1])
    asycols <- as.numeric(unlist(strsplit(grplst[2], ',')))
    for (i in 1:nrow(pkl)){
      if (pkl[i, blkcol] > 50000 & (pkl[i, blkcol] / min(as.numeric(unname(pkl[i, asycols])))) > 0.5){
        inblank <- c(inblank, i)
      }
    }  
  }
  return(unique(inblank))
}

imp_per_group <- function(df, groups, mdat){
  ##Output dataframe will have sample columns
  ##at the end of the dataframe. Blank should
  ##NOT be included at this step
  
  tdf <- as.data.frame(t(df))
  tdf$groups <- groups
  
  ##of non-zero/non-noise rows per group
  gil <- list()
  
  ##list of imputed values
  implist <- list()
  
  ##imputing for each group separately
  for (grp in unique(tdf$groups)){
    timp <- tdf[which(tdf$groups == grp),]
    imp <- t(timp)
    imp <- as.matrix(imp[-nrow(imp),])
    imp <- apply(imp, 2, as.numeric)
    imp[imp == 0] <- NA
    noise <- which(rowSums(is.na(imp)) >= (0.5 * ncol(imp)))
    good <- which(rowSums(is.na(imp)) < (0.5 * ncol(imp)))
    goodimp <- imp[-noise,]
    if (length(good) == 1){
      knnimp <- impute::impute.knn(t(as.matrix(goodimp)), k = 3, colmax = 1)$data
    } else knnimp <- t(impute::impute.knn(t(as.matrix(goodimp)), k = 3, colmax = 1)$data)
    gil[[grp]] <- good
    #names(gil)[length(names(gil))] <- grp
    implist[[grp]] <- round(knnimp)
  }
  
  ##making copy of df to remove rows from
  outdf <- df
  noiserws <- c()
  
  ##removing rows if they are zero/noise across ALL conditions
  for (rw in 1:nrow(df)){
    ## Is this row good in any condition? Record conditions if so.
    inn <- F
    conds <- c()
    for (vec in 1:length(gil)){
      if (rw %in% gil[[vec]]) {
        inn <- T
        conds <- c(conds, vec)
      }
    } 
    ## If zero/noise across all conds, add to noiserow.
    if (inn == F) {
      noiserws <- c(noiserws, rw)
      #outdf <- outdf[-rw,]
      #mdat <- mdat[-rw,]
    }else {
      ##if they are not zero/noise across ALL conditions,
      ##combine the zero/noise conditions with imputed conditions
      ##for each group it was good in
      goodcols <- c()
      goodvals <- c()
      ## Loop through conditions it was good in
      for (cnd in conds){
        ##get the columns of that group  
        cls <- which(tdf$group == names(gil)[cnd])
        goodcols <- c(goodcols, cls)
        ##get the imputed values for those groups
        impvals <- unname(implist[[cnd]][unname(which(gil[[cnd]] == rw)),])
        goodvals <- c(goodvals, impvals)
      }
      ##replace original values with these values
      for(cl in 1:length(goodcols)){
        outdf[rw, goodcols[cl]] <- goodvals[cl]
      }
      ## all other columns (noise ones) will be filled with 0
      zerovec <- c(seq(1, ncol(df)))
      zerovec <- zerovec[-goodcols]
      outdf[rw, zerovec] <- 0
    }
  }
  ## Remove zero/noise rows
  if (length(noiserws) != 0){
    outdf <- outdf[-noiserws,]
    mdat <- mdat[-noiserws,]
  }
  return(list(mdat, outdf))
}

vs_norm <- function(prepl, idcol, ccols, tcols, blkcol, sweights, sncol, lt, misscut,
                    snthresh, imp, telid, vsnn, isnm, snmatt, snsampt, blkgrp, qccols, introws, ph){
  
  #'prepl: peak list file
  #'idcol, sncol, qccols: column with Alignment IDs, col for SN, cols for qcs (if any, "" if none)
  #'c/tcols: columns of prepl containing samples values of metabolites in ctrl/treat conditions
  #'blkcol: column index of the blank
  #'sweights: vector of sample weights in the same order as samples.
  #'snthresh, snsampt: for removing lines with SN lower than this, for removing samples with sn lower than this
  #'misscut: for warning about samples with missing values more than this
  #'imp, vsnn, isnm, lt: booleans
  #'blkgrp: big string, to be parsed
  #'introws: numeric vector of IDs of interest
  #'snmatt, ph: additional files, sn matrix per sample, peak height
  #'
  #'First: S/N threshold. 2) noise threshold. 3: Blank filter. 4: sample weight 
  #'normalization. 5: Detect and remove replicate outliers via correlation 
  #'6: impute (if wanted) 7: vsn (if wanted)
  
  if (length(tcols) != 0){
    asycols <- c(ccols, tcols)
  }else asycols <- ccols
  #print(colnames(pl)[asycols])
  print(paste0('Total metabolites: ', nrow(prepl)))
  print(paste0('Total fragmented metabolites: ', length(which(prepl$MS.MS.spectrum != ''))))
  
  ##which rows have fragmented metabolites
  frags <- which(prepl$MS.MS.spectrum != '')
  
  ### CHECKS for %0 values and S/N values
  checkmissing(prepl[,asycols], misscut)
  checksn(asycols, snmatt, snsampt)
  
  ### 1) doing sn thresholding
  prepl[,sncol] <- as.numeric(prepl[,sncol])
  plr <- nrow(prepl)
  ##assuming internal std has survived SN threshold
  pl <- prepl
  snrow <- which(pl[,sncol] <= snthresh)
  #pl <- prepl[which(prepl[,sncol] > snthresh),]

  snfrags <- which(pl[snrow,]$MS.MS.spectrum != '')
  
  ##checking if IDs of interest were removed
  if (length(introws) > 0){
    badids <- pl[snrow, idcol]
    for (i in 1:length(introws)){
      if (!(introws[i] %in% badids)){
        print(paste0(introws[i], ' survived SN thresh!'))
        print(paste0('SN Threshod: ', snthresh))
        print(paste0('This metabolites SN: ', pl[which(pl[,idcol] == introws[i]), sncol]))
      }else{
        print(paste0(introws[i], ' was removed because of low SN'))
        print(paste0('SN Threshod: ', snthresh))
        print(paste0('This metabolites SN: ', prepl[which(prepl[,idcol] == introws[i]), sncol]))
      }
    }
  }
  
  ##getting internal standard row
  print(length(which(pl[,idcol] == telid)))
  telrow <- which(pl[,idcol] == telid)
  if (length(telrow) > 1) telrow <- telrow[1]
  #print(telrow)
  #print(pl[telrow,])
  
  ##getting sample values
  asy <- pl[,c(asycols, blkcol)]
  align <- pl[,1]
  #asy <- pl[,asycols]
  asy <- apply(asy, 2, as.numeric)
  #print('Correlation of raw values after s/n threshold')
  #print(cor(asy))
  asyr <-  nrow(asy)
  print(paste0('Num of metabolites removed from SN: ', length(snrow)))
  print(paste0("Num of fragmented metabolies removed: ", length(snfrags)))

  
  ### 2 Apply noise threshold: any metab with max height lower than 10k removed
  ### Should only do for non-qc normalized data
  phasy <- ph[,c(asycols, blkcol)]
  lowrow <- which(rowMaxs(as.matrix(phasy)) < 10000)
  badrow <- unique(c(snrow, lowrow))
  
  ##seeing if metabs of interest are removed.
  if (length(introws) > 0){
    for (i in 1:length(introws)){
      if (introws[i] %in% pl[lowrow,idcol]){
        print(paste0(introws[i], ' was removed because of low values'))
      }
    }
  }
  
  ##writing out correlation before blank
  cort <- cor(asy[-badrow,])
  print(paste0('Num of metabolites removed from low values: ', length(lowrow)))
  print(paste0('Num of fragmented metabolites removed: ', (length(which(pl[lowrow,]$MS.MS.spectrum != '')))))
  write.table(cort, file = paste0(args[1], '_cor_before_blank.tab'), sep = '\t', quote = F)

  ### 3) Blank filter: if average metabolite value in blank is >= 50% of
  ### min value in sample -- remove
  if (length(blkcol) == 1){
    inblank <- make_blank_single(asy)
  }else if (length(blkcol) > 1){
    inblank <- make_blank_mult(pl, blkcol, blkgrp)
  }
  
  ##seeing if metabs of interest are removed.
  if (length(introws) > 0){
    for (i in 1:length(introws)){
      if (introws[i] %in% pl[inblank,idcol]){
        print(paste0(introws[i], ' was in the blank'))
      }
    }
  }
  
  print(paste0('Num of metabolites in blank: ', length(inblank)))
  print(paste0('Num of fragmented metabolites removed: ', length(which(pl[inblank,]$MS.MS.spectrum != ''))))
  #print(inblank[which(inblank == telrow)])
  if (telrow %in% inblank) inblank <- inblank[-which(inblank == telrow)] ##need to keep internal stndard
  badrow <- c(badrow, inblank)
  badrow <- unique(badrow)
  #print('Correlation after blank filter')
  print(paste0('Num of metabolites removed from blank or low row or SN: ', length(badrow)))
  #print(colnames(asy))

  ### 4) Normalizing to sample weight
  for (i in 1:length(sweights)){
    asy[,i] <- asy[,i]/sweights[i]
  }
  
  ### 5) Normalizing to internal standard, if wanted
  if (isnm == 'yes'){
    intvals <- asy[telrow,]
    
    for (i in 1:length(intvals)){
      asy[,i] <- asy[,i]/intvals[i]
    }
    print('Done with intstd norm')
  }
  
  ##now need to add telrow to rows getting removed
  badrow <- c(badrow, telrow)
  badrow <- sort(badrow)
  if (length(blkcol) == 1){
    asy2 <- asy[-badrow, -ncol(asy)]
  } else {
    asy2 <- asy[-badrow, 1:(ncol(asy)-length(blkcol))]
  }
  
  ### 6) log transform (if wanted) added 11/17/20
  if (lt == 'yes'){
    asy2 <- log10(asy2 + 1)
  }
  
  ### 7) imputing missing values, ignoring blank
  if (imp == 'yes'){
    mdat <- pl[-badrow,]
    gps <- get_groups(asy2)
    knnlis <- imp_per_group(asy2, gps, mdat)
    mdat2 <- knnlis[[1]]
    knnasy <- knnlis[[2]]
    
    print(paste0('Num of metabolites removed from knn, noise or blank or SN: ',
                 nrow(prepl) - nrow(knnasy)))
    #print(length(badrow))
    #badrow <- c(badrow, missrow)
    #badrow <- unique(badrow)
    #print('Correlation after knn imputation')
    #print(cor(knnasy))
  }else{
    knnasy <- asy2
    gps <- get_groups(knnasy)
    #knnasy <- asy[-badrow,]
    mdat2 <- pl[-badrow,]
  }

  write.table(cor(knnasy), file = paste0(args[1], '_cor_after_blank_itsd_smpwght.tab'), sep = '\t', quote = F)
  
  ##8 removing non-correlated replicates
  
  ##getting list of colnames per group
  h <- list()
  tmp <-t(knnasy)
  tmp$groups <- gps
  for(gp in 1:length(unique(gps))){
    name <- unique(gps)[gp]
    h[[unique(gps)[gp]]] <- which(grepl(paste0('^',name), tmp$groups))
  }
  
  gdc <- c()
  for(i in 1:length(h)){
    cls <- h[[i]]
    corr <- cor(knnasy[,cls])
    
    ##which columns in h[[i]] have at least 3 correlations < 0.7
    fbad <- which(unname(colSums(corr < 0.7)) > 2)
    
    if (length(fbad) >= 1) {
      fgood <- h[[i]][-fbad]
    }
    else fgood <- h[[i]]
    gdc <- c(gdc, fgood)
  }
  
  knnasy <- knnasy[,gdc]
  
  newh <- list(h[[1]][h[[1]] %in% gdc], h[[2]][h[[2]] %in% gdc])
  newasy <- list(asycols[newh[[1]]], asycols[newh[[2]]])
  
  write.table(cor(knnasy), file = paste0(args[1], '_cor_low_cor_reps_removed.tab'), sep = '\t', quote = F)
  
  ##calculating fold change before normalization (Because vsn changes fold change)
  fldc <- log2(rowMeans(knnasy[,newh[[2]]]) / rowMeans(knnasy[,newh[[1]]]))
  
  ### 9) VSN (if wanted)
  ##ignoring rows that did not make it
  if (vsnn == 'yes'){
    vsnasy <- justvsn(knnasy)
    ##ordering columns
    vsnasy <- vsnasy[,order(colnames(vsnasy))]
    
    ##also doing vsn per experiment if all experients included
    if (any(startsWith(colnames(knnasy), 'Hydro')) &
        any(startsWith(colnames(knnasy), 'Sym')) &
        any(startsWith(colnames(knnasy), 'Tis'))){
      hydro <- justvsn(knnasy[,which(startsWith(colnames(knnasy), 'Hydro'))])
      sym <- justvsn(knnasy[,which(startsWith(colnames(knnasy), 'Sym'))])
      tis <- justvsn(knnasy[,which(startsWith(colnames(knnasy), 'Tis'))])
      separate_vsn <- cbind(hydro, sym, tis)
      ##ordering columns
      separate_vsn <- separate_vsn[,order(colnames(separate_vsn))]
    }
    
    #metadat <- mdat2[, -c(asycols)]
    metadat <- mdat2
    print('Correlation after vsn')
    write.table(cor(vsnasy), file = paste0(args[1], '_cor_after_vsn.tab'), sep = '\t', quote = F)
    if (exists('separate_vsn')) {
      write.table(cor(separate_vsn), file = paste0(args[1], '_cor_after_vsn_separate.tab'), sep = '\t', quote =F)
      return(list(newasy, data.frame(metadat[,1:32], vsnasy), data.frame(metadat[,1:32], separate_vsn)))
    }
    if (class(qccols) == 'numeric'){
      return(list(newasy, data.frame(metadat[,1:32], vsnasy, metadat[,qccols])))
    }else{
      return(list(newasy, data.frame(metadat[,1:32], vsnasy)))
    }
  }
  
  else {
    #metadat <- mdat2[, -c(asycols)]
    metadat <- mdat2
    if (class(qccols) == 'numeric'){
      return(list(newasy, data.frame(metadat[,1:32], knnasy[,order(colnames(knnasy))], metadat[,qccols], fldc)))
    }else{
      return(list(newasy, data.frame(metadat[,1:32], knnasy[,order(colnames(knnasy))], fldc)))
    }
  }
}

                #####################
                #### Coding time ####
                #####################


##parsing args
args <- commandArgs(T)
print(args)
if (length(args) != 20){
  stop('ARGS: 1) input peak area file, tsv 
       2) index of metab ID 
       3) indices of treated columns, separated by a comma 
       4) indices of control columns
       5) index of column with blank average(s) if multiple sepatate with comma
       6) vector of sample weights, separated by comma. First list
       weights of control columns, then treatment columns.
       7) column with s/n information 
       8) s/n cutoff 
       9) Peak Height file (to go with input area file)
       10) log transform? yes or no
       11) impute? yes OR no
       12) ID of internal standard
       13) Do you want vsn? yes OR no
       14) Do you want to normalize to an internal standard? yes OR no
       15) S/N matrix file 
       16) threshold for S/N value per individual sample
       17) blank groups. If multiple blanks, input as:
           [blank1-asycol1,asycol2][blank2-asycol4,asycol5]
       18) Qc columns, ("none" if there are none)
       19) IDs of metabolites to track, "none" for none
       20) Cutoff for missing values')
}


peakarea <- read.table(args[1], sep = '\t', header = T, stringsAsFactors = F, quote = "", fill = NA)
idcol <- as.numeric(args[2])
if (grepl(',', args[3], fixed = T)){
  treatcol <- as.numeric(unlist(strsplit(args[3], ',')))
}else treatcol <- 0
ctrlcol <- as.numeric(unlist(strsplit(args[4], ',')))

##are there multiple blanks?
if (grepl(',', args[5], fixed = T)) {
  blnkcol <- as.numeric(unlist(strsplit(args[5], ',', fixed = T)))
} else blnkcol <- as.numeric(args[5])

swt <- as.numeric(unlist(strsplit(args[6], ',')))
snc <- as.numeric(args[7])
sncut <- as.numeric(args[8])
peakheight <- read.table(args[9], sep = '\t', header = T, stringsAsFactors = F, quote = "", fill = NA)
logt <- as.character(args[10])
immp <- as.character(args[11])
tlid <- as.numeric(args[12])
varsn <- as.character(args[13])
isd <- as.character(args[14])
snmat <- read.table(args[15], sep = '\t', header = T, stringsAsFactors = F, quote = "", fill = NA)
snt <- as.numeric(args[16])
blankgroups <- args[17]

if (args[18] == 'none'){
  qccols <- ''
}else qccols <- as.numeric(unlist(strsplit(args[18], ',')))

if (args[19] == 'none'){
  introws <- c()
}else introws <- as.numeric(unlist(strsplit(args[19], ',')))
misscut <- as.numeric(args[20])
  
##normalizing
npeakarea_l <- vs_norm(prepl = peakarea, ccols = ctrlcol, tcols = treatcol,
                       blkcol = blnkcol, sweights = swt, sncol = snc, idcol = idcol,
                       snthresh = sncut, lt = logt, misscut = misscut, imp = immp, telid = tlid,
                       vsnn = varsn, isnm = isd, snmatt = snmat, snsampt = snt, blkgrp = blankgroups, 
                       qccols = qccols, introws = introws, ph = peakheight)

npeakarea <- as.data.frame(npeakarea_l[[2]])

##writing out normalized peak areas, of all surviving metabolites
outname <- c(args[1])
if (varsn == 'yes'){
  outname <- c(outname, '_VSN_')
}

if (isd == 'yes'){
  outname <- c(outname, '_IS_')
}

write.table(npeakarea, file = paste(c(outname, '_filtered_normalized_metabs.tab'), collapse = ''), 
            sep = '\t', row.names = F, quote = F)

if (length(npeakarea_l) ==3){
  write.table(npeakarea_l[[3]], file = paste(c(outname, '_vsn_separate.tab'), collapse = ''),
              sep = '\t', row.names = F, quote = F)
}

print('Done! :)')