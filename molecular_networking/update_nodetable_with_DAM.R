## Modifying MSFINDER mol networking node file to reflect if a metab is DA in each stress

## Get node table, matching msp folder and node file
## Match alignment ID to querys
## Get list of alignment IDs DE in each stress
## For each stress:
##  Make a dummy col in node file that is 1 if ID in stress, and 0 if not
## Output that

argg <- commandArgs(T)

if (length(argg) != 3){
  stop('ARGS: 1) Node table 2) mgf used to make MSP folder 
       4) Dir of metabolites DE in each stress')
}

get_id_from_mgf <- function(mgf){
  alignids <- c()
  for (i in 1:nrow(mgf)){
    if (grepl('SCANS', mgf$V1[i])){
      scan <- unlist(strsplit(mgf$V1[i], '='))[2]
      alignids <- c(alignids, scan)
    }
  }
  return(alignids)
}

reorder_query <- function(node){
  node$tmp <- node$Title
  node$tmp <- as.numeric(unlist(lapply(strsplit(node$tmp, '_'), '[[', 2)))
  node <- node[order(node$tmp),]
  node$tmp <- NULL
  return(node)
}

get_condname <- function(damfiles){
  bases <- basename(damfiles)
  conds <- unlist(lapply(strsplit(bases, '_'), tail, n = 1L))
  out <- substr(conds, 1, nchar(conds) - 4)
  return(out)
}

get_alignment_ids <- function(damfiles, conds){
  olist <- list()
  for (fil in 1:length(damfiles)){
    tmp <- read.table(damfiles[fil], sep = '\t', header = T, fill = NA)
    olist[[conds[fil]]] <- tmp$Alignment.ID
  }
  return(olist)
}

make_dummy_cols <- function(idlist, nodefile){
  tmpv <- matrix(0, nrow = nrow(nodefile), ncol = length(idlist))
  for (subl in 1:length(idlist)){
    print(names(idlist)[subl])
    dummycol <- ifelse(nodefile$Alignment.ID %in% idlist[[subl]], 1, 0)
    tmpv[,subl] <- dummycol
  }
  colnames(tmpv) <- names(idlist)
  tmpdam <- ifelse(rowSums(tmpv) == 0, 0, 1)
  out <- cbind.data.frame(nodefile, tmpv, tmpdam)
  return(out)
}

               ###################################################################
                                          ### CODING ###

mgf <- read.table(argg[2], sep = '\t', fill = NA, header = F)
nodetable <- read.table(argg[1], sep = '\t', header = T, fill = NA, quote = "")

alignmentids <- get_id_from_mgf(mgf)
nodetable <- reorder_query(nodetable)
nodetable$Alignment.ID <- alignmentids

damfiles <- list.files(argg[3], pattern = 'FDR', full.names = T)
condnames <- get_condname(damfiles)
idlist <- get_alignment_ids(damfiles, condnames)
newnode <- make_dummy_cols(idlist, nodetable)
write.table(newnode, file = paste0(argg[1], '_updated_DAM.tsv'), sep = '\t', row.names = F, quote = F)
print('Done!')
