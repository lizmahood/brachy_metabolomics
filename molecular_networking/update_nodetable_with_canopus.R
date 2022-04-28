## Read in mgf, get list of the alignment IDs.
## Read in node table (that came from that MGF!!), order querys
## cbind the alignment ID list with the node table
## Read in canopus, get df of just alignment ID and class
## Merge node table with small canopus table, keep all of node table
## Replace NA with None
## Write out new node table

argg <- commandArgs(T)

if (length(argg) != 3){
  stop('ARGS: 1) The mgf file that Mol Networking came from
       2) Node table (result of MSFINDER Mol Networking)
       3) Canopus summary file, made from same mgf as argg 1')
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

parse_id <- function(canopus){
  tl <- unlist(lapply(strsplit(canopus$name, '_'), '[[', 4))
  canopus$name <- substr(tl, 6, nchar(tl))
  return(canopus)
} 

mgfin <- read.table(argg[1], fill = NA, sep = '\t')
nodetable <- read.table(argg[2], sep = '\t', fill = NA, header = T)
canopus <- read.table(argg[3], sep = '\t', quote = "", fill = NA, header = T)

alignids <- get_id_from_mgf(mgfin)
nodetable <- reorder_query(nodetable)
nodetable$Alignment.ID <- alignids
#canopus <- parse_id(canopus)
smallcanopus <- canopus[,c('name', 'class', 'superclass')]
out <- merge(nodetable, smallcanopus, by.x = 'Alignment.ID', by.y = 'name', all.x = T)
out <- out[order(as.numeric(out$Alignment.ID)),]

out$class[which(is.na(out$class))] <- 'None'
out$superclass[which(is.na(out$superclass))] <- 'None'

write.table(out, file = paste0(argg[2], '_updated.tsv'), sep = '\t', row.names = F, quote = F)
