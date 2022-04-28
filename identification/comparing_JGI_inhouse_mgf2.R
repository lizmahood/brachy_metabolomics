## Linking up JGI's identified metabolites with ours

argg <- commandArgs(T)
if (length(argg) != 6){
  stop('ARGS: 1) file linking JGI ids to ours and with cosine score
       2) JGI GNPS cytoscape node metadata file
       3) JGI other id file (MSMS-blink-hits)
       4) Our normalized peak area file
       5) Directory to differentially abundant metabolites
       6) file of nodes in our MSFINDER network')
}

get_best_blink_id <- function(blinkid){
  if (class(blinkid$score) == 'character'){
    blinkid$score <- as.numeric(blinkid$score) 
  }
  tmp <- blinkid[order(blinkid$feature_id_query, -blinkid$score),]
  tmp <- tmp[-which(duplicated(tmp$feature_id_query)),]
  if (length(which(tmp$score <= 0.8)) > 0){
    tmp <- tmp[-which(tmp$score <= 0.8),]
  }
  tmp <- tmp[,c('name_ref', 'score', 'feature_id_query')]
  return(tmp)
}

get_gnps_id <- function(gnpsid){
  gnpsid <- gnpsid[,c('Compound_Name', 'MQScore', 'name')]
  if (length(which(gnpsid$MQScore <= 0.8)) > 0){
    gnpsid <- gnpsid[-which(gnpsid$MQScore <= 0.8),]
  }
  return(gnpsid)
}

get_condname <- function(damfiles){
  bases <- basename(damfiles)
  conds <- unlist(lapply(strsplit(bases, '_'), '[[', length(strsplit(bases, '_')) + 2))
  print(conds)
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

merge_all <- function(cosinelink, gnpsid, blinkid, peakheight, daml){
  ## filter matches between jgi and us, remove really bad ones
  cosinelink <- cosinelink[-which(cosinelink$cosine < 0.8),]
  
  ## which jgi ids are identified
  jgiblink <- merge(cosinelink, blinkid, by.x = 'JGI_ID', by.y = 'feature_id_query', all.x = T)
  jgiblinkgnps <- merge(jgiblink, gnpsid, by.x = 'JGI_ID', by.y = 'name', all.x = T)
  
  ## did any of the metabolites pass normalization/filters?
  tmp <- ifelse(jgiblinkgnps$our_id %in% peakheight$Alignment.ID, 'Yes', 'No')
  jgiblinkgnps$Survive_filter <- tmp
  
  ## Are any of our metabs differentially abundant?
  damtmp <- c()
  for (id in jgiblinkgnps$our_id){
    holder <- ''
    for (name in names(daml)){
      if (id %in% daml[[name]]){
        holder <- paste(holder, name, sep = '|')
      }
    }
    damtmp <- c(damtmp, holder)
  }
  
  jgiblinkgnps$DAM_conditions <- damtmp
  jgiblinkgnps$X <- NULL
  colnames(jgiblinkgnps)[3:7] <- c('cosine_jgi_us', 'Name_BLINK', 'Score_BLINK', 'Name_GNPS', 'Score_GNPS')
  return(jgiblinkgnps)
}


cosine_link <- read.table(argg[1], sep = '\t', header = T, fill = NA, quote = "")
gnps_id <- read.table(argg[2], sep = '\t', quote = "", fill = NA, header = T, stringsAsFactors = F)
blink_id <- read.table(argg[3], sep = '\t', header = T, quote = "", fill = NA, stringsAsFactors = F)
oursurviving <- read.table(argg[4], sep = '\t', header = T, fill = NA, quote = "")
network <- read.table(argg[6], header = T, fill = NA, quote = "", sep = '\t')

blinkid <- get_best_blink_id(blink_id)
gnpsid <- get_gnps_id(gnps_id)

damfiles <- list.files(argg[5], pattern = 'FDR', full.names = T)
condnames <- get_condname(damfiles)
idlist <- get_alignment_ids(damfiles, condnames)
idlist <- idlist[order(names(idlist))]

final <- merge_all(cosine_link, gnpsid, blinkid, oursurviving, idlist)
out <- merge(network, final, by.x = 'Alignment.ID', by.y = 'our_id')
out <- out[which(!(is.na(out$Score_GNPS)) | !(is.na(out$Score_BLINK))),]

write.table(final, file = paste0(argg[1], '_metabolite_annotations.tsv'), sep = '\t', row.names = F, quote = F)
write.table(out, file = paste0(argg[1], '_network_metablite_annotations.tsv'), sep = '\t', row.names = F, quote = F)
