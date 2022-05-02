## For getting identifications from MSDIAL, and combining them with 
## IDs from JGI and canopus annotations

argg <- commandArgs(T)

if (length(argg) != 6){
  stop('ARGS: 1) Fully processed and filtered peak area file, with identifications
       2) Canopus file 3) GNPS\' .graphml cytoscape node file 4) folder of DAMs
       5) file of our MSFINDER network nodes
       6) output folder and file name prefix')
}

get_msdial_id <- function(filtered_df){
  testy <- filtered_df

  ## got to change these columns to numeric
  testy$Dot.product[testy$Dot.product == 'null'] <- 0
  testy$Reverse.dot.product[testy$Reverse.dot.product == 'null'] <- 0
  
  testy$Dot.product <- as.numeric(testy$Dot.product)
  testy$Reverse.dot.product <- as.numeric(testy$Reverse.dot.product)
  
  ## Which ones have IDs?
  testy_id <- testy[which(testy$Dot.product > 80 & testy$Reverse.dot.product > 80),]
  
  ## Returning the ID'd ones
  return(testy_id[, c('Alignment.ID', 'Metabolite.name', 'Ontology',
                      'INCHIKEY', 'Dot.product', 'Reverse.dot.product')])
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

combine_all <- function(canopus, gnps_id, msdial_id){
  ## Changing some column names of msdial_id
  colnames(msdial_id)[c(2,3,4,5,6)] <- c('Name_MSDIAL', 'Ontology_MSDIAL',
                                         'INCHIKEY_MSDIAL', 'Score1_MSDIAL', 
                                         'Score2_MSDIAL')
  
  ## merging a bunch
  gnps_msdial <- merge(gnps_id, msdial_id, by.x = 'name', by.y = 'Alignment.ID',
                      all = T)
  
  gnps_msdial_canopus <- merge(gnps_msdial, canopus, by.x = 'name', by.y = 'name',
                              all.x = T)
  colnames(gnps_msdial_canopus)[1] <- 'Alignment.ID'
  
  return(gnps_msdial_canopus)  
}


## READING IN
filtered_df <- read.table(argg[1], sep = '\t', header = T, fill = NA, quote = "")
canopus <- read.table(argg[2], sep = '\t', header = T, fill = NA, quote = "")
gnps_id <- read.table(argg[3], sep = '\t', header = T, fill = NA, quote = "")
damfiles <- list.files(argg[4], pattern = 'FDR', full.names = T)
network <- read.table(argg[5], header = T, fill = NA, quote = "", sep = '\t')

## Getting only columns needed from id_df
id_df <- filtered_df[, c(1, 2, 3, 4, 12, 13, 26, 27, 31, 32)]

gnps_id$Analog.MQScore <- as.numeric(gnps_id$Analog.MQScore)
gnps_id <- gnps_id[which(gnps_id$MQScore > 0.8),]
gnps_id <- gnps_id[, c('name', 'Analog.MQScore', 'Analog.Compound_Name',
                       'Analog.INCHI', 'Data_Collector')]

msdial_id <- get_msdial_id(filtered_df)

condnames <- get_condname(damfiles)
idlist <- get_alignment_ids(damfiles, condnames)
idlist <- idlist[order(names(idlist))]

all_out <- combine_all(canopus, gnps_id, msdial_id)

damtmp <- c()
for (id in all_out$Alignment.ID){
  holder <- ''
  for (name in names(idlist)){
    if (id %in% idlist[[name]]){
      holder <- paste(holder, name, sep = '|')
    }
  }
  damtmp <- c(damtmp, holder)
}

all_out$DAM_conditions <- damtmp
all_out_nodes <- merge(all_out, network, by.x = 'Alignment.ID', 
                       by.y = 'Alignment.ID')

write.table(all_out, file = paste0(argg[6], '_gnps_msdial_canopus_id.tsv'),
            sep = '\t', row.names = F, quote = F)

write.table(all_out_nodes, file = paste0(argg[6], '_network_metablite_annotations.tsv'), 
            sep = '\t', row.names = F, quote = F)
