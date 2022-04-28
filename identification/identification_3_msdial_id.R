## For getting identifications from MSDIAL, and combining them with 
## IDs from JGI and canopus annotations

argg <- commandArgs(T)

if (length(argg) != 6){
  stop('ARGS: 1) Fully processed and filtered peak area file 2) New peak
       area file with MSDIAL identifications 3) Canopus file 4) JGI
       identifications file 5) folder of DAMs
       6) output folder and file name prefix')
}

get_msdial_id <- function(id_df, filtered_df){
  testy <- merge(id_df, filtered_df, by = c('Average.Rt.min.', 'Average.Mz',
                                            'MS1.isotopic.spectrum', 'MS.MS.spectrum'))
  testy <- testy[ , !names(testy) %in% c('Metabolite.name.y', 'Ontology.y', 'SMILES.y',
                                         'Dot.product.y', 'Reverse.dot.product.y')]
  
  ## got to change these columns to numeric
  testy$Dot.product.x[testy$Dot.product.x == 'null'] <- 0
  testy$Reverse.dot.product.x[testy$Reverse.dot.product.x == 'null'] <- 0
  
  testy$Dot.product.x <- as.numeric(testy$Dot.product.x)
  testy$Reverse.dot.product.x <- as.numeric(testy$Reverse.dot.product.x)
  
  ## Which ones have IDs?
  testy_id <- testy[which(testy$Dot.product.x > 80 & testy$Reverse.dot.product.x > 80),]
  
  ## Returning the ID'd ones
  colnames(testy_id)[5] <- 'NewAlignment.ID'
  colnames(testy_id)[11] <- 'Alignment.ID'
  
  return(testy_id[, c(5,6,7,8,9,10,11)])
}

combine_all <- function(canopus, jgi_id, msdial_id){
  ## Changing some column names of msdial_id
  colnames(msdial_id)[c(2,3,4,5,6)] <- c('Name_MSDIAL', 'Ontology_MSDIAL',
                                        'SMILES_MSDIAL', 'Score1_MSDIAL', 
                                        'Score2_MSDIAL')
  
  ## merging a bunch
  jgi_msdial <- merge(jgi_id, msdial_id, by.x = 'our_id', by.y = 'Alignment.ID',
                      all = T)
  
  jgi_msdial_canopus <- merge(jgi_msdial, canopus, by.x = 'our_id', by.y = 'name',
                              all.x = T)
  colnames(jgi_msdial_canopus)[1] <- 'Alignment.ID'
  jgi_msdial_canopus <- jgi_msdial_canopus[, -8]
  
  return(jgi_msdial_canopus)  
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

## Reading in
filtered_df <- read.table(argg[1], sep = '\t', header = T, fill = NA, quote = "")
id_df <- read.table(argg[2], sep = '\t', header = T, fill = NA, quote = "")
canopus <- read.table(argg[3], sep = '\t', header = T, fill = NA, quote = "")
jgi_id <- read.table(argg[4], sep = '\t', header = T, fill = NA, quote = "")
damfiles <- list.files(argg[5], pattern = 'FDR', full.names = T)

## Getting only columns needed from id_df
id_df <- id_df[, c(1, 2, 3, 4, 12, 14, 26, 27, 31, 32)]

## Now getting only identified metabs from JGI id
jgi_id <- jgi_id[which(jgi_id$Score_BLINK > 0.8 | jgi_id$Score_GNPS > 0.8),]

## Using functions
msdial_id <- get_msdial_id(id_df, filtered_df)

condnames <- get_condname(damfiles)
idlist <- get_alignment_ids(damfiles, condnames)
idlist <- idlist[order(names(idlist))]

damtmp <- c()
for (id in msdial_id$Alignment.ID){
  holder <- ''
  for (name in names(idlist)){
    if (id %in% idlist[[name]]){
      holder <- paste(holder, name, sep = '|')
    }
  }
  damtmp <- c(damtmp, holder)
}

msdial_id$DAM_conditions <- damtmp

all_out <- combine_all(canopus, jgi_id, msdial_id)

write.table(all_out, file = paste0(argg[6], '_jgi_msdial_canopus_id.tsv'),
            sep = '\t', row.names = F, quote = F)
