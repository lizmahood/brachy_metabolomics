## Read in the updated with superclass file, and the file of identified compounds in the network
## Merge them, make a is_identified column, output this

argg <- commandArgs(T)

if (length(argg) != 3){
  stop('ARGS: 1) Node file updated with canopus superclass 2) File of 
       identified network compounds 3) Name and prefix of output file')
}

canopus_node <- read.table(argg[1], sep = '\t', fill = NA, header = T, quote = "")
idd <- read.table(argg[2], sep = '\t', quote = "", fill = NA, header = T)

canopus_node$is_identified <- ifelse(canopus_node$Alignment.ID %in% idd$Alignment.ID, 1, 0)


## Getting the name of the identifications of the IDd compounds
small_id <- idd[, c('Alignment.ID', 'Analog.Compound_Name', 'Name_MSDIAL')]
canopus_idd <- merge(canopus_node, small_id, by = 'Alignment.ID')

write.table(canopus_idd, file = paste0(argg[3], 'updated_identifications.tsv'),
            row.names = F, quote = F, sep = '\t')
print('Done!')