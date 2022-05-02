library(matrixStats)

argg <- commandArgs(T)

if (length(argg) != 3){
  stop('ARGS: 1) Identification GNPS/MSDIAL and canopus file 
       2) canopus file (all canopus annotations)
       3) output file prefix')
}

ids <- read.table(argg[1], header = T, sep = '\t', fill = NA, quote = "")
canopus_full <- read.table(argg[2], header = T, sep = '\t', fill = NA, quote = "")

ids <- subset(ids, select = -c(molecularFormula, adduct, DAM_conditions,
                               Data_Collector, NewAlignment.ID, Ontology_MSDIAL))

## Column for if the compounds was identified by GNPS, MSDIAL, or both
ids$identification_source <- ifelse(ids$Analog.Compound_Name != ''& ids$Name_MSDIAL != '',
                                    'Both', ifelse(ids$Name_MSDIAL !='', 'MSDIAL',
                                           'GNPS'))


## Getting the highest score between the query and any database match
idsmat <- as.matrix(ids[, c(2,7,8)])
ids$match_score <- rowMaxs(idsmat, na.rm = T)
ids <- subset(ids, select = -c(Analog.MQScore, Score1_MSDIAL, Score2_MSDIAL))

## Column for name of compound (hoping MSDIAL doesn't disagree with canopus
## when they both ID it)
ids$match_name <- ifelse(ids$Analog.Compound_Name != '' & ids$Name_MSDIAL != '', 
                         paste0(ids$Analog.Compound_Name, '$$', ids$Name_MSDIAL),
                         ifelse(ids$Name_MSDIAL != '', ids$Name_MSDIAL,
                                ids$Analog.Compound_Name))
ids <- subset(ids, select = -c(Analog.Compound_Name, Name_MSDIAL))

## Merging with canopus
every_annot <- merge(ids, canopus_full, all = T, 
                     by.x = c('Alignment.ID', 'all.classifications','superclass', 
                              'class', 'subclass', 'level.5','most.specific.class'), 
                     by.y = c('name', 'all.classifications','superclass', 
                              'class', 'subclass', 'level.5','most.specific.class'))

## Removing not needed columns
every_annot <- subset(every_annot, select = -c(ClassyFy.Status, GNSP_INCHIKEY,
                                               Analog.INCHI, INCHIKEY_MSDIAL, molecularFormula,
                                               adduct))

## Fill NAs with 'None'
every_annot[is.na(every_annot)] <- ''

write.table(every_annot, file = paste0(argg[3], '_identifications_canopus_for_paper.tsv'),
            sep = '\t', row.names = F, quote = F)
