## For filtering CANOPUS output to include predictions only if they
## pass a confidence score threshold

argg <- commandArgs(T)

if (length(argg) != 3){
  stop('ARGS: 1) Canopus summary file 2) associated formula file 3) associated conf score file')
}

parse_id <- function(ids){
  print(head(ids))
  new <- unlist(lapply(strsplit(ids, 'scans'), '[[', 2))
  return(new)
}

canopus <- read.table(argg[1], sep = '\t', header = T, fill = NA, quote = "")
formula <- read.table(argg[2], sep = '\t', header = T, fill = NA, quote = "")
conf <- read.table(argg[3], sep = '\t', header = T, fill = NA, quote = "")

print('canopus')
canopus$name <- parse_id(canopus$name)
print(head(canopus$name))
print('hey! that worked')
print('formula')
formula$id <- parse_id(formula$id)

print(head(formula$id))
c_conf <- merge(canopus, conf, by.x = 'name', by.y = 'newnames')
f_conf <- merge(formula, conf, by.x = 'id', by.y = 'newnames')

c_conf <- c_conf[which(c_conf$prob >= 0.5),]
f_conf <- f_conf[which(f_conf$prob >= 0.5),]

c_conf$formula <- NULL; c_conf$class.y <- NULL; c_conf$prob <- NULL
f_conf$formula <- NULL; f_conf$class <- NULL; f_conf$prob <- NULL

write.table(c_conf, file = paste0(argg[1], '_filtered.tsv'), sep = '\t', quote = F, row.names = F)
write.table(f_conf, file = paste0(argg[2], '_filtered.tsv'), sep = '\t', quote = F, row.names = F)
