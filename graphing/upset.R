## get dir of all dam files, 
## list.files, pattern = 'FDR'
## also get a list of the condition names (basename on filename, substr and strsplit)
## save alignment IDs in files into a list, each element in list is name of condition
## make UpSet plot

library(UpSetR)
library(ComplexHeatmap)

argg <- commandArgs(T)

if (length(argg) != 4){
  stop('ARGS: 1) Directory containing all DAM files (for all conditions)
       2) Output name for upset plot 3) metabs OR trans 4) Do you want to output metabs/
       genes in certain sets? y OR n')
}

get_condname <- function(damfiles, typ){
  bases <- basename(damfiles)
  if (typ == 'metabs'){
    conds <- unlist(lapply(strsplit(bases, '_'), '[[', length(strsplit(bases, '_')) + 7))
    print(conds)
  }else conds <- bases
  out <- substr(conds, 1, nchar(conds) - 4)
  return(out)
}

get_alignment_ids <- function(damfiles, conds, typ){
  olist <- list()
  for (fil in 1:length(damfiles)){
    tmp <- read.table(damfiles[fil], sep = '\t', header = T, fill = NA)
    if (typ == 'metabs'){
      olist[[conds[fil]]] <- tmp$Alignment.ID
    }else if (typ == 'trans'){
      ntmp <- tmp[which(abs(tmp$log2FoldChange) >= 2),]
      if (length(which(is.na(ntmp$padj))) > 0){
        ntmp <- ntmp[-which(is.na(ntmp$padj)),]
      }
      if ('GeneID' %in% colnames(ntmp)) {
        olist[[conds[fil]]] <- ntmp$GeneID
      } else if ('X' %in% colnames(ntmp)){
        olist[[conds[fil]]] <- ntmp$X
      } else if (grepl('[a-z]', rownames(ntmp)[1])){
        olist[[conds[fil]]] <- rownames(ntmp)
      }
    }
  }
  return(olist)
}

make_colors <- function(conds){
  out <- c()
  for (cond in conds){
    if (grepl('Leaf|leaf', cond)){
      out <- c(out, 'forestgreen')
    }
    else out <- c(out, 'saddlebrown')
  }
  return(out)
}

typ <- argg[3]
if (typ == 'metabs'){
  damfiles <- list.files(argg[1], pattern = 'FDR', full.names = T)
}else damfiles <- list.files(argg[1], full.names = T)

## new 10.3.21 to plot leaves together and roots together
leafdam <- damfiles[grepl('Leaf|leaf', damfiles)]
rootdam <- damfiles[grepl('Root|root', damfiles)]

leafcondnames <- get_condname(leafdam, typ)
rootcondnames <- get_condname(rootdam, typ)
condnames <- c(leafcondnames, rootcondnames)

leafidlist <- get_alignment_ids(leafdam, leafcondnames, typ)
rootidlist <- get_alignment_ids(rootdam, rootcondnames, typ)

idlist <- c(leafidlist, rootidlist)
colorvec <- make_colors(condnames)
tisdata <- ifelse(grepl('Leaf|leaf', condnames), 'Leaf', 'Root')
metadata <- as.data.frame(cbind(condnames, tisdata))

##plotting
pdf(paste0(argg[2], '.pdf'), width = 8, height = 6)
#upset(fromList(idlist), sets = names(idlist), order.by = 'freq', keep.order = T, sets.bar.color = colorvec,
#      set.metadata = list(data = metadata, plots = list(list(type = "matrix_rows", 
#                          column = "tisdata", colors = c(Leaf = "forestgreen", Root = 'darkgoldenrod'), 
#                          alpha = 0.4))))

upset(fromList(idlist), sets = names(idlist), order.by = 'freq', keep.order = T, 
      number.angles = 30, sets.bar.color = colorvec, nintersects = NA,
      sets.x.label = 'DAMs per Condition', point.size = 1.5,
      mainbar.y.label = 'DAM Intersections', text.scale = c(1.3, 1, 1.3, 1, 1.3, 0.75),
      set.metadata = list(data = metadata, plots = list(list(type = "matrix_rows", 
                          column = "tisdata", colors = c(Leaf = "forestgreen", Root = 'darkgoldenrod'), 
                          alpha = 0.4))))
dev.off()

if (argg[4] == 'y'){
  if (class(idlist[[1]]) != 'character'){
    cidlist <- lapply(idlist, function(x) as.character(x))
  }else cidlist <- idlist
  combmat <- make_comb_mat(cidlist)
  all_combs <- comb_name(combmat)
  outn <- paste0(argg[1], '/', typ, '_set_')
  if ('0001000000' %in% all_combs) write.table(extract_comb(combmat, '0001000000'), quote = F, 
                                               sep = '\t', file = paste0(outn, 'heat_only.tsv'))
  if ('0000010000' %in% all_combs) write.table(extract_comb(combmat, '0000010000'), quote = F,
                                               sep = '\t', file = paste0(outn, 'heatcop_only.tsv'))
  if ('0100010000' %in% all_combs) write.table(extract_comb(combmat, '0100010000'), quote = F, 
                                               sep = '\t', file = paste0(outn, 'cop_and_heatcop.tsv'))
    
}


