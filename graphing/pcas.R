##Now, for PCAs
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(tidyr)
library(dplyr)
library(MASS)
library(reshape2)
library(cowplot)
library(pheatmap)
library(RColorBrewer)

get_groups <- function(df, typ){
  ct <- colnames(df)
  groups <- unlist(lapply(strsplit(ct, '_'), '[', 1))
  return(groups)
}

make_pca <- function(metabs, all){
  ##metabs is a data frame
  conds <- c(); tiss <- c()
  print(colnames(metabs))
  for (col in colnames(metabs)){
    if (grepl('Leaf|leaf', col)){
      tiss <- c(tiss, 'Leaf')
    }else if (grepl('Root', col)){
      tiss <- c(tiss, 'Root')
    }else if (grepl('Culm', col)){
      tiss <- c(tiss, 'Culm')
    }else if (grepl('Spike', col)){
      tiss <- c(tiss, 'Spike')
    }
    
    name <- strsplit(col, '_', fixed = T)[[1]][1]
    if (all == 'no'){
      conds <- c(conds, strsplit(name, '.', fixed = T)[[1]][1])
    }else if (all == 'yes'){
      conds <- c(conds, paste(strsplit(name, '.', fixed = T)[[1]][c(1,2)], collapse = '.'))
    }
  }
  
  tmetabs <- as.data.frame(t(metabs))
  tiss <- as.factor(tiss); conds <- as.factor(conds)
  tmetabs <- cbind(tmetabs, tiss, conds)
  pca1 <- PCA(tmetabs, graph = F, quali.sup = c(ncol(tmetabs), (ncol(tmetabs) -1)))
  return(list(pca1, tmetabs))
}

graph_pca <- function(pca1, tmetabs, odir, color_vec, tistog){
  tmetabs$pc1 <- pca1$ind$coord[, 1]
  tmetabs$pc2 <- pca1$ind$coord[, 2]
  varexp1 <- round(pca1$eig[1,2], 2)
  varexp2 <- round(pca1$eig[2,2], 2)
  
  pca.vars <- pca1$var$coord %>% data.frame
  pca.vars$vars <- rownames(pca.vars)
  print(head(pca.vars))
  pca.vars.m <- melt(pca.vars, id.vars = "vars")
  
  
  if (tistog == 'yes'){
    p <- ggplot(tmetabs, aes(pc1, pc2, color = conds, shape = tiss, size = 1, label = rownames(tmetabs))) +
      geom_point(alpha = 0.8) + geom_hline(yintercept = 0, lty = 2) + scale_color_manual(values = color_vec) +
      geom_vline(xintercept = 0, lty = 2) + theme_bw() + guides(size = 'none', label = 'none',
                                                                shape = guide_legend(override.aes = list(size = 3)),
                                                                color = guide_legend(override.aes = list(size = 3))) + 
      xlab(paste0('PC 1 (', varexp1, '%)')) + ylab(paste0('PC 2 (', varexp2, '%)')) + 
      scale_fill_manual(values = color_vec) + theme(legend.position = 'top', legend.box = 'vertical', aspect.ratio = 1)
  
    pdf(paste0(odir, '_pca.pdf'), height = 6, width = 6)
    print(p)
    dev.off()
  }else {
    if (grepl('ROOT', odir)) {
      shap = 15
    }else if (grepl('LEAF', odir)){
      shap = 17
    }
    p <- ggplot(tmetabs, aes(pc1, pc2, color = conds, shape = tiss, size = 1, label = rownames(tmetabs))) +
      geom_point(alpha = 0.8, shape = shap) + geom_hline(yintercept = 0, lty = 2) + 
      scale_color_manual(values = color_vec) +
      geom_vline(xintercept = 0, lty = 2) + theme_bw() + guides(size = 'none', label = 'none',
                                                                shape = guide_legend(override.aes = list(size = 3)),
                                                                color = guide_legend(override.aes = list(size = 3))) + 
      xlab(paste0('PC 1 (', varexp1, '%)')) + ylab(paste0('PC 2 (', varexp2, '%)')) + 
      scale_fill_manual(values = color_vec) + theme(legend.position = 'top', legend.box = 'vertical', aspect.ratio = 1)
    
    pdf(paste0(odir, '_pca.pdf'), height = 6, width = 6)
    print(p)
    dev.off()
  }
}

argg <- commandArgs(T)
if (length(argg) != 5){
  stop('ARGS: 1) path to input filtered PeakArea file 
       2) desired output directory for PCA
       3) All experiments? yes OR no
       4) Do you want tissues (Rt/ Lf) to be together? yes OR no
       5) What type of data? metab OR trans')
}

infil <- read.table(argg[1], header = T, sep = '\t', fill = NA, quote = "", stringsAsFactors = F)
odir <- argg[2]
all <- argg[3]
tistog <- argg[4]
typ <- argg[5]
if (typ == 'metab'){
  colstouse <- c(33:ncol(infil))
  metabs <- infil[,colstouse]
  groups <- get_groups(metabs, typ)
  ## new functionality - correlation heatmap 3.11.21
  sampleDists <- as.dist(1-cor(metabs))
  sampledistmat <- as.matrix(sampleDists)
  rownames(sampledistmat) <- NULL
  mycolors <- colorRampPalette( rev(brewer.pal(9, "YlGnBu")) )(255)
  pdf(paste0(odir, '_correlation_heatmap.pdf'), width = 11)
  pheatmap(sampledistmat, clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists, col = mycolors)
  dev.off()
  
}else if (typ == 'trans'){
  metabs <- infil[,2:ncol(infil)]
  colnames(metabs) <- sub('_', '.', colnames(metabs), fixed = T)
  colnames(metabs) <- sub('_', '.', colnames(metabs), fixed = T)
  groups <- get_groups(metabs, typ)
}


if (tistog == 'yes'){
  pcalist <- make_pca(metabs, all)
  pca1 <- pcalist[[1]]; tmetabs <- pcalist[[2]]
  ## new colors for publication 9/30/21
  cvec <- c('salmon', 'mediumslateblue','palegreen', 'goldenrod2',
            'lightblue', 'forestgreen', 'plum3', 'darkorchid4', 'pink',
            'azure4', 'black', 'peru')
  graph_pca(pca1, tmetabs, odir, cvec, tistog)
}else if (tistog == 'no'){
  rootmat <- metabs[,which(grepl('Root', colnames(metabs)))]
  leafmat <- metabs[,which(grepl('Leaf', colnames(metabs)))]
  pcarootlist <- make_pca(rootmat, all)
  pcaleaflist <- make_pca(leafmat, all)
  pcar1 <- pcarootlist[[1]]; trootmat <- pcarootlist[[2]]
  pcal1 <- pcaleaflist[[1]]; tleafmat <- pcaleaflist[[2]]
  lvec <- c('salmon', 'mediumslateblue', 'palegreen', 'goldenrod2',
            'lightblue','forestgreen', 'plum3', 'pink')
  rvec <- c('salmon', 'mediumslateblue','palegreen', 'goldenrod2',
             'lightblue', 'forestgreen','plum3')
  graph_pca(pcar1, trootmat, paste0(odir, 'ROOT'), rvec, tistog)
  graph_pca(pcal1, tleafmat, paste0(odir, 'LEAF'), lvec, tistog)
}

print('Done!')
