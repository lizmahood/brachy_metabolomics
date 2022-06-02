library(matrixStats)
library(vsn)
library(MetaboDiff)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(grid)

get_fc <- function(ctrlv, stressv){
  nctm <- mean(as.numeric(ctrlv))
  if (nctm <= 0.01) {nctm <- 1e-30}
  nstm <- mean(as.numeric(stressv))
  ## new 9.30.21
  if (nstm <= 0.01) {nstm <- 1e-30}
  logfc <- log2(nstm/nctm)
  #ifelse(nstm <= 0.01, logfc <- 0, logfc <- log2(nstm/nctm))
  return(logfc)
}

get_diff_exp <- function(pl, ctrl, stress, normpl, normctrl, normstress, topcut, bottomcut) {
  #'Ctrl and stress are both data frames with metabolite peak areas
  #'Ctrl has peak areas for control and stress has areas for infected
  #'pl is peaklist
  #'This function finds parametric and non-param p-values for
  #'each metabolite, then does FDR with Benjamani-Hochberg, with 
  #'adjusted p-value of 0.05
  #'
  #'Outputs: vector of rows containing differentially expressed
  #'metabolites for this stress.
  
  type <- data.frame('group' = c(rep('Control', ncol(normctrl)), rep('stress', ncol(normstress))))
  wpvals <- c()
  tpvals <- c()
  nfc <- c(); fc <- c(); cmean <- c(); tmean <- c()
  
  ##getting Wilcoxon and t-test p-values for each metabolite
  for (metab in 1:nrow(normctrl)){
    ## getting normalized fold change
    nfc <- c(nfc, get_fc(normctrl[metab,], normstress[metab,]))
    ## non-normalized fold change
    fc <- c(fc, get_fc(ctrl[metab,], stress[metab,]))
    cmean <- c(cmean, mean(as.numeric(ctrl[metab,])))
    tmean <- c(tmean, mean(as.numeric(stress[metab,])))
    
    vals <- as.numeric(unlist(c(normctrl[metab,], normstress[metab,])))
    if ((all(vals == 0.01)) | (all(vals < 0)) ){
      wpvals <- c(wpvals, 1)
      tpvals <- c(tpvals, 1)
    }else {
      ndf <- cbind(vals, type)
      ##getting wilcoxon p-values for this metab 
      wpvals <- c(wpvals, wilcox.test(vals ~ group, data = ndf)[[3]])
      
      tpval <- tryCatch({t.test(vals ~ group, data = ndf)[[3]]}, error = function(cond){return(1)})
      tpvals <- c(tpvals, tpval)
      
    }
  }
  ##finding Benjamani-Hochberg fdr corrected metabs
  diffmetab <- data.frame(pval = tpvals, wpval = wpvals,
                          adj_pval = p.adjust(tpvals, method = 'fdr'), 
                          nfold_change = nfc, org_fc = fc, 
                          ctrlmean = cmean, treatmean = tmean)

  ##which metabolites are Diff abund?
  low <- which((diffmetab$pval <= 0.05) & 
                 ((diffmetab$nfold_change >= topcut & diffmetab$org_fc > 2) | 
                    (diffmetab$nfold_change <= bottomcut & diffmetab$org_fc < -2)) &
                 (diffmetab$wpval <= 0.1 | diffmetab$wpval == min(diffmetab$wpval)) &
                 (diffmetab$ctrlmean >= 5000 | diffmetab$treatmean >= 5000))
  
  ##if there are any diff expressed metabolites before fdr, are there any after fdr?
  if (length(low) > 0){
    print('There are diff abund metabs before FDR!')
    print(length(low))
    
    ##setting up pvalue for adj p-values.
    pct <- 0.05
    adjlow <- which((diffmetab$adj_pval <= pct) & 
                      ((diffmetab$nfold_change >= topcut & diffmetab$org_fc > 2) | 
                         (diffmetab$nfold_change <= bottomcut & diffmetab$org_fc < -2)) &
                      (diffmetab$ctrlmean >= 5000 | diffmetab$treatmean >= 5000) &
                      diffmetab$wpval <= 0.1)
    
    if (length(adjlow) >= 1) {
      print('Got some FDR Diff abund metabolites! ')
      print(length(adjlow))
      out <- data.frame(cbind(normpl[adjlow,c(1,2,3,5)], normctrl[adjlow,],
                              normstress[adjlow,], diffmetab[adjlow,]))
      return(list(data.frame(cbind(normpl, diffmetab)), out))
      
    }else {
      print('No FDR Diff abund metabolites :(')
      #out <- data.frame(cbind(normpl[low,c(1,2,3)], normctrl[low,],
      #                        normstress[low,], diffmetab[low,]))
      return(data.frame(cbind(normpl, diffmetab)))
    }
    
  }else {print('No diff abund metabolites at all! :(((')}
  
  return(data.frame(normpl, diffmetab))
}

get_groups <- function(df){
  ct <- colnames(df)
  groups <- unlist(lapply(strsplit(ct, '_'), '[', 1))
  return(groups)
}

##Need to put in arguments
argg <- commandArgs(T)

if (length(argg) != 6){
  stop('ARGS: 1) Input VSN normalized peak area file
       2) Input peak area with KNN but no normalization
       3) output name for volcano plots 
       4) is this all experiments combined? yes OR no
       5) 99th percentile of normalized fold change when -2 < non-normalized FC > 2
       6) 0.01 percentile of above')
}

npeakarea <- read.table(argg[1], sep = '\t', fill = NA, quote = "", header = T)
pkareaorg <- read.table(argg[2], sep = '\t', fill = NA, quote = "", header = T)
odir <- argg[3]
allg <- argg[4]
topcut <- as.numeric(argg[5]); bottomcut <- as.numeric(argg[6])

newgroups <- get_groups(npeakarea[,-c(1:32)])
#leaves <- newgroups[which(grepl('Leaf', newgroups))]
#roots <- newgroups[which(grepl('Root', newgroups))]

#if (allg == 'yes'){
#  hleaves <- leaves[which(grepl('Hydro', leaves))]
#  hroots <- roots[which(grepl('Hydro', roots))]
#  sleaves <- leaves[which(grepl('Sym', leaves))]
#  sroots <- roots[which(grepl('Sym', roots))]
#  combined <- list(hleaves, hroots, sleaves, sroots)
#} else {
#  combined <- list(leaves, roots)
#}
newgroups <- newgroups[which(str_count(newgroups, fixed(".")) == 2)]

tissues <- unlist(lapply(strsplit(newgroups, '.', fixed = T), '[[', 3))
tissues <- unique(str_to_title(tissues))
combined <- list()
for (t in tissues){
  tmp <- newgroups[which(grepl(t, newgroups))]
  conds <- unique(lapply(strsplit(tmp, '.', fixed = T), '[[', 1))
  for (cd in conds){
    tmp2 <- tmp[which(grepl(cd, tmp))]
    combined[[paste0(t, cd, sep = '.')]] <- tmp2
  }
}


da_ids <- c()
for (tissue in combined){
  if (length(tissue) > 0){
    tname <- unique(tissue[grepl('Ctrl.', tissue)])
    ctrlcol <- which(grepl(tname, colnames(npeakarea)))
    print('Control: ')
    print(colnames(npeakarea[ctrlcol][1]))
    treat <- tissue[-which(grepl('Ctrl', tissue))]
    uniq <- unique(treat)
    for (cond in uniq){
      print(cond)
      ttreat <- which(grepl(paste0('^', cond), colnames(npeakarea)))
      if (length(ttreat) > 1 & length(ctrlcol) > 1){
        both <- get_diff_exp(pkareaorg, pkareaorg[,ctrlcol], pkareaorg[,ttreat], 
                             npeakarea, npeakarea[,ctrlcol], npeakarea[,ttreat],
                             topcut, bottomcut)
        if (inherits(both, "list")){
          allpval <- both[[1]]; diff <- both[[2]] 
          da_ids <- unique(c(da_ids, diff$Alignment.ID)) 
          write.table(diff, file = paste(c(argg[1], '_FDR_diff_', cond, '.tab'), collapse = ''), sep = '\t', row.names = F, quote = F)
          write.table(allpval, file = paste(c(argg[1], '_all_pval_fc_', cond, '.tab'), collapse = ''), sep = '\t', row.names = F, quote = F)
        }else{
          allpval <- both
          write.table(allpval, file = paste(c(argg[1], '_all_pval_fc_', cond, '.tab'), collapse = ''), sep = '\t', row.names = F, quote = F)
        }
        
        if (inherits(both, 'list')){
          ##making volcano plots
          if (grepl('Leaf', tissue)) {
            color <- 'darkgreen'
            tis <- 'Leaf'
          }else if (grepl('Root', tissue)) {
            color <- 'saddlebrown'
            tis <- 'Root'
          }
          
          allpval <- allpval[which(allpval$ctrlmean >= 5000 | allpval$treatmean >= 5000), ]
          
          ##adding additional columns
          ndat <- allpval %>% mutate(logpval = -log10(adj_pval)) %>%
            mutate(threshold = if_else((nfold_change >= topcut & logpval >= 1.30103 & org_fc > 2)
                                       |(nfold_change <= bottomcut & logpval >= 1.30103 & org_fc < -2),"A", "B"))
          ndat$nfold_change[ndat$nfold_change >= 2.5] <- 2.5
          ndat$nfold_change[ndat$nfold_change <= -2.5] <- -2.5
          
          grobup <- grobTree(textGrob(paste0("n = ", length(which(ndat$nfold_change >= topcut & ndat$logpval >= 1.30103 & ndat$org_fc >= 2))), 
                                      x=0.8,  y= 0.84, hjust=0, gp=gpar(col="blue", fontsize=15)))
          grobdown <- grobTree(textGrob(paste0("n = ", length(which(ndat$nfold_change <= bottomcut & ndat$logpval >= 1.30103 & ndat$org_fc <= -2))), 
                                        x=0.1,  y= 0.84, hjust=0, gp=gpar(col="blue", fontsize=15)))
          
          pdf(paste0(odir, '_', cond, tis, 'DAM_volcano_plot.pdf'))
          plot <- ggplot(ndat, aes(nfold_change,logpval, colour = threshold, size = threshold)) +
            geom_point(alpha = 0.5) + 
            geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5, size = 1) + 
            geom_vline(xintercept = topcut, linetype = 2, alpha = 0.5, size = 1) +
            geom_vline(xintercept = bottomcut, linetype = 2, alpha = 0.5, size = 1) +
            scale_colour_manual(values = c(color, 'black'), labels = cond) +
            scale_size_manual(values=c(3,1))+
            xlab("log2 fold change") + ylab("-log10 q-value")+labs(colour = '')+
            annotation_custom(grobup)+ annotation_custom(grobdown)+
            theme_bw()+guides(size = F, colour = F) + ggtitle(paste(tis, cond, sep = ' '))
          print(plot)
          dev.off()
        }
      }else {print('Not enough replicates to find Differentially Abundant metabolites!')}
    }
  }
}

if(length(da_ids) > 0){
  write.table(da_ids, file = paste0(argg[1], '_all_da_alignmentIDs.tab'), sep = '\t', row.names = F, quote = F)
}
print('Done! Wow!')
