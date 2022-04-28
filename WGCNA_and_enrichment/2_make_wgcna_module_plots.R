## Get module file
## Normalized expression file
## Canopus prediction file

library(ggplot2)
library(reshape2)
library(dplyr)
library(data.table)

argg <- commandArgs(T)
if (length(argg) != 6){
  stop('ARGS: 1) WGNCA module output file 2) Normalized Peak Area file
       3) CANOPUS class file 4) Any clusters you want to plot in different color?
       Separate with , or "none" 5) If 4 is not "none": How do you want to define 
       what peaks get what color? Provide condition names to look for peaks that
       are expressed vs. not expressed in. Syntax: [S.S.L,T.L] (use abbreviated names)
       6) output file directory and prefix')
}

get_groups <- function(df, typ){
  ct <- colnames(df)
  groups <- unlist(lapply(strsplit(ct, '_'), '[', 1))
  return(groups)
}

get_mean <- function(df, grps){
  outl <- list()
  for (grp in unique(grps)){
    tmp <- df[,which(grepl(grp, colnames(df)))]
    meantmp <- rowMeans(tmp)
    outl[[grp]] <- meantmp
  }
  return(do.call(cbind.data.frame, outl))
}

plot_line_colors <- function(mpeakmean, modvec, condvec){
  ## Breakup condvec (which is a string at input)
  ## add "color" col to mpeakmean, fill it with 0
  ## loop through mods in modvec
  ## get the mod from mpeakmean, remove it from mpeakmean
  ## find alignment IDs that are > 0 in each cond in condvec (a)
  ## find IDs that are 0 in each cond in condvec (b)
  ## for every instance of every ID in (a), put 1 in color col
  ## for every instance of every ID in (b), put 2 in color col
  ## reattach this back to mpeakmean
  
  ## splitting up condvec into list with what conds go with what module
  modgroups <- unlist(strsplit(condvec, '[', fixed = T))
  modgroups <- modgroups[modgroups != '']
  
  cond_list <- list()
  for (grp_n in 1:length(modgroups)){
    grp <- modgroups[grp_n]
    mod <- modvec[grp_n]
    grp <- gsub(']', '', grp, fixed = T)
    split_conds <- unlist(strsplit(grp, ','))
    cond_list[[mod]] <- split_conds
  }
  
  ## initialize color col
  mpeakmean$color <- 0
  
  for (mod in modvec){
    conds <- cond_list[[mod]]
    mpeakmean_mod <- mpeakmean[which(mpeakmean$cluster == mod),]
    mpeakmean <- mpeakmean[-which(mpeakmean$cluster == mod),]
    
    expr_ids <- unique(mpeakmean_mod[which(mpeakmean_mod$variable %in% conds &
                                      mpeakmean_mod$value > 0.010), 'Alignment.ID'])
    not_expr_ids <- unique(mpeakmean_mod[which(mpeakmean_mod$variable %in% conds &
                                                 mpeakmean_mod$value == 0.010), 'Alignment.ID'])
    mpeakmean_mod$color[mpeakmean_mod$Alignment.ID %in% expr_ids] <- 1
    mpeakmean_mod$color[mpeakmean_mod$Alignment.ID %in% not_expr_ids] <- 2
    mpeakmean <- rbind.data.frame(mpeakmean, mpeakmean_mod)
  }
  
  return(mpeakmean)
}


modules <- read.table(argg[1], sep = '\t', header = T, fill = NA, quote ="")
peakarea <- read.table(argg[2], sep = '\t', header = T, fill = NA, quote = "")
canopus <- read.table(argg[3], sep = '\t', header = T, fill = NA, quote = "")
modvec <- unlist(strsplit(argg[4], ','))
condstr <- argg[5]

## Finding CANOPUS class makeup of each module first
canopus_tomerge <- canopus[, c('name', 'class', 'superclass')]
module_canopus <- merge(modules, canopus_tomerge, all.x = T, by.x = 'X', by.y = 'name')

## Filling in NAs with 'None'
module_canopus[is.na(module_canopus$class), 'class'] <- 'None'
module_canopus[is.na(module_canopus$superclass), 'superclass'] <- 'None'

## Combining classes/superclasses with low numbers in a module to 'Other'
for (mod in unique(module_canopus$selectColors)) {
  print(mod)
  tmp <- module_canopus[module_canopus$selectColors == mod,]
  nums_classes <- as.data.frame(table(tmp$class))
  for (rw in 1:nrow(nums_classes)){
    if (nums_classes[rw, 'Freq'] < 2){
      small_class <- nums_classes[rw, 'Var1']
      module_canopus[module_canopus$selectColors == mod & module_canopus$class == small_class, 'class'] <- 'Other'
    }
  }
}

## Now making class names shorter
module_canopus$class <- gsub(' and derivatives', '', module_canopus$class)
module_canopus$class <- gsub(' and substituted derivatives', '', module_canopus$class)
module_canopus$class <- gsub(' and steroid derivatives', '', module_canopus$class)

## Removing class "None" for class barcharts
module_canopus_bars <- module_canopus[-which(module_canopus$class == 'None'),]

## Changing boxplot order so largest are first
module_tabs <- module_canopus_bars %>% group_by(selectColors, class) %>% summarize(class_counts = n())
molten <- setDT(module_tabs)
molten[, ord := sprintf("%02i", frank(molten, selectColors, -class_counts, ties.method = "first"))]

pdf(paste0(argg[6], '_wgcna_module_class_makeup.pdf'), width = 15, height = 15)
ggplot(molten, aes(ord, class_counts))+
  geom_col(width=0.7, fill="steelblue") + 
  facet_wrap(~selectColors, scales = 'free') + 
  scale_x_discrete(labels = molten[, setNames(as.character(class), ord)]) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

## Now for expression
peakinterest <- peakarea[which(peakarea$Alignment.ID %in% module_canopus$X),]
rownames(peakinterest) <- peakinterest$Alignment.ID
peakinterest <- peakinterest[, 33:ncol(peakinterest)]
groups <- get_groups(peakinterest)
peakmean <- get_mean(peakinterest, groups)

peakmean$Alignment.ID <- rownames(peakmean)
peakmean$cluster <- module_canopus$selectColors
mpeakmean <- melt(peakmean, id.vars = c('cluster', 'Alignment.ID'))

## Reordering and renaming factor levels
mpeakmean$variable <- factor(mpeakmean$variable, levels = c('Hydro.Ctrl.Root', 'Hydro.Cop.Root', 'Hydro.Heat.Root',
                                                            'Hydro.HeatCop.Root', 'Sym.Ctrl.Root', 'Sym.Spore.Root',
                                                            'Sym.SporeW.Root', 'Hydro.Ctrl.Leaf', 'Hydro.Cop.Leaf',
                                                            'Hydro.Heat.Leaf', 'Hydro.HeatCop.Leaf', 'Sym.Ctrl.Leaf',
                                                            'Sym.Spore.Leaf', 'Sym.SporeW.Leaf', 'Tis.Leaf', 'Tis.Culm',
                                                            'Tis.Spike'))

mpeakmean$variable <- recode_factor(mpeakmean$variable, Hydro.Ctrl.Root = 'H.C.R',
                                    Hydro.Cop.Root = 'H.NC.R', Hydro.Heat.Root = 'H.H.R',
                                    Hydro.HeatCop.Root = 'H.HNC.R', Sym.Ctrl.Root = 'S.C.R', 
                                    Sym.Spore.Root = 'S.S.R', Sym.SporeW.Root = 'S.SW.R',
                                    Hydro.Ctrl.Leaf = 'H.C.L', Hydro.Cop.Leaf = 'H.NC.L',
                                    Hydro.Heat.Leaf = 'H.H.L', Hydro.HeatCop.Leaf = 'H.HNC.L', 
                                    Sym.Ctrl.Leaf = 'S.C.L', Sym.Spore.Leaf = 'S.S.L', 
                                    Sym.SporeW.Leaf = 'S.SW.L', Tis.Leaf = 'T.L', 
                                    Tis.Culm = 'T.C', Tis.Spike = 'T.S')

mpeakmean <- plot_line_colors(mpeakmean, modvec, condstr)
cols <- c("0" = "black", "1" = "aquamarine3", "2" = "blue4")

## renaming names so that they have number of peaks in module


pdf(paste0(argg[6], '_wgcna_module_expression.pdf'), width = 15, height = 10)
ggplot(mpeakmean, aes(x = variable, y = value, group = Alignment.ID, color = as.factor(color))) + geom_line() + 
  geom_point() + facet_wrap(~cluster, scales = 'free') + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_colour_manual(values = cols) + guides(color = 'none')

dev.off()
