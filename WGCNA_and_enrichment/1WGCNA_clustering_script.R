## WGCNA because biclustering isn't working
library(WGCNA)
library(RColorBrewer)
library(dendextend)
library(pheatmap)

#setwd('E:/MS_Data/BrachyMetabolites/WGCNA/')
options(stringsAsFactors = FALSE)

argg <- commandArgs(T)
if (length(argg) != 2){
  stop('ARG: 1) Peak Area file 2) Desired path and prefix of output files')
}

datExpr <- read.table(argg[1], sep = '\t', header = T, fill = NA, quote = "")

## Removing non-fragmented metabs
datExpr <- datExpr[datExpr$MS.MS.spectrum != '',]

datExpr0 <- as.data.frame(t(datExpr[,33:ncol(datExpr)]))
ids <- datExpr$Alignment.ID

colnames(datExpr0) <- ids

## Making traits file
traits <- data.frame(matrix(nrow = nrow(datExpr0), ncol = 11))
colnames(traits) <- c('Leaf', 'Root', 'Spikelet', 'Culm', 'Hydroponics', 'Symbiosis',
                      'Soil', 'LowCu', 'Heat', 'LowP', 'AMS')
rownames(traits) <- rownames(datExpr0)
patternvec <- c('Leaf', 'Root', 'Spike', 'Culm', 'Hydro', 'Sym', 'Tis',
                'Cop', 'Heat', 'SporeW', 'Spore.')

for(col in 1:ncol(traits)){
  traits[,col] <- ifelse(grepl(patternvec[col], rownames(traits), fixed = T), 1, 0)
}

enableWGCNAThreads()

powers = c(c(1:10), seq(from = 12, to=2000, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sft$pick <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
softPower <- min(which(sft$pick > 0.8))  

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))

text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

adjacency = adjacency(datExpr0, power = softPower)
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM


geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

minModuleSize = 10
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")                    


## Ok, quite a few genes are grey -- don't care about these!
selected_genes <- which(dynamicColors != 'grey')
selectTOM = dissTOM[selected_genes, selected_genes]
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = dynamicColors[selected_genes]

newdynamicMods <- cutreeDynamic(dendro = selectTree, distM = selectTOM,
                                deepSplit = 2, pamRespectsDendro = FALSE,
                                minClusterSize = minModuleSize)

newdynamicColors <- labels2colors(newdynamicMods)
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
png(paste0(argg[2], '_wgcna_network_heatmap.png'))
TOMplot(1 - plotDiss, selectTree, newdynamicColors, main = "Network heatmap plot, selected genes",
        col = colorRampPalette( rev(brewer.pal(9, 'RdYlBu')))(255))
dev.off()
##writing out what metabs belong to what module
out <- cbind(colnames(datExpr0)[selected_genes], selectColors)
write.table(out, file = paste0(argg[2], 'wgcna_module_assignment.tsv'), sep = '\t', 
            row.names = F, quote = F)
collectGarbage()
## K not very satisfied with the above. Going back to plain old heirarchical 
## clustering. Only using the metabs not in the 'grey' modules
to_clust <- datExpr0[, selected_genes]
corr_mat <- cor(to_clust, method = 'spearman')

dist_mat <- as.dist(1 - corr_mat)
my_tree <- hclust(dist_mat, method = 'average')
pdf(paste0(argg[2], '_spearman_clustering_heatmap.pdf'))
pheatmap(corr_mat, clustering_method = 'average', show_colnames = F, 
         show_rownames = F)
dev.off()

## Module significance with traits
nSamples = nrow(datExpr0)
MEs0 = moduleEigengenes(datExpr0, dynamicColors, excludeGrey = T)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(datExpr0))

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 100, 3, 3));
# Display the correlation values within a heatmap plot
pdf(paste0(argg[2], '_module_importance_to_traits.pdf'), width = 7)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(100),
               textMatrix = textMatrix,
               setStdMargins = T,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
