## ----options, results="hide", include=FALSE, cache=FALSE, message=FALSE----
knitr::opts_chunk$set(fig.align="center", cache=TRUE,error=FALSE, #stop on error
fig.width=5, fig.height=5, autodep=TRUE, out.width="600px", out.height="600px",
results="markup", echo=TRUE, eval=TRUE)
#knitr::opts_knit$set(stop_on_error = 2L) #really make it stop
#knitr::dep_auto()
options(getClass.msg=FALSE)
graphics:::par(pch = 16, las = 1)
set.seed(12345) ## for reproducibility

library(slingshot, quietly = TRUE)

## ----data_sling------------------------------------------------------------
library(slingshot, quietly = FALSE)
data("slingshotExample")

dim(rd) # data representing cells in a reduced dimensional space
length(cl) # vector of cluster labels

## ----genefilt--------------------------------------------------------------
# filter genes down to potential cell-type markers
# at least M (15) reads in at least N (15) cells
geneFilter <- apply(assays(sim)$counts,1,function(x){
    sum(x >= 15) >= 15
})
sim <- sim[geneFilter, ]

## ----norm------------------------------------------------------------------
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sim)$norm <- FQnorm(assays(sim)$counts)

## ----pca, cache=TRUE-------------------------------------------------------
pca <- prcomp(t(log1p(assays(sim)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, col = topo.colors(100)[sim$Step], pch=16, asp = 1)

## ----dm, cache=TRUE--------------------------------------------------------
library(destiny, quietly = TRUE)
dm <- DiffusionMap(t(log1p(assays(sim)$norm)))
rd2 <- cbind(DC1 = dm$DC1, DC2 = dm$DC2)

plot(rd2, col = topo.colors(100)[sim$Step], pch=16, asp = 1)

## ----add_RDs, cache=TRUE---------------------------------------------------
reducedDims(sim) <- SimpleList(PCA = rd1, DiffMap = rd2)

## ----clustering_mclust-----------------------------------------------------
library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification
colData(sim)$GMM <- cl1

library(RColorBrewer)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

## ----clustering------------------------------------------------------------
cl2 <- kmeans(rd1, centers = 4)$cluster
colData(sim)$kmeans <- cl2

plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)

## ----sling_sce-------------------------------------------------------------
sce <- slingshot(sim, clusterLabels = 'GMM', reducedDim = 'PCA')

## ----plot_curve_1----------------------------------------------------------
summary(sce$slingPseudotime_1)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plot(reducedDims(sce)$PCA, col = colors[cut(sce$slingPseudotime_1,breaks=100)], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2)

## ----plot_curve_2----------------------------------------------------------
plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages')

## ----fitgam, cache = TRUE--------------------------------------------------
require(gam)
t <- sce$slingPseudotime_1

# for time, only look at the 1,000 most variable genes
Y <- log1p(assays(sim)$norm)
var1K <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:1000]
Y <- Y[var1K,]

# fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})

## ----heatmaps--------------------------------------------------------------
require(clusterExperiment)
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]
heatdata <- assays(sim)$norm[rownames(assays(sim)$norm) %in% topgenes, 
                        order(t, na.last = NA)]
heatclus <- sce$GMM[order(t, na.last = NA)]
ce <- ClusterExperiment(heatdata, heatclus, transformation = log1p)
plotHeatmap(ce, clusterSamplesData = "orderSamplesValue",visualizeData = 'transformed')

## ----sling_lines_unsup-----------------------------------------------------
lin1 <- getLineages(rd, cl, start.clus = '1')
lin1
plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
lines(lin1, lwd = 3)

## ----lines_sup_end---------------------------------------------------------
lin2 <- getLineages(rd, cl, start.clus= '1', end.clus = '3')

plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
lines(lin2, lwd = 3, show.constraints = TRUE)

## ----curves----------------------------------------------------------------
crv1 <- getCurves(lin1)
crv1
plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
lines(crv1, lwd = 3)

## ----session---------------------------------------------------------------
sessionInfo()

