setGeneric(
  name = "getLineages",
  signature = c('reducedDim','clusLabels'),
  def = function(reducedDim, clusLabels,
                 start.clus = NULL, end.clus = NULL,
                 dist.fun = NULL, omega = NULL, ...) {
    standardGeneric("getLineages")
  }
)

setGeneric(
  name = "getCurves",
  signature = 'sds',
  def = function(sds, 
                 clusLabels = sds@clusLabels,
                 lineages = sds@lineages,
                 shrink = TRUE, 
                 extend = 'y', 
                 reweight = TRUE,
                 drop.multi = TRUE, 
                 thresh = 0.001, maxit = 15, stretch = 2, 
                 smoother = 'smooth.spline', 
                 shrink.method = 'cosine', ...) {
    standardGeneric("getCurves")
  }
)

setGeneric(
  name = "slingshot",
  signature = c('reducedDim','clusLabels'),
  def = function(reducedDim, clusLabels,
                 start.clus = NULL, end.clus = NULL,
                 dist.fun = NULL, omega = NULL,
                 lineages = list(),
                 shrink = TRUE,
                 extend = 'y',
                 reweight = TRUE,
                 drop.multi = TRUE,
                 thresh = 0.001, maxit = 15, stretch = 2,
                 smoother = 'smooth.spline',
                 shrink.method = 'cosine', ...) {
    standardGeneric("slingshot")
  }
)
# accessor functions
setGeneric(name = "reducedDim",
           signature = "x",
           def = function(x) standardGeneric("reducedDim"))

setGeneric(name = "clusLabels",
           signature = "x",
           def = function(x) standardGeneric("clusLabels"))

setGeneric(name = "lineages",
           signature = "x",
           def = function(x) standardGeneric("lineages"))

setGeneric(name = "connectivity",
           signature = "x",
           def = function(x) standardGeneric("connectivity"))

setGeneric(name = "lineage.control",
           signature = "x",
           def = function(x) standardGeneric("lineage.control"))

setGeneric(name = "curves",
           signature = "x",
           def = function(x) standardGeneric("curves"))

setGeneric(name = "pseudotime",
           signature = "x",
           def = function(x) standardGeneric("pseudotime"))

setGeneric(name = "curveWeights",
           signature = "x",
           def = function(x) standardGeneric("curveWeights"))

setGeneric(name = "curve.control",
           signature = "x",
           def = function(x) standardGeneric("curve.control"))

# replacement functions
setGeneric(name = "reducedDim<-", 
           signature = "x",
           def = function(x, value) standardGeneric("reducedDim<-"))

setGeneric(name = "clusLabels<-", 
           signature = "x",
           def = function(x, value) standardGeneric("clusLabels<-"))

