setGeneric(
  name = "getLineages",
  signature = c('reducedDim','clus.labels'),
  def = function(reducedDim, clus.labels,
                 start.clus = NULL, end.clus = NULL,
                 dist.fun = NULL, omega = NULL, ...) {
    standardGeneric("getLineages")
  }
)

setGeneric(
  name = "getCurves",
  signature = 'sds',
  def = function(sds, 
                 clus.labels = sds@clus.labels,
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
  signature = c('reducedDim','clus.labels'),
  def = function(reducedDim, clus.labels,
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

setGeneric(name = "clus.labels",
           signature = "x",
           def = function(x) standardGeneric("clus.labels"))

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

setGeneric(name = "weights",
           signature = "x",
           def = function(x) standardGeneric("weights"))

setGeneric(name = "curve.control",
           signature = "x",
           def = function(x) standardGeneric("curve.control"))

# replacement functions
setGeneric(name = "reducedDim<-", 
           signature = "x",
           def = function(x, value) standardGeneric("reducedDim<-"))

setGeneric(name = "clus.labels<-", 
           signature = "x",
           def = function(x, value) standardGeneric("clus.labels<-"))

