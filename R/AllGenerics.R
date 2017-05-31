setGeneric(
  name = "getLineages",
  signature = c('reducedDim','clusterLabels'),
  def = function(reducedDim,
                 clusterLabels, ...) {
    standardGeneric("getLineages")
  }
)

setGeneric(
  name = "getCurves",
  signature = 'sds',
  def = function(sds, ...) {
    standardGeneric("getCurves")
  }
)

setGeneric(
  name = "slingshot",
  signature = c('reducedDim','clusterLabels'),
  def = function(reducedDim,
                 clusterLabels, ...) {
    standardGeneric("slingshot")
  }
)
# accessor functions
setGeneric(name = "reducedDim",
           signature = "x",
           def = function(x) standardGeneric("reducedDim"))

setGeneric(name = "clusterLabels",
           signature = "x",
           def = function(x) standardGeneric("clusterLabels"))

setGeneric(name = "lineages",
           signature = "x",
           def = function(x) standardGeneric("lineages"))

setGeneric(name = "connectivity",
           signature = "x",
           def = function(x) standardGeneric("connectivity"))

setGeneric(name = "lineageControl",
           signature = "x",
           def = function(x) standardGeneric("lineageControl"))

setGeneric(name = "curves",
           signature = "x",
           def = function(x) standardGeneric("curves"))

setGeneric(name = "pseudotime",
           signature = "x",
           def = function(x, ...) standardGeneric("pseudotime"))

setGeneric(name = "curveWeights",
           signature = "x",
           def = function(x) standardGeneric("curveWeights"))

setGeneric(name = "curveControl",
           signature = "x",
           def = function(x) standardGeneric("curveControl"))

# replacement functions
setGeneric(name = "reducedDim<-", 
           signature = "x",
           def = function(x, value) standardGeneric("reducedDim<-"))

setGeneric(name = "clusterLabels<-", 
           signature = "x",
           def = function(x, value) standardGeneric("clusterLabels<-"))

