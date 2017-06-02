#' @title Initialize an object of class \code{SlingshotDataSet}
#' @name newSlingshotDataSet
#' @docType methods
#'   
#' @description Constructs a \code{SlingshotDataSet} object. Additional helper 
#'   methods for manipulating \code{SlingshotDataSet} objects are  also 
#'   described below.
#'   
#' @param reducedDim matrix. An \code{n} by \code{p} numeric matrix or data
#'   frame giving the coordinates of the cells in a reduced dimensionality
#'   space.
#' @param clusterLabels character. A character vector of length \code{n}
#'   denoting each cell's cluster label.
#' @param lineages list. A list with each element a character vector of cluster 
#'   names representing a lineage as an ordered set of clusters.
#' @param connectivity matrix. A binary matrix describing the connectivity 
#'   between clusters induced by the minimum spanning tree.
#' @param lineageControl list. Additional parameters specifying how the minimum 
#'   spanning tree on clusters was constructed. \itemize{ 
#'   \item{\code{start.clus}}{character. The label of the root cluster.} 
#'   \item{\code{end.clus}}{character. Vector of cluster labels indicating the 
#'   terminal clusters.} \item{\code{start.given}}{logical. A logical value 
#'   indicating whether the initial state was pre-specified.} 
#'   \item{\code{end.given}}{logical. A vector of logical values indicating 
#'   whether each terminal state was pre-specified} \item{\code{dist}}{matrix. A
#'   numeric matrix of pairwise cluster distances.} }
#' @param curves list. A list of \code{principal.curve} objects produced by 
#'   \code{\link{getCurves}}.
#' @param curveControl list. Additional parameters specifying how the 
#'   simultaneous principal curves were constructed. \itemize{ 
#'   \item{\code{shrink}}{logical or numeric between 0 and 1. Determines whether
#'   and how much to shrink branching lineages toward their shared average 
#'   curve.} \item{\code{extend}}{character. Specifies the method for handling 
#'   root and leaf clusters of lineages when constructing the initial, 
#'   piece-wise linear curve. Accepted values are 'y' (default), 'n', and 'pc1'.
#'   See \code{\link{getCurves}} for details.} \item{\code{reweight}}{logical. 
#'   Indicates whether to reweight cells shared by multiple lineages during 
#'   curve-fitting. If \code{TRUE}, cells shared between lineages will have 
#'   lineage-specific weights determined by the ratio: (distance to nearest 
#'   curve) / (distance to specific curve).} \item{\code{drop.multi}}{logical. 
#'   Indicates whether to drop shared cells from lineages which do not fit them 
#'   well. If \code{TRUE}, shared cells with a distance to one lineage above the
#'   90th percentile and another lineage below the 50th percentile will be 
#'   dropped from the farther lineage.} \item{\code{shrink.method}}{character. 
#'   Denotes how to determine the amount of shrinkage for a branching lineage. 
#'   Accepted values are the same as for \code{kernel} in  the \code{density} 
#'   function (default is \code{"cosine"}), as well as \code{"tricube"} and 
#'   \code{"density"}. See \code{\link{getCurves}} for details.} \item{Other 
#'   parameters specified by \code{\link{principal.curve}}}. }
#'   
#' @return A \code{SlingshotDataSet} object with all specified values.
#'   
#' @examples
#' rd <- matrix(data=rnorm(100), ncol=2)
#' cl <- sample(letters[1:5], 50, replace = TRUE)
#' sds <- newSlingshotDataSet(rd, cl)
#' 
#' @import princurve
#' @import methods
#' @export
setGeneric(
  name = "newSlingshotDataSet",
  signature = c('reducedDim','clusterLabels'),
  def = function(reducedDim,  clusterLabels, ...) {
    standardGeneric("newSlingshotDataSet")
  }
)

#' @title Infer Lineage Structure from Clustered Samples
#' @name getLineages
#' @export
setGeneric(
  name = "getLineages",
  signature = c('reducedDim','clusterLabels'),
  def = function(reducedDim,
                 clusterLabels, ...) {
    standardGeneric("getLineages")
  }
)

#' @title Construct Smooth Lineage Curves
#' @name getCurves
#' @export
setGeneric(
  name = "getCurves",
  signature = 'sds',
  def = function(sds, ...) {
    standardGeneric("getCurves")
  }
)

#' @title Perform lineage inference with Slingshot
#' @name slingshot
#' @export
setGeneric(
  name = "slingshot",
  signature = c('reducedDim','clusterLabels'),
  def = function(reducedDim,
                 clusterLabels, ...) {
    standardGeneric("slingshot")
  }
)

# accessor functions
#' @title Returns the reduced dimensional representation of a dataset.
#'
#' @param x an object that describes a dataset or a model involving reduced
#'   dimensional data.
#' @return the matrix representing the reduced dimensional data.
#' @examples 
#' rd <- matrix(data=rnorm(100), ncol=2)
#' cl <- sample(letters[1:5], 50, replace = TRUE)
#' sds <- newSlingshotDataSet(rd, cl)
#' reducedDim(sds)
#' @export
setGeneric(name = "reducedDim",
           signature = "x",
           def = function(x) standardGeneric("reducedDim"))

#' @title Returns the cluster labels
#' @name clusterLabels
#'   
#' @param x an object that describes a dataset or a model involving cluster
#'   labels.
#' @return the vector of cluster labels.
#' @examples
#' rd <- matrix(data=rnorm(100), ncol=2)
#' cl <- sample(letters[1:5], 50, replace = TRUE)
#' sds <- newSlingshotDataSet(rd, cl)
#' clusterLabels(sds)
#' @export
setGeneric(name = "clusterLabels",
           signature = "x",
           def = function(x) standardGeneric("clusterLabels"))

#' @title Returns the lineages
#'   
#' @param x an object that describes a dataset or a model involving lineages
#' @return the list of lineages, represented by ordered sets of clusters.
#' @examples
#' data("slingshotExample")
#' sds <- getLineages(rd, cl)
#' lineages(sds)
#' @export
setGeneric(name = "lineages",
           signature = "x",
           def = function(x) standardGeneric("lineages"))

#' @title Returns the connectivity matrix
#'   
#' @param x an object that describes a dataset or a model involving a
#'   connectivity matrix.
#' @return the matrix of connections between clusters.
#' @examples
#' data("slingshotExample")
#' sds <- getLineages(rd, cl)
#' connectivity(sds)
#' @export
setGeneric(name = "connectivity",
           signature = "x",
           def = function(x) standardGeneric("connectivity"))

#' @title Returns the lineage control parameters
#'   
#' @param x an object that describes a dataset or a model involving lineages.
#' @return the list of additional lineage inference parameters.
#' @examples
#' data("slingshotExample")
#' sds <- getLineages(rd, cl, start.clus = '5')
#' lineageControl(sds)
#' @export
setGeneric(name = "lineageControl",
           signature = "x",
           def = function(x) standardGeneric("lineageControl"))

#' @title Returns the principal curves
#'   
#' @param x an object that describes a dataset or a model involving a set of
#'   principal curves.
#' @return the list of smooth lineage curves.
#' @examples
#' data("slingshotExample")
#' sds <- slingshot(rd, cl)
#' curves(sds)
#' @export
setGeneric(name = "curves",
           signature = "x",
           def = function(x) standardGeneric("curves"))

#' @title Returns the curve control parameters
#'   
#' @param x an object that describes a dataset or a model involving a set of
#'   principal curves.
#' @return the list of additional curve fitting parameters.
#' @examples
#' data("slingshotExample")
#' sds <- slingshot(rd, cl)
#' curveControl(sds)
#' @export
setGeneric(name = "curveControl",
           signature = "x",
           def = function(x) standardGeneric("curveControl"))

# replacement functions
#' @rdname reducedDim 
#' @return Updated object with new reduced dimensional matrix.
#' @export
setGeneric(name = "reducedDim<-", 
           signature = "x",
           def = function(x, value) standardGeneric("reducedDim<-"))

#' @rdname clusterLabels 
#' @return Updated object with new vector of cluster labels.
#' @export
setGeneric(name = "clusterLabels<-", 
           signature = "x",
           def = function(x, value) standardGeneric("clusterLabels<-"))

#' @title Get Slingshot pseudotime values
#' @name pseudotime
#' 
#' @description Extract the matrix of pseudotime values or cells' weights along
#'   each lineage.
#' 
#' @param x a \code{SlingshotDataSet} object.
#' @param na logical. If \code{TRUE} (default), cells that are not assigned to a
#'   lineage will have a pseudotime value of \code{NA}. Otherwise, their
#'   arclength along the curve will be returned.
#' @return an \code{n} by \code{L} matrix representing each cell's pseudotime
#'   along each lineage.
#' @examples
#' data("slingshotExample")
#' sds <- slingshot(rd, cl)
#' pseudotime(sds)
#' @export
setGeneric(name = "pseudotime",
           signature = "x",
           def = function(x, ...) standardGeneric("pseudotime"))

#' @rdname pseudotime 
#' @return an \code{n} by \code{L} matrix of cell weights along
#'   each lineage.
#' @examples
#' data("slingshotExample")
#' sds <- slingshot(rd, cl)
#' curveWeights(sds)
#' @export
setGeneric(name = "curveWeights",
           signature = "x",
           def = function(x) standardGeneric("curveWeights"))
