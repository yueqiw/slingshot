#' @title Class \code{SlingshotDataSet}
#' @aliases SlingshotDataSet-class slingshotDataSet-class
#' @name SlingshotDataSet-class
#' @docType class
#'   
#' @description The \code{SlingshotDataSet} class holds data relevant for
#'   performing lineage inference with the \code{slingshot} package, primarily a
#'   reduced dimensional representation of the data and a set of cluster labels.
#'   
#' @description All \code{slingshot} methods can take an object of the class 
#'   \code{SlingshotDataSet} as input and will output the same. Additionally, 
#'   simple helper methods for creating and manipulating objects of the class 
#'   \code{SlingshotDataSet} are described below.
#'   
#' @slot reducedDim matrix. An \code{n} by \code{p} numeric matrix or data frame giving the
#' coordinates of the cells in a reduced dimensionality space.
#' @slot clusterLabels character. A character vector of length \code{n} denoting each
#' cell's cluster label.
#' @slot lineages list. A list with each element a character vector of cluster names
#' representing a lineage as an ordered set of clusters.
#' @slot connectivity matrix. A binary matrix describing the connectivity between
#' clusters induced by the minimum spanning tree.
#' @slot lineageControl list. Additional parameters specifying how the minimum
#' spanning tree on clusters was constructed.
#' \itemize{
#' \item{\code{start.clus}}{character. The label of the root cluster.}
#' \item{\code{end.clus}}{character. Vector of cluster labels indicating the 
#' terminal clusters.}
#' \item{\code{start.given}}{logical. A logical value indicating whether the initial
#' state was pre-specified.}
#' \item{\code{end.given}}{logical. A vector of logical values indicating whether
#' each terminal state was pre-specified}
#' \item{\code{dist}}{matrix. A numeric matrix of pairwise cluster distances.}
#' }
#' @slot curves list. A list of \code{principal.curve} objects produced by
#' \code{\link{getCurves}}.
#' @slot curveControl list. Additional parameters specifying how the
#' simultaneous principal curves were constructed.
#' \itemize{
#' \item{\code{shrink}}{logical or numeric between 0 and 1. Determines whether
#' and how much to shrink branching lineages toward their shared average curve.}
#' \item{\code{extend}}{character. Specifies the method for handling root and
#' leaf clusters of lineages when constructing the initial, piece-wise linear
#' curve. Accepted values are 'y' (default), 'n', and 'pc1'. See
#' \code{\link{getCurves}} for details.}
#' \item{\code{reweight}}{logical. Indicates whether to reweight cells shared
#' by multiple lineages during curve-fitting. If \code{TRUE}, cells shared
#' between lineages will have lineage-specific weights determined by the ratio:
#' (distance to nearest curve) / (distance to specific curve).}
#' \item{\code{drop.multi}}{logical. Indicates whether to drop shared cells
#' from lineages which do not fit them well. If \code{TRUE}, shared cells 
#' with a distance to one lineage above the 90th percentile and another 
#' lineage below the 50th percentile will be dropped from the farther lineage.}
#' \item{\code{shrink.method}}{character. Denotes how to determine the
#' amount of shrinkage for a branching lineage. Accepted values are the same as 
#' for \code{kernel} in  the \code{density} function (default is 
#' \code{"cosine"}), as well as \code{"tricube"} and \code{"density"}. See 
#' \code{\link{getCurves}} for details.}
#' \item{Other parameters specified by \code{\link{principal.curve}}}.
#' }
#'
#' @import princurve
#' @import methods
#' @export
#'
setClass(
  Class = "SlingshotDataSet",
  slots = list(
    reducedDim = "matrix",
    clusterLabels = "character",
    lineages = "list",
    connectivity = "matrix",
    lineageControl = "list",
    curves = "list",
    pseudotime = "matrix", # make a function that pulls this from curves
    curveWeights = "matrix", #
    curveControl = "list"
  )
)

setValidity("SlingshotDataSet", function(object) {
  X <- object@reducedDim
  n <- nrow(X)
  p <- ncol(X)
  if(!is.numeric(X)) {
    return("Reduced dimensional coordinates must be numeric.")
  }
  if(nrow(X)==0){
    return('reducedDim has zero rows.')
  }
  if(ncol(X)==0){
    return('reducedDim has zero columns.')
  }
  if(length(object@clusterLabels) != n){
    return('nrow(reducedDim) must equal length(clusterLabels).')
  }
  # something requires row and column names. Princurve?
  if(is.null(rownames(object@reducedDim))){
    rownames(object@reducedDim) <- paste('Cell',seq_len(nrow(object@reducedDim)),sep='-')
  }
  if(is.null(colnames(object@reducedDim))){
    colnames(object@reducedDim) <- paste('Dim',seq_len(ncol(object@reducedDim)),sep='-')
  }
  
  # if lineages present
  if(length(object@lineages) > 0){
    L <- length(object@lineages)
    clus.names <- unique(object@clusterLabels)
    K <- length(clus.names)
    if(any(sapply(object@lineages,class) != 'character')){
      return("lineages must be a list of character vectors.")
    }
    if(!all(sapply(object@lineages, function(lin){ all(lin %in% clus.names) }))){
      return("lineages must be a list of character vectors composed of cluster names.")
    }
    if(!is.numeric(object@connectivity)) {
      return("Connectivity matrix must be numeric or logical.")
    }
    if(any(dim(object@connectivity) != K)){
      return("Connectivity matrix must be square with number of dimensions equal to number of clusters")
    }
    if(! is.null(object@lineageControl$start.clus)){
      if(!all(object@lineageControl$start.clus %in% clus.names)){
        return("Specified starting cluster not found in cluster labels")
      }
    }
    if(! is.null(object@lineageControl$end.clus)){
      if(!all(object@lineageControl$end.clus %in% clus.names)){
        return("Specified terminal cluster(s) not found in cluster labels")
      }
    }
    if(! is.null(object@lineageControl$dist.fun)){
      if(!is.function(object@lineageControl$dist.fun)){
        return("Pairwise cluster distance function must be a function.")
      }
    }
    if(! is.null(object@lineageControl$omega)){
      if(object@lineageControl$omega < 0 | (object@lineageControl$omega > 1 & object@lineageControl$omega != Inf)){
        return("Omega must be numeric element of [0,1] or Inf.")
      }
    }
  }
  
  # if curves present
  if(length(object@curves) > 0){
    if(length(object@lineages) > 0){
      L <- length(object@lineages)
      if(length(object@curves) != L){
        return("Number of curves does not match number of lineages")
      }
    }
    L <- length(object@curves)
    if(any(sapply(object@curves,class) != 'principal.curve')){
      return("curves must be a list of principal.curve objects.")
    }
    if(dim(object@pseudotime)[1] > 0){
      if(any(dim(object@pseudotime) != c(n,L))){
        return("Dimensions for pseudotime matrix are incorrect. Should be n (number of cells) by L (number of lineages).")
      }
    }
    if(dim(object@curveWeights)[1] > 0){
      if(any(dim(object@curveWeights) != c(n,L))){
        return("Dimensions for curveWeights matrix are incorrect. Should be n (number of cells) by L (number of lineages).")
      }
    }
    if(!is.null(object@curveControl$shrink)){
      if(object@curveControl$shrink < 0 | object@curveControl$shrink > 1){
        stop("shrink argument must be logical or numeric between 0 and 1.")
      }
    }
    if(!is.null(object@curveControl$extend)){
      if(! object@curveControl$extend %in% c('y','n','pc1')){
        stop("extend argument must be one of 'y', 'n', or 'pc1'.")
      }
    }
    if(!is.null(object@curveControl$reweight)){
      if(!is.logical(object@curveControl$reweight)){
        stop("reweight argument must be logical.")
      }
    }
    if(!is.null(object@curveControl$drop.multi)){
      if(!is.logical(object@curveControl$drop.multi)){
        stop("drop.multi argument must be logical.")
      }
    }
  }
  return(TRUE)
  })

#'@title Constructor function for the \code{SlingshotDataSet} class
#'@rdname SlingshotDataSet
#'  
#'@description The \code{SlingshotDataSet} constructor creates an object of the
#'  class \code{SlingshotDataSet}. Objects of this type will also be returned by
#'  \code{\link{getLineages}} and \code{\link{getCurves}}.
#'  
#'@description Additional functions extract or replace elements of a
#'  \code{SlingshotDataSet}.
#'  
#'@param reducedDim matrix. An \code{n} by \code{p} numeric matrix or data frame
#'  giving the coordinates of the cells in a reduced dimensionality space.
#'@param clusterLabels character. A character vector of length \code{n} denoting
#'  each cell's cluster label.
#'@param ... Additional elements to specify, see
#'  \code{\link{SlingshotDataSet-class}}.
#'  
#'@return A \code{SlingshotDataSet} object.
#'  
#'@seealso \code{\link{SlingshotDataSet-class}}
#'  
#'@examples
#'
#'rd <- matrix(data=rnorm(200), ncol=2)
#'cl <- sample(letters[1:5], 100, replace = TRUE)
#'sds <- SlingshotDataSet(rd, cl)
#'
#' @export
setGeneric(
  name = "SlingshotDataSet",
  def = function(reducedDim,  clusterLabels, ...) {
    standardGeneric("SlingshotDataSet")
  }
)
#' @rdname SlingshotDataSet
#' @export
setMethod(
  f = "SlingshotDataSet",
  signature = signature("data.frame","ANY"),
  definition = function(reducedDim, clusterLabels, ...){
    RD <- as.matrix(reducedDim)
    rownames(RD) <- rownames(reducedDim)
    SlingshotDataSet(RD, clusterLabels, ...)
  })
#' @rdname SlingshotDataSet
#' @export
setMethod(
  f = "SlingshotDataSet",
  signature = signature("matrix", "numeric"),
  definition = function(reducedDim, clusterLabels, ...){
    SlingshotDataSet(reducedDim, as.character(clusterLabels), ...)
  })
#' @rdname SlingshotDataSet
#' @export
setMethod(
  f = "SlingshotDataSet",
  signature = signature("matrix","factor"),
  definition = function(reducedDim, clusterLabels, ...){
    SlingshotDataSet(reducedDim, as.character(clusterLabels), ...)
  })
#' @rdname SlingshotDataSet
#' @export
setMethod(
  f = "SlingshotDataSet",
  signature = signature("matrix","character"),
  definition = function(reducedDim, clusterLabels,
                        lineages=list(),
                        connectivity=matrix(NA,0,0),
                        lineageControl=list(),
                        curves=list(),
                        pseudotime=matrix(NA,0,0),
                        curveWeights=matrix(NA,0,0),
                        curveControl=list()
  ){
    if(nrow(reducedDim) != length(clusterLabels)) {
      stop('nrow(reducedDim) must equal length(clusterLabels).')
    }
    # something requires row and column names. Princurve?
    if(is.null(rownames(reducedDim))){
      rownames(reducedDim) <- paste('Cell',seq_len(nrow(reducedDim)),sep='-')
    }
    if(is.null(colnames(reducedDim))){
      colnames(reducedDim) <- paste('Dim',seq_len(ncol(reducedDim)),sep='-')
    }
    if(is.null(names(clusterLabels))){
      names(clusterLabels) <- rownames(reducedDim)
    }
    out <- new("SlingshotDataSet",
               reducedDim=reducedDim,
               clusterLabels=clusterLabels,
               lineages=lineages,
               connectivity=connectivity,
               lineageControl=lineageControl,
               curves=curves,
               pseudotime=pseudotime,
               curveWeights=curveWeights,
               curveControl=curveControl
    )
    validObject(out)
    return(out)
    })

