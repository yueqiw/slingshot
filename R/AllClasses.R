#' @title Class \code{SlingshotDataSet}
#' @aliases slingshotDataSet
#'
#' @description The \code{SlingshotDataSet} class holds data relevant 
#' for performing lineage inference with the \code{slingshot} package, 
#' primarily a reduced dimensional representation of the data and a set
#' of cluster labels.
#'
#' @docType class
#' @aliases SlingshotDataSet SlingshotDataSet-class SlingshotDataSet SlingshotDataSet-class
#'
#' @description All \code{slingshot} methods can take an object of the class
#' \code{SlingshotDataSet} as input and will output the same. Additionally,
#' simple helper methods for creating and manipulating objects of the class
#' \code{SlingshotDataSet} are described below.
#'
#' @slot reducedDim matrix. An \code{n} by \code{p} numeric matrix or data frame giving the
#' coordinates of the cells in a reduced dimensionality space.
#' @slot clus.labels character. A character vector of length \code{n} denoting each
#' cell's cluster label.
#' @slot lineages list. A list with each element a character vector of cluster names
#' representing a lineage as an ordered set of clusters.
#' @slot connectivity matrix. A binary matrix describing the connectivity between
#' clusters induced by the minimum spanning tree.
#' @slot lineage.control list. Additional parameters specifying how the minimum
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
#' \link{\code{getCurves}}.
#' @slot pseudotime matrix. A matrix of size \code{n} by \code{L} (where 
#' \code{L} is the number of lineages) specifying each cell's pseudotime along
#' each lineage.
#' @slot weights matrix. An \code{n} by \code{L} matrix specifying the weight
#' of each cell's contribution to each lineage.
#' @slot curve.control list. Additional parameters specifying how the
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
#' @name SlingshotDataSet-class
#' @aliases SlingshotDataSet
#' @rdname SlingshotDataSet-class
#' @import princurve
#' @import methods
#' @export
#'
setClass(
  Class = "SlingshotDataSet",
  slots = list(
    reducedDim = "matrix",
    clus.labels = "character",
    lineages = "list",
    connectivity = "matrix",
    lineage.control = "list",
    curves = "list",
    pseudotime = "matrix",
    weights = "matrix",
    curve.control = "list"
  )
)

setValidity("SlingshotDataSet", function(object) {
  X <- object@reducedDim
  n <- nrow(X)
  p <- ncol(X)
  if(!is.numeric(X)) {
    return("Reduced dimensional coordinates must be numeric.")
  }
  if(length(object@clus.labels) != n){
    return('nrow(reducedDim) must equal length(clus.labels).')
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
    clus.names <- unique(object@clus.labels)
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
    if(! is.null(object@lineage.control$start.clus)){
      if(!all(object@lineage.control$start.clus %in% clus.names)){
        return("Specified starting cluster not found in cluster labels")
      }
    }
    if(! is.null(object@lineage.control$end.clus)){
      if(!all(object@lineage.control$end.clus %in% clus.names)){
        return("Specified terminal cluster(s) not found in cluster labels")
      }
    }
    if(! is.null(object@lineage.control$dist.fun)){
      if(!is.function(object@lineage.control$dist.fun)){
        return("Pairwise cluster distance function must be a function.")
      }
    }
    if(! is.null(object@lineage.control$omega)){
      if(object@lineage.control$omega < 0 | (object@lineage.control$omega > 1 & object@lineage.control$omega != Inf)){
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
    if(dim(object@weights)[1] > 0){
      if(any(dim(object@weights) != c(n,L))){
        return("Dimensions for weights matrix are incorrect. Should be n (number of cells) by L (number of lineages).")
      }
    }
    if(!is.null(object@curve.control$shrink)){
      if(object@curve.control$shrink < 0 | object@curve.control$shrink > 1){
        stop("shrink argument must be logical or numeric between 0 and 1.")
      }
    }
    if(!is.null(object@curve.control$extend)){
      if(! object@curve.control$extend %in% c('y','n','pc1')){
        stop("extend argument must be one of 'y', 'n', or 'pc1'.")
      }
    }
    if(!is.null(object@curve.control$reweight)){
      if(!is.logical(object@curve.control$reweight)){
        stop("reweight argument must be logical.")
      }
    }
    if(!is.null(object@curve.control$drop.multi)){
      if(!is.logical(object@curve.control$drop.multi)){
        stop("drop.multi argument must be logical.")
      }
    }
  }
  return(TRUE)
  })

#' @description The \code{SlingshotDataSet} constructor creates an object of
#' the class \code{SlingshotDataSet}. Objects of this type will also be 
#' returned by \code{\link{getLineages}} and \code{\link{getCurves}}.
#'
#' @param reducedDim matrix. An \code{n} by \code{p} numeric matrix or data frame giving the
#' coordinates of the cells in a reduced dimensionality space.
#' @param clus.labels character. A character vector of length \code{n} denoting each
#' cell's cluster label.
#'
#'@return A \code{SlingshotDataSet} object.
#'
#'@examples
#'
#'reducedDim <- matrix(data=rnorm(200), ncol=2)
#'clus.labels <- sample(letters[1:5], 100, replace = TRUE)
#'
#'sds <- SlingshotDataSet(reducedDim, clus.labels)
#'
#' @rdname SlingshotDataSet-class
#' @export
setGeneric(
  name = "SlingshotDataSet",
  def = function(reducedDim,  clus.labels, ...) {
    standardGeneric("SlingshotDataSet")
  }
)
#' @rdname SlingshotDataSet-class
#' @export
setMethod(
  f = "SlingshotDataSet",
  signature = signature("data.frame","ANY"),
  definition = function(reducedDim, clus.labels, ...){
    RD <- as.matrix(reducedDim)
    rownames(RD) <- rownames(reducedDim)
    SlingshotDataSet(RD, clus.labels, ...)
  })
#' @rdname SlingshotDataSet-class
setMethod(
  f = "SlingshotDataSet",
  signature = signature("matrix", "numeric"),
  definition = function(reducedDim, clus.labels, ...){
    SlingshotDataSet(reducedDim, as.character(clus.labels), ...)
  })
#' @rdname SlingshotDataSet-class
setMethod(
  f = "SlingshotDataSet",
  signature = signature("matrix","factor"),
  definition = function(reducedDim, clus.labels, ...){
    SlingshotDataSet(reducedDim, as.character(clus.labels), ...)
  })

setMethod(
  f = "SlingshotDataSet",
  signature = signature("matrix","character"),
  definition = function(reducedDim, clus.labels,
                        lineages=list(),
                        connectivity=matrix(NA,0,0),
                        lineage.control=list(),
                        curves=list(),
                        pseudotime=matrix(NA,0,0),
                        weights=matrix(NA,0,0),
                        curve.control=list()
  ){
    if(nrow(reducedDim) != length(clus.labels)) {
      stop('nrow(reducedDim) must equal length(clus.labels).')
    }
    # something requires row and column names. Princurve?
    if(is.null(rownames(reducedDim))){
      rownames(reducedDim) <- paste('Cell',seq_len(nrow(reducedDim)),sep='-')
    }
    if(is.null(colnames(reducedDim))){
      colnames(reducedDim) <- paste('Dim',seq_len(ncol(reducedDim)),sep='-')
    }
    if(is.null(names(clus.labels))){
      names(clus.labels) <- rownames(reducedDim)
    }
    out <- new("SlingshotDataSet",
               reducedDim=reducedDim,
               clus.labels=clus.labels,
               lineages=lineages,
               connectivity=connectivity,
               lineage.control=lineage.control,
               curves=curves,
               pseudotime=pseudotime,
               weights=weights,
               curve.control=curve.control
    )
    validObject(out)
    return(out)
    })

