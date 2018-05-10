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
#' @param ... additional components of a \code{SlingshotDataSet} to specify.
#'   This may include any of the following:
#' @param lineages list. A list with each element a character vector of cluster 
#'   names representing a lineage as an ordered set of clusters.
#' @param adjacency matrix. A binary matrix describing the connectivity 
#'   between clusters induced by the minimum spanning tree.
#' @param slingParams list. Additional parameters used by Slingshot. These may 
#'   specify how the minimum spanning tree on clusters was constructed: 
#'   \itemize{ 
#'   \item{\code{start.clus}}{character. The label of the root cluster.} 
#'   \item{\code{end.clus}}{character. Vector of cluster labels indicating the 
#'   terminal clusters.}
#'   \item{\code{start.given}}{logical. A logical value 
#'   indicating whether the initial state was pre-specified.} 
#'   \item{\code{end.given}}{logical. A vector of logical values indicating 
#'   whether each terminal state was pre-specified}
#'   \item{\code{dist}}{matrix. A
#'   numeric matrix of pairwise cluster distances.} }
#'   They may also specify how simultaneous principal curves were constructed:
#'   \itemize{ 
#'   \item{\code{shrink}}{logical or numeric between 0 and 1. Determines whether
#'   and how much to shrink branching lineages toward their shared average 
#'   curve.} 
#'   \item{\code{extend}}{character. Specifies the method for handling 
#'   root and leaf clusters of lineages when constructing the initial, 
#'   piece-wise linear curve. Accepted values are 'y' (default), 'n', and 'pc1'.
#'   See \code{\link{getCurves}} for details.} 
#'   \item{\code{reweight}}{logical. 
#'   Indicates whether to reweight cells shared by multiple lineages during 
#'   curve-fitting. If \code{TRUE}, cells shared between lineages will have 
#'   lineage-specific weights determined by the ratio: (distance to nearest 
#'   curve) / (distance to specific curve).} 
#'   \item{\code{drop.multi}}{logical. 
#'   Indicates whether to drop shared cells from lineages which do not fit them 
#'   well. If \code{TRUE}, shared cells with a distance to one lineage above the
#'   90th percentile and another lineage below the 50th percentile will be 
#'   dropped from the farther lineage.} 
#'   \item{\code{shrink.method}}{character. 
#'   Denotes how to determine the amount of shrinkage for a branching lineage. 
#'   Accepted values are the same as for \code{kernel} in  the \code{density} 
#'   function (default is \code{"cosine"}), as well as \code{"tricube"} and 
#'   \code{"density"}. See \code{\link{getCurves}} for details.} 
#'   \item{Other parameters specified by \code{\link{principal.curve}}}. }
#' @param curves list. A list of \code{principal.curve} objects produced by 
#'   \code{\link{getCurves}}.
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

#' @title Extract Slingshot output
#' @name SlingshotDataSet
#' @description This is a convenience function to extract a
#'   \code{SlingshotDataSet} from an object containing \code{\link{slingshot}}
#'   output.
#' @param data an object containing \code{slingshot} output.
#' @param ... additional arguments to pass to object-specific methods.
#' @return A \code{SlingshotDataSet} object containing the output of 
#' \code{slingshot}.
#' @export
setGeneric(
    name = "SlingshotDataSet",
    signature = c('data'),
    def = function(data, ...) {
        standardGeneric("SlingshotDataSet")
    }
)

#' @title Infer Lineage Structure from Clustered Samples
#' @name getLineages
#' @param ... Additional arguments to specify how lineages are constructed from
#'   clusters.
#' @export
setGeneric(
    name = "getLineages",
    signature = c('data','clusterLabels'),
    def = function(data,
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
    signature = c('data', 'clusterLabels'),
    def = function(data,
                   clusterLabels, ...) {
        standardGeneric("slingshot")
    }
)


#' @title Extract the Slingshot lineages
#' @name slingLineages
#' 
#' @description Extract lineages (represented by ordered sets of clusters) from
#'   a \code{SlingshotDataSet}.
#'   
#' @param x an object containing \code{\link{slingshot}} output.
#' @return the list of lineages, represented by ordered sets of clusters.
#' @examples
#' data("slingshotExample")
#' sds <- getLineages(rd, cl)
#' slingLineages(sds)
#' @export
setGeneric(name = "slingLineages",
           signature = "x",
           def = function(x) standardGeneric("slingLineages"))

#' @title Extract Slingshot adjacency matrix
#' @name slingAdjacency
#' @description Extract the adjacency matrix from an object containing
#'   \code{\link{slingshot}} output.
#' 
#' @param x an object containing \code{\link{slingshot}} output.
#' @return the matrix of connections between clusters, inferred by the MST.
#' @examples
#' data("slingshotExample")
#' sds <- getLineages(rd, cl)
#' slingAdjacency(sds)
#' @export
setGeneric(name = "slingAdjacency",
           signature = "x",
           def = function(x) standardGeneric("slingAdjacency"))

#' @title Methods for parameters used by Slingshot
#' @name slingParams
#' @description Extracts additional control parameters used by Slingshot in 
#' lineage inference and fitting simultaneous principal curves.
#'   
#' @param x an object containing \code{\link{slingshot}} output.
#' @return the list of additional parameters used by Slingshot.
#' @examples
#' data("slingshotExample")
#' sds <- slingshot(rd, cl, start.clus = '5')
#' slingParams(sds)
#' @export
setGeneric(name = "slingParams",
           signature = "x",
           def = function(x) standardGeneric("slingParams"))

#' @title Extract simultaneous principal curves
#' @name slingCurves
#' @description Extract the simultaneous principal curves from an object
#'   containing \code{\link{slingshot}} output.
#'   
#' @param x an object containing \code{\link{slingshot}} output.
#' @return the list of smooth lineage curves, each of which is a
#'   \code{\link{principal.curve}} object.
#' @examples
#' data("slingshotExample")
#' sds <- slingshot(rd, cl)
#' slingCurves(sds)
#' @export
setGeneric(name = "slingCurves",
           signature = "x",
           def = function(x) standardGeneric("slingCurves"))


#' @title Get Slingshot pseudotime values
#' @name slingPseudotime
#' 
#' @description Extract the matrix of pseudotime values or cells' weights along
#'   each lineage.
#' 
#' @param x an object containing \code{\link{slingshot}} output.
#' @param ... additional parameters to be passed to object-specific methods.
#' @return an \code{n} by \code{L} matrix representing each cell's pseudotime
#'   along each lineage.
#' @examples
#' data("slingshotExample")
#' sds <- slingshot(rd, cl)
#' slingPseudotime(sds)
#' @export
setGeneric(name = "slingPseudotime",
           signature = "x",
           def = function(x, ...) standardGeneric("slingPseudotime"))

#' @rdname slingPseudotime 
#' @return an \code{n} by \code{L} matrix of cell weights along
#'   each lineage.
#' @examples
#' data("slingshotExample")
#' sds <- slingshot(rd, cl)
#' slingCurveWeights(sds)
#' @export
setGeneric(name = "slingCurveWeights",
           signature = "x",
           def = function(x) standardGeneric("slingCurveWeights"))


# plotting
#' @title Plot Gene Expression by Pseudotime
#' @name plotGenePseudotime
#'   
#' @description Show the gene expression pattern for an individual gene along
#' lineages inferred by \code{\link{slingshot}}.
#' 
#' @export
setGeneric(
    name = "plotGenePseudotime",
    signature = c('gene','sds','exprs'),
    def = function(gene,
                   sds,
                   exprs,
                   ...) {
        standardGeneric("plotGenePseudotime")
    }
)

#' # plot3d generic, in case rgl is not loaded
#' #' @export
#' plot3d <- function(x, ...){
#'     UseMethod('plot3d')
#' }
