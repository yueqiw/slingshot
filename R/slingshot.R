#' @title Perform lineage inference with Slingshot
#' @aliases slingshot
#' 
#' @description Given a reduced-dimensional data matrix \code{n} by \code{p} and a vector of
#' cluster labels (potentially including \code{-1}'s for "unclustered"), this function
#' performs lineage inference using a cluster-based minimum spanning tree and 
#' constructing simulatenous principal curves for branching paths through the tree.
#' 
#' @description This wrapper function performs lineage inference in two steps: (1) 
#' identify lineage structure with a cluster-based minimum spanning tree with the
#' \code{\link{getLineages}} function and (2) construct smooth representations of 
#' each lineage using simultaneous principal curves from the function
#' \code{\link{getCurves}}.
#' 
#' @param reducedDim numeric matrix or \code{SlingshotDataSet} object containing low-
#' dimensional representation of single cells.
#' @param clusLabels character, a vector of length \code{n} denoting cluster labels,
#' optionally including \code{-1}'s for "unclustered." If \code{reducedDim} is a
#' \code{SlingshotDataSet}, cluster labels will be taken from it.
#' @param start.clus (optional) character, indicates the cluster(s) of origin. Lineages
#' will be represented by paths coming out of this cluster.
#' @param end.clus (optional) character, indicates the cluster(s) which will be
#'   forced leaf nodes. This introduces a constraint on the MST algorithm.
#' @param dist.fun (optional) function, method for calculating distances between
#'   clusters. Must take two matrices as input, corresponding to subsets of 
#'   \code{reducedDim}. If the minimum cluster size is larger than the
#'   number dimensions, the default is to use the joint covariance matrix to find
#'   squared distance between cluster centers. If not, the default is to use the
#'   diagonal of the joint covariance matrix.
#' @param omega (optional) numeric, granularity parameter between 0 and 1. \code{omega}
#'   determines the distance between every real cluster and the artificial
#'   cluster, \code{OMEGA}. It is parameterized as a fraction of the largest distance
#'   between two real clusters (hence, any value greater than 1 would result in a
#'   single tree being returned and would be equivalent to setting \code{omega = Inf},
#'   which is the default behavior).
#' @param lineages list, denotes which lineages each cluster is a part of and
#'   contains the \code{K x K} connectivity matrix constructed on the clusters by
#'   \code{\link{get_lineages}}.
#' @param thresh numeric, determines the convergence criterion. Percent change in
#'   the total distance from cells to their projections along curves must be less
#'   than \code{thresh}. Default is \code{0.001}, similar to
#'   \code{\link{principal.curve}}.
#' @param maxit numeric, maximum number of iterations, see
#'   \code{\link{principal.curve}}.
#' @param stretch numeric factor by which curves can be extrapolated beyond
#'   endpoints. Default is \code{2}, see \code{\link{principal.curve}}.
#' @param smoother, choice of scatter plot smoother. Same as
#'   \code{\link{principal.curve}}, but \code{"lowess"} option is replaced with
#'   \code{"loess"} for additional flexibility.
#' @param shrink logical or numeric between 0 and 1, determines whether and how 
#'   much to shrink branching lineages toward their average prior to the split.
#' @param extend character, how to handle root and leaf clusters of lineages when
#'   constructing the initial, piece-wise linear curve. Accepted values are
#'   \code{'y'} (default), \code{'n'}, and \code{'pc1'}. See 'Details' for more.
#' @param reweight logical, whether to allow cells shared between lineages to be
#'   reweighted during curve-fitting. If \code{TRUE}, cells shared between
#'   lineages will be weighted by: distance to nearest curve / distance to curve.
#' @param drop.multi logical, whether to drop shared cells from lineages which do 
#'   not fit them well. If \code{TRUE}, shared cells with a distance to one 
#'   lineage above the 90th percentile and another below the 50th will be dropped
#'   from the further lineage.
#' @param shrink.method character denoting how to determine the appropriate amount
#'   of shrinkage for a branching lineage. Accepted values are the same as for
#'   \code{kernel} in \code{\link{density}} (default is \code{"cosine"}), as well 
#'   as \code{"tricube"} and \code{"density"}. See 'Details' for more.
#' @param ... Additional parameters to pass to scatter plot smoothing function,
#'   \code{smoother}.
#' 
#' @details The \code{connectivity} matrix is learned by fitting a
#'   (possibly constrained) minimum-spanning tree on the clusters and the artificial
#'   cluster, OMEGA, which is a distance \code{omega} from every real cluster.
#'
#' @details Once the \code{connectivity} is known, lineages are identified in any tree
#'   with at least two clusters. For a given tree, if there is an annotated starting
#'   cluster, every possible path out of a starting cluster and ending in a leaf
#'   that isn't another starting cluster will be returned. If no starting cluster is
#'   annotated, every leaf will be considered as a potential starting cluster and
#'   whichever configuration produces the longest average lineage length (in terms
#'   of number of clusters included) will be returned.
#'
#' @details When there is only a single lineage, the curve-fitting algorithm is
#'   nearly identical to that of \code{\link{principal.curve}}. When there are 
#'   multiple lineages and \code{shrink=TRUE}, an additional step is added to the 
#'   iterative procedure, forcing curves to be similar in the neighborhood of 
#'   shared points (ie., before they branch).
#'   
#' @details The \code{extend} argument determines how to construct the piece-wise
#'   linear curve used to initiate the recursive algorithm. The initial curve is
#'   always based on the lines between cluster centers and if \code{extend = 'n'}, 
#'   this curve will terminate at the center of the endpoint clusters. Setting 
#'   \code{extend = 'y'} will allow the first and last segments to extend beyond
#'   the cluster center to the orthogonal projection of the furthest point. Setting
#'   \code{extend = 'pc1'} is similar to \code{'y'}, but uses the first principal
#'   component of the cluster to determine the direction of the curve beyond the
#'   cluster center. These options typically have little to no impact on the final
#'   curve, but can occasionally help with stability issues.
#'   
#' @details *** Explain shrink.method ***
#' 
#' @details *** acknowledge princurve ***
#' 
#' @return An object of class \code{\link{SlingshotDataSet}} containing the 
#' arguments provided to \code{getLineages} as well as the following output:
#' \itemize{
#' \item{\code{lineages}}{ a list of \code{L} items, where \code{L} is the number 
#'   of lineages identified. Each lineage is represented by a character vector 
#'   with the names of the clusters included in that lineage, in order.}
#' \item{\code{connectivity}}{ the inferred cluster connectivity matrix.}
#' \item{\code{lineage.control$start.given},\code{lineage.control$end.given}}
#'   { logical values indicating whether the starting and ending clusters were 
#'   specified a priori.} 
#' \item{\code{lineage.control$dist}}{ the pairwise cluster distance matrix.}}
#'
#' @examples
#' data("slingshotExample")
#' sds <- slingshot(reducedDim, clusLabels, start.clus = '5')
#' 
#' plot(reducedDim, col = clusLabels, asp = 1)
#' lines(sds, lwd = 3)
#' 
#' @export
#' 
setMethod(f = "slingshot",
          signature = signature(reducedDim = "matrix", clusLabels = "character"),
          definition = function(reducedDim, clusLabels,
                                start.clus = NULL, end.clus = NULL,
                                dist.fun = NULL, omega = NULL,
                                lineages = list(),
                                shrink = TRUE,
                                extend = 'y',
                                reweight = TRUE,
                                drop.multi = TRUE,
                                thresh = 0.001, maxit = 15, stretch = 2,
                                smoother = 'smooth.spline',
                                shrink.method = 'cosine', ...){
            sds <- getLineages(reducedDim, clusLabels,
                               start.clus = start.clus, end.clus = end.clus,
                               dist.fun = dist.fun, omega = omega)
            sds <- getCurves(sds,
                             shrink = shrink, extend = extend,
                             reweight = reweight, drop.multi = drop.multi,
                             thresh = thresh, maxit = maxit,
                             stretch = stretch, smoother = smoother,
                             shrink.method = shrink.method, ...)
            return(sds)
          }
)


setMethod(f = "slingshot",
          signature = signature(reducedDim = "SlingshotDataSet", clusLabels = "ANY"),
          definition = function(reducedDim,
                                clusLabels = reducedDim@clusLabels,
                                start.clus = NULL, end.clus = NULL,
                                dist.fun = NULL, omega = NULL,
                                lineages = list(),
                                shrink = TRUE,
                                extend = 'y',
                                reweight = TRUE,
                                drop.multi = TRUE,
                                thresh = 0.001, maxit = 15, stretch = 2,
                                smoother = 'smooth.spline',
                                shrink.method = 'cosine', ...){
            return(slingshot(reducedDim = reducedDim@reducedDim, 
                             clusLabels = reducedDim@clusLabels, 
                             start.clus = start.clus, end.clus = end.clus,
                             dist.fun = dist.fun, omega = omega,
                             shrink = shrink, extend = extend,
                             reweight = reweight, drop.multi = drop.multi,
                             thresh = thresh, maxit = maxit,
                             stretch = stretch, smoother = smoother,
                             shrink.method = shrink.method, ...))
          })

setMethod(f = "slingshot",
          signature = signature(reducedDim = "data.frame", clusLabels = "ANY"),
          definition = function(reducedDim, clusLabels, 
                                start.clus = NULL, end.clus = NULL,
                                dist.fun = NULL, omega = NULL,
                                lineages = list(),
                                shrink = TRUE,
                                extend = 'y',
                                reweight = TRUE,
                                drop.multi = TRUE,
                                thresh = 0.001, maxit = 15, stretch = 2,
                                smoother = 'smooth.spline',
                                shrink.method = 'cosine', ...){
            RD <- as.matrix(reducedDim)
            rownames(RD) <- rownames(reducedDim)
            return(slingshot(reducedDim = RD, 
                             clusLabels = clusLabels, 
                             start.clus = start.clus, end.clus = end.clus,
                             dist.fun = dist.fun, omega = omega,
                             shrink = shrink, extend = extend,
                             reweight = reweight, drop.multi = drop.multi,
                             thresh = thresh, maxit = maxit,
                             stretch = stretch, smoother = smoother,
                             shrink.method = shrink.method, ...))
          })

setMethod(f = "slingshot",
          signature = signature(reducedDim = "matrix", clusLabels = "numeric"),
          definition = function(reducedDim, clusLabels, 
                                start.clus = NULL, end.clus = NULL,
                                dist.fun = NULL, omega = NULL,
                                lineages = list(),
                                shrink = TRUE,
                                extend = 'y',
                                reweight = TRUE,
                                drop.multi = TRUE,
                                thresh = 0.001, maxit = 15, stretch = 2,
                                smoother = 'smooth.spline',
                                shrink.method = 'cosine', ...){
            return(slingshot(reducedDim = reducedDim, 
                             clusLabels = as.character(clusLabels), 
                             start.clus = start.clus, end.clus = end.clus,
                             dist.fun = dist.fun, omega = omega,
                             shrink = shrink, extend = extend,
                             reweight = reweight, drop.multi = drop.multi,
                             thresh = thresh, maxit = maxit,
                             stretch = stretch, smoother = smoother,
                             shrink.method = shrink.method, ...))
          })

setMethod(f = "slingshot",
          signature = signature(reducedDim = "matrix", clusLabels = "factor"),
          definition = function(reducedDim, clusLabels, 
                                start.clus = NULL, end.clus = NULL,
                                dist.fun = NULL, omega = NULL,
                                lineages = list(),
                                shrink = TRUE,
                                extend = 'y',
                                reweight = TRUE,
                                drop.multi = TRUE,
                                thresh = 0.001, maxit = 15, stretch = 2,
                                smoother = 'smooth.spline',
                                shrink.method = 'cosine', ...){
            return(slingshot(reducedDim = reducedDim, 
                             clusLabels = as.character(clusLabels), 
                             start.clus = start.clus, end.clus = end.clus,
                             dist.fun = dist.fun, omega = omega,
                             shrink = shrink, extend = extend,
                             reweight = reweight, drop.multi = drop.multi,
                             thresh = thresh, maxit = maxit,
                             stretch = stretch, smoother = smoother,
                             shrink.method = shrink.method, ...))
          })