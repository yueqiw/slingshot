#' @title Infer Lineage Structure from Clustered Samples
#' 
#' @description Given a reduced-dimension data matrix \code{n} by \code{p} and a vector of
#' cluster identities (potentially including -1's for "unclustered"), this function
#' infers a forest structure on the clusters and returns paths through the forest
#' that can be interpreted as lineages.
#' 
#' @param reducedDim numeric, the \code{n} by \code{p} matrix of samples in a reduced 
#' dimensionality space.
#' @param clus.labels character, a vector of length \code{n} denoting cluster labels,
#'   optionally including \code{-1}'s for "unclustered." If \code{reducedDim} is a
#'   \code{SlingshotDataSet}, cluster labels will be taken from it.
#' @param start.clus (optional) character, indicates the cluster(s) *from* which
#'   lineages will be drawn.
#' @param end.clus (optional) character, indicates the cluster(s) which will be
#'   forced leaf nodes in their trees.
#' @param dist.fun (optional) function, method for calculating distances between
#'   clusters. Must take two matrices as input, corresponding to points in
#'   reduced-dimensional space. If the minimum cluster size is larger than the
#'   number dimensions, the default is to use the joint covariance matrix to find
#'   squared distance between cluster centers. If not, the default is to use the
#'   diagonal of the joint covariance matrix.
#' @param omega (optional) numeric between 0 and 1. This granularity
#'   parameter determines the distance between every real cluster and the artificial
#'   cluster, OMEGA. It is parameterized as a fraction of the largest distance
#'   between two real clusters (hence, any value greater than 1 would result in a
#'   single tree being returned and would be equivalent to setting \code{omega = Inf},
#'   which is the default behavior).
#' 
#' @details The \code{connectivity} matrix is learned by fitting a
#'   (possibly constrained) minimum-spanning tree on the clusters and the artificial
#'   cluster, OMEGA, which is a fixed distance away from every real cluster. This 
#'   effectively limits the maximum branch length in the MST to twice the chosen 
#'   distance, meaning that the output may contain multiple trees.
#'
#' @details Once the \code{connectivity} is known, lineages are identified in any tree
#'   with at least two clusters. For a given tree, if there is an annotated starting
#'   cluster, every possible path out of a starting cluster and ending in a leaf
#'   that isn't another starting cluster will be returned. If no starting cluster is
#'   annotated, every leaf will be considered as a potential starting cluster and
#'   whichever configuration produces the longest average lineage length (in terms
#'   of number of clusters included) will be returned.
#'
#' @return An object of class \code{\link{SlingshotDataSet}} containing the 
#' arguments provided to \code{getLineages} as well as the following new
#' elements:
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
#' sds <- getLineages(reducedDim, clus.labels, start.clus = '5')
#' 
#' plot(reducedDim, col = clus.labels, asp = 1)
#' lines(sds, type = 'l', lwd = 3)
#' 
#' @export
#'
#' @importFrom igraph graph.adjacency
#' @importFrom igraph shortest_paths
#' @importFrom ape mst
#' 
setMethod(f = "getLineages",
          signature = signature(reducedDim = "matrix", clus.labels = "character"),
          definition = function(reducedDim, clus.labels,
                   start.clus = NULL, end.clus = NULL,
                   dist.fun = NULL, omega = NULL){
            
            X <- reducedDim
            # CHECKS
            clus.labels <- as.character(clus.labels)
            X <- as.matrix(X)
            if(nrow(X) != length(clus.labels)){
              stop('nrow(reducedDim) must equal length(clus.labels).')
            }
            if(any(is.na(X))){
              stop('reducedDim cannot contain missing values.')
            }
            if(!all(apply(X,2,is.numeric))){
              stop('reducedDim must only contain numeric values.')
            }
            if(is.null(rownames(X))){
              rownames(X) <- paste('cell',seq_len(nrow(X)),sep='-')
            }
            if(is.null(colnames(X))){
              colnames(X) <- paste('dim',seq_len(ncol(X)),sep='-')
            }
            
            # set up, remove unclustered cells (-1's)
            X.original <- X
            X <- X[clus.labels != -1, ,drop = FALSE]
            clus.labels <- clus.labels[clus.labels != -1]
            clusters <- unique(clus.labels)
            nclus <- length(clusters)
            if(!is.null(start.clus)){
              start.clus <- as.character(start.clus)
            }
            if(!is.null(end.clus)){
              end.clus <- as.character(end.clus)
            }
            
            
            ### get the connectivity matrix
            # get cluster centers
            centers <- t(sapply(clusters,function(clID){
              x.sub <- X[clus.labels == clID, ,drop = FALSE]
              return(colMeans(x.sub))
            }))
            
            # determine the distance function
            if(is.null(dist.fun)){
              min.clus.size <- min(table(clus.labels))
              if(min.clus.size <= ncol(X)){
                message('Using diagonal covariance matrix')
                dist.fun <- function(c1,c2) .dist_clusters_diag(c1,c2)
              }else{
                message('Using full covariance matrix')
                dist.fun <- function(c1,c2) .dist_clusters_full(c1,c2)
              }
            }
            
            ### get pairwise cluster distance matrix
            D <- as.matrix(sapply(clusters,function(clID1){
              sapply(clusters,function(clID2){
                clus1 <- X[clus.labels == clID1, ,drop = FALSE]
                clus2 <- X[clus.labels == clID2, ,drop = FALSE]
                return(dist.fun(clus1, clus2))
              })
            }))
            rownames(D) <- clusters
            colnames(D) <- clusters
            
            # if infinite, set omega to largest distance + 1
            if(is.null(omega)){
              omega <- max(D) + 1
            }else{
              if(omega >= 0 && omega <= 1){
                omega <- omega * max(D)
              }else{
                stop("omega must be between 0 and 1")
              }
            }
            D <- rbind(D, rep(omega, ncol(D)) )
            D <- cbind(D, c(rep(omega, ncol(D)), 0) )
            
            # draw MST on cluster centers + OMEGA (possibly excluding endpoint clusters)
            if(! is.null(end.clus)){
              end.idx <- which(clusters %in% end.clus)
              mstree <- ape::mst(D[-end.idx, -end.idx, drop = FALSE])
            }else{
              mstree <- ape::mst(D)
            }
            # (add in endpoint clusters)
            if(! is.null(end.clus)){
              forest <- D
              forest[forest != 0] <- 0
              forest[-end.idx, -end.idx] <- mstree
              for(clID in end.clus){
                cl.idx <- which(clusters == clID)
                dists <- D[! rownames(D) %in% end.clus, cl.idx, drop = FALSE]
                closest <- names(dists)[which.min(dists)] # get closest non-endpoint cluster
                closest.idx <- which.max(clusters == closest)
                forest[cl.idx, closest.idx] <- 1
                forest[closest.idx, cl.idx] <- 1
              }
            }else{
              forest <- mstree
            }
            forest <- forest[1:nclus, 1:nclus, drop = FALSE] # remove OMEGA
            rownames(forest) <- clusters
            colnames(forest) <- clusters
            
            ###############################
            ### use the "forest" to define lineages
            ###############################
            lineages <- list()
            
            # identify trees
            unused <- rownames(forest)
            trees <- list()
            ntree <- 0
            while(length(unused) > 0){
              ntree <- ntree + 1
              newtree <- .get_connections(unused[1], forest)
              trees[[ntree]] <- newtree
              unused <- unused[! unused %in% newtree]
            }
            trees <- trees[order(sapply(trees,length),decreasing = TRUE)]
            
            # identify lineages (paths through trees)
            for(tree in trees){
              if(length(tree) == 1){
                lineages[[length(lineages)+1]] <- tree
                next
              }
              tree.ind <- rownames(forest) %in% tree
              tree.graph <- forest[tree.ind, tree.ind, drop = FALSE]
              degree <- rowSums(tree.graph)
              g <- graph.adjacency(tree.graph, mode="undirected")
              
              # if you have starting cluster(s) in this tree, draw lineages to each leaf
              if(! is.null(start.clus)){
                if(sum(start.clus %in% tree) > 0){
                  starts <- start.clus[start.clus %in% tree]
                  ends <- rownames(tree.graph)[degree == 1 & ! rownames(tree.graph) %in% starts]
                  for(st in starts){
                    paths <- shortest_paths(g, from = st, to = ends, mode = 'out', output = 'vpath')$vpath
                    for(p in paths){
                      lineages[[length(lineages)+1]] <- names(p)
                    }
                  }
                }else{
                  # else, need a criteria for picking root
                  # highest average length (~parsimony, but this was just the easiest thing I came up with)
                  leaves <- rownames(tree.graph)[degree == 1]
                  avg.lineage.length <- sapply(leaves,function(l){
                    ends <- leaves[leaves != l]
                    paths <- shortest_paths(g, from = l, to = ends, mode = 'out', output = 'vpath')$vpath
                    mean(sapply(paths, length))
                  })
                  st <- names(avg.lineage.length)[which.max(avg.lineage.length)]
                  ends <- leaves[leaves != st]
                  paths <- shortest_paths(g, from = st, to = ends, mode = 'out', output = 'vpath')$vpath
                  for(p in paths){
                    lineages[[length(lineages)+1]] <- names(p)
                  }
                }
              }else{
                # else, need a criteria for picking root
                # highest average length (~parsimony, but this was just the easiest thing I came up with)
                leaves <- rownames(tree.graph)[degree == 1]
                avg.lineage.length <- sapply(leaves,function(l){
                  ends <- leaves[leaves != l]
                  paths <- shortest_paths(g, from = l, to = ends, mode = 'out', output = 'vpath')$vpath
                  mean(sapply(paths, length))
                })
                st <- names(avg.lineage.length)[which.max(avg.lineage.length)]
                ends <- leaves[leaves != st]
                paths <- shortest_paths(g, from = st, to = ends, mode = 'out', output = 'vpath')$vpath
                for(p in paths){
                  lineages[[length(lineages)+1]] <- names(p)
                }
              }
            }
            # sort by number of clusters included
            lineages <- lineages[order(sapply(lineages, length), decreasing = TRUE)]
            names(lineages) <- paste('Lineage',1:length(lineages),sep='')
            
            lineage.control <- list()
            first <- unique(sapply(lineages,function(l){ l[1] }))
            last <- unique(sapply(lineages,function(l){ l[length(l)] }))
            
            lineage.control$start.clus <- first
            lineage.control$end.clus <- last
            
            start.given <- first %in% start.clus
            end.given <- last %in% end.clus
            lineage.control$start.given <- start.given
            lineage.control$end.given <- end.given
            
            lineage.control$dist <- D[1:nclus,1:nclus, drop = FALSE]
            connectivity <- forest
            
            out <- SlingshotDataSet(reducedDim = X, clus.labels = clus.labels, lineages = lineages, connectivity = connectivity, lineage.control = lineage.control)
            
            validObject(out)
            return(out)
          }
)



setMethod(f = "getLineages",
          signature = signature(reducedDim = "SlingshotDataSet", clus.labels = "ANY"),
          definition = function(reducedDim,
                                clus.labels = reducedDim@clus.labels,
                                start.clus = NULL, end.clus = NULL,
                                dist.fun = NULL, omega = NULL){
            return(getLineages(reducedDim = reducedDim@reducedDim, 
                               clus.labels = reducedDim@clus.labels, 
                               start.clus = start.clus, end.clus = end.clus,
                               dist.fun = dist.fun, omega = omega))
          })

setMethod(f = "getLineages",
          signature = signature(reducedDim = "data.frame", clus.labels = "ANY"),
          definition = function(reducedDim, clus.labels, 
                                start.clus = NULL, end.clus = NULL,
                                dist.fun = NULL, omega = NULL){
            RD <- as.matrix(reducedDim)
            rownames(RD) <- rownames(reducedDim)
            return(getLineages(reducedDim = RD, 
                               clus.labels = clus.labels, 
                               start.clus = start.clus, end.clus = end.clus,
                               dist.fun = dist.fun, omega = omega))
          })

setMethod(f = "getLineages",
          signature = signature(reducedDim = "matrix", clus.labels = "numeric"),
          definition = function(reducedDim, clus.labels, 
                                start.clus = NULL, end.clus = NULL,
                                dist.fun = NULL, omega = NULL){
            return(getLineages(reducedDim = reducedDim, 
                               clus.labels = as.character(clus.labels), 
                               start.clus = start.clus, end.clus = end.clus,
                               dist.fun = dist.fun, omega = omega))
          })

setMethod(f = "getLineages",
          signature = signature(reducedDim = "matrix", clus.labels = "factor"),
          definition = function(reducedDim, clus.labels, 
                                start.clus = NULL, end.clus = NULL,
                                dist.fun = NULL, omega = NULL){
            return(getLineages(reducedDim = reducedDim, 
                               clus.labels = as.character(clus.labels), 
                               start.clus = start.clus, end.clus = end.clus,
                               dist.fun = dist.fun, omega = omega))
          })