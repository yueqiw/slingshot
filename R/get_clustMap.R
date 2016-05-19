#' @title Map Relationships Between Clusters
#' 
#' @description Given a reduced data matrix nxp and a vector of cluster identities (optionally including -1's for "unclustered"), this function infers a forest structure on the clusters.
#' 
#' @param X numeric, the nxp matrix of samples in a reduced dimensionality space
#' @param cluster character, a vector of length n denoting cluster labels
#' @param map.args (optional) list of additional arguments for controlling how mapping is performed.

#' @details The \code{forest} is learned by fitting a (possibly constrained) minimum-spanning tree on the clusters and the artificial cluster, OMEGA, which is a distance \code{omega} from every real cluster.
#'
#' Optional arguments to be passed via \code{map.args}:
#' @param omega (optional) numeric between 0 and 1 or Inf. This granularity parameter determines the distance between every real cluster and the artificial cluster, OMEGA. It is parameterized as a fraction of the largest distance between two real clusters (hence, any value greater than 1 would result in a single tree being returned and would be equivalent to setting omega = Inf, the default)
#' @param start.clus (optional) character, indicates the cluster(s) *from* which lineages will be drawn
#' @param end.clus (optional) character, indicates the cluster(s) which will be forced leaf nodes in their trees
#'
#' @return A KxK connectivity matrix, where K is the number of clusters, representing the learned graph structure.
#'
#' @examples
#' data("toy_data")
#' get_clustMap(X, cluster)
#' 
#' @export
#' 
#' @importFrom ape mst
#' @import igraph
#' 

setMethod(
  f = "get_clustMap",
  signature = signature(X = "matrix", cluster = "character"),
  definition = function(X, cluster, dist.fun = NULL, map.args = list()) {
    if(is.null(map.args$omega)){
      omega <- Inf
    }else{
      omega <- map.args$omega
    }
    start.clus <- map.args$start.clus
    end.clus <- map.args$end.clus
    # set up, remove unclustered cells (-1's)
    X.original <- X
    X <- X[cluster != -1,]
    cluster <- cluster[cluster != -1]
    clus.names <- unique(cluster)
    nclus <- length(clus.names)
    
    ###############################
    ### get the connectivity matrix that defines the "map"
    ###############################
    # get cluster centers
    centers <- t(sapply(clus.names,function(clID){
      x.sub <- X[cluster == clID, ,drop = FALSE]
      return(colMeans(x.sub))
    }))
    
    # determine whether to use full covariance matrix or componentwise for distances
    if(is.null(dist.fun)){
      min.clus.size <- min(table(cluster))
      if(min.clus.size <= ncol(X)){
        message('Using adjusted diagonal covariance matrix')
        dist.fun <- function(c1,c2) dist_clusters_partial(c1,c2,k=min.clus.size)
      }else{
        message('Using full covariance matrix')
        dist.fun <- function(c1,c2) dist_clusters_full(c1,c2)
      }
    }
    
    ### get pairwise cluster distance matrix
    D <- sapply(clus.names,function(clID1){
      sapply(clus.names,function(clID2){
        clus1 <- X[cluster == clID1, ,drop = FALSE]
        clus2 <- X[cluster == clID2, ,drop = FALSE]
        return(dist.fun(clus1, clus2))
      })
    })
    rownames(D) <- clus.names
    colnames(D) <- clus.names
    
    # if infinite, set omega to largest distance + 1
    if(omega == Inf){
      omega <- max(D) + 1
    }else{
      if(omega >= 0 & omega <= 1){
        omega <- omega * max(D)
      }else{
        stop("as it's currently set up, omega must be between 0 and 1")
      }
    }
    D <- rbind(D, rep(omega, ncol(D)) )
    D <- cbind(D, c(rep(omega, ncol(D)), 0) )
    
    # draw MST on cluster centers + OMEGA (possibly excluding endpoint clusters)
    if(! is.null(end.clus)){
      end.idx <- which(clus.names %in% end.clus)
      mstree <- ape::mst(D[-end.idx, -end.idx])
    }else{
      mstree <- ape::mst(D)
    }
    # (add in endpoint clusters)
    if(! is.null(end.clus)){
      map <- D
      map[map != 0] <- 0
      map[-end.idx, -end.idx] <- mstree
      for(clID in end.clus){
        cl.idx <- which(clus.names == clID)
        dists <- D[! rownames(D) %in% end.clus, cl.idx]
        closest <- names(dists)[which.min(dists)] # get closest non-endpoint cluster
        closest.idx <- which.max(clus.names == closest)
        map[cl.idx, closest.idx] <- 1
        map[closest.idx, cl.idx] <- 1
      }
    }else{
      map <- mstree
    }
    map <- map[1:nclus, 1:nclus] # remove OMEGA
    rownames(map) <- clus.names
    colnames(map) <- clus.names
    
    return(map)
  }
)

setMethod(
  f = "get_clustMap",
  signature = signature(X = "CellLineages"),
  definition = function(X, ...) {
    out <- X
    out@clustMap <- get_clustMap(reducedDim(X),colData(X)$cluster,...)
    return(out)
  }
)
