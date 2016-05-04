#' @title Infer Lineage Structure from Clustered Samples
#' 
#' @description Given a reduced data matrix nxp and a vector of cluster identities (optionally including -1's for "unclustered"), this function infers a forest structure on the clusters and returns paths through the forest that can be interpreted as lineages.
#' 
#' @param X numeric, the nxp matrix of samples in a reduced dimensionality space
#' @param clus.labels character, a vector of length n denoting cluster labels
#' @param omega (optional) numeric between 0 and 1 or Inf. This granularity parameter determines the distance between every real cluster and the artificial cluster, OMEGA. It is parameterized as a fraction of the largest distance between two real clusters (hence, any value greater than 1 would result in a single tree being returned and would be equivalent to setting omega = Inf, the default)
#' @param start.clus (optional) character, indicates the cluster(s) *from* which lineages will be drawn
#' @param end.clus (optional) character, indicates the cluster(s) which will be forced leaf nodes in their trees
#' 
#' @details The \code{forest} is learned by fitting a (possibly constrained) minimum-spanning tree on the clusters and the artificial cluster, OMEGA, which is a distance \code{omega} from every real cluster.
#'
#' Once the \code{forest} is known, lineages are identified in any tree with at least two clusters. For a given tree, if there is an annotated starting cluster, every possible path out of a starting cluster and ending in a leaf that isn't another starting cluster will be returned. If no starting cluster is annotated, every leaf will be considered as a potential starting cluster and whichever configuration produces the longest average lineage length (in terms of number of clusters included) will be returned.
#'
#' @return a list with L+2 items where L is the number of lineages identified. The first L items are character vectors with the names of the clusters included in that lineage. The last two items are \code{forest}, the connectivity matrix, and \code{C}, a clusters x lineages identity matrix.
#'
#' @examples
#' data("toy_data")
#' get_lineages(X, clus.labels)
#' get_lineages(X, clus.labels, start.clus = 'a')
#' 
#' @export
#' 
#' @importFrom ape mst
#' @import igraph
#' 

get_lineages <- function(X, clus.labels, omega = Inf, start.clus = NULL, end.clus = NULL, distout = FALSE){
  # set up, remove unclustered cells (-1's)
  X.original <- X
  X <- X[clus.labels != -1,]
  clus.labels <- clus.labels[clus.labels != -1]
  clusters <- unique(clus.labels)
  nclus <- length(clusters)
  
  ###############################
  ### get the connectivity matrix that defines the "forest"
  ###############################
  # get cluster centers
  centers <- t(sapply(clusters,function(clID){
    x.sub <- X[clus.labels == clID, ,drop = FALSE]
    return(colMeans(x.sub))
  }))
  
  # determine whether to use full covariance matrix or componentwise for distances
  min.clus.size <- min(table(clus.labels))
  if(min.clus.size <= ncol(X)){
    # componentwise
    message('Using adjusted diagonal covariance matrix')
    dist.fun <- function(clus1, clus2){
      mu1 <- colMeans(clus1)
      mu2 <- colMeans(clus2)
      diff <- mu1 - mu2
      s1 <- if(nrow(clus1) == 1) {matrix(0,ncol(clus1),ncol(clus1))} else {cov(clus1)}
      s1diag <- diag(s1)
      s2 <- if(nrow(clus2) == 1) {matrix(0,ncol(clus1),ncol(clus1))} else {cov(clus2)}
      s2diag <- diag(s2)
      jointCov <- s1 + s2
      jointCov[min.clus.size:ncol(X),] <- 0
      jointCov[,min.clus.size:ncol(X)] <- 0
      diag(jointCov) <- s1diag + s2diag
      if(all(jointCov == 0)) {jointCov <- diag(ncol(jointCov))}
      return(t(diff) %*% solve(jointCov) %*% diff)
    }
  }else{
    # full covariance
    dist.fun <- function(clus1, clus2){
      mu1 <- colMeans(clus1)
      mu2 <- colMeans(clus2)
      diff <- mu1 - mu2
      s1 <- cov(clus1)
      s2 <- cov(clus2)
      return(t(diff) %*% solve(s1 + s2) %*% diff)
    }
  }
  
  ### get pairwise cluster distance matrix
  D <- sapply(clusters,function(clID1){
    sapply(clusters,function(clID2){
      clus1 <- X[clus.labels == clID1, ,drop = FALSE]
      clus2 <- X[clus.labels == clID2, ,drop = FALSE]
      return(dist.fun(clus1, clus2))
    })
  })
  rownames(D) <- clusters
  colnames(D) <- clusters
  
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
    end.idx <- which(clusters %in% end.clus)
    mstree <- ape::mst(D[-end.idx, -end.idx])
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
      dists <- D[! rownames(D) %in% end.clus, cl.idx]
      closest <- names(dists)[which.min(dists)] # get closest non-endpoint cluster
      closest.idx <- which.max(clusters == closest)
      forest[cl.idx, closest.idx] <- 1
      forest[closest.idx, cl.idx] <- 1
    }
  }else{
    forest <- mstree
  }
  forest <- forest[1:nclus, 1:nclus] # remove OMEGA
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
  trees <- trees[order(sapply(trees,length),decreasing = T)]
  
  # identify lineages (paths through trees)
  for(tree in trees){
    if(length(tree) == 1){
      next # don't draw a lineage for a single-cluster tree
    }
    tree.ind <- rownames(forest) %in% tree
    tree.graph <- forest[tree.ind, tree.ind]
    degree <- rowSums(tree.graph)
    g <- igraph::graph.adjacency(tree.graph, mode="undirected")
    
    # if you have starting cluster(s) in this tree, draw lineages to each leaf
    if(! is.null(start.clus)){
      if(sum(start.clus %in% tree) > 0){
        starts <- start.clus[start.clus %in% tree]
        ends <- rownames(tree.graph)[degree == 1 & ! rownames(tree.graph) %in% starts]
        for(st in starts){
          paths <- igraph::shortest_paths(g, from = st, to = ends, mode = 'out', output = 'vpath')$vpath
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
          paths <- igraph::shortest_paths(g, from = l, to = ends, mode = 'out', output = 'vpath')$vpath
          mean(sapply(paths, length))
        })
        st <- names(avg.lineage.length)[which.max(avg.lineage.length)]
        ends <- leaves[leaves != st]
        paths <- igraph::shortest_paths(g, from = st, to = ends, mode = 'out', output = 'vpath')$vpath
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
  out <- lineages
  # include "forest" and clusters x lineages (C) matrices
  out$forest <- forest
  C <- sapply(lineages,function(lin){
    sapply(clusters,function(clID){
      as.numeric(clID %in% lin)
    })
  })
  rownames(C) <- clusters
  # should probably come up with a better name than C
  out$C <- C
  if(distout){
    out$dist <- D
  }
  return(out)
}
