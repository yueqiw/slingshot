#' @title Infer Lineage Structure from Clustered Samples
#' 
#' @description Given a reduced-dimension data matrix \code{n x p} and a vector of
#' cluster identities (potentially including -1's for "unclustered"), this function
#' infers a forest structure on the clusters and returns paths through the forest
#' that can be interpreted as lineages.
#' 
#' @param X numeric, the \code{n x p} matrix of samples in a reduced dimensionality 
#'   space.
#' @param clus.labels character, a vector of length n denoting cluster labels,
#'   potentially including -1's for "unclustered."
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
#' @param omega (optional) numeric between 0 and 1 or Inf. This granularity
#'   parameter determines the distance between every real cluster and the artificial
#'   cluster, OMEGA. It is parameterized as a fraction of the largest distance
#'   between two real clusters (hence, any value greater than 1 would result in a
#'   single tree being returned and would be equivalent to setting omega = Inf, the
#'   default)
#' @param distout (optional) logical, indicating whether the distance matrix between
#'   the clusters should be included with output.
#' 
#' @details The connectivity matrix, denoted \code{forest}, is learned by fitting a
#'   (possibly constrained) minimum-spanning tree on the clusters and the artificial
#'   cluster, OMEGA, which is a distance \code{omega} from every real cluster.
#'
#' @details Once the \code{forest} is known, lineages are identified in any tree
#'   with at least two clusters. For a given tree, if there is an annotated starting
#'   cluster, every possible path out of a starting cluster and ending in a leaf
#'   that isn't another starting cluster will be returned. If no starting cluster is
#'   annotated, every leaf will be considered as a potential starting cluster and
#'   whichever configuration produces the longest average lineage length (in terms
#'   of number of clusters included) will be returned.
#'
#' @return a list with at least \code{L} items where \code{L} is the number of
#'   lineages identified. Each lineage is represented by a character vector with the
#'   names of the clusters included in that lineage. Additional items may include:
#'   \itemize{\item{\code{forest}}{ the inferred connectivity matrix.} 
#'   \item{\code{C}}{ a \code{clusters x lineages} matrix indicating each cluster's
#'   inclusion in each lineage.} \item{\code{start.clus},\code{end.clus}}{ the 
#'   starting and ending cluster(s)} \item{\code{start.given},\code{end.given}}
#'   { logical values indicating whether the starting and ending clusters were 
#'   specified a priori} \item{\code{dist}}{ the pairwise cluster distance matrix.}}
#'
#' @examples
#' data("slingshot_example")
#' get_lineages(X, clus.labels)
#' lin <- get_lineages(X, clus.labels, start.clus = 'a')
#' plot_tree(X, clus.labels, lin)
#' 
#' @export
#'
#' @importFrom igraph graph.adjacency
#' @importFrom igraph shortest_paths
#' @importFrom ape mst
#' 

get_lineages <- function(X, clus.labels, start.clus = NULL, end.clus = NULL, dist.fun = NULL, omega = Inf, distout = FALSE){
  # CHECKS
  clus.labels <- as.character(clus.labels)
  X <- as.matrix(X)
  if(nrow(X) != length(clus.labels)){
    stop('nrow(X) must equal length(clus.labels)')
  }
  if(is.null(rownames(X))){
    rownames(X) <- paste('cell',seq_len(nrow(X)),sep='-')
  }
  if(is.null(colnames(X))){
    colnames(X) <- paste('dim',seq_len(ncol(X)),sep='-')
  }
  
  # set up, remove unclustered cells (-1's)
  X.original <- X
  X <- X[clus.labels != -1,]
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
  if(omega == Inf){
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
  trees <- trees[order(sapply(trees,length),decreasing = T)]
  
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
  out <- lineages
  names(out) <- paste('lineage',1:length(lineages),sep='')

  first <- unique(sapply(out,function(l){l[1]}))
  last <- unique(sapply(out,function(l){l[length(l)]}))
  start.given <- first %in% start.clus
  end.given <- last %in% end.clus
  out$start.clus <- first
  out$start.given <- start.given
  out$end.clus <- last
  out$end.given <- end.given
  
  # include "forest" and clusters x lineages (C) matrices
  out$forest <- forest
  C <- as.matrix(sapply(lineages,function(lin){
    sapply(clusters,function(clID){
      as.numeric(clID %in% lin)
    })
  }))
  rownames(C) <- clusters
  # should probably come up with a better name than C
  out$C <- C
  if(distout){
    out$dist <- D
  }
  return(out)
}


#' @title Construct Smooth Curves for Each Lineage
#' 
#' @description This function takes a reduced data matrix \code{n x p}, a vector of
#'   cluster identities (optionally including -1's for "unclustered"), and a set of
#'   lineages consisting of paths through a forest constructed on the clusters. It
#'   constructs smooth curves for each lineage and returns the points along these
#'   curves corresponding to the orthogonal projections of each data point, along
#'   with corresponding arclength (\code{pseudotime} or \code{lambda}) values.
#' 
#' @param X numeric, the \code{n x p} matrix of samples in a reduced dimensionality
#'   space.
#' @param clus.labels character, a vector of length n denoting cluster labels.
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
#'   constructing the initial, piece-wise linear curve. Accepted values are 'y'
#'   (default), 'n', and 'pc1'. See 'Details' for more.
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
#' @details When there is only a single lineage, the curve-fitting algorithm is
#'   identical to that of \code{\link{principal.curve}}. When there are multiple
#'   lineages and \code{shrink=TRUE}, an additional step is added to the iterative 
#'   procedure, forcing curves to be similar in the neighborhood of shared points
#'   (ie., before they branch).
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
#'
#' @return A list of length \code{L}, equal to the number of lineages. Each element
#'   is an object of class \code{principal.curve} containing the following objects: 
#'   \itemize{ \item{\code{s}}{ a matrix of points along the curve corresponding to
#'   the projections of points in \code{X} onto the curve, ordered by pseudotime.}
#'   \item{\code{lambda}}{ a vector of pseudotime values in the same order as 
#'   \code{s}, representing each point's arclength along the curve.}
#'   \item{\code{dist}}{ the total squared distance between points used in the 
#'   construction of the curve and their projections onto the curve.}
#'   \item{\code{pseudotime}}{ a vector of pseudotime values of length \code{n},
#'   containing \code{NA} values for cells not represented by this lineage} }
#'
#' @examples
#' data("slingshot_example")
#' lin <- get_lineages(X, clus.labels, start.clus = 'a')
#' crv <- get_curves(X, clus.labels, lin)
#' plot_curves(X, clus.labels, crv)
#' 
#' @export
#' 

get_curves <- function(X, clus.labels, lineages, shrink = TRUE, extend = 'y', reweight = TRUE, drop.multi = TRUE, thresh = 0.001, maxit = 15, stretch = 2, smoother = 'smooth.spline', shrink.method = 'cosine', ...){
  # CHECKS
  X <- as.matrix(X)
  clus.labels <- as.character(clus.labels)
  shrink <- as.numeric(shrink)
  if(shrink < 0 | shrink > 1){
    stop("shrink must be logical or numeric between 0 and 1")
  }
  if(nrow(X) != length(clus.labels)){
    stop('nrow(X) must equal length(clus.labels)')
  }
  if(is.null(rownames(X))){
    rownames(X) <- paste('cell',seq_len(nrow(X)),sep='-')
  }
  if(is.null(colnames(X))){
    colnames(X) <- paste('dim',seq_len(ncol(X)),sep='-')
  }
  # DEFINE SMOOTHER FUNCTION
  smootherFcn <- switch(smoother, loess = function(lambda, xj, w = NULL, ...){
    loess(xj ~ lambda, weights = w, ...)$fitted
  }, smooth.spline = function(lambda, xj, w = NULL, ..., df = 5){
    fit <- smooth.spline(lambda, xj, w = w, ..., df = df, keep.data = FALSE)
    predict(fit, x = lambda)$y
  })
  
  # remove unclustered cells
  X.original <- X
  clus.labels.original <- clus.labels
  X <- X[clus.labels != -1,]
  clus.labels <- clus.labels[clus.labels != -1]
  # SETUP
  L <- length(grep("lineage",names(lineages))) # number of lineages
  clusters <- unique(clus.labels)
  d <- dim(X); n <- d[1]; p <- d[2]
  nclus <- length(clusters)
  centers <- t(sapply(clusters,function(clID){
    x.sub <- X[clus.labels == clID, ,drop = FALSE]
    return(colMeans(x.sub))
  }))
  rownames(centers) <- clusters
  W <- sapply(seq_len(L),function(l){as.numeric(clus.labels %in% lineages[[l]])}) # weighting matrix
  rownames(W) <- rownames(X); colnames(W) <- names(lineages)[seq_len(L)]
  W.orig <- W
  D <- W; D[,] <- NA
  
  # determine curve hierarchy
  C <- as.matrix(sapply(lineages[seq_len(L)], function(lin) {
    sapply(clusters, function(clID) {
      as.numeric(clID %in% lin)
    })
  }))
  rownames(C) <- clusters
  segmnts <- unique(C[rowSums(C)>1,,drop = FALSE])
  segmnts <- segmnts[order(rowSums(segmnts),decreasing = FALSE),,drop = FALSE]
  avg.order <- list()
  for(i in seq_len(nrow(segmnts))){
    idx <- segmnts[i,] == 1
    avg.order[[i]] <- colnames(segmnts)[idx]
    new.col <- rowMeans(segmnts[,idx, drop = FALSE])
    segmnts <- cbind(segmnts[, !idx, drop = FALSE],new.col)
    colnames(segmnts)[ncol(segmnts)] <- paste('average',i,sep='')
  }
  
  # initial curves are piecewise linear paths through the tree
  pcurves <- list()
  for(l in seq_len(L)){
    idx <- W[,l] > 0
    clus.sub <- clus.labels[idx]
    line.initial <- centers[clusters %in% lineages[[l]], , drop = FALSE]
    line.initial <- line.initial[match(lineages[[l]],rownames(line.initial)),  ,drop = FALSE]
    K <- nrow(line.initial)
    # special case: single-cluster lineage
    if(K == 1){
      pca <- prcomp(X[idx, ,drop = FALSE])
      ctr <- line.initial
      line.initial <- rbind(ctr - 10*pca$sdev[1]*pca$rotation[,1], ctr, ctr + 10*pca$sdev[1]*pca$rotation[,1])
      curve <- .get_lam(X[idx, ,drop = FALSE], s = line.initial, stretch = 9999)
      # do this twice because all points should have projections on all 
      # lineages, but only those points on the lineage should extend it.
      pcurve <- .get_lam(X, s = curve$s[curve$tag,], stretch=0)
      pcurve$lambda <- pcurve$lambda - min(pcurve$lambda, na.rm=TRUE) # force pseudotime to start at 0
      pcurve$w <- W[,l]
      pcurves[[l]] <- pcurve
      D[,l] <- pcurve$dist
      next
    }
    
    if(extend == 'y'){
      curve <- .get_lam(X[idx, ,drop = FALSE], s = line.initial, stretch = 9999)
    }
    if(extend == 'n'){
      curve <- .get_lam(X[idx, ,drop = FALSE], s = line.initial, stretch = 0)
    }
    if(extend == 'pc1'){
      pc1.1 <- prcomp(X[clus.labels == lineages[[l]][1],])
      pc1.1 <- pc1.1$rotation[,1] * pc1.1$sdev[1]^2
      leg1 <- line.initial[2,] - line.initial[1,]
      # pick the direction most "in line with" the first branch
      if(sum(pc1.1*leg1) > 0){ # dot prod < 0 => cos(theta) < 0 => larger angle
        pc1.1 <- -pc1.1 
      }
      pc1.2 <- prcomp(X[clus.labels == lineages[[l]][K],])
      pc1.2 <- pc1.2$rotation[,1] * pc1.2$sdev[1]^2
      leg2 <- line.initial[K-1,] - line.initial[K,]
      if(sum(pc1.2*leg2) > 0){ # dot prod < 0 => cos(theta) < 0 => larger angle
        pc1.2 <- -pc1.2 
      }
      line.initial <- rbind(line.initial[1] + pc1.1, line.initial)
      line.initial <- rbind(line.initial, line.initial[K] + pc1.2)
      curve <- .get_lam(X[idx, ,drop = FALSE], s = line.initial, stretch = 9999)
    }
    
    pcurve <- .get_lam(X, s = curve$s[curve$tag,], stretch=0)
    pcurve$lambda <- pcurve$lambda - min(pcurve$lambda, na.rm=TRUE) # force pseudotime to start at 0
    pcurve$w <- W[,l]
    pcurves[[l]] <- pcurve
    D[,l] <- pcurve$dist
  }

  # track distances between curves and data points to determine convergence
  dist.new <- sum(D[W>0], na.rm=T)
  
  it <- 0
  hasConverged <- FALSE
  while (!hasConverged && it < maxit){
    it <- it + 1
    dist.old <- dist.new
    
    if(reweight){
      W[,] <- t(sapply(seq_len(nrow(W)),function(i){
        ds <- D[i,]
        out <- min(ds)/ds
        return(out)
      }))
      W[W > 1] <- 1
      W[W.orig==0] <- 0
    }
    if(drop.multi){
      Z <- D; Z[,] <- NA
      for(l in seq_len(L)){
        idx <- W[,l] > 0
        Z[idx,l] <- rank(D[idx,l]) / sum(idx)
      }
      W[,] <- t(sapply(seq_len(nrow(W)),function(i){
        zs <- Z[i,]
        out <- W[i,]
        if(max(zs,na.rm=TRUE) > .9 && min(zs,na.rm=TRUE) <= .5){
          out[!is.na(zs) & zs > .9] <- 0
        }
        return(out)
      }))
    }
    
    # predict each dimension as a function of lambda (pseudotime)
    for(l in seq_len(L)){
      pcurve <- pcurves[[l]]
      s <- pcurve$s
      ord <- order(pcurve$lambda)
      for(jj in seq_len(p)){
        s[, jj] <- smootherFcn(pcurve$lambda, X[,jj], w = pcurve$w, ...)[ord]
      }
      new.pcurve <- .get_lam(X, s = s, stretch = stretch)
      new.pcurve$lambda <- new.pcurve$lambda - min(new.pcurve$lambda, na.rm = TRUE)
      new.pcurve$w <- W[,l]
      pcurves[[l]] <- new.pcurve
    }
    D[,] <- sapply(pcurves, function(p){ p$dist })
    
    # shrink together lineages near shared clusters
    if(shrink > 0){
      if(max(rowSums(C)) > 1){
        
        segmnts <- unique(C[rowSums(C)>1,,drop=FALSE])
        segmnts <- segmnts[order(rowSums(segmnts),decreasing = FALSE),,drop = FALSE]
        seg.mix <- segmnts
        avg.lines <- list()
        pct.shrink <- list()
        
        # determine average curves and amount of shrinkage
        for(i in seq_along(avg.order)){
          ns <- avg.order[[i]]
          to.avg <- lapply(ns,function(n){
            if(grepl('lineage',n)){
              l.ind <- as.numeric(gsub('lineage','',n))
              return(pcurves[[l.ind]])
            }
            if(grepl('average',n)){
              a.ind <- as.numeric(gsub('average','',n))
              return(avg.lines[[a.ind]])
            }
          })
          avg <- .avg_curves(to.avg, X, stretch = stretch)
          avg.lines[[i]] <- avg
          common.ind <- rowMeans(sapply(to.avg,function(crv){crv$w > 0})) == 1
          pct.shrink[[i]] <- lapply(to.avg,function(crv){
            .percent_shrinkage(crv, common.ind, method = shrink.method)
          })
        }
        # do the shrinking in reverse order
        for(j in rev(seq_along(avg.lines))){
          ns <- avg.order[[j]]
          avg <- avg.lines[[j]]
          to.shrink <- lapply(ns,function(n){
            if(grepl('lineage',n)){
              l.ind <- as.numeric(gsub('lineage','',n))
              return(pcurves[[l.ind]])
            }
            if(grepl('average',n)){
              a.ind <- as.numeric(gsub('average','',n))
              return(avg.lines[[a.ind]])
            }
          })
          shrunk <- lapply(seq_along(ns),function(jj){
            crv <- to.shrink[[jj]]
            return(.shrink_to_avg(crv, avg, pct.shrink[[j]][[jj]], X, stretch = stretch)) # sorry
          })
          for(jj in seq_along(ns)){
            n <- ns[jj]
            if(grepl('lineage',n)){
              l.ind <- as.numeric(gsub('lineage','',n))
              pcurves[[l.ind]] <- shrunk[[jj]]
            }
            if(grepl('average',n)){
              a.ind <- as.numeric(gsub('average','',n))
              avg.lines[[a.ind]] <- shrunk[[jj]]
            }
          }
        }
        
      }
    }
    D[,] <- sapply(pcurves, function(p){ p$dist })
    
    dist.new <- sum(D[W>0], na.rm=TRUE)
    hasConverged <- (abs((dist.old - dist.new)/dist.old) <= thresh)
  }
  
  if(reweight){
    W[,] <- t(sapply(seq_len(nrow(W)),function(i){
      ds <- D[i,]
      out <- min(ds)/ds
      return(out)
    }))
    W[W > 1] <- 1
    W[W.orig==0] <- 0
  }
  if(drop.multi){
    Z <- D; Z[,] <- NA
    for(l in seq_len(L)){
      idx <- W[,l] > 0
      Z[idx,l] <- rank(D[idx,l]) / sum(idx)
    }
    W[,] <- t(sapply(seq_len(nrow(W)),function(i){
      zs <- Z[i,]
      out <- W[i,]
      if(max(zs,na.rm=TRUE) > .9 && min(zs,na.rm=TRUE) <= .5){
        out[!is.na(zs) & zs > .9] <- 0
      }
      return(out)
    }))
  }
  
  for(l in seq_len(L)){
    pcurve <- pcurves[[l]]
    
    pcurve$pseudotime <- pcurve$lambda
    pcurve$w <- W[,l]
    pcurve$pseudotime[pcurve$w==0] <- NA
    class(pcurve) <- 'principal.curve'
    
    pcurves[[l]] <- pcurve
  }
  names(pcurves) <- paste('curve',1:length(pcurves),sep='')
  return(pcurves)
}













