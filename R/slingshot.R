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
  # set up, remove unclustered cells (-1's)
  X.original <- X
  X <- X[clus.labels != -1,]
  clus.labels <- clus.labels[clus.labels != -1]
  clusters <- unique(clus.labels)
  nclus <- length(clusters)
  
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
#'   contains the matrix defining the forest structure drawn on the clusters by
#'   \code{\link{get_lineages}}.
#' @param thresh (optional) see documentation for \code{\link{principal.curve}}.
#' @param maxit (optional) see documentation for \code{\link{principal.curve}}.
#' @param stretch (optional) see documentation for \code{\link{principal.curve}}.
#' @param shrink logical or numeric between 0 and 1, determines whether and how 
#'   much to shrink branching lineages toward their average prior to the split.
#' @param extend character, how to handle root and leaf clusters of lineages when
#'   constructing the initial, piece-wise linear curve. Accepted values are 'y'
#'   (default), 'n', and 'pc1'. See 'Details' for more. 
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
#' @importFrom princurve get.lam
#' 

get_curves <- function(X, clus.labels, lineages, thresh = 0.0001, maxit = 100, stretch = 2, shrink = .5, extend = 'y'){
  shrink <- as.numeric(shrink)
  if(shrink < 0 | shrink > 1){
    stop("shrink must be logical or numeric between 0 and 1")
  }
  smoother <- "smooth.spline"
  smootherFcn <- function(lambda, xj, ..., df = 5) {
    o <- order(lambda)
    lambda <- lambda[o]
    xj <- xj[o]
    fit <- smooth.spline(x = lambda, y = xj, ..., df = df, keep.data = FALSE)
    predict(fit, x = lambda)$y
  }
  # remove unclustered cells
  X.original <- X
  clus.labels.original <- clus.labels
  X <- X[clus.labels != -1,]
  clus.labels <- clus.labels[clus.labels != -1]
  # setup
  L <- length(grep("lineage",names(lineages))) # number of lineages
  clusters <- unique(clus.labels)
  C <- sapply(lineages[1:L],function(lin){
    sapply(clusters,function(clID){
      as.numeric(clID %in% lin)
    })
  })
  rownames(C) <- clusters
  clust.sizes <- table(clus.labels)
  
  d <- dim(X)
  n <- d[1]
  p <- d[2]
  nclus <- length(clusters) # number of clusters
  centers <- t(sapply(clusters,function(clID){
    x.sub <- X[clus.labels == clID, ,drop = FALSE]
    return(colMeans(x.sub))
  }))
  rownames(centers) <- clusters
  
  # initial curves are piecewise linear paths through the tree
  pcurves <- lapply(1:L,function(l){
    x.sub <- X[clus.labels %in% lineages[[l]], ,drop = FALSE]
    clus.sub <- clus.labels[clus.labels %in% lineages[[l]]]
    line.initial <- centers[clusters %in% lineages[[l]], , drop = FALSE]
    line.initial <- line.initial[match(lineages[[l]],rownames(line.initial)),]
    K <- nrow(line.initial)
    
    if(extend == 'y'){
      s <- .project_points_to_lineage(line.initial, x.sub)
      group1idx <- apply(s, 1, function(x){
        identical(x, line.initial[1, ])
      }) & (clus.sub == lineages[[l]][1])
      group2idx <- apply(s, 1, function(x){
        identical(x, line.initial[K, ])
      }) & (clus.sub == lineages[[l]][K])
      proj1 <- .project_points_to_ray(line.initial[2,], line.initial[1,], x.sub[group1idx, ])
      proj2 <- .project_points_to_ray(line.initial[K-1,], line.initial[K,], x.sub[group2idx, ])
      line.initial <- rbind(proj1[nrow(proj1), ], line.initial)
      line.initial <- rbind(line.initial, proj2[nrow(proj2), ])
      s <- .project_points_to_lineage(line.initial, x.sub)
    }
    if(extend == 'n'){
      s <- .project_points_to_lineage(line.initial, x.sub)
    }
    if(extend == 'pc1'){
      s <- .project_points_to_lineage(line.initial, x.sub)
      group1idx <- apply(s, 1, function(x){
        identical(x, line.initial[1, ])
      }) & (clus.sub == lineages[[l]][1])
      group2idx <- apply(s, 1, function(x){
        identical(x, line.initial[K, ])
      }) & (clus.sub == lineages[[l]][K])
      pc1.1 <- prcomp(X[clus.labels == lineages[[l]][1]])$rotation[,1]
      leg1 <- line.initial[2,] - line.initial[1,]
      # pick the direction most "in line with" the first branch
      if(sum(pc1.1*leg1) > 0){ # dot prod < 0 => cos(theta) < 0 => larger angle
        pc1.1 <- -pc1.1 
      }
      pc1.2 <- prcomp(X[clus.labels == lineages[[l]][K]])$rotation[,1]
      leg2 <- line.initial[K-1,] - line.initial[K,]
      if(sum(pc1.2*leg2) > 0){ # dot prod < 0 => cos(theta) < 0 => larger angle
        pc1.2 <- -pc1.2 
      }
      proj1 <- .project_points_to_ray(line.initial[1,], line.initial[1,]+pc1.1, x.sub[group1idx, ])
      proj2 <- .project_points_to_ray(line.initial[K,], line.initial[K,]+pc1.2, x.sub[group2idx, ])
      line.initial <- rbind(proj1[nrow(proj1), ], line.initial)
      line.initial <- rbind(line.initial, proj2[nrow(proj2), ])
      s <- .project_points_to_lineage(line.initial, x.sub)
    }
    
    # get total squared distance to lineage
    dist <- sum(.dist_points_to_lineage(line.initial, x.sub)^2)
    lambda <- apply(s,1,function(sp){
      K <- nrow(line.initial)
      dists <- sapply(1:(K-1), function(k){
        .dist_point_to_segment(line.initial[k,],line.initial[k+1,],sp)
      })
      seg <- which.min(dists)
      partial <- rbind(line.initial[1:seg,],sp)
      return(.lineage_length(partial))
    })
    tag <- order(lambda)
    #ord <- order(lambda)
    #s <- s[ord,]
    #lambda <- lambda[ord]
    start <- list(s = s, lambda = lambda, tag = tag, dist = dist)
    return(start)
  })
  
  # track distances between curves and data points to determine convergence
  dist.new <- sum(sapply(pcurves, function(pcv){ pcv$dist }))
  
  it <- 0
  hasConverged <- FALSE
  while (!hasConverged && it < maxit){
    it <- it + 1
    dist.old <- dist.new
    # predict each dimension as a function of lambda (pseudotime)
    for(l in 1:L){
      pcurve <- pcurves[[l]]
      s <- pcurve$s
      x.sub <- X[clus.labels %in% lineages[[l]],]
      for(jj in 1:p){
        s[, jj] <- smootherFcn(pcurve$lambda, x.sub[,jj])
      }
      new.pcurve <- get.lam(x.sub, s = s, stretch = stretch)
      new.pcurve$lambda <- new.pcurve$lambda - min(new.pcurve$lambda, na.rm = TRUE) # start at 0 instead of mean 0
      pcurves[[l]] <- new.pcurve
    }
    
    # shrink together lineages near shared clusters
    if(shrink > 0){
      if(max(rowSums(C)) > 1){
        
        segmnts <- unique(C[rowSums(C)>1,])
        segmnts <- segmnts[order(rowSums(segmnts),decreasing = FALSE),,drop = FALSE]
        seg.mix <- segmnts
        avg.lines <- list()
        
        bws <- sapply(seq_len(nrow(segmnts)),function(ii){
          seg <- segmnts[ii,]
          sapply(1:L,function(l){
            if(seg[l]==1){
              cls <- rownames(C)[apply(C,1,function(x){all(x==seg)})]
              ind <- clus.labels %in% lineages[[l]]
              x <- pcurves[[l]]$lambda
              out <- tryCatch(bw.SJ(x[clus.labels[ind] %in% cls]), error=function(e) 0)
              return(out)
            }else{
              return(0)
            }
          })
        })
        bw.med <- median(bws[bws > 0])
        den.lines <- lapply(1:L,function(l) density(pcurves[[l]]$lambda, bw=bw.med))
        
        # Calculate two curves' avg, shrink, repeat.
        for(i in seq_len(nrow(segmnts))){
          seg <- segmnts[i,]
          lines <- which(seg==1)
          ind <- seg.mix[i,] == 1
          ns <- colnames(seg.mix)[ind]
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
          avg <- .avg_curves(to.avg)
          avg.lines[[i]] <- avg
          new.col <- rowMeans(seg.mix[,ind, drop = FALSE])
          seg.mix <- cbind(seg.mix[, !ind, drop = FALSE],new.col)
          colnames(seg.mix)[ncol(seg.mix)] <- paste('average',i,sep='')
          
          cls <- rownames(C)[apply(C,1,function(x){all(x[ind]==1)})]
          pct <- lapply(lines,function(l){
            pcurve <- pcurves[[l]]
            ind <- clus.labels %in% lineages[[l]]
            pst <- pcurve$lambda
            return(.percent_shrinkage(pst, den.lines[[l]], clus.labels[ind] %in% cls, bw.med))
          })
          names(pct) <- lines
          pcurves.shrink <- lapply(lines,function(l){
            pcurve <- pcurves[[l]]
            pct.i <- pct[[which(names(pct) == l)]] * shrink
            s <- sapply(1:p,function(jj){
              lam <- pcurve$lambda
              avg.jj <- avg$s[match(lam,avg$lambda),jj]
              orig.jj <- pcurve$s[,jj]
              out <- avg.jj * pct.i + orig.jj * (1-pct.i)
              out[is.na(out)] <- orig.jj[is.na(out)]
              return(out)
            })
            pcurve$s <- s
            pcurve$tag <- order(pcurve$lambda)
            return(pcurve)
          })
          pcurves[lines] <- pcurves.shrink
        }
        
      }
    }
    
    dist.new <- sum(sapply(pcurves, function(pcv){ pcv$dist }))
    hasConverged <- (abs((dist.old - dist.new)/dist.old) <= thresh) || (it >= maxit)
  }
  # lines are set, but because shrinking happens second, the points defining
  # the lines are not the projections of the data points yet
  for(l in 1:L){
    pcurve <- pcurves[[l]]
    x.sub <- X[clus.labels %in% lineages[[l]],]
    new.pcurve <- get.lam(x.sub, s = pcurve$s, tag = pcurve$tag, stretch = stretch)
    new.pcurve$lambda <- new.pcurve$lambda - min(new.pcurve$lambda, na.rm = TRUE) # start at 0 instead of mean 0
    new.pcurve$pseudotime <- new.pcurve$lambda
    names(new.pcurve$pseudotime) <- rownames(x.sub)
    new.pcurve$pseudotime <- new.pcurve$pseudotime[match(rownames(X.original), names(new.pcurve$pseudotime))]
    names(new.pcurve$pseudotime) <- rownames(X.original)
    rownames(new.pcurve$s) <- rownames(x.sub)
    names(new.pcurve$lambda) <- rownames(x.sub)
    ord <- new.pcurve$tag
    new.pcurve$s <- new.pcurve$s[ord,]
    new.pcurve$lambda <- new.pcurve$lambda[ord]
    new.pcurve$tag <- NULL
    pcurves[[l]] <- new.pcurve
  }
  names(pcurves) <- paste('curve',1:length(pcurves),sep='')
  return(pcurves)
}













