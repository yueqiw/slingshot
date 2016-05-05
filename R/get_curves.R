#' @title Construct Smooth Curves for Each Lineage
#' 
#' @description This function takes a reduced data matrix nxp, a vector of cluster identities (optionally including -1's for "unclustered"), and a set of lineages consisting of paths through a forest constructed on the clusters. It constructs smooth curves for each lineage and returns the points along these curves corresponding to the orthogonal projections of each data point, along with lambda (pseudotime) values.
#' 
#' @param X numeric, the nxp matrix of samples in a reduced dimensionality space
#' @param clus.labels character, a vector of length n denoting cluster labels
#' @param lineages list, denotes which lineages each cluster is a part of and contains the matrix defining the forest structure drawn on the clusters by \code{get_lineages}.
#' 
#' @details TODO
#'
#' @return TODO
#'
#' @examples
#' data("toy_data")
#' lineages <- get_lineages(X, clus.labels, start.clus = 'a')
#' get_curves(X, clus.labels, lineages)
#' 
#' @export
#' 

get_curves <- function(X, clus.labels, lineages, shrink = FALSE){
  thresh = 0.0001
  maxit = 100
  stretch = 2
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
  C <- lineages$C
  d <- dim(X)
  n <- d[1]
  p <- d[2]
  L <- ncol(C) # number of lineages
  clusters <- unique(clus.labels)
  nclus <- length(clusters) # number of clusters
  centers <- t(sapply(clusters,function(clID){
    x.sub <- X[clus.labels == clID, ,drop = FALSE]
    return(colMeans(x.sub))
  }))
  rownames(centers) <- clusters
  
  # initial curves are piecewise linear paths through the tree
  pcurves <- lapply(1:L,function(l){
    x.sub <- X[clus.labels %in% lineages[[l]], ,drop = FALSE]
    line.centers <- centers[clusters %in% lineages[[l]], , drop = FALSE]
    line.centers <- line.centers[match(lineages[[l]],rownames(line.centers)),]
    K <- nrow(line.centers)
    s <- .project_points_to_lineage(line.centers, x.sub)
    # adjust - extend lineage past endpoint cluster centers (prevents clumping at ends)
    group1idx <- apply(s,1,function(x){identical(x,line.centers[1,])})
    group2idx <- apply(s,1,function(x){identical(x,line.centers[K,])})
    s[group1idx,] <- .project_points_to_line(line.centers[1,],line.centers[2,], x.sub[group1idx,])
    s[group2idx,] <- .project_points_to_line(line.centers[K-1,],line.centers[K,], x.sub[group2idx,])
    # adjust line.centers to reflect extended lineage
    group1.dist2center <- apply(s[group1idx,],1,function(p){.dist_point_to_segment(line.centers[1,],line.centers[2,],p)})
    group2.dist2center <- apply(s[group2idx,],1,function(p){.dist_point_to_segment(line.centers[K-1,],line.centers[K,],p)})
    line.centers <- rbind(s[group1idx,][which.max(group1.dist2center),], line.centers)
    line.centers <- rbind(line.centers, s[group2idx,][which.max(group2.dist2center),])
    # get total squared distance to lineage
    dist <- sum(.dist_points_to_lineage(line.centers, x.sub)^2)
    lambda <- apply(s,1,function(sp){
      K <- nrow(line.centers)
      dists <- sapply(1:(K-1), function(k){
        .dist_point_to_segment(line.centers[k,],line.centers[k+1,],sp)
      })
      seg <- which.min(dists)
      partial <- rbind(line.centers[1:seg,],sp)
      return(.lineage_length(partial))
    })
    tag <- order(lambda)
    start <- list(s = s, tag = tag, lambda = lambda, dist = dist)
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
      for (jj in 1:p) {
        s[, jj] <- smootherFcn(pcurve$lambda, x.sub[,jj])
      }
      pcurve <- get.lam(x.sub, s = s, stretch = stretch)
      pcurve$lambda <- pcurve$lambda - min(pcurve$lambda, na.rm = TRUE) # start at 0 instead of mean 0
      pcurves[[l]] <- pcurve
    }
    # shrink together lineages near shared clusters
    for(c in 1:nclus){
      clID <- clusters[c]
      if(sum(C[c,]) > 1){
        lines <- which(C[c,]==1)
        avg <- .avg_curves(pcurves[lines])
        pct <- lapply(lines,function(l){
          pcurve <- pcurves[[l]]
          ind <- clus.labels %in% lineages[[l]]
          x <- pcurve$lambda
          #d1 <- density(x)
          #d2 <- density(x[clus.labels[ind] == clID], bw = d1$bw)
          d2 <- density(x[clus.labels[ind] == clID])
          d1 <- density(x, bw = d2$bw)
          scale <- sum(clus.labels[ind] == clID)/length(x)
          pct.l <- sapply(x,function(x){
            (approx(d2$x,d2$y,xout = x, yleft = 0, yright = 0)$y * scale) / approx(d1$x,d1$y,xout = x, yleft = 0, yright = 0)$y
          })
          return(pct.l)
        })
        names(pct) <- lines
        pcurves.shrink <- lapply(lines,function(l){
          pcurve <- pcurves[[l]]
          pct.i <- pct[[which(names(pct) == l)]]
          s <- t(sapply(1:length(pcurve$lambda),function(i){
            lam <- pcurve$lambda[i]
            sapply(1:p,function(jj){
              if(lam %in% avg$lambda){
                avg.jj <- avg$avg[avg$lambda == lam,jj]
                orig.jj <- pcurve$s[i,jj]
                return(avg.jj * pct.i[i] + orig.jj * (1-pct.i[i]))
              }else{
                return(pcurve$s[i,jj])
              }
            })
          }))
          pcurve$s <- s
          return(pcurve)
        })
        pcurves[lines] <- pcurves.shrink
      }
    }
    
    dist.new <- sum(sapply(pcurves, function(pcv){ pcv$dist }))
    hasConverged <- (abs((dist.old - dist.new)/dist.old) <= thresh) || (it >= maxit)
  }
  
  # lines are set, but because shrinking happens second, the points defining
  # the lines are not the projections of the data points yet
  for(l in 1:L){
    x.sub <- X[clus.labels %in% lineages[[l]],]
    line <- pcurves[[l]]$s[pcurves[[l]]$tag,]
    s <- .project_points_to_lineage(line,x.sub)
    rownames(s) <- rownames(x.sub)
    lambda <- apply(s,1,function(sp){
      K <- nrow(line)
      dists <- sapply(1:(K-1), function(k){
        .dist_point_to_segment(line[k,],line[k+1,],sp)
      })
      seg <- which.min(dists)
      if(seg == 1){
        partial <- rbind(line[1,],sp)
      }else{
        partial <- rbind(line[1:(seg-1),],sp)
      }
      return(.lineage_length(partial))
    })
    names(lambda) <- rownames(x.sub)
    tag <- order(lambda)
    names(tag) <- rownames(x.sub)
    pcurves[[l]]$s <- s[tag,]
    pcurves[[l]]$lambda <- lambda[tag]
    pcurves[[l]]$tag <- tag[tag]
  }
  return(pcurves)
}

