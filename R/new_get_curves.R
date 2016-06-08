
new_get_curves <- function(X, clus.labels, lineages, thresh = 0.0001, maxit = 100, stretch = 2, trace = FALSE, shrink = TRUE){
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
    #ord <- order(lambda)
    #s <- s[ord,]
    #lambda <- lambda[ord]
    start <- list(s = s, lambda = lambda, tag = tag, dist = dist)
    return(start)
  })
  
  #for(i in 1:length(pcurves)){lines(pcurves[[i]]$s[pcurves[[i]]$tag,],lwd=2)}
  
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
        s[, jj] <- smootherFcn(pcurve$lambda, x.sub[,jj]) # !!! ordering problem !!!
      }
      new.pcurve <- get.lam(x.sub, s = s, stretch = stretch)
      new.pcurve$lambda <- new.pcurve$lambda - min(new.pcurve$lambda, na.rm = TRUE) # start at 0 instead of mean 0
      pcurves[[l]] <- new.pcurve
    }
    
    # shrink together lineages near shared clusters
    if(max(rowSums(C)) > 1){
      segmnts <- unique(C[rowSums(C)>1,])
      segmnts <- segmnts[order(rowSums(segmnts),decreasing = TRUE),,drop = FALSE]
      avg.lines <- lapply(1:nrow(segmnts),function(i){
        ind <- which(segmnts[i,] == 1)
        .avg_curves(pcurves[ind])
      })
      bws <- sapply(1:nrow(segmnts),function(ii){
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
      
      for(i in 1:nrow(segmnts)){
        seg <- segmnts[i,]
        lines <- which(seg==1)
        cls <- rownames(C)[apply(C,1,function(x){all(x==seg)})]
        avg <- avg.lines[[i]]
        pct <- lapply(lines,function(l){
          pcurve <- pcurves[[l]]
          ind <- clus.labels %in% lineages[[l]]
          pst <- pcurve$lambda
          return(.percent_shrinkage(pst, den.lines[[l]], clus.labels[ind] %in% cls, bw.med))
        })
        names(pct) <- lines
        pcurves.shrink <- lapply(lines,function(l){
          pcurve <- pcurves[[l]]
          pct.i <- pct[[which(names(pct) == l)]]
          s <- sapply(1:p,function(jj){
            lam <- pcurve$lambda
            avg.jj <- avg$line[match(lam,avg$lambda),jj]
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
    
    dist.new <- sum(sapply(pcurves, function(pcv){ pcv$dist }))
    hasConverged <- (abs((dist.old - dist.new)/dist.old) <= thresh) || (it >= maxit)
    hasConverged
    #for(i in 1:length(pcurves)){lines(pcurves[[i]]$s[pcurves[[i]]$tag,],lwd=2)}
  }
  # lines are set, but because shrinking happens second, the points defining
  # the lines are not the projections of the data points yet
  for(l in 1:L){
    pcurve <- pcurves[[l]]
    x.sub <- X[clus.labels %in% lineages[[l]],]
    line <- pcurve$s[pcurve$tag,]
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
    pcurve$s <- s
    pcurve$lambda <- lambda
    pcurve$tag <- tag
    pcurves[[l]] <- pcurve
  }
  return(pcurves)
}

