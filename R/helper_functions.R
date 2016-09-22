#########################
### Helper functions for Slingshot
#########################
.project_point_to_segment <- function(A,B,p){
  AB <- B-A
  AB_squared <- sum(AB*AB)
  if(AB_squared==0){
    return(A)
  }
  Ap <- p-A
  t <- sum(Ap*AB)/AB_squared
  if(t < 0){
    return(A)
  }
  if(t > 1){
    return(B)
  }
  return(A + t*AB)
}
.project_points_to_segment <- function(A,B,pts,order=TRUE){
  if(class(pts)=='numeric'){
    pts <- matrix(pts,ncol=length(pts))
  }
  n <- nrow(pts)
  bigA <- t(matrix(A,nrow=length(A),ncol=n))
  bigB <- t(matrix(B,nrow=length(A),ncol=n))
  AB <- B-A
  bigAB <- t(matrix(AB,nrow=length(A),ncol=n))
  AB_squared <- sum(AB*AB)
  if(AB_squared==0){
    return(bigA)
  }
  bigAp <- pts-bigA
  t <- diag(bigAp %*% t(bigAB))/AB_squared
  final <- t(sapply(t,function(ti){
    if(ti < 0){
      return(A)
    }
    if(ti > 1){
      return(B)
    }
    return(A + ti*AB)
  }))
  rownames(final) <- rownames(pts)
  colnames(final) <- colnames(pts)
  if(order){
    return(final[order(t),])
  }
  return(final)
}
.project_points_to_line <- function(A,B,pts,order=TRUE){
  if(class(pts)=='numeric'){
    pts <- matrix(pts,ncol=length(pts))
  }
  n <- nrow(pts)
  bigA <- t(matrix(A,nrow=length(A),ncol=n))
  bigB <- t(matrix(B,nrow=length(A),ncol=n))
  AB <- B-A
  bigAB <- t(matrix(AB,nrow=length(A),ncol=n))
  AB_squared <- sum(AB*AB)
  if(AB_squared==0){
    return(bigA)
  }
  bigAp <- pts-bigA
  t <- diag(bigAp %*% t(bigAB))/AB_squared
  final <- t(sapply(t,function(ti){
    A + ti*AB
  }))
  rownames(final) <- rownames(pts)
  colnames(final) <- colnames(pts)
  if(order){
    return(final[order(t),])
  }
  return(final)
}
.project_points_to_ray <- function(A,B,pts,order=FALSE){
  if(class(pts)=='numeric'){
    pts <- matrix(pts,ncol=length(pts))
  }
  n <- nrow(pts)
  bigA <- t(matrix(A,nrow=length(A),ncol=n))
  bigB <- t(matrix(B,nrow=length(A),ncol=n))
  AB <- B-A
  bigAB <- t(matrix(AB,nrow=length(A),ncol=n))
  AB_squared <- sum(AB*AB)
  if(AB_squared==0){
    return(bigA)
  }
  bigAp <- pts-bigA
  t <- diag(bigAp %*% t(bigAB))/AB_squared
  final <- t(sapply(t,function(ti){
    if(ti < 0){
      return(A)
    }
    return(A + ti*AB)
  }))
  rownames(final) <- rownames(pts)
  colnames(final) <- colnames(pts)
  if(order){
    return(final[order(t),])
  }
  return(final)
}
.dist_point_to_segment <- function(A,B,p){
  AB <- B-A
  AB_squared <- sum(AB*AB)
  if(AB_squared==0){
    return(sqrt(sum((A-p)^2)))
  }
  Ap <- p-A
  t <- sum(Ap*AB)/AB_squared
  if(t < 0){
    return(sqrt(sum((A-p)^2)))
  }
  if(t > 1){
    return(sqrt(sum((B-p)^2)))
  }
  q <- (A + t*AB)
  return(sqrt(sum((q-p)^2)))
}
.project_points_to_lineage <- function(lineage,pts, extend.ends=FALSE){
  n <- nrow(pts)
  K <- nrow(lineage)
  if(K == 2){
    projs <- sapply(1:n, function(i){
      p <- pts[i,]
      return(.project_point_to_segment(lineage[1,],lineage[2,], p))
    })
  }else{
    dists <- t(apply(pts,1,function(p){
      sapply(1:(K-1),function(k){
        .dist_point_to_segment(lineage[k,],lineage[k+1,],p)
      })
    }))
    projs <- sapply(1:n, function(i){
      k <- which.min(dists[i,])
      p <- pts[i,]
      return(.project_point_to_segment(lineage[k,],lineage[k+1,], p))
    })
  }
  projs <- t(projs)
  if(extend.ends){
    group1idx <- apply(projs,1,function(x){identical(x,lineage[1,])})
    group2idx <- apply(projs,1,function(x){identical(x,lineage[K,])})
    projs[group1idx,] <- .project_points_to_line(lineage[1,],lineage[2,], pts[group1idx,])
    projs[group2idx,] <- .project_points_to_line(lineage[K-1,],lineage[K,], pts[group2idx,])
  }
  return(projs)
}
.dist_points_to_lineage <- function(lineage,pts){
  d <- apply(pts,1,function(p){
    K <- nrow(lineage)
    min(sapply(1:(K-1),function(k){
      .dist_point_to_segment(lineage[k,],lineage[k+1,],p)
    }))
  })
  return(d)
}
.lineage_length <- function(lineage){
  if(class(lineage)=="numeric"){
    return(0)
  }
  K <- nrow(lineage)
  d <- 0
  for(i in 1:(K-1)){
    d <- d + sqrt(sum((lineage[i,]-lineage[i+1,])^2))
  }
  return(d)
}
.get_connections <- function(clus, forest, parent = NULL){
  children.idx <- forest[,clus] == 1
  children <- rownames(forest)[children.idx]
  if(is.null(parent)){
    out <- clus
    for(child in children){
      out <- c(out, Recall(child, forest, clus))
    }
  }else{
    children <- children[children != parent]
    out <- clus
    for(child in children){
      out <- c(out, Recall(child, forest, clus))
    }
  }
  return(out)
}
.scale01 <- function(x){
  return((x - min(x,na.rm=T))/max(x - min(x,na.rm=T)))
}
.avg_curves <- function(pcurves){
  p <- ncol(pcurves[[1]]$s)
  lambdas.all <- lapply(pcurves, function(pcv){pcv$lambda})
  lambdas.all <- unique(unlist(lambdas.all))
  max.shared.lambda <- min(sapply(pcurves, function(pcv){max(pcv$lambda)}))
  lambdas.all <- sort(lambdas.all[lambdas.all <= max.shared.lambda])
  pcurves.dense <- lapply(pcurves,function(pcv){
    sapply(1:p,function(jj){
      interpolated <- approx(pcv$lambda, pcv$s[,jj], xout = lambdas.all)$y
      return(interpolated)
    })
  })
  avg <- sapply(1:p,function(jj){
    dim.all <- sapply(1:length(pcurves.dense),function(i){ pcurves.dense[[i]][,jj] })
    return(rowMeans(dim.all))
  })
  return(list(s=avg,lambda=lambdas.all))
}
.dist_clusters_full <- function(c1,c2){
  mu1 <- colMeans(c1)
  mu2 <- colMeans(c2)
  diff <- mu1 - mu2
  s1 <- cov(c1)
  s2 <- cov(c2)
  return(as.numeric(t(diff) %*% solve(s1 + s2) %*% diff))
}
.dist_clusters_diag <- function(c1,c2){
  mu1 <- colMeans(c1)
  mu2 <- colMeans(c2)
  diff <- mu1 - mu2
  if(nrow(c1)==1){
    s1 <-  diag(ncol(c1))
  }else{
    s1 <- diag(diag(cov(c1)))
  }
  if(nrow(c2)==1){
    s2 <-  diag(ncol(c2))
  }else{
    s2 <- diag(diag(cov(c2)))
  }
  return(as.numeric(t(diff) %*% solve(s1 + s2) %*% diff))
}
.cumMin <- function(x,time){
  sapply(seq_along(x),function(i){ min(x[time <= time[i]]) })
}
.percent_shrinkage <- function(pst, lineage.density, share.idx, bw){
  #d1 <- density(pst)
  #d2 <- density(pst[share.idx], bw = bw.med)
  d2 <- density(pst[share.idx], bw = bw)
  d1 <- lineage.density
  scale <- mean(share.idx)
  pct.l <- (approx(d2$x,d2$y,xout = pst, yleft = 0, yright = 0)$y * scale) / approx(d1$x,d1$y,xout = pst, yleft = 0, yright = 0)$y
  pct.l[is.na(pct.l)] <- 0
  pct.l <- .cumMin(pct.l, pst)
  return(pct.l)
}
.shrink_to_avg <- function(pcurve, avg.curve, pct){
  lam <- pcurve$lambda
  avg.curve$avg <- avg.curve$avg[avg.curve$lambda %in% lam,]
  avg.curve$lambda <- avg.curve$lambda[avg.curve$lambda %in% lam]
  s <- sapply(1:p,function(jj){
    avg.jj <- avg.curve$avg[,jj]
    orig.jj <- pcurve$s[,jj]
    return(avg.jj * pct + orig.jj * (1-pct))
  })
  pcurve$s <- s
  return(pcurve)
}
.sq_segment_lengths <- function(from, to){
  if(any(dim(from) != dim(to))){
    stop('input matrices must have same dimensions')
  }
  sqdists <- sapply(1:nrow(from),function(i){
    sum((from[i,]-to[i,])^2)
  })
  return(sqdists)
}
.get_lam <- function(x, s, tag, stretch = 2){
  storage.mode(x) <- "double"
  storage.mode(s) <- "double"
  storage.mode(stretch) <- "double"
  if (!missing(tag)) 
    s <- s[tag, ]
  np <- dim(x)
  if (length(np) != 2) 
    stop("get.lam needs a matrix input")
  n <- np[1]
  p <- np[2]
  tt <- .Fortran("getlam", n, p, x, s = x, lambda = double(n), 
                 tag = integer(n), dist = double(n), as.integer(nrow(s)), 
                 s, stretch, double(p), double(p), PACKAGE = "princurve")[c("s", 
                                                                            "tag", "lambda", "dist")]
  #tt$dist <- sum(tt$dist)
  class(tt) <- "principal.curve"
  tt
}
