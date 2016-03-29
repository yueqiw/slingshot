#########################
### Helper functions for Slingshot
#########################
project_point_to_segment <- function(A,B,p){
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
project_points_to_segment <- function(A,B,pts){
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
  return(final[order(t),])
}
dist_point_to_segment <- function(A,B,p){
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
project_points_to_lineage <- function(lineage,pts){
  n <- nrow(pts)
  K <- nrow(lineage)
  if(K == 2){
    projs <- sapply(1:n, function(i){
      p <- pts[i,]
      return(project_point_to_segment(lineage[1,],lineage[2,], p))
    })
  }else{
    dists <- t(apply(pts,1,function(p){
      sapply(1:(K-1),function(k){
        dist_point_to_segment(lineage[k,],lineage[k+1,],p)
      })
    }))
    projs <- sapply(1:n, function(i){
      k <- which.min(dists[i,])
      p <- pts[i,]
      return(project_point_to_segment(lineage[k,],lineage[k+1,], p))
    })
  }
  return(t(projs))
}
dist_points_to_lineage <- function(lineage,pts){
  d <- apply(pts,1,function(p){
    K <- nrow(lineage)
    min(sapply(1:(K-1),function(k){
      dist_point_to_segment(lineage[k,],lineage[k+1,],p)
    }))
  })
  return(d)
}
lineage_length <- function(lineage){
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
get_connections <- function(clus, forest, parent = NULL){
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
scale01 <- function(x){
  return((x - min(x,na.rm=T))/max(x - min(x,na.rm=T)))
}
avg.curves <- function(pcurves){
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
  ind <- !duplicated(lambdas.all)
  return(list(avg=avg[ind,],lambda=lambdas.all[ind]))
}
