#########################
### Helper functions for Slingshot
#########################
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
.scaleAB <- function(x,a=0,b=1){
  ((x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE)))*(b-a)+a
}
.avg_curves <- function(pcurves, X, stretch = 2){
  p <- ncol(pcurves[[1]]$s)
  lambdas.all <- lapply(pcurves, function(pcv){pcv$lambda})
  lambdas.all <- unique(unlist(lambdas.all))
  max.shared.lambda <- min(sapply(pcurves, function(pcv){max(pcv$lambda)}))
  lambdas.all <- sort(lambdas.all[lambdas.all <= max.shared.lambda])
  pcurves.dense <- lapply(pcurves,function(pcv){
    sapply(seq_len(p),function(jj){
      interpolated <- approx(pcv$lambda, pcv$s[,jj], xout = lambdas.all)$y
      return(interpolated)
    })
  })
  avg <- sapply(seq_len(p),function(jj){
    dim.all <- sapply(1:length(pcurves.dense),function(i){ pcurves.dense[[i]][,jj] })
    return(rowMeans(dim.all))
  })
  avg.curve <- .get_lam(X, avg, stretch=stretch)
  avg.curve$w <- rowMeans(sapply(pcurves, function(p){ p$w }))
  return(avg.curve)
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
.percent_shrinkage <- function(crv, share.idx, method = 'cosine'){
  pst <- crv$lambda
  if(method %in% eval(formals(density.default)$kernel)){
    dens <- density(0, bw=1, kernel = method)
    surv <- list(x = dens$x, y = (sum(dens$y) - cumsum(dens$y))/sum(dens$y))
    surv$x <- .scaleAB(surv$x, a = min(pst,na.rm=TRUE), b = max(pst[share.idx],na.rm=TRUE))
    pct.l <- approx(surv$x, surv$y, pst, rule = 2)$y
  }
  if(method == 'tricube'){
    tc <- function(x){ ifelse(abs(x) <= 1, (70/81)*((1-abs(x)^3)^3), 0) }
    dens <- list(x = seq(-3,3,length.out = 512))
    dens$y <- tc(dens$x)
    surv <- list(x = dens$x, y = (sum(dens$y) - cumsum(dens$y))/sum(dens$y))
    surv$x <- .scaleAB(surv$x, a = min(pst,na.rm=TRUE), b = max(pst[share.idx],na.rm=TRUE))
    pct.l <- approx(surv$x, surv$y, pst, rule = 2)$y
  }
  if(method == 'density'){
    bw1 <- bw.SJ(pst)
    bw2 <- bw.SJ(pst[share.idx])
    bw <- (bw1 + bw2) / 2
    d2 <- density(pst[share.idx], bw = bw, weights = crv$w[share.idx]/sum(crv$w[share.idx]))
    d1 <- density(pst, bw = bw, weights = crv$w/sum(crv$w))
    scale <- sum(crv$w[share.idx]) / sum(crv$w)
    pct.l <- (approx(d2$x,d2$y,xout = pst, yleft = 0, yright = 0)$y * scale) / approx(d1$x,d1$y,xout = pst, yleft = 0, yright = 0)$y
    pct.l[is.na(pct.l)] <- 0
    pct.l <- .cumMin(pct.l, pst)
  }
  return(pct.l)
}
.shrink_to_avg <- function(pcurve, avg.curve, pct, X, stretch = 2){
  p <- ncol(pcurve$s)
  lam <- pcurve$lambda
  s <- sapply(1:p,function(jj){
    orig.jj <- pcurve$s[,jj]
    avg.jj <- approx(x = avg.curve$lambda, y = avg.curve$s[,jj], xout = lam, rule = 2)$y
    return(avg.jj * pct + orig.jj * (1-pct))
  })
  w <- pcurve$w
  pcurve <- .get_lam(X, s, pcurve$tag, stretch = stretch)
  pcurve$w <- w
  return(pcurve)
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
