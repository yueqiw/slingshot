#' @describeIn newSlingshotDataSet returns a \code{SlingshotDataSet} object.
#' @export
setMethod(
  f = "newSlingshotDataSet",
  signature = signature("data.frame","ANY"),
  definition = function(reducedDim, clusterLabels, ...){
    RD <- as.matrix(reducedDim)
    rownames(RD) <- rownames(reducedDim)
    newSlingshotDataSet(RD, clusterLabels, ...)
  })
#' @describeIn newSlingshotDataSet returns a \code{SlingshotDataSet} object.
#' @export
setMethod(
  f = "newSlingshotDataSet",
  signature = signature("matrix", "numeric"),
  definition = function(reducedDim, clusterLabels, ...){
    newSlingshotDataSet(reducedDim, as.character(clusterLabels), ...)
  })
#' @describeIn newSlingshotDataSet returns a \code{SlingshotDataSet} object.
#' @export
setMethod(
  f = "newSlingshotDataSet",
  signature = signature("matrix","factor"),
  definition = function(reducedDim, clusterLabels, ...){
    newSlingshotDataSet(reducedDim, as.character(clusterLabels), ...)
  })
#' @describeIn newSlingshotDataSet returns a \code{SlingshotDataSet} object.
#' @export
setMethod(
  f = "newSlingshotDataSet",
  signature = signature("matrix","ANY"),
  definition = function(reducedDim, clusterLabels, ...){
    if(missing(clusterLabels)){
      message('Unclustered data detected.')
      clusterLabels <- rep('1', nrow(reducedDim))
    }
    newSlingshotDataSet(reducedDim, as.character(clusterLabels), ...)
  })
#' @describeIn newSlingshotDataSet returns a \code{SlingshotDataSet} object.
#' @export
setMethod(
  f = "newSlingshotDataSet",
  signature = signature("matrix","character"),
  definition = function(reducedDim, clusterLabels,
                        lineages=list(),
                        connectivity=matrix(NA,0,0),
                        lineageControl=list(),
                        curves=list(),
                        pseudotime=matrix(NA,0,0),
                        curveWeights=matrix(NA,0,0),
                        curveControl=list()
  ){
    if(nrow(reducedDim) != length(clusterLabels)) {
      stop('nrow(reducedDim) must equal length(clusterLabels).')
    }
    # something requires row and column names. Princurve?
    if(is.null(rownames(reducedDim))){
      rownames(reducedDim) <- paste('Cell',seq_len(nrow(reducedDim)),sep='-')
    }
    if(is.null(colnames(reducedDim))){
      colnames(reducedDim) <- paste('Dim',seq_len(ncol(reducedDim)),sep='-')
    }
    if(is.null(names(clusterLabels))){
      names(clusterLabels) <- rownames(reducedDim)
    }
    out <- new("SlingshotDataSet",
               reducedDim=reducedDim,
               clusterLabels=clusterLabels,
               lineages=lineages,
               connectivity=connectivity,
               lineageControl=lineageControl,
               curves=curves,
               pseudotime=pseudotime,
               curveWeights=curveWeights,
               curveControl=curveControl
    )
    validObject(out)
    return(out)
  })

#' @describeIn SlingshotDataSet a short summary of \code{SlingshotDataSet}
#'   object.
#'   
#' @param object a \code{SlingshotDataSet} object.
#' @export
setMethod(
  f = "show",
  signature = "SlingshotDataSet",
  definition = function(object) {
    cat("class:", class(object), "\n\n")
    df <- data.frame(Samples = dim(reducedDim(object))[1], 
                     Dimensions = dim(reducedDim(object))[2])
    print(df, row.names = FALSE)
    cat('\nlineages:', length(lineages(object)), "\n")
    for(i in seq_len(length(lineages(object)))){
      cat('Lineage',i,": ", paste(lineages(object)[[i]],' '), "\n", sep='')
    }
    cat('\ncurves:', length(curves(object)), "\n")
    for(i in seq_len(length(curves(object)))){
      cat('Curve',i,": ", "Length: ", signif(max(curves(object)[[i]]$lambda), 
                                             digits = 5), "\tSamples: ", 
          round(sum(curves(object)[[i]]$w), digits = 2), "\n", sep='')
    }
  }
)
# accessor methods
#' @describeIn SlingshotDataSet returns the matrix representing the reduced
#'   dimensional dataset.
#' @param x a \code{SlingshotDataSet} object.
#' @export
setMethod(
  f = "reducedDim",
  signature = "SlingshotDataSet",
  definition = function(x) x@reducedDim
)
#' @describeIn SlingshotDataSet returns the vector of cluster labels.
#' @export
setMethod(
  f = "clusterLabels",
  signature = "SlingshotDataSet",
  definition = function(x) x@clusterLabels
)
#' @describeIn SlingshotDataSet returns the connectivity matrix between
#'   clusters.
#' @export
setMethod(
  f = "connectivity",
  signature = "SlingshotDataSet",
  definition = function(x) x@connectivity
)
#' @describeIn SlingshotDataSet returns the list of lineages, represented by
#'   ordered sets of clusters.
#' @export
setMethod(
  f = "lineages",
  signature = "SlingshotDataSet",
  definition = function(x) x@lineages
)
#' @describeIn SlingshotDataSet returns the list of additional lineage inference
#'   parameters.
#' @export
setMethod(
  f = "lineageControl",
  signature = "SlingshotDataSet",
  definition = function(x) x@lineageControl
)
#' @describeIn SlingshotDataSet returns the list of smooth lineage curves.
#' @export
setMethod(
  f = "curves",
  signature = "SlingshotDataSet",
  definition = function(x) x@curves
)
#' @describeIn SlingshotDataSet returns the list of additional curve fitting
#'   parameters.
#' @export
setMethod(
  f = "curveControl",
  signature = "SlingshotDataSet",
  definition = function(x) x@curveControl
)

# replacement methods
#' @describeIn SlingshotDataSet Updated object with new reduced dimensional
#'   matrix.
#' @export
setReplaceMethod(
  f = "reducedDim", 
  signature = "SlingshotDataSet",
  definition = function(x, value) initialize(x, reducedDim = value))

#' @describeIn SlingshotDataSet Updated object with new vector of cluster
#'   labels.
#' @export
setReplaceMethod(
  f = "clusterLabels", 
  signature = "SlingshotDataSet",
  definition = function(x, value) initialize(x, clusterLabels = value))
setMethod(f = "[", 
          signature = c("SlingshotDataSet", "ANY", "ANY", "ANY"),
          function(x, i, j, ..., drop=FALSE)
          {
            rd <- reducedDim(x)[i,j]
            cl <- clusterLabels(x)[i]
            initialize(x, reducedDim = rd,
                       clusterLabels  = cl,
                       lineages = list(),
                       connectivity = matrix(NA,0,0),
                       lineageControl = lineageControl(x),
                       curves = list(),
                       pseudotime = matrix(NA,0,0),
                       curveWeights = matrix(NA,0,0),
                       curveControl = curveControl(x))
          })


#' @describeIn SlingshotDataSet returns the matrix of pseudotime values.
#' @export
setMethod(
  f = "pseudotime",
  signature = "SlingshotDataSet",
  definition = function(x, na = TRUE){
    if(length(curves(x))==0){
      stop('No curves detected.')
    }
    pst <- sapply(curves(x), function(pc) {
      t <- pc$lambda
      if(na){
        t[pc$w == 0] <- NA
      }
      return(t)
    })
    rownames(pst) <- rownames(reducedDim(x))
    colnames(pst) <- names(curves(x))
    return(pst)
  }
)

#' @describeIn SlingshotDataSet returns the matrix of cell weights along each
#'   lineage.
#' @export
setMethod(
  f = "curveWeights",
  signature = "SlingshotDataSet",
  definition = function(x){
    if(length(curves(x))==0){
      stop('No curves detected.')
    }
    weights <- sapply(curves(x), function(pc) { pc$w })
    rownames(weights) <- rownames(reducedDim(x))
    colnames(weights) <- names(curves(x))
    return(weights)
  }
)


# internal functions
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
    dim.all <- sapply(1:length(pcurves.dense),function(i){
      pcurves.dense[[i]][,jj]
    })
    return(rowMeans(dim.all))
  })
  avg.curve <- .get_lam(X, avg, stretch=stretch)
  avg.curve$w <- rowMeans(sapply(pcurves, function(p){ p$w }))
  return(avg.curve)
}
# export?
.dist_clusters_full <- function(c1,c2){
  mu1 <- colMeans(c1)
  mu2 <- colMeans(c2)
  diff <- mu1 - mu2
  s1 <- cov(c1)
  s2 <- cov(c2)
  return(as.numeric(t(diff) %*% solve(s1 + s2) %*% diff))
}
# export?
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
    box.vals <- boxplot(pst, plot = FALSE)$stats
    surv$x <- .scaleAB(surv$x, a = box.vals[1], b = box.vals[5])
    pct.l <- approx(surv$x, surv$y, pst, rule = 2)$y
  }
  if(method == 'tricube'){
    tc <- function(x){ ifelse(abs(x) <= 1, (70/81)*((1-abs(x)^3)^3), 0) }
    dens <- list(x = seq(-3,3,length.out = 512))
    dens$y <- tc(dens$x)
    surv <- list(x = dens$x, y = (sum(dens$y) - cumsum(dens$y))/sum(dens$y))
    box.vals <- boxplot(pst, plot = FALSE)$stats
    surv$x <- .scaleAB(surv$x, a = box.vals[1], b = box.vals[5])
    pct.l <- approx(surv$x, surv$y, pst, rule = 2)$y
  }
  if(method == 'density'){
    bw1 <- bw.SJ(pst)
    bw2 <- bw.SJ(pst[share.idx])
    bw <- (bw1 + bw2) / 2
    d2 <- density(pst[share.idx], bw = bw, 
                  weights = crv$w[share.idx]/sum(crv$w[share.idx]))
    d1 <- density(pst, bw = bw, weights = crv$w/sum(crv$w))
    scale <- sum(crv$w[share.idx]) / sum(crv$w)
    pct.l <- (approx(d2$x,d2$y,xout = pst, yleft = 0, yright = 0)$y * scale) / 
      approx(d1$x,d1$y,xout = pst, yleft = 0, yright = 0)$y
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
    avg.jj <- approx(x = avg.curve$lambda, y = avg.curve$s[,jj], xout = lam, 
                     rule = 2)$y
    return(avg.jj * pct + orig.jj * (1-pct))
  })
  w <- pcurve$w
  pcurve <- .get_lam(X, s, pcurve$tag, stretch = stretch)
  pcurve$w <- w
  return(pcurve)
}
# export?
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
                 s, stretch, double(p), double(p), 
                 PACKAGE = "princurve")[c("s","tag", "lambda", "dist")]
  #tt$dist <- sum(tt$dist)
  class(tt) <- "principal.curve"
  tt
}

