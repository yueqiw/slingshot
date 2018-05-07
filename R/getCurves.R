#' @rdname getCurves
#'   
#' @description This function takes a reduced data matrix \code{n} by \code{p}, 
#'   a vector of cluster identities (optionally including \code{-1}'s for 
#'   "unclustered"), and a set of lineages consisting of paths through a forest 
#'   constructed on the clusters. It constructs smooth curves for each lineage 
#'   and returns the points along these curves corresponding to the orthogonal 
#'   projections of each data point, along with corresponding arclength 
#'   (\code{pseudotime} or \code{lambda}) values.
#'   
#' @param sds The \code{SlingshotDataSet} for which to construct simultaneous 
#'   principal curves. This should already have lineages identified by 
#'   \code{\link{getLineages}}.
#' @param shrink logical or numeric between 0 and 1, determines whether and how 
#'   much to shrink branching lineages toward their average prior to the split.
#' @param extend character, how to handle root and leaf clusters of lineages 
#'   when constructing the initial, piece-wise linear curve. Accepted values are
#'   \code{'y'} (default), \code{'n'}, and \code{'pc1'}. See 'Details' for more.
#' @param reweight logical, whether to allow cells shared between lineages to be
#'   reweighted during curve-fitting. If \code{TRUE}, cells shared between 
#'   lineages will be weighted by: distance to nearest curve / distance to
#'   curve.
#' @param drop.multi logical, whether to drop shared cells from lineages which 
#'   do not fit them well. If \code{TRUE}, shared cells with a distance to one 
#'   lineage above the 90th percentile and another below the 50th will be
#'   dropped from the further lineage.
#' @param thresh numeric, determines the convergence criterion. Percent change 
#'   in the total distance from cells to their projections along curves must be 
#'   less than \code{thresh}. Default is \code{0.001}, similar to 
#'   \code{\link{principal.curve}}.
#' @param maxit numeric, maximum number of iterations, see 
#'   \code{\link{principal.curve}}.
#' @param stretch numeric factor by which curves can be extrapolated beyond 
#'   endpoints. Default is \code{2}, see \code{\link{principal.curve}}.
#' @param smoother, choice of scatter plot smoother. Same as 
#'   \code{\link{principal.curve}}, but \code{"lowess"} option is replaced with 
#'   \code{"loess"} for additional flexibility.
#' @param shrink.method character denoting how to determine the appropriate 
#'   amount of shrinkage for a branching lineage. Accepted values are the same
#'   as for \code{kernel} in \code{\link{density}} (default is \code{"cosine"}),
#'   as well as \code{"tricube"} and \code{"density"}. See 'Details' for more.
#' @param allow.breaks logical, determines whether curves that branch very close
#'   to the origin should be allowed to have different starting points.
#' @param ... Additional parameters to pass to scatter plot smoothing function, 
#'   \code{smoother}.
#'   
#' @details When there is only a single lineage, the curve-fitting algorithm is 
#'   nearly identical to that of \code{\link{principal.curve}}. When there are 
#'   multiple lineages and \code{shrink == TRUE}, an additional step is added to
#'   the iterative procedure, forcing curves to be similar in the neighborhood 
#'   of shared points (ie., before they branch).
#'   
#' @details The \code{extend} argument determines how to construct the
#'   piece-wise linear curve used to initiate the recursive algorithm. The
#'   initial curve is always based on the lines between cluster centers and if
#'   \code{extend = 'n'}, this curve will terminate at the center of the
#'   endpoint clusters. Setting \code{extend = 'y'} will allow the first and
#'   last segments to extend beyond the cluster center to the orthogonal
#'   projection of the furthest point. Setting \code{extend = 'pc1'} is similar
#'   to \code{'y'}, but uses the first principal component of the cluster to
#'   determine the direction of the curve beyond the cluster center. These
#'   options typically have little to no impact on the final curve, but can
#'   occasionally help with stability issues.
#'   
#' @details When \code{shink == TRUE}, we compute a shrinkage curve,
#'   \eqn{w_l(t)}, for each lineage, a non-increasing function of pseudotime
#'   that determines how much that lineage should be shrunk toward a shared
#'   average curve. We set \eqn{w_l(0) = 1}, so that the curves will perfectly
#'   overlap the average curve at pseudotime \code{0}. The weighting curve
#'   decreases from \code{1} to \code{0} over the non-outlying pseudotime values
#'   of shared cells (where outliers are defined by the \code{1.5*IQR} rule).
#'   The exact shape of the curve in this region is controlled by
#'   \code{shrink.method}, and can follow the shape of any standard kernel
#'   function's cumulative density curve (or more precisely, survival curve,
#'   since we require a decreasing function). Different choices of
#'   \code{shrink.method} seem to have little impact on the final curves, in 
#'   most cases.
#'   
#' @return An updated \code{\link{SlingshotDataSet}} object containing the 
#'   oringinal input, arguments provided to \code{getCurves} as well as the 
#'   following new elements: \itemize{ \item{curves}{A list of
#'   \code{\link{principal.curve}} objects.} \item{slingParams}{Additional
#'   parameters used for fitting simultaneous principal curves.}}
#'   
#' @references Hastie, T., and Stuetzle, W. (1989). "Principal Curves."
#'   \emph{Journal of the American Statistical Association}, 84:502-516.
#'   
#' @seealso \code{\link{slingshot}}
#'   
#' @examples
#' data("slingshotExample")
#' sds <- getLineages(rd, cl, start.clus = '1')
#' sds <- getCurves(sds)
#' 
#' plot(rd, col = cl, asp = 1)
#' lines(sds, type = 'c', lwd = 3)
#' 
#' @export
#' 
setMethod(f = "getCurves",
          signature = signature(sds = "SlingshotDataSet"),
          definition = function(sds,
                                shrink = TRUE, 
                                extend = 'y', 
                                reweight = TRUE,
                                drop.multi = TRUE, 
                                thresh = 0.001, maxit = 15, stretch = 2, 
                                smoother = 'smooth.spline', 
                                shrink.method = 'cosine',
                                allow.breaks = TRUE, ...){
              
              X <- reducedDim(sds)
              clusterLabels <- clusterLabels(sds)
              lineages <- slingLineages(sds)
              
              sds@slingParams$shrink <- shrink
              sds@slingParams$extend <- extend
              sds@slingParams$reweight <- reweight
              sds@slingParams$drop.multi <- drop.multi
              sds@slingParams$shrink.method <- shrink.method
              
              shrink <- as.numeric(shrink)
              
              # CHECKS
              if(shrink < 0 | shrink > 1){
                  stop("'shrink' parameter must be logical or numeric between 0 
                       and 1")
              }
              
              if(nrow(X)==0){
                  stop('reducedDim has zero rows.')
              }
              if(ncol(X)==0){
                  stop('reducedDim has zero columns.')
              }
              if(nrow(X) != nrow(clusterLabels)){
                  stop('nrow(reducedDim) must equal nrow(clusterLabels).')
              }
              if(any(is.na(X))){
                  stop('reducedDim cannot contain missing values.')
              }
              if(!all(apply(X,2,is.numeric))){
                  stop('reducedDim must only contain numeric values.')
              }
              if (is.null(rownames(X)) &
                  is.null(rownames(clusterLabels))) {
                  rownames(X) <- paste('Cell', seq_len(nrow(X)), sep = '-')
                  rownames(clusterLabels) <-
                      paste('Cell', seq_len(nrow(X)), sep = '-')
              }
              if(is.null(colnames(X))){
                  colnames(X) <- paste('Dim',seq_len(ncol(X)),sep='-')
              }
              if(is.null(colnames(clusterLabels))) {
                  colnames(clusterLabels) <- seq_len(ncol(clusterLabels))
              }
              if(any(rownames(X)=='')){
                  miss.ind <- which(rownames(X) == '')
                  rownames(X)[miss.ind] <- paste('Cell',miss.ind,sep='-')
              }
              if(any(colnames(X)=='')){
                  miss.ind <- which(colnames(X) == '')
                  colnames(X)[miss.ind] <- paste('Dim',miss.ind,sep='-')
              }
              if(is.null(rownames(clusterLabels)) & 
                 !is.null(rownames(X))){
                  rownames(clusterLabels) <- rownames(X)
              }
              if(is.null(rownames(X)) & 
                 !is.null(rownames(clusterLabels))){
                  rownames(X) <- rownames(clusterLabels)
              }
              if(any(rowSums(clusterLabels)>1)){
                  rs <- matrix(rowSums(clusterLabels), 
                               ncol = ncol(clusterLabels))
                  clusterLabels <- clusterLabels / rs
              }
              if(any(colSums(clusterLabels)==0)){
                  clusterLabels <- clusterLabels[, colSums(clusterLabels)!=0, 
                                                 drop = FALSE]
              }
              
              # DEFINE SMOOTHER FUNCTION
              smootherFcn <- switch(smoother, loess = function(lambda, xj, 
                                                               w = NULL, ...){
                  loess(xj ~ lambda, weights = w, ...)$fitted
              }, smooth.spline = function(lambda, xj, w = NULL, ..., df = 5, 
                                          tol = 1e-4){
                  # fit <- smooth.spline(lambda, xj, w = w, ..., df = df, 
                  #                      tol = tol, keep.data = FALSE)
                  fit <- tryCatch({
                      smooth.spline(lambda, xj, w = w, ..., df = df, 
                                    tol = tol, keep.data = FALSE)
                  }, error = function(e){
                      smooth.spline(lambda, xj, w = w, ..., df = df, 
                                    tol = tol, keep.data = FALSE, spar = 1)
                  })
                  predict(fit, x = lambda)$y
              })
              
              # remove unclustered cells
              X.original <- X
              clusterLabels.original <- clusterLabels
              clusterLabels <- clusterLabels[, colnames(clusterLabels) != -1,
                                             drop = FALSE]
              # SETUP
              L <- length(grep("Lineage",names(lineages))) # number of lineages
              clusters <- colnames(clusterLabels)
              d <- dim(X); n <- d[1]; p <- d[2]
              nclus <- length(clusters)
              centers <- t(sapply(clusters,function(clID){
                  w <- clusterLabels[, clID, drop = FALSE]
                  return(apply(X, 2, weighted.mean, w = w))
              }))
              if(p == 1){
                  centers <- t(centers)
              }
              rownames(centers) <- clusters
              W <- sapply(seq_len(L),function(l){
                  rowSums(clusterLabels[, lineages[[l]], drop = FALSE])
              }) # weighting matrix
              rownames(W) <- rownames(X)
              colnames(W) <- names(lineages)[seq_len(L)]
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
              segmnts <- segmnts[order(rowSums(segmnts),decreasing = FALSE), ,
                                 drop = FALSE]
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
                  line.initial <- centers[clusters %in% lineages[[l]], , 
                                          drop = FALSE]
                  line.initial <- line.initial[match(lineages[[l]],
                                                     rownames(line.initial)),  ,
                                               drop = FALSE]
                  K <- nrow(line.initial)
                  # special case: single-cluster lineage
                  if(K == 1){
                      pca <- prcomp(X[idx, ,drop = FALSE])
                      ctr <- line.initial
                      line.initial <- rbind(ctr - 10*pca$sdev[1] * 
                                                pca$rotation[,1], ctr, 
                                            ctr + 10*pca$sdev[1] *
                                                pca$rotation[,1])
                      curve <- .get_lam(X[idx, ,drop = FALSE], s = line.initial,
                                        stretch = 9999)
                      # do this twice because all points should have projections
                      # on all lineages, but only those points on the lineage
                      # should extend it
                      pcurve <- .get_lam(X, s = curve$s[curve$tag,], stretch=0)
                      pcurve$dist <- abs(pcurve$dist) 
                      # ^ force non-negative distances
                      pcurve$lambda <- pcurve$lambda - min(pcurve$lambda, 
                                                           na.rm=TRUE)
                      # ^ force pseudotime to start at 0
                      pcurve$w <- W[,l]
                      pcurves[[l]] <- pcurve
                      D[,l] <- abs(pcurve$dist)
                      next
                  }
                  
                  if(extend == 'y'){
                      curve <- .get_lam(X[idx, ,drop = FALSE], s = line.initial,
                                        stretch = 9999)
                      curve$dist <- abs(curve$dist)
                  }
                  if(extend == 'n'){
                      curve <- .get_lam(X[idx, ,drop = FALSE], s = line.initial,
                                        stretch = 0)
                      curve$dist <- abs(curve$dist)
                  }
                  if(extend == 'pc1'){
                      cl1.idx <- clusterLabels[ , lineages[[l]][1] , 
                                                drop = FALSE] > 0
                      pc1.1 <- prcomp(X[cl1.idx, ])
                      pc1.1 <- pc1.1$rotation[,1] * pc1.1$sdev[1]^2
                      leg1 <- line.initial[2,] - line.initial[1,]
                      # pick the direction most "in line with" the first branch
                      # dot prod < 0 => cos(theta) < 0 => larger angle
                      if(sum(pc1.1*leg1) > 0){ 
                          pc1.1 <- -pc1.1 
                      }
                      cl2.idx <- clusterLabels[ , lineages[[l]][K] , 
                                                drop = FALSE] > 0
                      pc1.2 <- prcomp(X[cl2.idx, ])
                      pc1.2 <- pc1.2$rotation[,1] * pc1.2$sdev[1]^2
                      leg2 <- line.initial[K-1,] - line.initial[K,]
                      # dot prod < 0 => cos(theta) < 0 => larger angle
                      if(sum(pc1.2*leg2) > 0){ 
                          pc1.2 <- -pc1.2 
                      }
                      line.initial <- rbind(line.initial[1] + pc1.1, 
                                            line.initial)
                      line.initial <- rbind(line.initial, 
                                            line.initial[K] + pc1.2)
                      curve <- .get_lam(X[idx, ,drop = FALSE], s = line.initial,
                                        stretch = 9999)
                      curve$dist <- abs(curve$dist)
                  }
                  
                  pcurve <- .get_lam(X, s = curve$s[curve$tag, ,drop=FALSE], 
                                     stretch=0)
                  # force non-negative distances
                  pcurve$dist <- abs(pcurve$dist)
                  # force pseudotime to start at 0
                  pcurve$lambda <- pcurve$lambda - min(pcurve$lambda, 
                                                       na.rm=TRUE) 
                  pcurve$w <- W[,l]
                  pcurves[[l]] <- pcurve
                  D[,l] <- abs(pcurve$dist)
              }
              
              # track distances between curves and data points to determine 
              # convergence
              dist.new <- sum(abs(D[W>0]), na.rm=TRUE)
              
              it <- 0
              hasConverged <- FALSE
              while (!hasConverged && it < maxit){
                  it <- it + 1
                  dist.old <- dist.new
                  
                  if(reweight){
                      W[,] <- t(sapply(seq_len(nrow(W)),function(i){
                          ds <- D[i,]
                          out <- min(ds)/ds
                          out[is.nan(out)] <- 1 # handle 0/0
                          return(out)
                      }))
                      W[W > 1] <- 1
                      W[W < 0] <- 0
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
                          if(max(zs,na.rm=TRUE) > .9 && 
                             min(zs,na.rm=TRUE) <= .5){
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
                          s[, jj] <- smootherFcn(pcurve$lambda, X[,jj], 
                                                 w = pcurve$w, ...)[ord]
                      }
                      new.pcurve <- .get_lam(X, s = s, stretch = stretch)
                      new.pcurve$dist <- abs(new.pcurve$dist)
                      new.pcurve$lambda <- new.pcurve$lambda - 
                          min(new.pcurve$lambda, na.rm = TRUE)
                      new.pcurve$w <- W[,l]
                      pcurves[[l]] <- new.pcurve
                  }
                  D[,] <- sapply(pcurves, function(p){ p$dist })
                  
                  # shrink together lineages near shared clusters
                  if(shrink > 0){
                      if(max(rowSums(C)) > 1){
                          
                          segmnts <- unique(C[rowSums(C)>1,,drop=FALSE])
                          segmnts <- segmnts[order(rowSums(segmnts),
                                                   decreasing = FALSE),
                                             , drop = FALSE]
                          seg.mix <- segmnts
                          avg.lines <- list()
                          pct.shrink <- list()
                          
                          # determine average curves and amount of shrinkage
                          for(i in seq_along(avg.order)){
                              ns <- avg.order[[i]]
                              to.avg <- lapply(ns,function(n){
                                  if(grepl('Lineage',n)){
                                      l.ind <- as.numeric(gsub('Lineage','',n))
                                      return(pcurves[[l.ind]])
                                  }
                                  if(grepl('average',n)){
                                      a.ind <- as.numeric(gsub('average','',n))
                                      return(avg.lines[[a.ind]])
                                  }
                              })
                              avg <- .avg_curves(to.avg, X, stretch = stretch)
                              avg.lines[[i]] <- avg
                              common.ind <- rowMeans(sapply(to.avg,
                                                            function(crv){
                                  crv$w > 0})) == 1
                              pct.shrink[[i]] <- lapply(to.avg,function(crv){
                                  .percent_shrinkage(crv, common.ind, 
                                                     method = shrink.method)
                              })
                              # check for degenerate case (if one curve won't be
                              # shrunk, then the other curve shouldn't be,
                              # either)
                              new.avg.order <- avg.order
                              all.zero <- sapply(pct.shrink[[i]], function(pij){
                                  return(all(pij == 0))
                              })
                              if(any(all.zero)){
                                  if(allow.breaks){
                                      new.avg.order[[i]] <- NULL
                                      message('Curves for', ns[1], 'and', ns[2],
                                              'appear to be going in opposite 
                                              directions. No longer forcing them
                                              to share an initial point. To 
                                              manually override this, set 
                                              allow.breaks = FALSE.')
                                  }
                                  pct.shrink[[i]] <- lapply(pct.shrink[[i]], 
                                                            function(pij){
                                      pij[] <- 0
                                      return(pij)
                                  })
                              }
                          }
                          # do the shrinking in reverse order
                          for(j in rev(seq_along(avg.lines))){
                              ns <- avg.order[[j]]
                              avg <- avg.lines[[j]]
                              to.shrink <- lapply(ns,function(n){
                                  if(grepl('Lineage',n)){
                                      l.ind <- as.numeric(gsub('Lineage','',n))
                                      return(pcurves[[l.ind]])
                                  }
                                  if(grepl('average',n)){
                                      a.ind <- as.numeric(gsub('average','',n))
                                      return(avg.lines[[a.ind]])
                                  }
                              })
                              shrunk <- lapply(seq_along(ns),function(jj){
                                  crv <- to.shrink[[jj]]
                                  return(.shrink_to_avg(crv, avg, 
                                                        pct.shrink[[j]][[jj]], 
                                                        X, stretch = stretch))
                              })
                              for(jj in seq_along(ns)){
                                  n <- ns[jj]
                                  if(grepl('Lineage',n)){
                                      l.ind <- as.numeric(gsub('Lineage','',n))
                                      pcurves[[l.ind]] <- shrunk[[jj]]
                                  }
                                  if(grepl('average',n)){
                                      a.ind <- as.numeric(gsub('average','',n))
                                      avg.lines[[a.ind]] <- shrunk[[jj]]
                                  }
                              }
                          }
                          avg.order <- new.avg.order
                      }
                  }
                  D[,] <- sapply(pcurves, function(p){ p$dist })
                  
                  dist.new <- sum(D[W>0], na.rm=TRUE)
                  hasConverged <- (abs((dist.old - 
                                            dist.new)/dist.old) <= thresh)
              }
              
              if(reweight){
                  W[,] <- t(sapply(seq_len(nrow(W)),function(i){
                      ds <- D[i,]
                      out <- min(ds)/ds
                      out[is.nan(out)] <- 1 # handle 0/0
                      return(out)
                  }))
                  W[W > 1] <- 1
                  W[W < 0] <- 0
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
                  class(pcurves[[l]]) <- 'principal.curve'
              }
              names(pcurves) <- paste('curve',1:length(pcurves),sep='')
              
              for(l in seq_len(L)){
                  pcurve$pseudotime <- pcurve$lambda
                  pcurve$w <- W[,l]
                  pcurve$pseudotime[pcurve$w==0] <- NA
              }
              
              #pseudotime <- sapply(pcurves, function(pc) { pc$lambda })
              #weights <- sapply(seq_len(L), function(l) { W[,l] })
              #rownames(weights) <- rownames(X)
              #colnames(weights) <- names(pcurves)
              #pseudotime[weights == 0] <- NA
              #rownames(pseudotime) <- rownames(X)
              #colnames(pseudotime) <- names(pcurves)
              
              sds@curves <- pcurves
              #sds@pseudotime <- pseudotime
              #sds@curveWeights <- weights
              
              validObject(sds)
              return(sds)
          })

