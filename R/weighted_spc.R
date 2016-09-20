# for testing:
require(slingshot)
data("slingshot_example")
lineages <- get_lineages(X,clus.labels, start.clus = 'a')
thresh = 1e-04; maxit = 100; stretch = 2; shrink = T; extend = "y"

simple_gc <- function (X, clus.labels, lineages, thresh = 1e-04, maxit = 100, stretch = 2, shrink = 0.5, extend = "y"){
  # CHECKS
  shrink <- as.numeric(shrink)
  if (shrink < 0 | shrink > 1) {
    stop("shrink must be logical or numeric between 0 and 1")
  }
  if(length(clus.labels) != nrow(X)){
    stop("cluster vector length must match the number of rows in X")
  }
  # SETUP
  smoother <- "smooth.spline"
  smootherFcn <- function(lambda, xj, w, ..., df = 5) {
    o <- order(lambda)
    lambda <- lambda[o]
    xj <- xj[o]
    fit <- smooth.spline(x = lambda, y = xj, w = w, ..., df = df, keep.data = FALSE)
    predict(fit, x = lambda)$y
  }
  X.original <- X
  clus.labels.original <- clus.labels
  X <- X[clus.labels != -1, ]
  clus.labels <- clus.labels[clus.labels != -1]
  L <- length(grep("lineage", names(lineages)))
  clusters <- unique(clus.labels)
  C <- sapply(lineages[1:L], function(lin) {
    sapply(clusters, function(clID) {
      as.numeric(clID %in% lin)
    })
  })
  rownames(C) <- clusters
  clust.sizes <- table(clus.labels)
  d <- dim(X); n <- d[1]; p <- d[2]
  nclus <- length(clusters)
  centers <- t(sapply(clusters, function(clID) {
    x.sub <- X[clus.labels == clID, , drop = FALSE]
    return(colMeans(x.sub))
  }))
  rownames(centers) <- clusters
  lin.ID <- sapply(seq_len(L), function(l){ clus.labels %in% lineages[[l]] })
  # INITIAL GUESSES
  pcurves <- lapply(seq_len(L), function(l){
    line.initial <- centers[clusters %in% lineages[[l]], , drop = FALSE]
    line.initial <- line.initial[match(lineages[[l]], rownames(line.initial)),]
    K <- nrow(line.initial)
    s <- matrix(NA, nrow=n, ncol=p)
    if (extend == "y") {
      s[lin.ID[,l]] <- .project_points_to_lineage(line.initial, X[lin.ID[,l],])
      group1idx <- sapply(seq_len(n), function(i) {
        identical(s[i,], line.initial[1, ])
      }) & (clus.labels == lineages[[l]][1])
      group2idx <- sapply(seq_len(n), function(i) {
        identical(s[i,], line.initial[K, ])
      }) & (clus.labels == lineages[[l]][K])
      #
      proj1 <- .project_points_to_ray(line.initial[2, ], 
                                      line.initial[1, ], x.sub[group1idx, ])
      proj2 <- .project_points_to_ray(line.initial[K - 
                                                     1, ], line.initial[K, ], x.sub[group2idx, ])
      line.initial <- rbind(proj1[nrow(proj1), ], line.initial)
      line.initial <- rbind(line.initial, proj2[nrow(proj2), 
                                                ])
      s <- .project_points_to_lineage(line.initial, x.sub)
    }
    if (extend == "n") {
      s <- .project_points_to_lineage(line.initial, x.sub)
    }
    if (extend == "pc1") {
      s <- .project_points_to_lineage(line.initial, x.sub)
      group1idx <- apply(s, 1, function(x) {
        identical(x, line.initial[1, ])
      }) & (clus.sub == lineages[[l]][1])
      group2idx <- apply(s, 1, function(x) {
        identical(x, line.initial[K, ])
      }) & (clus.sub == lineages[[l]][K])
      pc1.1 <- prcomp(X[clus.labels == lineages[[l]][1]])$rotation[, 
                                                                   1]
      leg1 <- line.initial[2, ] - line.initial[1, ]
      if (sum(pc1.1 * leg1) > 0) {
        pc1.1 <- -pc1.1
      }
      pc1.2 <- prcomp(X[clus.labels == lineages[[l]][K]])$rotation[, 
                                                                   1]
      leg2 <- line.initial[K - 1, ] - line.initial[K, ]
      if (sum(pc1.2 * leg2) > 0) {
        pc1.2 <- -pc1.2
      }
      proj1 <- .project_points_to_ray(line.initial[1, ], 
                                      line.initial[1, ] + pc1.1, x.sub[group1idx, ])
      proj2 <- .project_points_to_ray(line.initial[K, ], 
                                      line.initial[K, ] + pc1.2, x.sub[group2idx, ])
      line.initial <- rbind(proj1[nrow(proj1), ], line.initial)
      line.initial <- rbind(line.initial, proj2[nrow(proj2), 
                                                ])
      s <- .project_points_to_lineage(line.initial, x.sub)
    }
    dist.vec <- .dist_points_to_lineage(line.initial, x.sub)
    dist <- sum(dist.vec^2)
    lambda <- apply(s, 1, function(sp) {
      K <- nrow(line.initial)
      dists <- sapply(1:(K - 1), function(k) {
        .dist_point_to_segment(line.initial[k, ], line.initial[k + 
                                                                 1, ], sp)
      })
      seg <- which.min(dists)
      partial <- rbind(line.initial[1:seg, ], sp)
      return(.lineage_length(partial))
    })
    tag <- order(lambda)
    start <- list(s = s, lambda = lambda, tag = tag, dist = dist, dist.vec = dist.vec)
    return(start)
  })
  dist.new <- sum(sapply(pcurves, function(pcv) {
    pcv$dist
  }))
  
  dist.mat <- sapply(pcurves, function(pcv){
    out <- pcv$dist.vec[match(rownames(X),names(pcv$dist.vec))]^2
    names(out) <- rownames(X)
    return(out)
  })
  wts <- t(apply(dist.mat,1,function(dis){
    out <- dis*sum(!is.na(dis)) / sum(dis, na.rm = T)
    out[out > 1] <- 1
    return(out^2)
  }))
  
  it <- 0
  hasConverged <- FALSE
  while (!hasConverged && it < maxit) {
    it <- it + 1
    dist.old <- dist.new
    for (l in 1:L) {
      pcurve <- pcurves[[l]]
      s <- pcurve$s
      x.sub <- X[clus.labels %in% lineages[[l]], ]
      w <- wts[clus.labels %in% lineages[[l]],l]
      for (jj in 1:p) {
        s[, jj] <- smootherFcn(pcurve$lambda, x.sub[, jj], w=w)
      }
      new.pcurve <- get.lam(x.sub, s = s, stretch = stretch)
      new.pcurve$lambda <- new.pcurve$lambda - min(new.pcurve$lambda, 
                                                   na.rm = TRUE)
      new.pcurve$dist.vec <- .dist_points_to_lineage(new.pcurve$s[order(new.pcurve$lambda),], x.sub)
      pcurves[[l]] <- new.pcurve
    }
    dist.mat <- sapply(pcurves, function(pcv){
      out <- pcv$dist.vec[match(rownames(X),names(pcv$dist.vec))]^2
      names(out) <- rownames(X)
      return(out)
    })
    wts <- t(apply(dist.mat,1,function(dis){
      out <- dis*sum(!is.na(dis)) / sum(dis, na.rm = T)
      out[out > 1] <- 1
      return(out^2)
    }))
    
    if (shrink > 0) {
      if (max(rowSums(C)) > 1) {
        segmnts <- unique(C[rowSums(C) > 1, ])
        segmnts <- segmnts[order(rowSums(segmnts), decreasing = FALSE), 
                           , drop = FALSE]
        seg.mix <- segmnts
        avg.lines <- list()
        bws <- sapply(seq_len(nrow(segmnts)), function(ii) {
          seg <- segmnts[ii, ]
          sapply(1:L, function(l) {
            if (seg[l] == 1) {
              cls <- rownames(C)[apply(C, 1, function(x) {
                all(x == seg)
              })]
              ind <- clus.labels %in% lineages[[l]]
              x <- pcurves[[l]]$lambda
              out <- tryCatch(bw.SJ(x[clus.labels[ind] %in% 
                                        cls]), error = function(e) 0)
              return(out)
            }
            else {
              return(0)
            }
          })
        })
        bw.med <- median(bws[bws > 0])
        den.lines <- lapply(1:L, function(l) density(pcurves[[l]]$lambda, bw = bw.med, weights = wts[,l]))
        for (i in seq_len(nrow(segmnts))) {
          seg <- segmnts[i, ]
          lines <- which(seg == 1)
          ind <- seg.mix[i, ] == 1
          ns <- colnames(seg.mix)[ind]
          to.avg <- lapply(ns, function(n) {
            if (grepl("lineage", n)) {
              l.ind <- as.numeric(gsub("lineage", "", 
                                       n))
              return(pcurves[[l.ind]])
            }
            if (grepl("average", n)) {
              a.ind <- as.numeric(gsub("average", "", 
                                       n))
              return(avg.lines[[a.ind]])
            }
          })
          avg <- .avg_curves(to.avg)
          avg.lines[[i]] <- avg
          new.col <- rowMeans(seg.mix[, ind, drop = FALSE])
          seg.mix <- cbind(seg.mix[, !ind, drop = FALSE], 
                           new.col)
          colnames(seg.mix)[ncol(seg.mix)] <- paste("average", 
                                                    i, sep = "")
          cls <- rownames(C)[apply(C, 1, function(x) {
            all(x[ind] == 1)
          })]
          pct <- lapply(lines, function(l) {
            pcurve <- pcurves[[l]]
            ind <- clus.labels %in% lineages[[l]]
            pst <- pcurve$lambda
            return(.percent_shrinkage(pst, den.lines[[l]], 
                                      clus.labels[ind] %in% cls, bw.med))
          })
          names(pct) <- lines
          pcurves.shrink <- lapply(lines, function(l) {
            pcurve <- pcurves[[l]]
            pct.i <- pct[[which(names(pct) == l)]] * 
              shrink
            s <- sapply(1:p, function(jj) {
              lam <- pcurve$lambda
              avg.jj <- avg$s[match(lam, avg$lambda), 
                              jj]
              orig.jj <- pcurve$s[, jj]
              out <- avg.jj * pct.i + orig.jj * (1 - 
                                                   pct.i)
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
    dist.new <- sum(sapply(pcurves, function(pcv) {
      pcv$dist
    }))
    hasConverged <- (abs((dist.old - dist.new)/dist.old) <= 
                       thresh) || (it >= maxit)
  }
  for (l in 1:L) {
    pcurve <- pcurves[[l]]
    x.sub <- X[clus.labels %in% lineages[[l]], ]
    new.pcurve <- get.lam(x.sub, s = pcurve$s, tag = pcurve$tag, 
                          stretch = stretch)
    new.pcurve$lambda <- new.pcurve$lambda - min(new.pcurve$lambda, 
                                                 na.rm = TRUE)
    new.pcurve$pseudotime <- new.pcurve$lambda
    names(new.pcurve$pseudotime) <- rownames(x.sub)
    new.pcurve$pseudotime <- new.pcurve$pseudotime[match(rownames(X.original), 
                                                         names(new.pcurve$pseudotime))]
    names(new.pcurve$pseudotime) <- rownames(X.original)
    rownames(new.pcurve$s) <- rownames(x.sub)
    names(new.pcurve$lambda) <- rownames(x.sub)
    ord <- new.pcurve$tag
    new.pcurve$s <- new.pcurve$s[ord, ]
    new.pcurve$lambda <- new.pcurve$lambda[ord]
    new.pcurve$tag <- NULL
    pcurves[[l]] <- new.pcurve
  }
  names(pcurves) <- paste("curve", 1:length(pcurves), sep = "")
  return(pcurves)
}