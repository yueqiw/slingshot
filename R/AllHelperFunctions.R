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
    definition = function(reducedDim, clusterLabels, ...){
        if(nrow(reducedDim) != length(clusterLabels)) {
            stop('nrow(reducedDim) must equal length(clusterLabels).')
        }
        # something requires row and column names. Princurve?
        if(is.null(rownames(reducedDim))){
            rownames(reducedDim) <- paste('Cell',seq_len(nrow(reducedDim)),
                                          sep='-')
        }
        if(is.null(colnames(reducedDim))){
            colnames(reducedDim) <- paste('Dim',seq_len(ncol(reducedDim)),
                                          sep='-')
        }
        if(is.null(names(clusterLabels))){
            names(clusterLabels) <- rownames(reducedDim)
        }
        clusW <- table(rownames(reducedDim), clusterLabels)
        clusW <- clusW[match(rownames(reducedDim),rownames(clusW)), ,drop=FALSE]
        class(clusW) <- 'matrix'
        return(newSlingshotDataSet(reducedDim, clusW, ...))
    })

#' @describeIn newSlingshotDataSet returns a \code{SlingshotDataSet} object.
#' @export
setMethod(
    f = "newSlingshotDataSet",
    signature = signature("matrix","matrix"),
    definition = function(reducedDim, clusterLabels, 
                          lineages=list(),
                          adjacency=matrix(NA,0,0),
                          curves=list(),
                          slingParams=list()
    ){
        if(nrow(reducedDim) != nrow(clusterLabels)) {
            stop('nrow(reducedDim) must equal nrow(clusterLabels).')
        }
        # something requires row and column names. Princurve?
        if(is.null(rownames(reducedDim))){
            rownames(reducedDim) <- paste('Cell',seq_len(nrow(reducedDim)),
                                          sep='-')
        }
        if(is.null(colnames(reducedDim))){
            colnames(reducedDim) <- paste('Dim',seq_len(ncol(reducedDim)),
                                          sep='-')
        }
        if(is.null(rownames(clusterLabels))){
            rownames(clusterLabels) <- rownames(reducedDim)
        }
        if(is.null(colnames(clusterLabels))){
            colnames(clusterLabels) <- seq_len(ncol(clusterLabels))
        }
        out <- new("SlingshotDataSet",
                   reducedDim = reducedDim,
                   clusterLabels = clusterLabels,
                   lineages = lineages,
                   adjacency = adjacency,
                   curves = curves,
                   slingParams = slingParams
        )
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
        df <- data.frame(Samples = nrow(reducedDim(object)), 
                         Dimensions = ncol(reducedDim(object)))
        print(df, row.names = FALSE)
        cat('\nlineages:', length(slingLineages(object)), "\n")
        for(i in seq_len(length(slingLineages(object)))){
            cat('Lineage',i,": ", paste(slingLineages(object)[[i]],' '), "\n", 
                sep='')
        }
        cat('\ncurves:', length(slingCurves(object)), "\n")
        for(i in seq_len(length(slingCurves(object)))){
            cat('Curve',i,": ", "Length: ", 
                signif(max(slingCurves(object)[[i]]$lambda), digits = 5), 
                "\tSamples: ", round(sum(slingCurves(object)[[i]]$w), 
                                     digits = 2), 
                "\n", sep='')
        }
    }
)
# accessor methods
#' @describeIn SlingshotDataSet returns the matrix representing the reduced
#'   dimensional dataset.
#' @param x a \code{SlingshotDataSet} object.
#' @importFrom SingleCellExperiment reducedDim
#' @export
setMethod(
    f = "reducedDim",
    signature = "SlingshotDataSet",
    definition = function(x) x@reducedDim
)
#' @describeIn SlingshotDataSet returns the matrix representing the reduced
#'   dimensional dataset.
#' @importFrom SingleCellExperiment reducedDims
#' @export
setMethod(
    f = "reducedDims",
    signature = "SlingshotDataSet",
    definition = function(x) x@reducedDim
)
#' @describeIn SlingshotDataSet extracts cluster labels from a
#'   \code{SlingshotDataSet} object.
#' @description Extract cluster labels, either a character vector or matrix of 
#'   weights.
#' @return A matrix of cluster weights for each cell or a vector of cluster
#'   assignments.
#' @examples
#' rd <- matrix(data=rnorm(100), ncol=2)
#' cl <- sample(letters[seq_len(5)], 50, replace = TRUE)
#' sds <- newSlingshotDataSet(rd, cl)
#' clusterLabels(sds)
#' @export
setMethod(
    f = "clusterLabels",
    signature = signature(x="SlingshotDataSet"),
    definition = function(x){
        cl <- x@clusterLabels
        # if(all(cl %in% 0:1) || all(cl %in% c(TRUE,FALSE))){
        #     cl <- colnames(cl)[apply(cl,1,which.max)]
        # }
        return(cl)
    }
)

#' @describeIn slingAdjacency returns the adjacency matrix between
#'   clusters from a \code{\link{SlingshotDataSet}} object.
#' @export
setMethod(
    f = "slingAdjacency",
    signature = "SlingshotDataSet",
    definition = function(x) x@adjacency
)
#' @describeIn slingAdjacency returns the adjacency matrix between
#'   clusters from a \code{\link{SingleCellExperiment}} object.
#' @export
setMethod(
    f = "slingAdjacency",
    signature = "SingleCellExperiment",
    definition = function(x) x@int_metadata$slingshot@adjacency
)

#' @describeIn slingLineages returns the list of lineages, represented by
#'   ordered sets of clusters from a \code{\link{SlingshotDataSet}} object.
#' @export
setMethod(
    f = "slingLineages",
    signature = "SlingshotDataSet",
    definition = function(x) x@lineages
)
#' @describeIn slingLineages returns the list of lineages, represented by 
#'   ordered sets of clusters from a \code{\link{SingleCellExperiment}} object.
#' @export
setMethod(
    f = "slingLineages",
    signature = "SingleCellExperiment",
    definition = function(x) x@int_metadata$slingshot@lineages
)

#' @describeIn slingCurves returns the list of smooth lineage curves from a
#'   \code{\link{SlingshotDataSet}} object.
#' @export
setMethod(
    f = "slingCurves",
    signature = "SlingshotDataSet",
    definition = function(x) x@curves
)
#' @describeIn slingCurves returns the list of smooth lineage curves from a
#'   \code{\link{SingleCellExperiment}} object.
#' @export
setMethod(
    f = "slingCurves",
    signature = "SingleCellExperiment",
    definition = function(x) x@int_metadata$slingshot@curves
)

#' @describeIn slingParams returns the list of additional parameters used
#' by Slingshot from a \code{\link{SlingshotDataSet}} object.
#' @export
setMethod(
    f = "slingParams",
    signature = "SlingshotDataSet",
    definition = function(x) x@slingParams
)
#' @describeIn slingParams returns the list of additional parameters used
#' by Slingshot from a \code{\link{SingleCellExperiment}} object.
#' @export
setMethod(
    f = "slingParams",
    signature = "SingleCellExperiment",
    definition = function(x) x@int_metadata$slingshot@slingParams
)


# replacement methods
# #' @describeIn SlingshotDataSet Updated object with new reduced dimensional
# #'   matrix.
# #' @param value matrix, the new reduced dimensional dataset.
# #' 
# #' @details 
# #' Warning: this will remove any existing lineages or curves from the 
# #' \code{SlingshotDataSet} object.
# #' @importFrom SingleCellExperiment reducedDim<-
# #' @export
# setReplaceMethod(
#     f = "reducedDim",
#     signature = "SlingshotDataSet",
#     definition = function(x, value) initialize(x, reducedDim = value,
#                                          clusterLabels = clusterLabels(x)))
# 
#' # replacement methods
# #' @describeIn SlingshotDataSet Updated object with new reduced dimensional
# #'   matrix.
# #' @param value matrix, the new reduced dimensional dataset.
# #' 
# #' @details 
# #' Warning: this will remove any existing lineages or curves from the 
# #' \code{SlingshotDataSet} object.
# #' @importFrom SingleCellExperiment reducedDims<-
# #' @export
# setReplaceMethod(
#     f = "reducedDims",
#     signature = "SlingshotDataSet",
#     definition = function(x, value) initialize(x, reducedDim = value,
#                                          clusterLabels = clusterLabels(x)))
# 
# #' @describeIn SlingshotDataSet Updated object with new vector of cluster
# #'   labels.
# #' @param value character, the new vector of cluster labels.
# #' 
# #' @details 
# #' Warning: this will remove any existing lineages or curves from the 
# #' \code{SlingshotDataSet} object.
# #' @export
# setReplaceMethod(
#     f = "clusterLabels",
#     signature = c(object = "SlingshotDataSet", value = "ANY"),
#     definition = function(object, value) initialize(x,
#                                                reducedDim = reducedDim(x),
#                                                clusterLabels = value))

#' @describeIn SlingshotDataSet Subset dataset and cluster labels.
#' @param i indices to be applied to rows (cells) of the reduced dimensional
#' matrix and cluster labels.
#' @param j indices to be applied to the columns (dimensions) of the reduced
#' dimensional matrix.
#' @details 
#' Warning: this will remove any existing lineages or curves from the 
#' \code{SlingshotDataSet} object.
#' @export
setMethod(f = "[", 
          signature = c("SlingshotDataSet", "ANY", "ANY", "ANY"),
          function(x, i, j)
          {
              rd <- reducedDim(x)[i,j, drop=FALSE]
              cl <- clusterLabels(x)[i]
              initialize(x, reducedDim = rd,
                         clusterLabels  = cl,
                         lineages = list(),
                         adjacency = matrix(NA,0,0),
                         curves = list(),
                         slingParams = slingParams(x))
          })


#' @describeIn slingPseudotime returns the matrix of pseudotime values from a
#'   \code{\link{SlingshotDataSet}} object.
#' @param na logical. If \code{TRUE} (default), cells that are not assigned to a
#'   lineage will have a pseudotime value of \code{NA}. Otherwise, their
#'   arclength along each curve will be returned.
#' @export
setMethod(
    f = "slingPseudotime",
    signature = "SlingshotDataSet",
    definition = function(x, na = TRUE){
        if(length(slingCurves(x))==0){
            stop('No curves detected.')
        }
        pst <- vapply(slingCurves(x), function(pc) {
            t <- pc$lambda
            if(na){
                t[pc$w == 0] <- NA
            }
            return(t)
        }, rep(0,nrow(reducedDim(x))))
        rownames(pst) <- rownames(reducedDim(x))
        colnames(pst) <- names(slingCurves(x))
        return(pst)
    }
)
#' @describeIn slingPseudotime returns the matrix of pseudotime values from a
#'   \code{\link{SingleCellExperiment}} object.
#' @export
setMethod(
    f = "slingPseudotime",
    signature = "SingleCellExperiment",
    definition = function(x, na = TRUE){
        return(slingPseudotime(x@int_metadata$slingshot))
    }
)

#' @describeIn slingPseudotime returns the matrix of cell weights along each
#'   lineage from a \code{\link{SlingshotDataSet}} object.
#' @export
setMethod(
    f = "slingCurveWeights",
    signature = "SlingshotDataSet",
    definition = function(x){
        if(length(slingCurves(x))==0){
            stop('No curves detected.')
        }
        weights <- vapply(slingCurves(x), function(pc) { pc$w },
            rep(0, nrow(reducedDim(x))))
        rownames(weights) <- rownames(reducedDim(x))
        colnames(weights) <- names(slingCurves(x))
        return(weights)
    }
)
#' @describeIn slingPseudotime returns the matrix of cell weights along each
#'   lineage from a \code{\link{SingleCellExperiment}} object.
#' @export
setMethod(
    f = "slingCurveWeights",
    signature = "SingleCellExperiment",
    definition = function(x){
        return(slingCurveWeights(x@int_metadata$slingshot))
    }
)

#' @rdname SlingshotDataSet
#' @import SingleCellExperiment
#' @export
setMethod(
    f = "SlingshotDataSet",
    signature = "SingleCellExperiment",
    definition = function(data){
        if(! "slingshot" %in% names(data@int_metadata)){
            stop('No slingshot results found.')
        }
        return(data@int_metadata$slingshot)
    }
)



##########################
### Internal functions ###
##########################
#' @import stats
#' @import graphics
`.slingParams<-` <- function(x, value) {
    x@slingParams <- value
    x
}
`.slingCurves<-` <- function(x, value) {
    x@curves <- value
    x
}
# to avoid confusion between the clusterLabels argument and function
.getClusterLabels <- function(x){
    x@clusterLabels
}
.scaleAB <- function(x,a=0,b=1){
    ((x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE)))*(b-a)+a
}
.avg_curves <- function(pcurves, X, stretch = 2){
    n <- nrow(pcurves[[1]]$s)
    p <- ncol(pcurves[[1]]$s)
    lambdas.all <- lapply(pcurves, function(pcv){pcv$lambda})
    lambdas.all <- unique(unlist(lambdas.all))
    max.shared.lambda <- min(vapply(pcurves, function(pcv){max(pcv$lambda)},0))
    lambdas.all <- sort(lambdas.all[lambdas.all <= max.shared.lambda])
    pcurves.dense <- lapply(pcurves,function(pcv){
        vapply(seq_len(p),function(jj){
            interpolated <- approx(pcv$lambda, pcv$s[,jj], xout = lambdas.all)$y
            return(interpolated)
        }, rep(0,length(lambdas.all)))
    })
    avg <- vapply(seq_len(p),function(jj){
        dim.all <- vapply(seq_along(pcurves.dense),function(i){
            pcurves.dense[[i]][,jj]
        }, rep(0,length(lambdas.all)))
        return(rowMeans(dim.all))
    }, rep(0,length(lambdas.all)))
    avg.curve <- .get_lam(X, avg, stretch=stretch)
    avg.curve$w <- rowMeans(vapply(pcurves, function(p){ p$w }, rep(0,n)))
    return(avg.curve)
}
# export?
#' @import matrixStats
.dist_clusters_full <- function(X, w1, w2){
    if(length(w1) != nrow(X) | length(w2) != nrow(X)){
        stop("Reduced dimensional matrix and weights vector contain different
             numbers of points.")
    }
    mu1 <- colWeightedMeans(X, w = w1)
    mu2 <- colWeightedMeans(X, w = w2)
    diff <- mu1 - mu2
    s1 <- cov.wt(X, wt = w1)$cov
    s2 <- cov.wt(X, wt = w2)$cov
    return(as.numeric(t(diff) %*% solve(s1 + s2) %*% diff))
}
# export?
.dist_clusters_diag <- function(X, w1, w2){
    if(length(w1) != nrow(X) | length(w2) != nrow(X)){
        stop("Reduced dimensional matrix and weights vector contain different
             numbers of points.")
    }
    mu1 <- colWeightedMeans(X, w = w1)
    mu2 <- colWeightedMeans(X, w = w2)
    diff <- mu1 - mu2
    if(sum(w1>0)==1){
        s1 <-  diag(ncol(X))
    }else{
        s1 <- diag(diag(cov.wt(X, wt = w1)$cov))
    }
    if(sum(w2>0)==1){
        s2 <-  diag(ncol(X))
    }else{
        s2 <- diag(diag(cov.wt(X, wt = w2)$cov))
    }
    return(as.numeric(t(diff) %*% solve(s1 + s2) %*% diff))
}
.cumMin <- function(x,time){
    vapply(seq_along(x),function(i){ min(x[time <= time[i]]) }, 0)
}
.percent_shrinkage <- function(crv, share.idx, method = 'cosine'){
    pst <- crv$lambda
    if(method %in% eval(formals(density.default)$kernel)){
        dens <- density(0, bw=1, kernel = method)
        surv <- list(x = dens$x, y = (sum(dens$y) - cumsum(dens$y))/sum(dens$y))
        box.vals <- boxplot(pst[share.idx], plot = FALSE)$stats
        surv$x <- .scaleAB(surv$x, a = box.vals[1], b = box.vals[5])
        if(box.vals[1]==box.vals[5]){
            pct.l <- rep(0, length(pst))
        }else{
            pct.l <- approx(surv$x, surv$y, pst, rule = 2)$y
        }
    }
    if(method == 'tricube'){
        tc <- function(x){ ifelse(abs(x) <= 1, (70/81)*((1-abs(x)^3)^3), 0) }
        dens <- list(x = seq(-3,3,length.out = 512))
        dens$y <- tc(dens$x)
        surv <- list(x = dens$x, y = (sum(dens$y) - cumsum(dens$y))/sum(dens$y))
        box.vals <- boxplot(pst[share.idx], plot = FALSE)$stats
        surv$x <- .scaleAB(surv$x, a = box.vals[1], b = box.vals[5])
        if(box.vals[1]==box.vals[5]){
            pct.l <- rep(0, length(pst))
        }else{
            pct.l <- approx(surv$x, surv$y, pst, rule = 2)$y
        }
    }
    if(method == 'density'){
        bw1 <- bw.SJ(pst)
        bw2 <- bw.SJ(pst[share.idx])
        bw <- (bw1 + bw2) / 2
        d2 <- density(pst[share.idx], bw = bw, 
                      weights = crv$w[share.idx]/sum(crv$w[share.idx]))
        d1 <- density(pst, bw = bw, weights = crv$w/sum(crv$w))
        scale <- sum(crv$w[share.idx]) / sum(crv$w)
        pct.l <- (approx(d2$x,d2$y,xout = pst, yleft = 0, 
                         yright = 0)$y * scale) / 
            approx(d1$x,d1$y,xout = pst, yleft = 0, yright = 0)$y
        pct.l[is.na(pct.l)] <- 0
        pct.l <- .cumMin(pct.l, pst)
    }
    return(pct.l)
}
.shrink_to_avg <- function(pcurve, avg.curve, pct, X, stretch = 2){
    n <- nrow(pcurve$s)
    p <- ncol(pcurve$s)
    lam <- pcurve$lambda
    s <- vapply(seq_len(p),function(jj){
        orig.jj <- pcurve$s[,jj]
        avg.jj <- approx(x = avg.curve$lambda, y = avg.curve$s[,jj], xout = lam,
                         rule = 2)$y
        return(avg.jj * pct + orig.jj * (1-pct))
    }, rep(0,n))
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



################
### Datasets ###
################

#' @title Parameters for simulating a single lineage
#' @name splatterParameters
#' @aliases splatterParameters params
#'
#' @description A \code{SplatParams} object containing the parameters
#'   needed to perform the simulation shown in the Slingshot vignette. Values
#'   are designed to produce a single, non-branching trajectory. Code for
#'   producing these parameters can be found in the Slingshot vignette.
#'   
#' @format A \code{\link[splatter:SplatParams]{SplatParams}} object with some
#'   path-related parameters set manually.
#' @source Most parameters were learned by
#'   \code{\link[splatter:splatEstimate]{splatEstimate}} and based on the
#'   \code{\link[HSMMSingleCell:HSMM]{HSMM}} dataset.
"params"


#' @title Bifurcating lineages data
#' @name slingshotExample
#' @aliases slingshotExample rd cl
#'   
#' @description A simulated low-dimensional representation of two bifurcating
#'   lineages (\code{rd}) and a vector of cluster labels generated by $k$-means
#'   with $K = 5$ (\code{cl}).
#'   
#' @format A matrix of coordinates in two dimensions, representing $140$ cells,
#'   and a numeric vector $140$ corresponding cluster labels.
#' @source Simulated data provided with the \code{slingshot} package.
"rd"

#' @rdname slingshotExample
"cl"