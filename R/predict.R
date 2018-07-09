#' @rdname predict.SlingshotDataSet
#' @title Predict from a Slingshot model
#' 
#' @description Map new observations onto simultaneous principal curves fitted
#'   by \code{slingshot}.
#'   
#' @param object a \code{\link{SlingshotDataSet}} containing simultaneous 
#'   principal curves to use for prediction.
#' @param newdata a matrix or data frame of new points in the same 
#'   reduced-dimensional space as the original input to \code{slingshot} (or
#'   \code{getLineages}).
#'   
#' @details This function is a method for the generic function \code{predict}
#'   with \code{signature(object = "SlingshotDataSet")}. If no \code{newdata}
#'   argument is provided, it will return the original results, given by
#'   \code{object}.
#' 
#' @return A \code{SlingshotDataSet} object based on the input \code{newdata}.
#'   New cells are treated as "unclustered" and the \code{lineages} and
#'   \code{adjacency} slots are intentionally left blank, to distinguish these
#'   results from the original \code{slingshot} output. The \code{curves} slot
#'   represents the projections of each new cell onto the existing curves. As
#'   with standard \code{slingshot} output, the lineage-specific pseudotimes and
#'   assignment weights can be accessed via the functions
#'   \code{\link{slingPseudotime}} and \code{\link{slingCurveWeights}}.
#'   
#' @seealso \code{\link{slingshot}}, \code{\link{SlingshotDataSet}}
#'   
#' @examples
#' data("slingshotExample")
#' sds <- slingshot(rd, cl, start.clus = '1')
#' 
#' x <- cbind(runif(100, min = -5, max = 10), runif(100, min = -4, max = 4))
#' predict(sds, x)
#' 
#' @export
#' 
setMethod(f = "predict",
    signature = signature(object = "SlingshotDataSet"),
    definition = function(object, newdata = NULL){
        # setup
        sds <- object
        if(is.null(newdata)){
            return(sds)
        }else{
            x <- as.matrix(newdata)
        }
        n0 <- nrow(reducedDim(sds))
        curves <- slingCurves(sds)
        L <- length(curves)
        
        # checks
        if(ncol(x) != ncol(reducedDim(sds))){
            stop('New data does not match original number of dimensions.\n',
                'Original: ',ncol(reducedDim(sds)),' columns\n',
                'New data: ',ncol(x),' columns')
        }
        if(nrow(x)==0){
            stop('newdata has zero rows.')
        }
        if(ncol(x)==0){
            stop('newdata has zero columns.')
        }
        if(any(is.na(x))){
            stop('newdata cannot contain missing values.')
        }
        if(!all(apply(x,2,is.numeric))){
            stop('newdata must only contain numeric values.')
        }
        if(is.null(rownames(x))){
            rownames(x) <- paste('newCell', seq_len(nrow(x)), sep = '-')
        }
        if(is.null(colnames(x))){
            colnames(x) <- colnames(reducedDim(sds))
        }
        if(any(colnames(x) != colnames(reducedDim(sds)))){
            colnames(x) <- colnames(reducedDim(sds))
        }
        if(any(rownames(x)=='')){
            miss.ind <- which(rownames(x) == '')
            rownames(x)[miss.ind] <- paste('newCell',miss.ind,sep='-')
        }

        D.orig <- vapply(curves, function(crv){ crv$dist_ind }, rep(0, n0))
        W.orig <- vapply(curves, function(crv){ crv$w }, rep(0, n0))
        
        ordD.orig <- order(D.orig)
        W.prob.orig <- W.orig/rowSums(W.orig)
        W.prob.orig[is.nan(W.prob.orig)] <- 0
        WrnkD.orig <- cumsum(W.prob.orig[ordD.orig]) / sum(W.prob.orig)
        Z.orig <- D.orig
        Z.orig[ordD.orig] <- WrnkD.orig
        
        dists.tofit <- as.numeric(D.orig[W.orig > 0])
        eps <- .5*min(dists.tofit[dists.tofit>0])
        logdists.tofit <- log(dists.tofit + eps)
        quants.tofit <- as.numeric(Z.orig[W.orig > 0])
        
        suppressWarnings({
            fit <- glm(quants.tofit ~ logdists.tofit, family = 'binomial')
        })
        
        crv.proj <- lapply(curves, function(crv){
            project_to_curve(x, crv$s[crv$ord, , drop = FALSE], stretch = 0)
        })
        
        D.proj <- vapply(crv.proj, function(crv){ crv$dist_ind }, 
            rep(0,nrow(x)))
        Z.proj <- D.proj
        Z.proj[,] <- predict(fit, 
            newdata = data.frame(logdists.tofit = as.numeric(log(D.proj+eps))), 
            type = 'response')
        
        Z.prime.proj <- 1-Z.proj^2
        W.proj <- Z.prime.proj / rowMaxs(Z.prime.proj,na.rm = TRUE)
        W.proj[is.nan(W.proj)] <- 1 # handle 0/0
        W.proj[is.na(W.proj)] <- 0
        W.proj[W.proj > 1] <- 1
        W.proj[W.proj < 0] <- 0
        
        # add if z < .5
        idx <- Z.proj < .5
        W.proj[idx] <- 1
        # drop if z > .9 and w < .1
        ridx <- rowMaxs(Z.proj, na.rm = TRUE) > .9 & 
            rowMins(W.proj, na.rm = TRUE) < .1
        W0 <- W.proj[ridx, ]
        Z0 <- Z.proj[ridx, ]
        W0[!is.na(Z0) & Z0 > .9 & W0 < .1] <- 0
        W.proj[ridx, ] <- W0
        
        for(l in seq_len(L)){
            crv.proj[[l]]$w <- W.proj[,l]
        }
        
        cl.mat <- matrix(0, nrow = nrow(x), ncol = ncol(clusterLabels(sds)))
        colnames(cl.mat) <- colnames(clusterLabels(sds))
        rownames(cl.mat) <- rownames(x)
        
        out <- newSlingshotDataSet(reducedDim = x, 
            clusterLabels = cl.mat, 
            curves = crv.proj,
            slingParams = slingParams(sds))
        
        return(out)
    })
