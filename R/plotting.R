#' @title Plot Slingshot output 
#' @name plot-SlingshotDataSet
#' @aliases plot-SlingshotDataSet plot,SlingshotDataSet,ANY-method
#'   
#' @description Tools for visualizing lineages inferred by \code{slingshot}.
#'   
#' @param x a \code{SlingshotDataSet} with results to be plotted.
#' @param type character, the type of output to be plotted, can be one of 
#'   \code{"lineages"}, \code{"curves"}, or \code{"both"} (by partial matching),
#'   see Details for more.
#' @param show.constraints logical, whether or not the user-specified initial 
#'   and terminal clusters should be specially denoted by green and red dots,
#'   respectively.
#' @param add logical, indicates whether the output should be added to an
#'   existing plot.
#' @param dims numeric, which dimensions to plot (default is \code{1:2}).
#' @param asp numeric, the y/x aspect ratio, see \code{\link{plot.window}}.
#' @param cex numeric, amount by which points should be magnified, see
#'   \code{\link{par}}.
#' @param lwd numeric, the line width, see \code{\link{par}}.
#' @param col character or numeric, color(s) for lines, see \code{\link{par}}.
#' @param ... additional parameters to be passed to \code{\link{lines}}.
#'   
#' @details If \code{type == 'lineages'}, straight line connectors between
#'   cluster centers will be plotted. If \code{type == 'curves'}, simultaneous
#'   principal curves will be plotted.
#'   
#' @details When \code{type} is not specified, the function will first check the
#'   \code{curves} slot and plot the curves, if present. Otherwise,
#'   \code{lineages} will be plotted, if present.
#'   
#' @return returns \code{NULL}.
#'   
#' @examples
#' data("slingshotExample")
#' sds <- slingshot(rd, cl, start.clus = "1")
#' plot(sds, type = 'b')
#' 
#' # add to existing plot
#' plot(rd, col = 'grey50')
#' lines(sds, lwd = 3)
#' 
#' @import graphics
#' @import grDevices
#' @export
setMethod(
    f = "plot",
    signature = signature(x = "SlingshotDataSet"),
    definition = function(x, type = NULL,
                          show.constraints = FALSE,
                          add = FALSE,
                          dims = seq_len(2),
                          asp = 1,
                          cex = 2,
                          lwd = 2,
                          col = 1,
                          ...) {
        col <- rep(col, length(slingLineages(x)))
        curves <- FALSE
        lineages <- FALSE
        if(is.null(type)){
            if(length(slingCurves(x)) > 0){
                type <- 'curves'
            }else if(length(slingLineages(x)) > 0){
                type <- 'lineages'
            }else{
                stop('No lineages or curves detected.')
            }
        }else{
            type <- c('curves','lineages','both')[pmatch(type,
                c('curves','lineages','both'))]
            if(is.na(type)){
                stop('Unrecognized type argument.')
            }
        }
        
        if(type %in% c('lineages','both')){
            lineages <- TRUE
        }
        if(type %in% c('curves','both')){
            curves <- TRUE
        }
        
        if(lineages & (length(slingLineages(x))==0)){
            stop('No lineages detected.')
        }
        if(curves & (length(slingCurves(x))==0)){
            stop('No curves detected.')
        }
        
        if(lineages){
            X <- reducedDim(x)
            clusterLabels <- clusterLabels(x)
            connectivity <- slingAdjacency(x)
            clusters <- rownames(connectivity)
            nclus <- nrow(connectivity)
            centers <- t(vapply(clusters,function(clID){
                w <- clusterLabels[,clID]
                return(apply(X, 2, weighted.mean, w = w))
            }, rep(0,ncol(X))))
            rownames(centers) <- clusters
            X <- X[rowSums(clusterLabels) > 0, , drop = FALSE]
            clusterLabels <- clusterLabels[rowSums(clusterLabels) > 0, , 
                                           drop = FALSE]
            linC <- slingParams(x)
        }
        
        if(!add){
            xs <- NULL
            ys <- NULL
            if(lineages){
                xs <- c(xs, centers[,dims[1]])
                ys <- c(ys, centers[,dims[2]])
            }
            if(curves){
                xs <- c(xs, as.numeric(vapply(slingCurves(x), 
                    function(c){ c$s[,dims[1]] }, rep(0,nrow(reducedDim(x))))))
                ys <- c(ys, as.numeric(vapply(slingCurves(x), 
                    function(c){ c$s[,dims[2]] }, rep(0,nrow(reducedDim(x))))))
            }
            plot(x = NULL, y = NULL, asp = asp,
                 xlim = range(xs), ylim = range(ys),
                 xlab = colnames(reducedDim(x))[dims[1]],
                 ylab = colnames(reducedDim(x))[dims[2]])
        }
        
        if(lineages){
            for(i in seq_len(nclus-1)){
                for(j in seq(i+1,nclus)){
                    if(connectivity[i,j]==1){
                        lines(centers[c(i,j),], lwd = lwd, col = col, ...)
                    }
                }
            }
            points(centers, cex = cex, pch = 16, col = col)
            if(show.constraints){
                if(any(linC$start.given)){
                    points(centers[clusters %in% 
                                       linC$start.clus[linC$start.given],,
                                   drop=FALSE], cex = cex / 2, 
                           col = 'green3', pch = 16)
                }
                if(any(linC$end.given)){
                    points(centers[clusters %in% linC$end.clus[linC$end.given],,
                                   drop=FALSE], cex = cex / 2, 
                           col = 'red2', pch = 16)
                }
            }
            
        }
        if(curves){
            for(ii in seq_along(slingCurves(x))){
                c <- slingCurves(x)[[ii]]
                lines(c$s[c$ord,dims], lwd = lwd, col = col[ii], ...)
            }
        }
        invisible(NULL)
    }
)

#' @rdname plot-SlingshotDataSet
#' @export
setMethod(
    f = "lines",
    signature = "SlingshotDataSet",
    definition = function(x,
                          type = NULL,
                          dims = seq_len(2),
                          ...) {
        plot(x, type = type, add = TRUE, dims = dims, ...)
        invisible(NULL)
    }
)

## Individual gene plots
#' @rdname plotGenePseudotime
#' @title Plot Gene Expression over Pseudotime
#'
#' @param gene the gene to be plotted. If \code{exprs} is provided, this may be
#'   either the gene name or its row index in \code{exprs}. Otherwise, this is
#'   assumed to be a vector of scaled expression values.
#' @param exprs the genes-by-samples matrix of scaled expression values (log
#'   counts or normalized log counts).
#' @param loess logical, whether to include a loess fit in each plot (default is
#' \code{TRUE}).
#' @param loessCI logical, whether to include a confidence band around the loess
#' curve (default is \code{TRUE}).
#' @param ... additional parameters to be passed to \code{\link{plot}}.
#'   
#' @return returns \code{NULL}.
#'   
#' @examples
#' data("slingshotExample")
#' sds <- slingshot(rd, cl, start.clus = "1")
#' ex <- matrix(c(rchisq(100,1),rchisq(20,3),rchisq(20,6)),nrow=1)
#' rownames(ex) <- 'Gene-1'
#' plotGenePseudotime(sds, 'Gene-1', ex)
#' 
#' @export
setMethod(
    f = "plotGenePseudotime",
    signature = signature(data = "SlingshotDataSet"),
    definition = function(data, gene, exprs,
        loess = TRUE, loessCI = TRUE, ...) {
        if(length(gene) > 1 & is.numeric(gene)){
            y <- gene
        }
        if(length(gene) == 1){
            y <- exprs[gene, ,drop=FALSE][1,]
        }
        pst <- slingPseudotime(data)
        w <- slingCurveWeights(data)
        L <- length(slingLineages(data))
        
        par(mfrow = c(L,1))
        for(l in seq_len(L)){
            plot(pst[,l], y, xlab = 'Pseudotime', ylab = 'Expression', 
                main=paste(gene, ', Lineage ',l, sep=''), ...)
            if(loess | loessCI){
                l <- loess(y ~ pst[,l], weights = w[,l])
            }
            if(loessCI){
                pl <- predict(l, se=TRUE)
                polygon(c(l$x[order(l$x)],rev(l$x[order(l$x)])), 
                    c((pl$fit+qt(0.975,pl$df)*pl$se)[order(l$x)], 
                        rev((pl$fit-qt(0.975,pl$df)*pl$se)[order(l$x)])),
                    border = NA, col = rgb(0,0,0,.3))
            }
            if(loess){
                lines(l$x[order(l$x)], l$fitted[order(l$x)], lwd=2)
            }
        }
        par(mfrow = c(1,1))
        invisible(NULL)
    }
)

#' @rdname plotGenePseudotime
#' @importFrom SummarizedExperiment assays
#' @export
setMethod(
    f = "plotGenePseudotime",
    signature = signature(data = "SingleCellExperiment"),
    definition = function(data, gene, exprs,
        loess = TRUE, loessCI = TRUE, ...) {
        if(missing(exprs)){
            EX <- assays(data)[[1]]
        }else{
            if(length(exprs)==1){
                EX <- assays(data)[[exprs]]
            }else{
                EX <- exprs
            }
        }
        plotGenePseudotime(SlingshotDataSet(data), gene, exprs = EX,
            loess = loess, loessCI = loessCI, ...)
    }
)


#' @name plot3d-SlingshotDataSet
#' @title Plot Slingshot output in 3D
#' 
#' @description Tools for visualizing lineages inferred by \code{slingshot}.
#'   
#' @param x a \code{SlingshotDataSet} with results to be plotted.
#' @param type character, the type of output to be plotted, can be one of 
#'   \code{"lineages"}, \code{curves}, or \code{both} (by partial matching), see
#'   Details for more.
#' @param add logical, indicates whether the output should be added to an
#'   existing plot.
#' @param dims numeric, which dimensions to plot (default is \code{1:3}).
#' @param aspect either a logical indicating whether to adjust the aspect ratio 
#'   or a new ratio, see \code{\link[rgl:plot3d]{plot3d}}.
#' @param ... additional parameters to be passed to \code{lines3d}.
#'   
#' @details If \code{type == 'lineages'}, straight line connectors between
#'   cluster centers will be plotted. If \code{type == 'curves'}, simultaneous
#'   princiapl curves will be plotted.
#'   
#' @details When \code{type} is not specified, the function will first check the
#'   \code{curves} slot and plot the curves, if present. Otherwise,
#'   \code{lineages} will be plotted, if present.
#'   
#' @return returns \code{NULL}.
#'   
#' @examples
#' \dontrun{
#' library(rgl)
#' data("slingshotExample")
#' rd <- cbind(rd, rnorm(nrow(rd)))
#' sds <- slingshot(rd, cl, start.clus = "1")
#' plot3d(sds, type = 'b')
#' 
#' # add to existing plot
#' plot3d(rd, col = 'grey50', aspect = 'iso')
#' plot3d(sds, lwd = 3, add = TRUE)
#' }
#' @importFrom rgl plot3d
#' @export
plot3d.SlingshotDataSet <- function(x,
                                    type = NULL,
                                    add = FALSE,
                                    dims = seq_len(3),
                                    aspect = 'iso',
                                    ...){
    if (!requireNamespace("rgl", quietly = TRUE)) {
        stop("Package 'rgl' is required for 3D plotting.",
             call. = FALSE)
    }
    n <- nrow(reducedDim(x))
    curves <- FALSE
    lineages <- FALSE
    if(is.null(type)){
        if(length(slingCurves(x)) > 0){
            type <- 'curves'
        }else if(length(slingLineages(x)) > 0){
            type <- 'lineages'
        }else{
            stop('No lineages or curves detected.')
        }
    }else{
        type <- c('curves','lineages','both')[pmatch(type,c('curves','lineages',
                                                            'both'))]
        if(is.na(type)){
            stop('Unrecognized type argument.')
        }
    }
    
    if(type %in% c('lineages','both')){
        lineages <- TRUE
    }
    if(type %in% c('curves','both')){
        curves <- TRUE
    }
    
    if(lineages & (length(slingLineages(x))==0)){
        stop('No lineages detected.')
    }
    if(curves & (length(slingCurves(x))==0)){
        stop('No curves detected.')
    }
    
    if(lineages){
        X <- reducedDim(x)
        clusterLabels <- clusterLabels(x)
        connectivity <- slingAdjacency(x)
        clusters <- rownames(connectivity)
        nclus <- nrow(connectivity)
        centers <- t(vapply(clusters,function(clID){
            w <- clusterLabels[,clID]
            return(apply(X, 2, weighted.mean, w = w))
        }, rep(0,ncol(X))))
        rownames(centers) <- clusters
        X <- X[rowSums(clusterLabels) > 0, , drop = FALSE]
        clusterLabels <- clusterLabels[rowSums(clusterLabels) > 0, , 
                                       drop = FALSE]
    }
    
    if(!add){
        xs <- NULL
        ys <- NULL
        zs <- NULL
        if(lineages){
            xs <- c(xs, centers[,dims[1]])
            ys <- c(ys, centers[,dims[2]])
            zs <- c(zs, centers[,dims[3]])
        }
        if(curves){
            xs <- c(xs, as.numeric(vapply(slingCurves(x), function(c){ 
                c$s[,dims[1]] }, rep(0,n))))
            ys <- c(ys, as.numeric(vapply(slingCurves(x), function(c){ 
                c$s[,dims[2]] }, rep(0,n))))
            zs <- c(zs, as.numeric(vapply(slingCurves(x), function(c){ 
                c$s[,dims[3]] }, rep(0,n))))
        }
        rgl::plot3d(x = NULL, y = NULL, z = NULL, aspect = aspect,
                    xlim = range(xs), ylim = range(ys), zlim = range(zs),
                    xlab = colnames(reducedDim(x))[dims[1]],
                    ylab = colnames(reducedDim(x))[dims[2]],
                    zlab = colnames(reducedDim(x))[dims[3]])
    }
    
    if(lineages){
        for(i in seq_len(nclus-1)){
            for(j in seq(i+1,nclus)){
                if(connectivity[i,j]==1){
                    rgl::lines3d(x = centers[c(i,j),dims[1]], 
                                 y = centers[c(i,j),dims[2]],
                                 z = centers[c(i,j),dims[3]], ...)
                }
            }
        }
    }
    if(curves){
        for(c in slingCurves(x)){ rgl::lines3d(c$s[c$ord,dims], ...) }
    }
    invisible(NULL)
}

# #' @rdname SlingshotDataSet-plot3d
# #' @export
# setMethod(
#   f = "plot3d",
#   signature = "SlingshotDataSet",
#   definition = function(x,
#                         type = NULL,
#                         add = FALSE,
#                         dims = 1:3,
#                         aspect = 'iso',
#                         ...) {
#     plot3d.SlingshotDataSet(x, type = type, add = add, dims = dims, 
#                             aspect = aspect, ...)
#   }
# )

# #' @rdname SlingshotDataSet-plot3d
# #' @examples
# #' %\dontrun{
# #' data("slingshotExample")
# #' rd <- cbind(rd, rnorm(nrow(rd)))
# #' sds <- slingshot(rd, cl, start.clus = "1")
# #' plot3d(rd, col = 'grey50', aspect = 'iso')
# #' lines3d(sds, lwd = 3)
# #' }
# #' @export
# lines3d.SlingshotDataSet <- function(x,
#                                      type = NULL,
#                                      dims = 1:3,
#                                      ...) {
#   plot3d(x, type = type, add = TRUE, dims = dims, ...)
#   invisible(NULL)
# }
# #' @rdname SlingshotDataSet-plot3d
# #' @export
# setMethod(
#   f = "lines3d",
#   signature = "SlingshotDataSet",
#   definition = function(x,
#                         type = NULL,
#                         dims = 1:3,
#                         ...) {
#     lines3d.SlingshotDataSet(x, type = type, add = TRUE, dims = dims, ...)
#   }
# )


#' @title Pairs plot of Slingshot output
#' @name pairs-SlingshotDataSet
#'   
#' @description A tool for quickly visualizing lineages inferred by
#'   \code{slingshot}.
#'   
#' @param x a \code{SlingshotDataSet} with results to be plotted.
#' @param type character, the type of output to be plotted, can be one of 
#'   \code{"lineages"}, \code{curves}, or \code{both} (by partial matching), see
#'   Details for more.
#' @param show.constraints logical, whether or not the user-specified initial 
#'   and terminal clusters should be specially denoted by green and red dots,
#'   respectively.
#' @param col character, color vector for points.
#' @param pch integer or character specifying the plotting symbol, see 
#'   \code{\link{par}}.
#' @param cex numeric, amount by which points should be magnified, see
#'   \code{\link{par}}.
#' @param lwd numeric, the line width, see \code{\link{par}}.
#' @param ... additional parameters for \code{plot} or \code{axis}, see 
#'   \code{\link{pairs}}.
#' @param labels character, the names of the variables, see \code{\link{pairs}}.
#' @param horInd see \code{\link{pairs}}.
#' @param verInd see \code{\link{pairs}}.
#' @param lower.panel see \code{\link{pairs}}.
#' @param upper.panel see \code{\link{pairs}}.
#' @param diag.panel see \code{\link{pairs}}.
#' @param text.panel see \code{\link{pairs}}.
#' @param label.pos see \code{\link{pairs}}.
#' @param line.main see \code{\link{pairs}}.
#' @param cex.labels see \code{\link{pairs}}.
#' @param font.labels see \code{\link{pairs}}.
#' @param row1attop see \code{\link{pairs}}.
#' @param gap see \code{\link{pairs}}.
#'   
#' @details If \code{type == 'lineages'}, straight line connectors between
#'   cluster centers will be plotted. If \code{type == 'curves'}, simultaneous
#'   princiapl curves will be plotted.
#'   
#' @details When \code{type} is not specified, the function will first check the
#'   \code{curves} slot and plot the curves, if present. Otherwise,
#'   \code{lineages} will be plotted, if present.
#'   
#' @return returns \code{NULL}.
#'   
#' @examples
#' data("slingshotExample")
#' sds <- slingshot(rd, cl, start.clus = "1")
#' pairs(sds, type = 'curves')
#' 
#' @export
pairs.SlingshotDataSet <-
    function (x, type = NULL, show.constraints = FALSE, col = NULL, 
              pch = 16, cex=1, lwd=2, ...,
              labels, horInd = seq_len(nc), verInd = seq_len(nc),
              lower.panel = FALSE, upper.panel = TRUE,
              diag.panel = NULL, text.panel = textPanel,
              label.pos = 0.5 + has.diag/3, line.main = 3,
              cex.labels = NULL, font.labels = 1,
              row1attop = TRUE, gap = 1)
    {
        #####
        lp.sling <- lower.panel
        up.sling <- upper.panel
        panel <- points
        if(!up.sling){
            upper.panel <- NULL
        }else{
            upper.panel <- panel
        }
        if(!lower.panel){
            lower.panel <- NULL
        }else{
            lower.panel <- panel
        }
        log = ""
        sds <- x
        x <- reducedDim(sds)
        curves <- FALSE
        lineages <- FALSE
        if(is.null(type)){
            if(length(slingCurves(sds)) > 0){
                type <- 'curves'
            }else if(length(slingLineages(sds)) > 0){
                type <- 'lineages'
            }else{
                stop('No lineages or curves detected.')
            }
        }else{
            type <- c('curves','lineages','both')[pmatch(type,
                                                         c('curves','lineages',
                                                           'both'))]
            if(is.na(type)){
                stop('Unrecognized type argument.')
            }
        }
        if(type %in% c('lineages','both')){
            lineages <- TRUE
        }
        if(type %in% c('curves','both')){
            curves <- TRUE
        }
        if(lineages & (length(slingLineages(sds))==0)){
            stop('No lineages detected.')
        }
        if(curves & (length(slingCurves(sds))==0)){
            stop('No curves detected.')
        }
        if(lineages){
            forest <- slingAdjacency(sds)
            clusters <- rownames(forest)
            nclus <- nrow(forest)
            centers <- t(vapply(clusters,function(clID){
                w <- clusterLabels[,clID]
                return(apply(x, 2, weighted.mean, w = w))
            }, rep(0,ncol(reducedDim(sds)))))
            rownames(centers) <- clusters
            linC <- slingParams(sds)
        }
        range.max <- max(apply(x,2,function(xi){
            r <- range(xi, na.rm = TRUE)
            return(abs(r[2] - r[1]))
        }))
        plot.ranges <- apply(x,2,function(xi){
            mid <- (max(xi,na.rm = TRUE) + min(xi,na.rm = TRUE))/2
            return(c(mid - range.max/2, mid + range.max/2))
        })
        if(is.null(col)){
            if(requireNamespace("RColorBrewer", quietly = TRUE)) {
                cc <- c(RColorBrewer::brewer.pal(9, "Set1")[-c(1,3,6)], 
                        RColorBrewer::brewer.pal(7, "Set2")[-2], 
                        RColorBrewer::brewer.pal(6, "Dark2")[-5], 
                        RColorBrewer::brewer.pal(8, "Set3")[-c(1,2)])
            } else {
                cc <- seq_len(100)
            }
            col <- cc[apply(clusterLabels(sds),1,which.max)]
        }
        #####
        if(doText <- missing(text.panel) || is.function(text.panel))
            textPanel <-
            function(x = 0.5, y = 0.5, txt, cex, font)
                text(x, y, txt, cex = cex, font = font)
        
        localAxis <- function(side, x, y, xpd, bg, col=NULL, lwd=NULL, main, 
                              oma, ...) {
            ## Explicitly ignore any color argument passed in as
            ## it was most likely meant for the data points and
            ## not for the axis.
            xpd <- NA
            if(side %% 2L == 1L && xl[j]) xpd <- FALSE
            if(side %% 2L == 0L && yl[i]) xpd <- FALSE
            if(side %% 2L == 1L) Axis(x, side = side, xpd = xpd, ...)
            else Axis(y, side = side, xpd = xpd, ...)
        }
        
        localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
        localLowerPanel <- function(..., main, oma, font.main, cex.main)
            lower.panel(...)
        localUpperPanel <- function(..., main, oma, font.main, cex.main)
            upper.panel(...)
        localDiagPanel <- function(..., main, oma, font.main, cex.main)
            diag.panel(...)
        
        dots <- list(...); nmdots <- names(dots)
        if (!is.matrix(x)) {
            x <- as.data.frame(x)
            for(i in seq_along(names(x))) {
                if(is.factor(x[[i]]) || is.logical(x[[i]]))
                    x[[i]] <- as.numeric(x[[i]])
                if(!is.numeric(unclass(x[[i]])))
                    stop("non-numeric argument to 'pairs'")
            }
        } else if (!is.numeric(x)) stop("non-numeric argument to 'pairs'")
        panel <- match.fun(panel)
        if((has.lower <- !is.null(lower.panel)) && !missing(lower.panel))
            lower.panel <- match.fun(lower.panel)
        if((has.upper <- !is.null(upper.panel)) && !missing(upper.panel))
            upper.panel <- match.fun(upper.panel)
        if((has.diag  <- !is.null( diag.panel)) && !missing( diag.panel))
            diag.panel <- match.fun( diag.panel)
        
        if(row1attop) {
            tmp <- lower.panel; lower.panel <- upper.panel; upper.panel <- tmp
            tmp <- has.lower; has.lower <- has.upper; has.upper <- tmp
        }
        
        nc <- ncol(x)
        if (nc < 2L) stop("only one column in the argument to 'pairs'")
        if(!all(horInd >= 1L && horInd <= nc))
            stop("invalid argument 'horInd'")
        if(!all(verInd >= 1L && verInd <= nc))
            stop("invalid argument 'verInd'")
        if(doText) {
            if (missing(labels)) {
                labels <- colnames(x)
                if (is.null(labels)) labels <- paste("var", 1L:nc)
            }
            else if(is.null(labels)) doText <- FALSE
        }
        oma <- if("oma" %in% nmdots) dots$oma
        main <- if("main" %in% nmdots) dots$main
        if (is.null(oma))
            oma <- c(4, 4, if(!is.null(main)) 6 else 4, 4)
        opar <- par(mfrow = c(length(horInd), length(verInd)),
                    mar = rep.int(gap/2, 4), oma = oma)
        on.exit(par(opar))
        dev.hold(); on.exit(dev.flush(), add = TRUE)
        
        xl <- yl <- logical(nc)
        if (is.numeric(log)) xl[log] <- yl[log] <- TRUE
        else {xl[] <- grepl("x", log); yl[] <- grepl("y", log)}
        for (i in if(row1attop) verInd else rev(verInd))
            for (j in horInd) {
                l <- paste0(ifelse(xl[j], "x", ""), ifelse(yl[i], "y", ""))
                localPlot(x[, j], x[, i], xlab = "", ylab = "",
                          axes = FALSE, type = "n", ..., log = l,
                          xlim = plot.ranges[,j], ylim = plot.ranges[,i])
                if(i == j || (i < j && has.lower) || (i > j && has.upper) ) {
                    box()
                    if(i == 1  && (!(j %% 2L) || !has.upper || !has.lower ))
                        localAxis(1L + 2L*row1attop, x[, j], x[, i], ...)
                    if(i == nc && (  j %% 2L  || !has.upper || !has.lower ))
                        localAxis(3L - 2L*row1attop, x[, j], x[, i], ...)
                    if(j == 1  && (!(i %% 2L) || !has.upper || !has.lower ))
                        localAxis(2L, x[, j], x[, i], ...)
                    if(j == nc && (  i %% 2L  || !has.upper || !has.lower ))
                        localAxis(4L, x[, j], x[, i], ...)
                    mfg <- par("mfg")
                    if(i == j) {
                        if (has.diag) localDiagPanel(as.vector(x[, i]), ...)
                        if (doText) {
                            par(usr = c(0, 1, 0, 1))
                            if(is.null(cex.labels)) {
                                l.wid <- strwidth(labels, "user")
                                cex.labels <- max(0.8, min(2, .9 / max(l.wid)))
                            }
                            xlp <- if(xl[i]) 10^0.5 else 0.5
                            ylp <- if(yl[j]) 10^label.pos else label.pos
                            text.panel(xlp, ylp, labels[i],
                                       cex = cex.labels, font = font.labels)
                        }
                    } else if(i < j){
                        if(up.sling){
                            points(as.vector(x[, j]), as.vector(x[, i]), 
                                   col = col, cex = cex, pch=pch, ...)
                            if(lineages){
                                for(ii in seq_len(nclus-1)){
                                    for(jj in seq(ii+1,nclus)){
                                        if(forest[ii,jj]==1){
                                            seg.col <- 1
                                            lines(centers[c(ii,jj),j], 
                                                  centers[c(ii,jj),i], 
                                                  lwd = lwd, col = seg.col, ...)
                                        }
                                    }
                                }
                                points(centers[,j],centers[,i], pch = pch, 
                                       cex=2*cex)
                                if(show.constraints){
                                    if(any(linC$start.given)){
                                        st.ind <- clusters %in% 
                                            linC$start.clus[linC$start.given]
                                        points(centers[st.ind,j],
                                               centers[st.ind,i], cex = cex, 
                                               col = 'green3', 
                                               pch = pch)
                                    }
                                    if(any(linC$end.given)){
                                        en.ind <- clusters %in% 
                                            linC$end.clus[linC$end.given]
                                        points(centers[en.ind,j],
                                               centers[en.ind,i], cex = cex, 
                                               col = 'red2', pch = pch)
                                    }
                                }
                            }
                            if(curves){
                                for(c in slingCurves(sds)){
                                    lines(c$s[c$ord,c(j,i)], lwd = lwd, 
                                          col=1, ...)
                                }
                            }
                        }
                    }
                    else{
                        if(lp.sling){
                            points(as.vector(x[, j]), as.vector(x[, i]), 
                                   col = col, cex = cex, pch=pch, ...)
                            if(lineages){
                                for(ii in seq_len(nclus-1)){
                                    for(jj in seq(ii+1,nclus)){
                                        if(forest[ii,jj]==1){
                                            if(clusters[ii] %in% 
                                               linC$start.clus | 
                                               clusters[jj] %in% 
                                               linC$start.clus){
                                                seg.col <- 'green3'
                                            }else if(clusters[ii] %in% 
                                                     linC$end.clus[
                                                         linC$end.given] |
                                                     clusters[jj] %in% 
                                                     linC$end.clus[
                                                         linC$end.given]){
                                                seg.col <- 'red2'
                                            }else{
                                                seg.col <- 1
                                            }
                                            lines(centers[c(ii,jj),j], 
                                                  centers[c(ii,jj),i], 
                                                  lwd = lwd, col = seg.col,...)
                                        }
                                    }
                                }
                                points(centers[,j],centers[,i], pch = pch, 
                                       cex = 2*cex)
                            }
                            if(curves){
                                for(c in slingCurves(sds)){ 
                                    lines(c$s[c$ord,c(j,i)],lwd = lwd, 
                                          col=1, ...) 
                                }
                            }
                        }
                    }
                    if (any(par("mfg") != mfg))
                        stop("the 'panel' function made a new plot")
                } else par(new = FALSE)
                
            }
        if (!is.null(main)) {
            font.main <- if("font.main" %in% nmdots){
                dots$font.main
            }else par("font.main")
            cex.main <- if("cex.main" %in% 
                           nmdots) dots$cex.main else par("cex.main")
            mtext(main, 3, line.main, outer=TRUE, at = 0.5, cex = cex.main, 
                  font = font.main)
        }
        invisible(NULL)
    }


