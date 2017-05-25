## plot
#' @title Plot Slingshot output
#' 
#' @description Tools for visualizing lineages inferred by \code{slingshot}.
#' 
#' @param x a \code{SlingshotDataSet} with results to be plotted.
#' @param type character, the type of output to be plotted, can be one of
#' \code{"lineages"}, \code{curves}, or \code{both} (by partial matching), see
#' Details for more.
#' @param add logical, indicates whether the output should be added to an existing plot.
#' @param dims numeric, which dimensions to plot (default is \code{1:2}).
#' @param asp numeric, the y/x aspect ratio, see \code{\link{plot.window}}.
#' @param ... additional parameters to be passed to \code{\link{lines}}.
#' 
#' @details If \code{type == 'lineages'}, straight line connectors between cluster
#' centers will be plotted. If \code{type == 'curves'}, simultaneous princiapl curves
#' will be plotted.
#'
#' @details When \code{type} is not specified, the function will first check the 
#' \code{curves} slot and plot the curves, if present. Otherwise, \code{lineages} 
#' will be plotted, if present.
#'
#' @return returns \code{NULL}.
#'
#' @examples
#' data("slingshotExample")
#' sds <- slingshot(reducedDim, clusLabels, start.clus = "5")
#' plot(sds, type = 'b')
#' 
#' # add to existing plot
#' plot(reducedDim, col = 'grey50')
#' lines(sds, lwd = 3)
#' 
#' @export
setMethod(
  f = "plot",
  signature = "SlingshotDataSet",
  definition = function(x,
                        type = NULL,
                        add = FALSE,
                        dims = 1:2,
                        asp = 1,
                        ...) {
    curves <- FALSE
    lineages <- FALSE
    if(is.null(type)){
      if(length(x@curves) > 0){
        type <- 'curves'
      }else if(length(x@lineages) > 0){
        type <- 'lineages'
      }else{
        stop('No lineages or curves detected.')
      }
    }else{
      type <- c('curves','lineages','both')[pmatch(type,c('curves','lineages','both'))]
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
    
    if(lineages & (length(x@lineages)==0)){
      stop('No lineages detected.')
    }
    if(curves & (length(x@curves)==0)){
      stop('No curves detected.')
    }
    
    if(lineages){
      X <- x@reducedDim
      clusLabels <- x@clusLabels
      connectivity <- x@connectivity
      clusters <- rownames(connectivity)
      nclus <- nrow(connectivity)
      centers <- t(sapply(clusters,function(clID){
        x.sub <- X[clusLabels == clID,]
        return(colMeans(x.sub))
      }))
      rownames(centers) <- clusters
      X <- X[clusLabels %in% clusters,]
      clusLabels <- clusLabels[clusLabels %in% clusters]
    }
    
    if(!add){
      xs <- NULL
      ys <- NULL
      if(lineages){
        xs <- c(xs, centers[,dims[1]])
        ys <- c(ys, centers[,dims[2]])
      }
      if(curves){
        xs <- c(xs, as.numeric(sapply(x@curves, function(c){ c$s[,dims[1]] })))
        ys <- c(ys, as.numeric(sapply(x@curves, function(c){ c$s[,dims[2]] })))        
      }
      plot(x = NULL, y = NULL, asp = asp,
           xlim = range(xs), ylim = range(ys),
           xlab = colnames(x@reducedDim)[dims[1]],
           ylab = colnames(x@reducedDim)[dims[2]])
    }
    
    if(lineages){
      for(i in 1:(nclus-1)){
        for(j in (i+1):nclus){
          if(connectivity[i,j]==1){
            lines(x = centers[c(i,j),dims[1]], y = centers[c(i,j),dims[2]], ...)
          }
        }
      }
    }
    if(curves){
      for(c in x@curves){ lines(c$s[c$tag,dims], ...) }
    }
    invisible(NULL)
  }
)

## lines --- HOW DO I COMBINE THE DOC FILES FOR THIS AND PLOT?
#' @rdname SlingshotDataSet-methods
#' @export
setMethod(
  f = "lines",
  signature = "SlingshotDataSet",
  definition = function(x,
                        type = NULL,
                        dims = 1:2,
                        ...) {
    plot(x, type = type, add = TRUE, dims = dims, ...)
    invisible(NULL)
  }
)

## plot3d
#' @title Plot Slingshot output in 3D
#' 
#' @description Tools for visualizing lineages inferred by \code{slingshot}.
#' 
#' @param x a \code{SlingshotDataSet} with results to be plotted.
#' @param type character, the type of output to be plotted, can be one of
#' \code{"lineages"}, \code{curves}, or \code{both} (by partial matching), see
#' Details for more.
#' @param add logical, indicates whether the output should be added to an existing plot.
#' @param dims numeric, which dimensions to plot (default is \code{1:3}).
#' @param aspect either a logical indicating whether to adjust the aspect ratio 
#' or a new ratio, see \code{\link{plot3d}}.
#' @param ... additional parameters to be passed to \code{\link{lines3d}}.
#' 
#' @details If \code{type == 'lineages'}, straight line connectors between cluster
#' centers will be plotted. If \code{type == 'curves'}, simultaneous princiapl curves
#' will be plotted.
#'
#' @details When \code{type} is not specified, the function will first check the 
#' \code{curves} slot and plot the curves, if present. Otherwise, \code{lineages} 
#' will be plotted, if present.
#'
#' @return returns \code{NULL}.
#'
#' @examples
#' data("slingshotExample")
#' reducedDim <- cbind(reducedDim, rnorm(nrow(reducedDim)))
#' sds <- slingshot(reducedDim, clusLabels, start.clus = "5")
#' plot3d(sds, type = 'b')
#' 
#' # add to existing plot
#' plot3d(reducedDim, col = 'grey50', aspect = 'iso')
#' lines3d(sds, lwd = 3)
#' 
#' @importFrom rgl plot3d
#' 
#' @export
# plot3d.SlingshotDataSet <- function(x,
#                                     type = NULL,
#                                     add = FALSE,
#                                     dims = 1:3,
#                                     aspect = 'iso',
#                                     ...){
#   require(rgl)
#   curves <- FALSE
#   lineages <- FALSE
#   if(is.null(type)){
#     if(length(x@curves) > 0){
#       type <- 'curves'
#     }else if(length(x@lineages) > 0){
#       type <- 'lineages'
#     }else{
#       stop('No lineages or curves detected.')
#     }
#   }else{
#     type <- c('curves','lineages','both')[pmatch(type,c('curves','lineages','both'))]
#     if(is.na(type)){
#       stop('Unrecognized type argument.')
#     }
#   }
#   
#   if(type %in% c('lineages','both')){
#     lineages <- TRUE
#   }
#   if(type %in% c('curves','both')){
#     curves <- TRUE
#   }
#   
#   if(lineages & (length(x@lineages)==0)){
#     stop('No lineages detected.')
#   }
#   if(curves & (length(x@curves)==0)){
#     stop('No curves detected.')
#   }
#   
#   if(lineages){
#     X <- x@reducedDim
#     clusLabels <- x@clusLabels
#     connectivity <- x@connectivity
#     clusters <- rownames(connectivity)
#     nclus <- nrow(connectivity)
#     centers <- t(sapply(clusters,function(clID){
#       x.sub <- X[clusLabels == clID,]
#       return(colMeans(x.sub))
#     }))
#     rownames(centers) <- clusters
#     X <- X[clusLabels %in% clusters,]
#     clusLabels <- clusLabels[clusLabels %in% clusters]
#   }
#   
#   if(!add){
#     xs <- NULL
#     ys <- NULL
#     zs <- NULL
#     if(lineages){
#       xs <- c(xs, centers[,dims[1]])
#       ys <- c(ys, centers[,dims[2]])
#       zs <- c(zs, centers[,dims[3]])
#     }
#     if(curves){
#       xs <- c(xs, as.numeric(sapply(x@curves, function(c){ c$s[,dims[1]] })))
#       ys <- c(ys, as.numeric(sapply(x@curves, function(c){ c$s[,dims[2]] })))        
#       zs <- c(zs, as.numeric(sapply(x@curves, function(c){ c$s[,dims[3]] })))        
#     }
#     plot3d(x = NULL, y = NULL, z = NULL, aspect = aspect,
#            xlim = range(xs), ylim = range(ys), zlim = range(zs), 
#            xlab = colnames(x@reducedDim)[dims[1]],
#            ylab = colnames(x@reducedDim)[dims[2]],
#            zlab = colnames(x@reducedDim)[dims[3]])
#   }
#   
#   if(lineages){
#     for(i in 1:(nclus-1)){
#       for(j in (i+1):nclus){
#         if(connectivity[i,j]==1){
#           lines3d(x = centers[c(i,j),dims[1]], y = centers[c(i,j),dims[2]], z = centers[c(i,j),dims[3]], ...)
#         }
#       }
#     }
#   }
#   if(curves){
#     for(c in x@curves){ lines3d(c$s[c$tag,dims], ...) }
#   }
#   invisible(NULL)
#   
# }


## lines3d
#' @rdname SlingshotDataSet-methods
#' @export
# lines3d.SlingshotDataSet <- function(x,
#                                      type = NULL,
#                                      dims = 1:3,
#                                      ...) {
#   plot3d(x, type = type, add = TRUE, dims = dims, ...)
#   invisible(NULL)
# }






# setMethod(
#   f = "plot3d",
#   signature = "SlingshotDataSet",
#   definition = function(x,
#                         type = NULL,
#                         add = FALSE,
#                         dims = 1:3,
#                         aspect = 'iso',
#                         ...) {
#     require(rgl)
#     curves <- FALSE
#     lineages <- FALSE
#     if(is.null(type)){
#       if(length(x@curves) > 0){
#         type <- 'curves'
#       }else if(length(x@lineages) > 0){
#         type <- 'lineages'
#       }else{
#         stop('No lineages or curves detected.')
#       }
#     }else{
#       type <- c('curves','lineages','both')[pmatch(type,c('curves','lineages','both'))]
#       if(is.na(type)){
#         stop('Unrecognized type argument.')
#       }
#     }
#     
#     if(type %in% c('lineages','both')){
#       lineages <- TRUE
#     }
#     if(type %in% c('curves','both')){
#       curves <- TRUE
#     }
#     
#     if(lineages & (length(x@lineages)==0)){
#       stop('No lineages detected.')
#     }
#     if(curves & (length(x@curves)==0)){
#       stop('No curves detected.')
#     }
#     
#     if(lineages){
#       X <- x@reducedDim
#       clusLabels <- x@clusLabels
#       connectivity <- x@connectivity
#       clusters <- rownames(connectivity)
#       nclus <- nrow(connectivity)
#       centers <- t(sapply(clusters,function(clID){
#         x.sub <- X[clusLabels == clID,]
#         return(colMeans(x.sub))
#       }))
#       rownames(centers) <- clusters
#       X <- X[clusLabels %in% clusters,]
#       clusLabels <- clusLabels[clusLabels %in% clusters]
#     }
#     
#     if(!add){
#       xs <- NULL
#       ys <- NULL
#       zs <- NULL
#       if(lineages){
#         xs <- c(xs, centers[,dims[1]])
#         ys <- c(ys, centers[,dims[2]])
#         zs <- c(zs, centers[,dims[3]])
#       }
#       if(curves){
#         xs <- c(xs, as.numeric(sapply(x@curves, function(c){ c$s[,dims[1]] })))
#         ys <- c(ys, as.numeric(sapply(x@curves, function(c){ c$s[,dims[2]] })))        
#         zs <- c(zs, as.numeric(sapply(x@curves, function(c){ c$s[,dims[3]] })))        
#       }
#       plot3d(x = NULL, y = NULL, z = NULL, aspect = aspect,
#              xlim = range(xs), ylim = range(ys), zlim = range(zs), 
#              xlab = colnames(x@reducedDim)[dims[1]],
#              ylab = colnames(x@reducedDim)[dims[2]],
#              zlab = colnames(x@reducedDim)[dims[3]])
#     }
#     
#     if(lineages){
#       for(i in 1:(nclus-1)){
#         for(j in (i+1):nclus){
#           if(connectivity[i,j]==1){
#             lines3d(x = centers[c(i,j),dims[1]], y = centers[c(i,j),dims[2]], z = centers[c(i,j),dims[3]], ...)
#           }
#         }
#       }
#     }
#     if(curves){
#       for(c in x@curves){ lines3d(c$s[c$tag,dims], ...) }
#     }
#     invisible(NULL)
#   }
# )










# setMethod(
#   f = "lines3d",
#   signature = "SlingshotDataSet",
#   definition = function(x,
#                         type = NULL,
#                         dims = 1:3,
#                         ...) {
#     plot3d(x, type = type, add = TRUE, dims = dims, ...)
#     invisible(NULL)
#   }
# )
