#' @title Plot cluster connectivity
#' 
#' @description \code{plot_tree} visualizes the output from \code{\link{get_lineages}}, 
#' showing cells (colored by cluster) in a reduced-dimensional space with connecting 
#' lines between cluster centers representing the inferred structure.
#' 
#' @param X numeric, the \code{n x p} matrix of samples in a reduced dimensionality 
#'  space.
#' @param clus.labels character, a vector of length n denoting cluster labels.
#' @param lineages list, the out put of \code{\link{get_lineages}}, this list denotes 
#'  which clusters belong to each lineage and contains the inferred connectivity 
#'  between clusters.
#' @param threeD logical indicates whether to make a 3D plot with the \code{rgl} 
#'  package.
#' @param dim total number of dimensions to be shown in a series of two-dimensional 
#'  plots, similar to \code{pairs} plots (only applicable if \code{threeD} is false).
#' @param col.clus (optional) vector of colors to use for denoting clusters.
#' @param labels logical indicating whether to include labels on cluster centers.
#' 
#' @details Plots cells as points in a reduced-dimensional space, colored by cluster.
#'  If a \code{lineages} argument is given, the plot will include straight lines
#'  between cluster centers for connections between clusters.
#'
#' @details The default behavior is to produce a series of 2-D plots representing
#'  all pairwise combinations of dimensions in \code{X}, or the first \code{dim}
#'  dimensions. If \code{threeD=TRUE}, a three-dimensional plot will be created,
#'  using the \code{rgl} package.
#'
#' @return returns \code{NULL}.
#'
#' @examples
#' data("slingshot_example")
#' lin <- get_lineages(X, clus.labels, start.clus = 'a')
#' plot_tree(X, clus.labels, lin)
#' 
#' @export
#'
#' @import RColorBrewer
#' @import rgl
#' 

plot_tree <- function(X, clus.labels, lineages = NULL, threeD = FALSE, dim = NA, col.clus = NULL, labels = TRUE){
  forest <- lineages$forest
  clusters <- rownames(forest)
  nclus <- nrow(forest)
  centers <- t(sapply(clusters,function(clID){
    x.sub <- X[clus.labels == clID,]
    return(colMeans(x.sub))
  }))
  rownames(centers) <- clusters
  X <- X[clus.labels %in% clusters,]
  clus.labels <- clus.labels[clus.labels %in% clusters]
  if(is.null(col.clus)){
    cc <- c(brewer.pal(9, "Set1")[-c(1,3)], brewer.pal(7, "Set2")[-2], brewer.pal(6, "Dark2")[-5], brewer.pal(8, "Set3")[-c(1,2)])
    center.col <- cc[1:nclus]
  }else{
    center.col <- rep_len(col.clus, length.out = nclus)
  }
  clus.col <- vapply(clus.labels,function(clID) {center.col[which(clusters==clID)]}, character(1))
  if(threeD){
    plot3d(X[,1:3],col=clus.col,size=5,box=FALSE,aspect = 'iso')
    for(i in 1:(nclus-1)){
      for(j in (i+1):nclus){
        if(forest[i,j]==1){
          if(clusters[i] %in% lineages$start.clus | clusters[j] %in% lineages$start.clus){
            seg.col <- brewer.pal(4,'Set1')[3]
          }else if(clusters[i] %in% lineages$end.clus[lineages$end.given] | clusters[j] %in% lineages$end.clus[lineages$end.given]){
            seg.col <- brewer.pal(4,'Set1')[1]
          }else{
            seg.col <- 1
          }
          lines3d(centers[c(i,j),1], centers[c(i,j),2], centers[c(i,j),3], lwd=3, col = seg.col)
        }
      }
    }
    plot3d(centers, size = 8, add = TRUE, col = center.col)
    plot3d(centers, size=9, add = TRUE)
    if(labels){
      text.col <- rep(1,nclus)
      #text.col[rownames(centers) %in% lineages$start.clus] <- brewer.pal(4,'Set1')[3]
      #text.col[rownames(centers) %in% lineages$end.clus] <- brewer.pal(4,'Set1')[1]
      text3d(centers, texts = rownames(centers), adj = c(1,1), add = TRUE, color = text.col)
    }
  }else{
    if(is.na(dim)){
      dim <- ncol(X)
    }
    par(mfrow=c(dim-1,dim-1),mar=c(4,4,.2,.2))
    for(ii in 1:(dim-1)){
      for(jj in 2:dim){
        if(ii<jj){
          y <- X[,ii]; x <- X[,jj]
          if(abs(ii-jj)>1){
            plot(x,y,col=clus.col,pch=16,asp=1,ylab='',xlab='',axes=F); box()
          }else{
            plot(x,y,col=clus.col,pch=16,asp=1,ylab=colnames(X)[ii],xlab=colnames(X)[jj])
          }
          for(i in 1:(nclus-1)){
            for(j in (i+1):nclus){
              if(forest[i,j]==1){
                if(clusters[i] %in% lineages$start.clus | clusters[j] %in% lineages$start.clus){
                  seg.col <- brewer.pal(4,'Set1')[3]
                }else if(clusters[i] %in% lineages$end.clus[lineages$end.given] | clusters[j] %in% lineages$end.clus[lineages$end.given]){
                  seg.col <- brewer.pal(4,'Set1')[1]
                }else{
                  seg.col <- 1
                }
                lines(centers[c(i,j),jj], centers[c(i,j),ii], lwd=1.5, col = seg.col)
              }
            }
          }
          points(centers[,jj],centers[,ii])
          points(centers[,jj],centers[,ii],pch=16,col=center.col)
          #text(centers, labels = rownames(centers))
        }else{
          plot.new()
        }
      }
    }
    par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))
  }
  invisible(NULL)
}


#' @title Plot smooth curves
#' 
#' @description \code{plot_curves} visualizes the output from \code{\link{get_curves}}, 
#'  showing cells (colored by cluster) in a reduced-dimensional space with smooth 
#'  curves representing each inferred lineage.
#' 
#' @param X numeric, the \code{n x p} matrix of samples in a reduced dimensionality 
#'  space.
#' @param clus.labels character, a vector of length n denoting cluster labels.
#' @param curves list, the output of \code{\link{get_curves}}, this list includes 
#'  matrices for each curve in the reduced-dimensional space.
#' @param threeD logical indicates whether to make a 3D plot with the \code{rgl} 
#'  package.
#' @param dim total number of dimensions to be shown in a series of two-dimensional 
#'  plots, similar to \code{pairs} plots (only applicable if \code{threeD} is false).
#' @param col.clus (optional) vector of colors to use for denoting clusters.
#' @param col.lin (optional) vector of colors to use for denoting lineages.
#' 
#' @details Plots cells as points in a reduced-dimensional space, colored by cluster.
#'  If a \code{curves} argument is given, the plot will include curves representing
#'  each lineage.
#'
#' @details The default behavior is to produce a series of 2-D plots representing
#'  all pairwise combinations of dimensions in \code{X}, or the first \code{dim}
#'  dimensions. If \code{threeD=TRUE}, a three-dimensional plot will be created,
#'  using the \code{rgl} package.
#'
#' @return returns \code{NULL}.
#'
#' @examples
#' data("slingshot_example")
#' lin <- get_lineages(X, clus.labels, start.clus = 'a')
#' crv <- get_curves(X, clus.labels, lin)
#' plot_curves(X, clus.labels, crv, threeD = FALSE)
#' 
#' @export
#'
#' @import RColorBrewer
#' @import rgl
#' 

plot_curves <- function(X, clus.labels, curves = NULL, threeD = FALSE, dim = NA, col.clus = NULL, col.lin = NULL){
  X <- X[clus.labels != '-1',]
  clus.labels <- clus.labels[clus.labels != '-1']
  clusters <- unique(clus.labels)
  nclus <- length(clusters)
  if(is.null(col.clus)){
    cc <- c(brewer.pal(9, "Set1")[-c(1,3)], brewer.pal(7, "Set2")[-2], brewer.pal(6, "Dark2")[-5], brewer.pal(8, "Set3")[-c(1,2)])
    center.col <- cc[1:nclus]
  }else{
    center.col <- rep_len(col.clus, length.out = nclus)
  }
  if(! is.null(col.lin)){
    lin.col <- rep_len(col.lin, length.out = length(curves))
  }else{
    lin.col <- rep(1, length(curves))
  }
  clus.col <- sapply(clus.labels,function(clID){center.col[which(clusters==clID)]})
  if(threeD){
    plot3d(X[,1:3],col=clus.col,size=5,box=FALSE,aspect='iso')
    for(i in 1:length(curves)){lines3d(curves[[i]]$s[curves[[i]]$tag,],lwd=3, col=lin.col[i])}
  }else{
    if(is.na(dim)){
      dim <- ncol(X)
    }
    par(mfrow=c(dim-1,dim-1),mar=c(4,4,.2,.2))
    for(ii in 1:(dim-1)){
      for(jj in 2:dim){
        if(ii<jj){
          y <- X[,ii]; x <- X[,jj]
          if(abs(ii-jj)>1){
            plot(x,y,col=clus.col,pch=16,asp=1,ylab='',xlab='',axes=F); box()
          }else{
            plot(x,y,col=clus.col,pch=16,asp=1,ylab=colnames(X)[ii],xlab=colnames(X)[jj])
          }
          for(i in 1:length(curves)){lines(curves[[i]]$s[curves[[i]]$tag,c(jj,ii)],lwd=2, col = lin.col[i])}
        }else{
          plot.new()
        }
      }
    }
    par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))
  }
  invisible(NULL)
}

