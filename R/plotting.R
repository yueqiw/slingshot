require(RColorBrewer)
require(rgl)


plot_tree <- function(X, clus.labels, lineages, threeD = FALSE, dim = NA, col = NA, labels = TRUE){
  forest <- lineages$forest
  clusters <- rownames(forest)
  nclus <- nrow(forest)
  centers <- t(sapply(clusters,function(clID){
    x.sub <- X[clus.labels == clID,]
    return(colMeans(x.sub))
  }))
  rownames(centers) <- clusters
  if(is.na(col)){
    cc <- c(brewer.pal(9, "Set1")[-c(1,3)], brewer.pal(7, "Set2")[-2], brewer.pal(6, "Dark2")[-5], brewer.pal(8, "Set3")[-c(1,2)])
    center.col <- cc[1:nclus]
  }else{
    center.col <- rep_len(col, length.out = nclus)
  }
  clus.col <- sapply(clus.labels,function(clID){center.col[which(clusters==clID)]})
  if(threeD){
    plot3d(X[,1:3],col=clus.col,size=5,box=FALSE,aspect = 'iso')
    for(i in 1:(nclus-1)){
      for(j in (i+1):nclus){
        if(forest[i,j]==1){
          if(clusters[i] %in% lineages$start.clus | clusters[j] %in% lineages$start.clus){
            seg.col <- brewer.pal(4,'Set1')[3]
          }else if(clusters[i] %in% lineages$end.clus | clusters[j] %in% lineages$end.clus){
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
      text.col[rownames(centers) %in% lineages$start.clus] <- brewer.pal(4,'Set1')[3]
      text.col[rownames(centers) %in% lineages$end.clus] <- brewer.pal(4,'Set1')[1]
      text3d(centers, texts = rownames(centers), adj = c(1,1), add = TRUE, color = text.col)
    }
  }else{
    if(is.na(dim)){
      dim <- ncol(X)
    }
    par(mfrow=c(dim-1,dim-1),mar=c(4,4,.5,.5))
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
}


plot_curves <- function(X,clus.labels,curves, threeD = TRUE, dim = NA, col = NA){
  clusters <- unique(clus.labels)
  nclus <- length(clusters)
  if(is.na(col)){
    cc <- c(brewer.pal(9, "Set1")[-c(1,3)], brewer.pal(7, "Set2")[-2], brewer.pal(6, "Dark2")[-5], brewer.pal(8, "Set3")[-c(1,2)])
    center.col <- cc[1:nclus]
  }else{
    center.col <- rep_len(col, length.out = nclus)
  }
  clus.col <- sapply(clus.labels,function(clID){center.col[which(clusters==clID)]})
  if(threeD){
    plot3d(X[,1:3],col=cc[as.factor(clus.labels)],size=5,box=FALSE,aspect='iso')
    for(i in 1:length(curves)){lines3d(curves[[i]]$s,lwd=3)}
  }else{
    if(is.na(dim)){
      dim <- ncol(X)
    }
    par(mfrow=c(dim-1,dim-1),mar=c(4,4,.5,.5))
    for(ii in 1:(dim-1)){
      for(jj in 2:dim){
        if(ii<jj){
          y <- X[,ii]; x <- X[,jj]
          if(abs(ii-jj)>1){
            plot(x,y,col=clus.col,pch=16,asp=1,ylab='',xlab='',axes=F); box()
          }else{
            plot(x,y,col=clus.col,pch=16,asp=1,ylab=colnames(X)[ii],xlab=colnames(X)[jj])
          }
          for(i in 1:length(curves)){lines(curves[[i]]$s[,c(jj,ii)],lwd=2)}
        }else{
          plot.new()
        }
      }
    }
    par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))
  }
}

