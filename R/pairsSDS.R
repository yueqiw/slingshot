pairs.SlingshotDataSet <-
  function (x, labels, col = NULL, cex=1, lwd=2, ...,
            horInd = 1:nc, verInd = 1:nc,
            lower.panel = FALSE, upper.panel = TRUE,
            diag.panel = NULL, text.panel = textPanel,
            label.pos = 0.5 + has.diag/3, line.main = 3,
            cex.labels = NULL, font.labels = 1,
            row1attop = TRUE, gap = 1, show.constraints = FALSE,
            type = NULL, pch =16)
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
      if(length(curves(sds)) > 0){
        type <- 'curves'
      }else if(length(lineages(sds)) > 0){
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
    if(lineages & (length(lineages(sds))==0)){
      stop('No lineages detected.')
    }
    if(curves & (length(curves(sds))==0)){
      stop('No curves detected.')
    }
    if(lineages){
      forest <- connectivity(sds)
      clusters <- rownames(forest)
      nclus <- nrow(forest)
      centers <- t(sapply(clusters,function(clID){
        x.sub <- x[clusterLabels(sds) == clID,]
        return(colMeans(x.sub))
      }))
      rownames(centers) <- clusters
      linC <- lineageControl(sds)
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
      cc <- c(brewer.pal(9, "Set1")[-c(1,3,6)], brewer.pal(7, "Set2")[-2], brewer.pal(6, "Dark2")[-5], brewer.pal(8, "Set3")[-c(1,2)])
      col <- cc[as.factor(clusterLabels(sds))]
    }
    #####
    if(doText <- missing(text.panel) || is.function(text.panel))
      textPanel <-
      function(x = 0.5, y = 0.5, txt, cex, font)
        text(x, y, txt, cex = cex, font = font)
    
    localAxis <- function(side, x, y, xpd, bg, col=NULL, lwd=NULL, main, oma, ...) {
      ## Explicitly ignore any color argument passed in as
      ## it was most likely meant for the data points and
      ## not for the axis.
      xpd <- NA
      if(side %% 2L == 1L && xl[j]) xpd <- FALSE
      if(side %% 2L == 0L && yl[i]) xpd <- FALSE
      if(side %% 2L == 1L) Axis(x, side = side, xpd = xpd, ...)
      else Axis(y, side = side, xpd = xpd, ...)
    }
    
    slingshotPanel <- function(..., main, oma, font.main, cex.main){
      
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
              points(as.vector(x[, j]), as.vector(x[, i]), col = col, cex = cex, pch=pch, ...)
              if(lineages){
                for(ii in 1:(nclus-1)){
                  for(jj in (i+1):nclus){
                    if(forest[ii,jj]==1){
                      # if(show.constraints & (clusters[ii] %in% linC$start.clus | clusters[jj] %in% linC$start.clus)){
                      #   seg.col <- brewer.pal(4,'Set1')[3]
                      # }else if(show.constraints & (clusters[ii] %in% linC$end.clus[linC$end.given] | clusters[jj] %in% linC$end.clus[linC$end.given])){
                      #   seg.col <- brewer.pal(4,'Set1')[1]
                      # }else{
                      #   seg.col <- 1
                      # }
                      seg.col <- 1
                      lines(centers[c(ii,jj),j], centers[c(ii,jj),i], lwd = lwd, col = seg.col, ...)
                    }
                  }
                }
                points(centers[,j],centers[,i], pch = pch, cex=2*cex)
                if(show.constraints){
                  if(any(linC$start.given)){
                    st.ind <- clusters %in% linC$start.clus[linC$start.given]
                    points(centers[st.ind,j],centers[st.ind,i], cex = cex, col = brewer.pal(4,'Set1')[3], pch = pch)
                  }
                  if(any(linC$end.given)){
                    en.ind <- clusters %in% linC$end.clus[linC$end.given]
                    points(centers[en.ind,j],centers[en.ind,i], cex = cex, col = brewer.pal(4,'Set1')[1], pch = pch)
                  }
                }
              }
              if(curves){
                for(c in curves(sds)){ lines(c$s[c$tag,c(j,i)], lwd = lwd, col=1, ...) }
              }
            }
          }
          else{
            if(lp.sling){
              points(as.vector(x[, j]), as.vector(x[, i]), col = col, cex = cex, pch=pch, ...)
              if(lineages){
                for(ii in 1:(nclus-1)){
                  for(jj in (i+1):nclus){
                    if(forest[ii,jj]==1){
                      if(clusters[ii] %in% linC$start.clus | clusters[jj] %in% linC$start.clus){
                        seg.col <- brewer.pal(4,'Set1')[3]
                      }else if(clusters[ii] %in% linC$end.clus[linC$end.given] | clusters[jj] %in% linC$end.clus[linC$end.given]){
                        seg.col <- brewer.pal(4,'Set1')[1]
                      }else{
                        seg.col <- 1
                      }
                      lines(centers[c(ii,jj),j], centers[c(ii,jj),i], lwd = lwd, col = seg.col, ...)
                    }
                  }
                }
                points(centers[,j],centers[,i], pch = pch, cex=2*cex)
              }
              if(curves){
                for(c in curves(sds)){ lines(c$s[c$tag,c(j,i)], lwd = lwd, col=1, ...) }
              }
            }
          }
          if (any(par("mfg") != mfg))
            stop("the 'panel' function made a new plot")
        } else par(new = FALSE)
        
      }
    if (!is.null(main)) {
      font.main <- if("font.main" %in% nmdots) dots$font.main else par("font.main")
      cex.main <- if("cex.main" %in% nmdots) dots$cex.main else par("cex.main")
      mtext(main, 3, line.main, outer=TRUE, at = 0.5, cex = cex.main, font = font.main)
    }
    invisible(NULL)
  }