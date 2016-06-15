merge_lineages <- function(X, clus.labels, lineages, curves = NA, ...){
  if(is.na(curves)){
    curves <- get_curves(X, clus.labels, lineages, ...)
  }
  
  

}




data("toy_data")
source('~/Projects/slingshot/R/helper_functions.R')
source('~/Projects/slingshot/R/slingshot.R')
source('~/Projects/slingshot/R/plotting.R')
lines <- get_lineages(X, clus.labels, start.clus = 'a')
plot_tree(X,clus.labels,lines,threeD=F,dim=2)

L <- ncol(lines$C)

test_lineage_split <- function(line1, line2, X, clus.labels, lineages, curves = NA, ...){
  L <- length(grep("lineage",names(lineages))) # number of lineages
  suppressWarnings(if(is.na(curves)){
    curves <- get_curves(X, clus.labels, lineages, ...)
  })
  line1.clus <- lineages[[line1]]
  line2.clus <- lineages[[line2]]
  clus2combine <- c(setdiff(line1.clus,line2.clus), setdiff(line2.clus,line1.clus))
  sizes <- sapply(clus2combine,function(clID){sum(clus.labels == clID)})
  name <- clus2combine[which.max(sizes)]
  name <- paste('new.',name,sep='')
  clus.combine <- clus.labels
  clus.combine[clus.combine %in% clus2combine] <- name
  
  lineages.combine <- list()
  lineage.new <- line1.clus[line1.clus %in% line2.clus]
  lineage.new <- c(lineage.new, name)
  lineages.combine[[1]] <- lineage.new
  names(lineages.combine)[1] <- 'lineage.new'
  if(L > 2){
    lineages.combine[2:(L-1)] <- lineages[1:L][-c(line1,line2)]
    names(lineages.combine)[2:(L-1)] <- names(lineages[1:L][-c(line1,line2)])
  }
  curves.combine <- get_curves(X, clus.combine, lineages.combine, ...)
  
  lam1 <- get.lam(X[clus.labels %in% clus2combine,], curves.combine[[1]]$s)
  lam2.1 <- get.lam(X[clus.labels %in% clus2combine & clus.labels %in% line1.clus,], curves[[line1]]$s)
  lam2.2 <- get.lam(X[clus.labels %in% clus2combine & clus.labels %in% line2.clus,], curves[[line2]]$s)

  d1 <- .sq_segment_lengths(X[clus.labels %in% clus2combine,], lam1$s)
  d2.1 <- .sq_segment_lengths(X[clus.labels %in% clus2combine & clus.labels %in% line1.clus,], lam2.1$s)
  d2.2 <- .sq_segment_lengths(X[clus.labels %in% clus2combine & clus.labels %in% line2.clus,], lam2.2$s)
  d2 <- c(d2.1,d2.2)
  
  n <- sum(clus.labels %in% clus2combine)
  s1 <- sum(d1)/(n-1)
  s2 <- sum(d2)/(n-1)
  stat <- 2*(sum(dchisq(d1/s1^2,1)) - sum(d2/s2^2,1))
  #pval <- pchisq(stat,df=2,lower.tail = F) # how many degrees of freedom?
  return(stat)
}




matrixJitter <- function(mat){
  apply(mat,2,jitter,amount=0)
}

pvals <- sapply(1:100,function(ii){
  clus.labels <- rep(letters[1:7],each=50)
  x <- c(rep(-4,50),rep(-1,50),rep(0,100),rep(c(3,6,9),each=50)) + rnorm(350)
  y <- c(rep(-4,50),rep(c(-1,4,8),each=50),rep(-2,100),rep(0,50)) + rnorm(350)
  X <- cbind(x,y); rm(x,y)
  X <- matrixJitter(matrixJitter(X))
  lineages <- get_lineages(X,clus.labels,start.clus = 'a')
  if(all(lineages[[2]] == c('a','b','c','d')) & all(lineages[[1]] == c('a','b','e','f','g'))){
    fake_lineages <- list()
    fake_lineages[[1]] <- c('a','b','e','g')
    fake_lineages[[2]] <- c('a','b','e','f')
    fake_lineages[[3]] <- c('a','b','c','d')
    names(fake_lineages) <- c('lineage1','lineage2','lineage3')
    p1 <- test_lineage_split(1,2,X,clus.labels,lineages)
    p2 <- test_lineage_split(1,2,X,clus.labels,fake_lineages)
    return(c(p1,p2))
  }else{
    return(c(NA,NA))
  }
})












