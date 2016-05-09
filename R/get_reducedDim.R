

get_reducedDim <- function(expr, method="pca", k=2, ...){
  if(method=="pca"){
    return(prcomp(t(expr),...)$x[,1:k])
  }
  if(method=="laplacian"){
    require(scater)
    require(embeddr)
    sce <- newSCESet(exprsData = as.data.frame(expr))
    sce <- embeddr(sce,p=k,...)
    return(sce@reducedDimension)
  }
  if(method=='ica'){
    require(fastICA)
    return(fastICA(t(expr), n.comp = k))
  }
}