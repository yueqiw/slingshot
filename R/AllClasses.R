#' @title Class CellLineages
#'
#' @description \code{CellLineages} is a class that extends
#' \code{SummarizedExperiment} and is used to store 
#'
#' @docType class
#' @aliases CellLineages CellLineages-class cellLineages
#'
#' @description In addition to the slots of the \code{SummarizedExperiment}
#' class, the \code{CellLineages} object has the additional slots described
#' in the Slots section.
#'
#' @description There are several methods implemented for this class. The most
#' important methods (e.g., \code{\link{clusterMany}}, \code{\link{combineMany}},
#' ...) have their own help page. Simple helper methods are described in the
#' Methods section. For a comprehensive list of methods specific to this class
#' see the Reference Manual.
#'
#' @slot transformation function. Function to transform the data by when methods
#' that assume normal-like data (e.g. log)
#' @slot clusterMatrix matrix. A matrix giving the integer-valued cluster ids
#' for each sample. The rows of the matrix correspond to clusterings and columns
#' to samples. The integer values are assigned in the order that the clusters
#' were found, if found by setting sequential=TRUE in clusterSingle. "-1" indicates
#' the sample was not clustered.
#' @slot primaryIndex numeric. An index that specifies the primary set of
#' labels.
#' @slot clusterInfo list. A list with info about the clustering.
#' If created from \code{\link{clusterSingle}}, clusterInfo will include the
#' parameter used for the call, and the call itself. If \code{sequential = TRUE}
#' it will also include the following components.
#' \itemize{
#' \item{\code{clusterInfo}}{if sequential=TRUE and clusters were successfully
#' found, a matrix of information regarding the algorithm behavior for each
#' cluster (the starting and stopping K for each cluster, and the number of
#' iterations for each cluster).}
#' \item{\code{whyStop}}{if sequential=TRUE and clusters were successfully
#' found, a character string explaining what triggered the algorithm to stop.}
#' }
#' @slot clusterTypes character vector with the origin of each column of
#' clusterMatrix.
#' @slot dendro_samples dendrogram. A dendrogram containing the cluster
#' relationship (leaves are samples; see \code{\link{makeDendrogram}} for
#' details).
#' @slot dendro_clusters dendrogram. A dendrogram containing the cluster
#' relationship (leaves are clusters; see \code{\link{makeDendrogram}} for
#' details).
#' @slot dendro_index numeric. An integer giving the cluster that was used to
#'   make the dendrograms. NA_real_ value if no dendrograms are saved.
#' @slot coClustering matrix. A matrix with the cluster co-occurrence
#' information; this can either be based on subsampling or on co-clustering
#' across parameter sets (see \code{clusterMany}). The matrix is a square matrix
#' with number of rows/columns equal to the number of samples.
#' @slot clusterLegend a list, one per cluster in \code{clusterMatrix}. Each
#' element of the list is a matrix with nrows equal to the number of different
#' clusters in the clustering, and consisting of at least two columns with the
#' following column names: "clusterId" and "color".
#' @slot orderSamples a numeric vector (of integers) defining the order of
#' samples to be used for plotting of samples. Usually set internally by other
#' functions.
#'
#' @name CellLineages-class
#' @aliases CellLineages
#' @rdname CellLineages-class
#' @import SummarizedExperiment
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
#'
setClass(
  Class = "CellLineages",
  contains = "SummarizedExperiment",
  slots = list(
    cluster = "character",
    reducedDim = "matrix",
    clusterMap = "matrix",
    lineages = "list",
    pseudotime = "matrix",
    curves = "list"
  )
)


setValidity("CellLineages", function(object){
  #browser()
  # reducedDim
  if(!is.null(object@reducedDim)){
    if(!(nrow(object@reducedDim) == ncol(object))){
      return("If present, `reducedDim` must have as many rows as cells.")
    }
  }
  
  cluster <- colData(object)$cluster
  if(!is.null(cluster)){
    clusterNames <- unique(cluster[! cluster %in% c(-1,NA)])
    K <- length(clusterNames)
    
    # clusterMap
    if(!is.null(object@clusterMap) &&
       (nrow(object@clusterMap) != ncol(object@clusterMap)
        | ncol(object@clusterMap) != K)) {
      return("`clusterMap` must be a cluster by cluster matrix.")
    }
    
    # lineages
    if(!is.null(lineages)){
      if(!all(lapply(lineages,inherits,'character'))){
        return("`lineages` must be a list of character vectors.")
      }
      if(!all(unlist(lineages)%in%clusterNames)){
        return("`lineages` must only contain cluster names")
      }
    }
    
  }
  
  
  
  
  if(length(assays(object)) < 1) {
    return("There must be at least one assay slot.")
  }
  if(!is.numeric(assay(object))) {
    return("The data must be numeric.")
  }
  if(any(is.na(assay(object)))) {
    return("NA values are not allowed.")
  }
  tX <- try(transform(object),silent=TRUE)
  if(inherits(tX, "try-error")){
    stop(paste("User-supplied `transformation` produces error on the input data
               matrix:\n",x))
  }
  if(any(is.na(tX))) {
    return("NA values after transforming data matrix are not allowed.")
  }
  
  if(!all(is.na((object@clusterMatrix))) &
     !(NROW(object@clusterMatrix) == NCOL(object))) {
    return("If present, `clusterMatrix` must have as many row as cells.")
  }
  if(!is.numeric(object@clusterMatrix)) {
    return("`clusterMatrix` must be a numeric matrix.")
  }
  
  if(NCOL(object@clusterMatrix)!= length(object@clusterTypes)) {
    return("length of clusterTypes must be same as NCOL of the clusterMatrix")
  }
  
  if(NCOL(object@clusterMatrix)!= length(object@clusterInfo)) {
    return("length of clusterInfo must be same as NCOL of the clusterMatrix")
  }
  
  ##Check dendrograms
  if(!is.null(object@dendro_samples)){
    if(nobs(object@dendro_samples) != NCOL(object)) {
      return("dendro_samples must have the same number of leaves as the number of samples")
    }
  }
  else{
    if(!is.null(object@dendro_clusters)) return("dendro_samples should not be null if dendro_clusters is non-null")
  }
  if(!is.null(object@dendro_clusters)){
    if(is.na(object@dendro_index)) return("if dendrogram slots are filled, must have corresponding dendro_index defined.")
    dcluster<-clusterMatrix(object)[,object@dendro_index]
    if(nobs(object@dendro_clusters) != max(dcluster)) {
      return("dendro_clusters must have the same number of leaves as the number of (non-negative) clusters")
    }
  }
  else{
    if(!is.null(object@dendro_samples)) return("dendro_clusters should not be null if dendro_samples is non-null")
  }
  if(!is.null(object@coClustering) &&
     (NROW(object@coClustering) != NCOL(object@coClustering)
      | NCOL(object@coClustering) != NCOL(object))) {
    return("`coClustering` must be a sample by sample matrix.")
  }
  if(!all(is.na(object@clusterMatrix))){ #what does this mean, how can they be all NA?
    #check primary index
    if(length(object@primaryIndex) != 1) {
      if(length(object@primaryIndex) == 0) return("If more than one set of clusterings, a primary cluster must
                                                  be specified.")
      if(length(object@primaryIndex) > 0) return("Only a single primary index may be specified")
    }
    if(object@primaryIndex > NCOL(object@clusterMatrix) |
       object@primaryIndex < 1) {
      return("`primaryIndex` out of bounds.")
    }
    #check clusterTypes
    if(NCOL(object@clusterMatrix) != length(object@clusterTypes)) {
      return("`clusterTypes` must be the same length as NCOL of
             `clusterMatrix`.")
    }
    #check internally stored as integers
    testConsecIntegers<-apply(object@clusterMatrix,2,function(x){
      whCl<-which(!x %in% c(-1,-2))
      uniqVals<-unique(x[whCl])
      return(all(sort(uniqVals)==1:length(uniqVals)))
    })
    #browser()
    if(!all(testConsecIntegers)) return("the cluster ids in clusterMatrix must be stored internally as consecutive integer values")
    
    ####
    #test that colnames of clusterMatrix appropriately aligns with everything else
    ####
    if(is.null(colnames(object@clusterMatrix))) return("clusterMatrix must have column names")
    if(any(duplicated(colnames(object@clusterMatrix)))) return("clusterMatrix must have unique column names")
    if(!is.null(names(object@clusterTypes))) return("clusterTypes should not have names")
    if(!is.null(names(object@clusterInfo))) return("clusterInfo should not have names")
    if(!is.null(names(object@clusterLegend))) return("clusterLegend should not have names")
    ####
    #test that @clusterLegend is proper form
    ####
    if(length(object@clusterLegend) != NCOL(object@clusterMatrix)) {
      return("`clusterLegend` must be list of same length as NCOL of
             `clusterMatrix`")
    }
    testIsMatrix <- sapply(object@clusterLegend,
                           function(x) {!is.null(dim(x))})
    if(!all(testIsMatrix)) {
      return("Each element of `clusterLegend` list must be a matrix")
    }
    testColorRows <- sapply(object@clusterLegend, function(x){nrow(x)})
    testClusterMat <- apply(object@clusterMatrix, 2, function(x) {
      length(unique(x))})
    if(!all(testColorRows == testClusterMat)) {
      return("each element of `clusterLegend` must be matrix with number of
             rows equal to the number of clusters (including -1 or -2 values)
             in `clusterMatrix`")
    }
    testColorCols1 <- sapply(object@clusterLegend, function(x) {
      "color" %in% colnames(x)})
    testColorCols2 <- sapply(object@clusterLegend, function(x) {
      "clusterIds" %in% colnames(x)})
    testColorCols3 <- sapply(object@clusterLegend, function(x) {
      "name" %in% colnames(x)})
    if(!all(testColorCols1) || !all(testColorCols2) || !all(testColorCols3)) {
      return("each element of `clusterLegend` must be matrix with at least 3
             columns, and at least 3 columns have names `clusterIds`,
             `color` and `name`")
    }
    #     testUniqueName <- sapply(object@clusterLegend, function(x) {
    #       any(duplicated(x[,"name"]))})
    #     if(any(testUniqueName)) return("the column")
    testColorCols1 <- sapply(object@clusterLegend, function(x){is.character(x)})
    if(!all(testColorCols1)) {
      return("each element of `clusterLegend` must be matrix of character
             values")
    }
    testColorCols1 <- sapply(1:length(object@clusterLegend), function(ii){
      col<-object@clusterLegend[[ii]]
      x<-object@clusterMatrix[,ii]
      y<-as.numeric(col[,"clusterIds"])
      all(y %in% x)
    })
    if(!all(testColorCols1)) {
      return("each element of `clusterLegend` must be matrix with column
             `clusterIds` matching the corresponding integer valued
             clusterMatrix values")
    }
    }
  if(length(object@orderSamples)!=NCOL(assay(object))) {
    return("`orderSamples` must be of same length as number of samples
           (NCOL(assay(object)))")
  }
  if(any(!object@orderSamples %in% 1:NCOL(assay(object)))) {
    return("`orderSamples` must be values between 1 and the number of samples.")
  }
  return(TRUE)
})
