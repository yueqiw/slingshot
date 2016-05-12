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
    clusterMap = "matrix",
    reducedDim = "matrix",
    lineages = "list"
  )
)