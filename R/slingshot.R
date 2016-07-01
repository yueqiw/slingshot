#' @title Identify transcriptional trajectories and order cells
#' 
#' @description Performs the \code{slingshot} lineage inference procedure for continuous development processes.
#' 
#' @param X numeric, the nxp matrix of samples
#' 
#' @details TODO
#'
#' @return TODO
#'
#' @examples
#' data("toy_data")
#' 
# @export
#' 

slingshot <- function(X, clus.labels = NULL, start.clus = NULL, end.clus = NULL, dist.fun = NULL, omega = Inf, distout = TRUE, thresh = 0.0001, maxit = 100, stretch = 2, shrink = TRUE){
  lineages <- get_lineages(X, clus.labels, start.clus, end.clus, dist.fun, omega, distout)
  curves <- get_curves(X, clus.labels, lineages, thresh = 0.0001, maxit = 100, stretch = 2, shrink = TRUE)
  return(curves)
}
