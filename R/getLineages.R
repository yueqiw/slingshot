#' @rdname getLineages
#'   
#' @description Given a reduced-dimension data matrix \code{n} by \code{p} and a
#'   vector of cluster identities (potentially including -1's for
#'   "unclustered"), this function infers a forest structure on the clusters and
#'   returns paths through the forest that can be interpreted as lineages.
#'   
#' @param data a data object containing the matrix of coordinates to be used for
#'   lineage inference. Supported types include \code{matrix}, 
#'   \code{\link{SingleCellExperiment}}, and \code{\link{SlingshotDataSet}}.
#' @param clusterLabels character, a vector of length \code{n} denoting cluster
#'   labels, optionally including \code{-1}'s for "unclustered." If
#'   \code{reducedDim} is a \code{SlingshotDataSet}, cluster labels will be
#'   taken from it.
#' @param reducedDim (optional) identifier to be used if \code{reducedDim(data)}
#'   contains multiple elements. Otherwise, the first element will be used by
#'   default.
#' @param start.clus (optional) character, indicates the cluster(s) *from* which
#'   lineages will be drawn.
#' @param end.clus (optional) character, indicates the cluster(s) which will be 
#'   forced leaf nodes in their trees.
#' @param dist.fun (optional) function, method for calculating distances between
#'   clusters. Must take two matrices as input, corresponding to points in 
#'   reduced-dimensional space. If the minimum cluster size is larger than the 
#'   number dimensions, the default is to use the joint covariance matrix to
#'   find squared distance between cluster centers. If not, the default is to
#'   use the diagonal of the joint covariance matrix.
#' @param omega (optional) numeric, this granularity parameter determines the 
#'   distance between every real cluster and the artificial cluster, OMEGA. It 
#'   is parameterized such that this distance is \code{omega / 2}, making 
#'   \code{omega} the maximum distance between two connected clusters. By 
#'   default, \code{omega = Inf}.
#'   
#' @details The \code{connectivity} matrix is learned by fitting a (possibly
#'   constrained) minimum-spanning tree on the clusters and the artificial 
#'   cluster, OMEGA, which is a fixed distance away from every real cluster.
#'   This effectively limits the maximum branch length in the MST to twice the
#'   chosen distance, meaning that the output may contain multiple trees.
#'   
#' @details Once the \code{connectivity} is known, lineages are identified in
#'   any tree with at least two clusters. For a given tree, if there is an
#'   annotated starting cluster, every possible path out of a starting cluster
#'   and ending in a leaf that isn't another starting cluster will be returned.
#'   If no starting cluster is annotated, every leaf will be considered as a
#'   potential starting cluster and whichever configuration produces the longest
#'   average lineage length (in terms of number of clusters included) will be
#'   returned.
#'   
#' @return An object of class \code{\link{SlingshotDataSet}} containing the 
#'   arguments provided to \code{getLineages} as well as the following new 
#'   elements: \itemize{ \item{\code{lineages}}{ a list of \code{L} items, where
#'   \code{L} is the number of lineages identified. Each lineage is represented
#'   by a character vector with the names of the clusters included in that
#'   lineage, in order.} \item{\code{connectivity}}{ the inferred cluster
#'   connectivity matrix.} 
#'   \item{\code{slingParams$start.given},\code{slingParams$end.given}} {
#'   logical values indicating whether the starting and ending clusters were 
#'   specified a priori.} \item{\code{slingParams$dist}}{ the pairwise
#'   cluster distance matrix.}}
#'   
#' @examples
#' data("slingshotExample")
#' sds <- getLineages(rd, cl, start.clus = '1')
#' 
#' plot(rd, col = cl, asp = 1)
#' lines(sds, type = 'l', lwd = 3)
#' 
#' @export
#' 
#' @importFrom igraph graph.adjacency
#' @importFrom igraph shortest_paths
#' @importFrom ape mst
#' @import matrixStats
#'   
setMethod(f = "getLineages",
    signature = signature(data = "matrix", 
        clusterLabels = "matrix"),
    definition = function(data, clusterLabels, reducedDim = NULL,
        start.clus = NULL, end.clus = NULL,
        dist.fun = NULL, omega = NULL){
        
        X <- as.matrix(data)
        clusterLabels <- as.matrix(clusterLabels)
        # CHECKS
        if(nrow(X)==0){
            stop('reducedDim has zero rows.')
        }
        if(ncol(X)==0){
            stop('reducedDim has zero columns.')
        }
        if(nrow(X) != nrow(clusterLabels)){
            stop('nrow(data) must equal nrow(clusterLabels).')
        }
        if(any(is.na(X))){
            stop('reducedDim cannot contain missing values.')
        }
        if(!all(apply(X,2,is.numeric))){
            stop('reducedDim must only contain numeric values.')
        }
        if (is.null(rownames(X)) &
                is.null(rownames(clusterLabels))) {
            rownames(X) <- paste('Cell', seq_len(nrow(X)), sep = '-')
            rownames(clusterLabels) <-
                paste('Cell', seq_len(nrow(X)), sep = '-')
        }
        if(is.null(colnames(X))){
            colnames(X) <- paste('Dim',seq_len(ncol(X)),sep='-')
        }
        if(is.null(colnames(clusterLabels))) {
            colnames(clusterLabels) <- seq_len(ncol(clusterLabels))
        }
        if(any(colnames(clusterLabels) == "")){
            colnames(clusterLabels)[colnames(clusterLabels)==""] <-
                which(colnames(clusterLabels)=="")
        }
        if(any(rownames(X)=='')){
            miss.ind <- which(rownames(X) == '')
            rownames(X)[miss.ind] <- paste('Cell',miss.ind,sep='-')
        }
        if(any(colnames(X)=='')){
            miss.ind <- which(colnames(X) == '')
            colnames(X)[miss.ind] <- paste('Dim',miss.ind,sep='-')
        }
        if(is.null(rownames(clusterLabels)) & 
                !is.null(rownames(X))){
            rownames(clusterLabels) <- rownames(X)
        }
        if(is.null(rownames(X)) & 
                !is.null(rownames(clusterLabels))){
            rownames(X) <- rownames(clusterLabels)
        }
        if(any(rowSums(clusterLabels)>1)){
            rs <- rowSums(clusterLabels)
            clusterLabels <- clusterLabels / rs
        }
        if(any(colSums(clusterLabels)==0)){
            clusterLabels <- clusterLabels[, colSums(clusterLabels)!=0, 
                drop = FALSE]
        }
        
        # set up, remove unclustered cells (-1's)
        X.original <- X
        clusterLabels <- clusterLabels[, colnames(clusterLabels) != -1, 
            drop = FALSE]
        clusters <- colnames(clusterLabels)
        nclus <- length(clusters)
        if(!is.null(start.clus)){
            start.clus <- as.character(start.clus)
        }
        if(!is.null(end.clus)){
            end.clus <- as.character(end.clus)
        }
        
        ### get the connectivity matrix
        # get cluster centers
        centers <- t(vapply(clusters,function(clID){
            w <- clusterLabels[,clID]
            return(colWeightedMeans(X, w = w))
        }, rep(0,ncol(X))))

        # determine the distance function
        if(is.null(dist.fun)){
            min.clus.size <- min(colSums(clusterLabels))
            if(min.clus.size <= ncol(X)){
                message('Using diagonal covariance matrix')
                dist.fun <- function(X,w1,w2) .dist_clusters_diag(X,w1,w2)
            }else{
                message('Using full covariance matrix')
                dist.fun <- function(X,w1,w2) .dist_clusters_full(X,w1,w2)
            }
        }
        
        ### get pairwise cluster distance matrix
        D <- as.matrix(vapply(clusters,function(clID1){
            vapply(clusters,function(clID2){
                w1 <- clusterLabels[,clID1]
                w2 <- clusterLabels[,clID2]
                return(dist.fun(X, w1, w2))
            },0)
        },rep(0,nclus)))
        rownames(D) <- clusters
        colnames(D) <- clusters
        
        # if infinite, set omega to largest distance + 1
        if(is.null(omega)){
            omega <- max(D) + 1
        }else{
            if(omega > 0){
                if(omega == Inf){
                    omega <- max(D) + 1
                }else{
                    omega <- omega / 2
                }
            }else{
                stop("omega must be a positive number.")
            }
        }
        D <- rbind(D, rep(omega, ncol(D)) )
        D <- cbind(D, c(rep(omega, ncol(D)), 0) )
        
        # draw MST on cluster centers + OMEGA
        # (possibly excluding endpoint clusters)
        if(! is.null(end.clus)){
            end.idx <- which(clusters %in% end.clus)
            mstree <- ape::mst(D[-end.idx, -end.idx, drop = FALSE])
        }else{
            mstree <- ape::mst(D)
        }
        # (add in endpoint clusters)
        if(! is.null(end.clus)){
            forest <- D
            forest[forest != 0] <- 0
            forest[-end.idx, -end.idx] <- mstree
            for(clID in end.clus){
                cl.idx <- which(clusters == clID)
                dists <- D[! rownames(D) %in% end.clus, cl.idx]
                # get closest non-endpoint cluster
                closest <- names(dists)[which.min(dists)] 
                closest.idx <- which.max(clusters == closest)
                forest[cl.idx, closest.idx] <- 1
                forest[closest.idx, cl.idx] <- 1
            }
        }else{
            forest <- mstree
        }
        # remove OMEGA
        forest <- forest[seq_len(nclus), seq_len(nclus), drop = FALSE] 
        rownames(forest) <- clusters
        colnames(forest) <- clusters
        
        ###############################
        ### use the "forest" to define lineages
        ###############################
        lineages <- list()
        
        # identify sub-trees
        subtrees <- subtrees.update <- forest
        diag(subtrees) <- 1
        while(sum(subtrees.update) > 0){
            subtrees.new <- apply(subtrees,2,function(col){
                rowSums(subtrees[,as.logical(col), drop=FALSE]) > 0
            })
            subtrees.update <- subtrees.new - subtrees
            subtrees <- subtrees.new
        }
        subtrees <- unique(subtrees)
        trees <- lapply(seq_len(nrow(subtrees)),function(ri){
            colnames(forest)[subtrees[ri,]]
        })
        trees <- trees[order(vapply(trees,length,0),decreasing = TRUE)]
        ntree <- length(trees)
        
        # identify lineages (paths through trees)
        for(tree in trees){
            if(length(tree) == 1){
                lineages[[length(lineages)+1]] <- tree
                next
            }
            tree.ind <- rownames(forest) %in% tree
            tree.graph <- forest[tree.ind, tree.ind, drop = FALSE]
            degree <- rowSums(tree.graph)
            g <- graph.adjacency(tree.graph, mode="undirected")
            
            # if you have starting cluster(s) in this tree, draw lineages
            # to each leaf
            if(! is.null(start.clus)){
                if(sum(start.clus %in% tree) > 0){
                    starts <- start.clus[start.clus %in% tree]
                    ends <- rownames(tree.graph)[
                        degree == 1 & ! rownames(tree.graph) %in% starts]
                    for(st in starts){
                        paths <- shortest_paths(g, from = st, to = ends, 
                            mode = 'out', 
                            output = 'vpath')$vpath
                        for(p in paths){
                            lineages[[length(lineages)+1]] <- names(p)
                        }
                    }
                }else{
                    # else, need a criteria for picking root
                    # highest average length (~parsimony)
                    leaves <- rownames(tree.graph)[degree == 1]
                    avg.lineage.length <- vapply(leaves,function(l){
                        ends <- leaves[leaves != l]
                        paths <- shortest_paths(g, from = l, to = ends, 
                            mode = 'out', 
                            output = 'vpath')$vpath
                        mean(vapply(paths, length, 0))
                    }, 0)
                    st <- names(avg.lineage.length)[
                        which.max(avg.lineage.length)]
                    ends <- leaves[leaves != st]
                    paths <- shortest_paths(g, from = st, to = ends, 
                        mode = 'out',
                        output = 'vpath')$vpath
                    for(p in paths){
                        lineages[[length(lineages)+1]] <- names(p)
                    }
                }
            }else{
                # else, need a criteria for picking root
                # highest average length (~parsimony)
                leaves <- rownames(tree.graph)[degree == 1]
                avg.lineage.length <- vapply(leaves,function(l){
                    ends <- leaves[leaves != l]
                    paths <- shortest_paths(g, from = l, to = ends, 
                        mode = 'out',
                        output = 'vpath')$vpath
                    mean(vapply(paths, length, 0))
                }, 0)
                st <- names(avg.lineage.length)[
                    which.max(avg.lineage.length)]
                ends <- leaves[leaves != st]
                paths <- shortest_paths(g, from = st, to = ends, 
                    mode = 'out',
                    output = 'vpath')$vpath
                for(p in paths){
                    lineages[[length(lineages)+1]] <- names(p)
                }
            }
        }
        # sort by number of clusters included
        lineages <- lineages[order(vapply(lineages, length, 0), 
            decreasing = TRUE)]
        names(lineages) <- paste('Lineage',seq_along(lineages),sep='')
        
        lineageControl <- list()
        first <- unique(vapply(lineages,function(l){ l[1] },''))
        last <- unique(vapply(lineages,function(l){ l[length(l)] },''))
        
        lineageControl$start.clus <- first
        lineageControl$end.clus <- last
        
        start.given <- first %in% start.clus
        end.given <- last %in% end.clus
        lineageControl$start.given <- start.given
        lineageControl$end.given <- end.given
        
        lineageControl$dist <- D[seq_len(nclus),seq_len(nclus), 
            drop = FALSE]
        
        out <- newSlingshotDataSet(reducedDim = X, 
            clusterLabels = clusterLabels, 
            lineages = lineages, 
            adjacency = forest, 
            slingParams = lineageControl)
        
        validObject(out)
        return(out)
    }
)

#' @rdname getLineages
#' @export
setMethod(f = "getLineages",
    signature = signature(data = "matrix", 
        clusterLabels = "character"),
    definition = function(data, clusterLabels, reducedDim = NULL,
        start.clus = NULL, end.clus = NULL,
        dist.fun = NULL, omega = NULL){
        
        # CHECKS
        clusterLabels <- as.character(clusterLabels)
        X <- as.matrix(data)
        if(nrow(X)==0){
            stop('reducedDim has zero rows.')
        }
        if(ncol(X)==0){
            stop('reducedDim has zero columns.')
        }
        if(nrow(X) != length(clusterLabels)){
            stop('nrow(data) must equal length(clusterLabels).')
        }
        if(any(is.na(clusterLabels))){
            message("Cluster labels of 'NA' being treated as unclustered.")
            clusterLabels[is.na(clusterLabels)] <- '-1'
        }
        
        # convert clusterLabels into cluster weights matrix
        clusters <- unique(clusterLabels)
        clusWeight <- vapply(clusters,function(clID){
            as.numeric(clusterLabels == clID)
        },rep(0,nrow(X)))
        colnames(clusWeight) <- clusters
        return(getLineages(data = data, clusterLabels = clusWeight,
            reducedDim = reducedDim,
            start.clus = start.clus, end.clus = end.clus,
            dist.fun = dist.fun, omega = omega))
    }
)

#' @rdname getLineages
#' @export
setMethod(f = "getLineages",
    signature = signature(data = "matrix", clusterLabels = "ANY"),
    definition = function(data, clusterLabels, reducedDim = NULL,
        start.clus = NULL, end.clus = NULL,
        dist.fun = NULL, omega = NULL){
        if(missing(clusterLabels)){
            message('No cluster labels provided. Continuing with ',
                'one cluster.')
            clusterLabels <- rep('1', nrow(data))
        }
        if(! any(c(length(clusterLabels), nrow(clusterLabels)) == 
                nrow(data))){
            stop("clusterLabels must have length or number of rows equal',
                'to nrow(data).")
        }
        return(getLineages(data = data, clusterLabels = clusterLabels,
            reducedDim = reducedDim,
            start.clus = start.clus, end.clus = end.clus,
            dist.fun = dist.fun, omega = omega))
    })

#' @rdname getLineages
#' @export
setMethod(f = "getLineages",
    signature = signature(data = "SlingshotDataSet", 
        clusterLabels = "ANY"),
    definition = function(data, clusterLabels,
        reducedDim = NULL,
        start.clus = NULL, end.clus = NULL,
        dist.fun = NULL, omega = NULL){
        return(getLineages(data = reducedDim(data), 
            clusterLabels = .getClusterLabels(data),
            reducedDim = reducedDim,
            start.clus = start.clus, end.clus = end.clus,
            dist.fun = dist.fun, omega = omega))
    })

#' @rdname getLineages
#' @export
setMethod(f = "getLineages",
    signature = signature(data = "data.frame", 
        clusterLabels = "ANY"),
    definition = function(data, clusterLabels, reducedDim = NULL,
        start.clus = NULL, end.clus = NULL,
        dist.fun = NULL, omega = NULL){
        RD <- as.matrix(data)
        rownames(RD) <- rownames(data)
        return(getLineages(data = RD, clusterLabels = clusterLabels,
            reducedDim = reducedDim,
            start.clus = start.clus, end.clus = end.clus,
            dist.fun = dist.fun, omega = omega))
    })

#' @rdname getLineages
#' @export
setMethod(f = "getLineages",
    signature = signature(data = "matrix", 
        clusterLabels = "numeric"),
    definition = function(data, clusterLabels, reducedDim = NULL,
        start.clus = NULL, end.clus = NULL,
        dist.fun = NULL, omega = NULL){
        return(getLineages(data = data, 
            clusterLabels = as.character(clusterLabels),
            reducedDim = reducedDim,
            start.clus = start.clus, end.clus = end.clus,
            dist.fun = dist.fun, omega = omega))
    })

#' @rdname getLineages
#' @export
setMethod(f = "getLineages",
    signature = signature(data = "matrix", 
        clusterLabels = "factor"),
    definition = function(data, clusterLabels, reducedDim = NULL,
        start.clus = NULL, end.clus = NULL,
        dist.fun = NULL, omega = NULL){
        return(getLineages(data = data, 
            clusterLabels = as.character(clusterLabels), 
            reducedDim = reducedDim,
            start.clus = start.clus, end.clus = end.clus,
            dist.fun = dist.fun, omega = omega))
    })