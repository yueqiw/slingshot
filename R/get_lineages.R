#' @title Infer Lineage Structure from Clustered Samples with Known Structure
#' 
#' @description Given a connectivity matrix, this function returns paths through the graph that can be interpreted as lineages.
#' 
#' @param map matrix, the binary, KxK connectivity matrix
#' @param start.clus (optional) character, indicates the cluster(s) *from* which lineages will be drawn
#' @param end.clus (optional) character, indicates the cluster(s) which will be forced leaf nodes in their trees
#' 
#' @details Lineages are identified in any connected subgraph with at least two clusters. For a given subgraph, if there is an annotated starting cluster, every possible path out of a starting cluster and ending in a leaf that isn't another starting cluster will be returned. If no starting cluster is annotated, every leaf will be considered as a potential starting cluster and whichever configuration produces the longest average lineage length (in terms of number of clusters included) will be returned.
#'
#' @return a list with L+2 items where L is the number of lineages identified. The first L items are character vectors with the names of the clusters included in that lineage. The last two items are \code{forest}, the connectivity matrix, and \code{C}, a clusters x lineages identity matrix.
#'
#' @examples
#' data("toy_data")
#' map <- get_clustMap(X, clus.labels)
#' get_lineages(map, start.clus = 'a')
#' 
#' @export
#' 
#' @import igraph
#' 

get_lineages <- function(map, start.clus = NULL, end.clus = NULL){
  clusters <- rownames(map)
  lineages <- list()
  
  # identify trees
  unused <- rownames(map)
  trees <- list()
  ntree <- 0
  while(length(unused) > 0){
    ntree <- ntree + 1
    newtree <- get_connections(unused[1], map)
    trees[[ntree]] <- newtree
    unused <- unused[! unused %in% newtree]
  }
  trees <- trees[order(sapply(trees,length),decreasing = T)]
  
  # identify lineages (paths through trees)
  for(tree in trees){
    if(length(tree) == 1){
      next # don't draw a lineage for a single-cluster tree
    }
    tree.ind <- rownames(map) %in% tree
    tree.graph <- map[tree.ind, tree.ind]
    degree <- rowSums(tree.graph)
    g <- igraph::graph.adjacency(tree.graph, mode="undirected")
    
    # if you have starting cluster(s) in this tree, draw lineages to each leaf
    if(! is.null(start.clus)){
      if(sum(start.clus %in% tree) > 0){
        starts <- start.clus[start.clus %in% tree]
        ends <- rownames(tree.graph)[degree == 1 & ! rownames(tree.graph) %in% starts]
        for(st in starts){
          paths <- igraph::shortest_paths(g, from = st, to = ends, mode = 'out', output = 'vpath')$vpath
          for(p in paths){
            lineages[[length(lineages)+1]] <- names(p)
          }
        }
      }else{
        # else, need a criteria for picking root
        # highest average length (~parsimony, but this was just the easiest thing I came up with)
        leaves <- rownames(tree.graph)[degree == 1]
        avg.lineage.length <- sapply(leaves,function(l){
          ends <- leaves[leaves != l]
          paths <- igraph::shortest_paths(g, from = l, to = ends, mode = 'out', output = 'vpath')$vpath
          mean(sapply(paths, length))
        })
        st <- names(avg.lineage.length)[which.max(avg.lineage.length)]
        ends <- leaves[leaves != st]
        paths <- igraph::shortest_paths(g, from = st, to = ends, mode = 'out', output = 'vpath')$vpath
        for(p in paths){
          lineages[[length(lineages)+1]] <- names(p)
        }
      }
    }else{
      # else, need a criteria for picking root
      # highest average length (~parsimony, but this was just the easiest thing I came up with)
      leaves <- rownames(tree.graph)[degree == 1]
      avg.lineage.length <- sapply(leaves,function(l){
        ends <- leaves[leaves != l]
        paths <- shortest_paths(g, from = l, to = ends, mode = 'out', output = 'vpath')$vpath
        mean(sapply(paths, length))
      })
      st <- names(avg.lineage.length)[which.max(avg.lineage.length)]
      ends <- leaves[leaves != st]
      paths <- shortest_paths(g, from = st, to = ends, mode = 'out', output = 'vpath')$vpath
      for(p in paths){
        lineages[[length(lineages)+1]] <- names(p)
      }
    }
  }
  # sort by number of clusters included
  lineages <- lineages[order(sapply(lineages, length), decreasing = TRUE)]
  out <- lineages
  # include "map" and clusters x lineages (C) matrices
  out$map <- map
  C <- sapply(lineages,function(lin){
    sapply(clusters,function(clID){
      as.numeric(clID %in% lin)
    })
  })
  rownames(C) <- clusters
  # should probably come up with a better name than C
  out$C <- C
  return(out)
}
