#' find_unique_feature
#'
#' The function finds the clustered unique feature genes
#' @param so a SiFINeT object
#' @param t1 feature gene selection parameter, lower threshold for 1st order connectivity in absolute network
#' @param t2 feature gene selection parameter, lower threshold for 2nd order connectivity in absolute network
#' @param t3 feature gene selection parameter, lower threshold for 3rd order connectivity in absolute network
#' @param t1p unique feature gene selection parameter, lower threshold for 1st order connectivity in positive sub-network
#' @param t2p unique feature gene selection parameter, lower threshold for 2nd order connectivity in positive sub-network
#' @param t3p unique feature gene selection parameter, lower threshold for 3rd order connectivity in positive sub-network
#' @param resolution resolution for louvain clustering of unique feature genes
#' @param min_set_size minimum size for a unique feature gene cluster to be a separate unique feature gene set
#' @return SiFINeT object with fg_id (candidate feature gene index), uni_fg_id (candidate unique feature gene index), conn2 (connectivities in positive sub-network),
#' uni_cluster (cluster of candidate unique feature genes), selected_cluster (selected unique feature gene clusters), and unique feature genes in featureset updated.
#' @details SiFINeT first find genes with high 1st (>= t1), 2nd (>= t2) and 3rd (>= t3) order connectivities in absolute network (conn) to be candidate feature genes.
#' Then a positive sub-network is created where only candidate feature gene nodes (fg_id) and edges with positive coexpression patterns (coexp >= thres) are included.
#' Feature genes genes with high 1st (>= t1p), 2nd (>= t2p) and at least moderate 3rd (>= t3p) order connectivities in positive sub-network (conn2) are chosen to be candidate unique feature genes.
#' Note that when the network is not too sparse, t3p should usually be smaller than t2p for the detection of unique feature genes in transition cell types.
#' The candidate unique feature genes are then separated into groups by louvain clustering (with resolution defined by the resolution parameter),
#' and among them large groups (number of genes greater than min_set_size) are chosen to be unique feature gene sets that represent different cell types.
#' @import igraph
#' @export
#'
find_unique_feature <- function(so,
                                t1 = 5, t2 = 0.4, t3 = 0.3,
                                t1p = 5, t2p = 0.7, t3p = 0.5,
                                resolution = 1, min_set_size = 5) {
  so@fg_id <- which((so@conn$C1 >= t1) & (so@conn$C2 >= t2) & (so@conn$C3 >= t3))
  coex_f <- so@coexp[so@kset[so@fg_id], so@kset[so@fg_id]] - so@est_ms$mean
  diag(coex_f) <- 0
  conn2 <- cal_conn(coex_f, so@thres, abso = FALSE)
  conn2$name <- so@gene.name[so@kset[so@fg_id]]
  for (i in 1:3) {
    conn2[[i]][is.na(conn2[[i]])] <- 0
  }
  so@conn2 <- conn2

  so@uni_fg_id <- which((conn2$C1 >= t1p) & (conn2$C2 >= t2p) & (conn2$C3 >= t3p))
  coex_uf <- coex_f[so@uni_fg_id, so@uni_fg_id]

  g <- graph_from_adjacency_matrix(coex_uf >= so@thres, mode = "undirected")
  cl_g <- cluster_louvain(g, resolution = resolution)
  so@uni_cluster <- cl_g$membership
  tab <- table(so@uni_cluster)
  so@selected_cluster <- as.numeric(names(tab))[tab >= min_set_size]

  for (i in 1:length(so@selected_cluster)) {
    temp <- so@conn2$name[so@uni_fg_id[so@uni_cluster == so@selected_cluster[i]]]
    so@featureset$unique[[i]] <- temp[order(temp)]
  }
  return(so)
}
