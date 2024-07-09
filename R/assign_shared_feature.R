#' assign_shared_feature
#'
#' The function assigns non-unique candidate feature genes as shared feature genes into unique feature gene sets
#' @param so a SiFINeT object
#' @param min_edge_prop minimum proportion of edges between a gene and a unique feature gene set for the new gene to be assigned to the set
#' @return SiFINeT object with shared feature genes in featureset updated.
#' @details Candidate feature genes that are not chosen as unique feature genes would be reconsidered as shared feature genes.
#' A non-unique candidate feature gene would be assigned to a unique feature gene group if it is connected to more than min_edge_prop of the unique genes in the group.
#' @export
#'
assign_shared_feature <- function(so, min_edge_prop = 0.4) {
  coex_f <- so@coexp[so@kset[so@fg_id], so@kset[so@fg_id]] - so@est_ms$mean
  diag(coex_f) <- 0
  edge_mat <- coex_f >= so@thres
  shared_cand <- setdiff(1:length(so@fg_id), so@uni_fg_id)

  if (length(shared_cand) > 0) {
    idx <- list()
    for (i in 1:length(so@selected_cluster)) {
      idx[[i]] <- so@uni_fg_id[so@uni_cluster == so@selected_cluster[i]]
    }

    m_conn <- matrix(0, length(shared_cand), length(so@selected_cluster))
    for (i in 1:length(shared_cand)) {
      for (j in 1:length(so@selected_cluster)) {
        m_conn[i, j] <- mean(edge_mat[shared_cand[i], idx[[j]]])
      }
    }

    m_conn_a <- m_conn >= min_edge_prop
    for (i in 1:length(so@selected_cluster)) {
      temp <- so@conn2$name[shared_cand[(m_conn_a[, i] == T)]]
      so@featureset$shared[[i]] <- temp[order(temp)]
    }
  } else {
    for (i in 1:length(so@selected_cluster)) {
      so@featureset$shared[[i]] <- vector()
    }
  }

  return(so)
}
