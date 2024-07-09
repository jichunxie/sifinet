#' enrich_feature_set
#'
#' The function chooses genes that are not found to be feature genes as enriched feature genes and assigns them into unique+shared feature gene sets
#' @param so a SiFINeT object
#' @param min_edge_prop minimum proportion of edges between a gene and a unique+shared feature gene set for the new feature to be assigned to the set
#' @return SiFINeT object with enriched feature genes in geneset updated.
#' @details Genes that are not selected as feature genes would be added in the enriched section of the feature gene set if they are connected with
#' more than min_edge_prop of the unique and shared feature genes in each of the feature gene group.
#' @export
#'
enrich_feature_set <- function(so, min_edge_prop = 0.9) {
  featureset_all <- unique(c(unlist(so@featureset$unique), unlist(so@featureset$shared)))
  for (i in 1:length(so@selected_cluster)) {
    coex_temp <- so@coexp[so@kset[match(c(so@featureset$unique[[i]], so@featureset$shared[[i]]), so@gene.name[so@kset])], so@kset] - so@est_ms$mean
    temp <- colSums(coex_temp >= so@thres) >= min_edge_prop * nrow(coex_temp)
    temp2 <- setdiff(so@gene.name[so@kset[temp]], featureset_all)
    so@featureset$enriched[[i]] <- temp2[order(temp2)]
  }
  return(so)
}
