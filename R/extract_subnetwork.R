#' extract_subnetwork
#' 
#' The function extract a subnetwork from the co-expression network
#' @param so a SiFINeT object
#' @param target_gene_name the names of the target genes in the output network
#' @param target_gene_id the indices of the target genes in the output network, not used when target_gene_name is not Null
#' @param positive whether only positive (default) co-expressions or all co-expressions are considered in assigning edges
#' @return an adjacency matrix of the output subnetwork
#' @export
#' 
extract_subnetwork <- function(so, target_gene_name = NULL, target_gene_id = NULL, positive = TRUE){
  if (is.null(target_gene_name)){
    if (is.null(target_gene_id)){
      idx <- so@kset
    } else {
      idx <- intersect(so@kset, target_gene_id)
    }
  } else {
    match_gene <- match(target_gene_name, so@gene.name)
    idx <- intersect(so@kset, match_gene)
  }
  if (positive){
    out <- (so@coexp[idx, idx] - so@est_ms$mean) > so@thres
  } else {
    out <- abs(so@coexp[idx, idx] - so@est_ms$mean) > so@thres
  }
  rownames(out) <- so@gene.name[idx]
  colnames(out) <- so@gene.name[idx]
  return(out)
}