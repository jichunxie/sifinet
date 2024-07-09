#' cal_connectivity
#'
#' The function calculates the 1st, 2nd and 3rd order connectivities for all genes
#' @param so a SiFINeT object
#' @param m number of neighbors sampled each time for the calculation of 3rd order connectivity
#' @param niter number of samples created for the calculation of 3rd order connectivity
#' @return SiFINeT object with conn (absolute network connectivities) updated.
#' @details For gene i, First order connectivity is defined as the number of edges connected to gene i (degree of the gene node i in the network).
#' Second order connectivity is defined as the proportion of edges between the neighbors of gene i, calculated as number of observed edges between the neighbors of gene i divided by the number of possible edges between the neighbors.
#' Third order connectivity is defined as a weighted proportion of edges between neighbors and neighbors of neighbors of gene i. Third order connectivity is calculated as the mean of edge proportions across weighted samples.
#' Each gene is weighted by the number of edges it has with the neighbors of gene i. Then SiFINeT repeatedly samples m genes for niter times. For each sample, the edge proportion (number of observed edges / number of possible edges) is calculated. And the mean edge proportion across the sample is the 3rd order connectivity for gene i.
#' @export
#'
cal_connectivity <- function(so, m = 10, niter = 100) {
  conn <- cal_conn(so@coexp[so@kset, so@kset] - so@est_ms$mean, so@thres,
    abso = TRUE, m = m, niter = niter
  )
  for (i in 1:3) {
    conn[[i]][is.na(conn[[i]])] <- 0
  }
  conn$name <- so@gene.name[so@kset]
  so@conn <- conn
  return(so)
}
