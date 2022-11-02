#' create_network
#' 
#' The function estimates the null distribution of coexpression patterns 
#' and generates coexpression network
#' @param so a SiFINeT object
#' @param alpha the Type I error rate used for FDR control procedure
#' @param manual whether to manually set threshold for edge assignment
#' @param least_edge_prop the minimum proportion of edges. Only used when manual = TRUE
#' @return SiFINeT object with est_ms (estimated mean and sd) and thres (network edge threshold) updated.
#' @details Theoretically the distribution of coexpression patterns would converge to standard Gaussian if either one of the gene pair is not feature gene. 
#' However in genomics analysis, empirical null could be much more variable than theoretical null. 
#' SiFINeT uses estimated null mean and standard deviation to find the threshold for network edges. An edge is assigned to a pair of gene if the absolute value of coexpression pattern between the 2 genes is greater than the threshold
#' Assuming the distribution to be Gaussian, with the estimated null mean and standard deviation, SiFINeT uses SQUAC to control the false discovery rate (FDR) for coexpression patterns. 
#' In case the signal is not strong enough and the coexpression network is too sparse, SiFINeT also accept user-defined lower bound for the least proportion of edges. 
#' Usually a coexpression network with edge proportion between 0.5% - 10% would have better performance for the detection of feature gene sets.
#' @references Jiashun Jin and Tony T. Cai. “Estimating the Null and the Proportion of Non-Null Effects in Large-Scale Multiple Comparisons”. In: Journal of the American Statistical Association 102 (478 2004), pp. 495–506. doi: 10.1198/016214507000000167.
#' @references Jichun Xie and Ruosha Li. “False discovery rate control for high dimensional networks of quantile associations conditioning on covariates”. In: J R Stat Soc Series B Stat Methodol (2018). doi: 10.1111/rssb.12288.
#' @export
#' 
create_network <- function(so, alpha = 0.05, manual = FALSE, least_edge_prop = 0.01){
  coex_vec <- so@coexp[upper.tri(so@coexp, diag = FALSE)]
  so@est_ms <- EstNull(coex_vec)
  so@thres <- norm_FDR_SQAUC(coex_vec, so@est_ms$mean, so@est_ms$std, alpha)
  if (manual){
    so@thres <- min(so@thres, quantile(abs(coex_vec - so@est_ms$mean), (1 - least_edge_prop)))
  }
  return(so)
}