#' feature_coexp
#' 
#' The function calculates coexpression patterns between genes
#' @param so a SiFINeT object
#' @return SiFINeT object with coexp (gene coexpression matrix) updated.
#' @details The coexpression pattern of a pair of genes is a normalized co-occurrence of high (or equivalently low) expression level of the 2 genes in the classified count matrix.
#' The normalization is based on the estimated quantiles of the low-high separation instead of the quantiles used for quantile regressions. 
#' Theoretically, the distribution of coexpression patterns should asymptotically follow standard Gaussian distribution if at least one of the 2 genes is not differentially expressed feature gene. 
#' @export
#' 
feature_coexp <- function(so){
  if (so@sparse == FALSE){
    so@coexp <- cal_coexp(so@data.thres[[so@data.name]])
  } else {
    so@coexp <- cal_coexp_sp(so@data.thres[[so@data.name]])
  }
  diag(so@coexp) <- 0
  return(so)
}
