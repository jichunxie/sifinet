#' norm_FDR_SQAUC
#' The function controls the false discovery rate (FDR) of coexpression patterns
#' using SQAUC method
#' Jichun Xie and Ruosha Li. "False discovery rate control for high dimensional 
#' networks of quantile associations conditioning on covariates". 
#' In: J R Stat Soc Series B Stat Methodol (2018). doi: 10.1111/rssb.12288.
#' @param value the vector of coexpression patterns
#' @param sam_mean the estimated sample mean
#' @param sam_sd the estimated sample sd
#' @param alpha the type I error rate
#' @param n the number of cells
#' @param p the number of genes
#' @return lower bound threshold for genes to be significantly coexpressed
#' @export

norm_FDR_SQAUC <- function(value, sam_mean, sam_sd, alpha, n, p){
  tp = 2*sqrt(log(max(n,p)))*sam_sd
  d <- length(value)
  value_a <- abs(value - sam_mean)
  value.ord <- order(value_a, decreasing=TRUE)
  
  value_a_s <- value_a[value.ord]
  P22s <- 2*(1-pnorm(value_a_s, mean = 0, sd = sam_sd))
  FDR2h <- P22s*d/(1:d)
  R2 <- max(which(FDR2h<=alpha))
  if (is.infinite(R2)){
    return(tp)
  } else if (value_a_s[R2] > tp){
    return(tp)
  } else {
    return(value_a_s[R2])
  }
}