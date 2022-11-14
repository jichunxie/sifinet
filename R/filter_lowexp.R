#' filter_lowexp
#' 
#' The function filters out genes with low expression rate and 
#' high positive coexpression with genes of same expression level
#' @param so a SiFINeT object
#' @param t1 threshold for number of total edges connecting the feature node. Lower t1 leads to stricter filtering.
#' @param t2 threshold for the proportion of positive edges. Lower t2 leads to stricter filtering.
#' @param t3 threshold for the proportion of edges with features of same expression level. Lower t3 leads to stricter filtering.
#' @return SiFINeT object with kset (kept index set) updated.
#' @details When using only mean expression level as independent variable in quantile regression, 
#' it is observed that genes with low expression level tend to have large positive $S_{ij}$ with genes that have same median expression level.
#' To reduce the coexpression noise caused by low expression level, it is preferred to filter out genes which have large amount and high proportion of positive coexpressions with genes sharing same median expression level.
#' @export
#' 
filter_lowexp <- function(so, t1 = 10, t2 = 0.9, t3 = 0.9){
  r_set <- c()
  for (target in 0:2){
    screen_set <- which(so@q5 == target)
    if (length(screen_set) >= 2){
      abs_sum <- colSums(abs(so@coexp[, screen_set] - so@est_ms$mean) > so@thres)
      pos_sum_w <- colSums(so@coexp[screen_set, screen_set] > (so@thres + so@est_ms$mean))
      pos_sum <- colSums(so@coexp[, screen_set] > (so@thres + so@est_ms$mean))
      r_set <- c(r_set, screen_set[which((abs_sum >= t1) & (pos_sum / abs_sum >= t2) & (pos_sum_w / pos_sum >= t3))])
    }
  }
  
  so@kset <- setdiff(1:length(so@q5), r_set)
  return(so)
}