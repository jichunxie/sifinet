#' quantile_thres
#'
#' The function classifies count data into binary low (0) - high (1) data, based on whether the count number is greater than a threshold.
#' @param so a SiFINeT object
#' @return SiFINeT object with data.thres (categorized count matrix) updated.
#' @details The threshold used for classification is defined by quantile regression on each gene using Frisch–Newton interior point method ("fn" option for method variable in quantreg package, rq function).
#' By default if no meta data is provided, the quantile regression would be applied on the mean expression level of each cell.
#' The quantile to be estimated in the quantile regression is set to be the estimated 50% quantile of the non-zero part of the expression level for each gene.
#' If the expression level of a gene is low (with median 0), then the threshold is set to be 0.
#' @references Koenker, R. and S. Portnoy (1997) The Gaussian Hare and the Laplacean Tortoise: Computability of Squared-error vs Absolute Error Estimators, (with discussion). Statistical Science, 12, 279-300.
#' @references Roger Koenker (2022). quantreg: Quantile Regression. R package version 5.94. https://CRAN.R-project.org/package=quantreg
#' @import quantreg
#' @import Matrix
#' @export
#'
quantile_thres <- function(so) {
  n <- so@n
  p <- so@p
  Z <- so@meta.data
  if (ncol(Z) == 0) {
    Z <- rowMeans(so@data[[so@data.name]])
  }
  if (so@sparse == FALSE) {
    dt <- matrix(0, n, p)
    for (j in 1:p) {
      temp <- so@data[[so@data.name]][, j]
      v5 <- quantile(temp, 0.5)
      if ((v5 == 0) | (v5 >= quantile(temp, 0.99))) {
        dt[, j] <- (temp > 0) * 1
      } else {
        quant <- sum(temp <= v5) / n
        fit <- rq(temp ~ Z, quant, method = "fn")
        Q <- fit$fitted.values
        dt[, j] <- (temp > Q) * 1
      }
    }
  } else {
    so@data[[so@data.name]] <- as(so@data[[so@data.name]], "dgCMatrix")
    out <- list()
    for (j in 1:p) {
      temp <- as.vector(so@data[[so@data.name]][, j])
      v5 <- quantile(temp, 0.5)
      if ((v5 == 0) | (v5 >= quantile(temp, 0.99))) {
        out[[j]] <- which(temp > 0) - 1
      } else {
        quant <- sum(temp <= v5) / n
        fit <- rq(temp ~ Z, quant, method = "fn")
        Q <- fit$fitted.values
        out[[j]] <- which(temp > Q) - 1
      }
    }
    temp <- sapply(out, length)
    dt <- new("dgCMatrix",
      i = as.integer(unlist(out)),
      p = as.integer(c(0, cumsum(temp))),
      x = rep(1, length(as.integer(unlist(out)))),
      Dim = c(n, p)
    )
  }
  so@data.thres <- list(dt)
  names(so@data.thres) <- so@data.name
  return(so)
}
