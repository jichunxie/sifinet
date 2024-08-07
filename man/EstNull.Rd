% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{EstNull}
\alias{EstNull}
\title{EstNull
This function is a Rcpp version of Wenguang Sun and Tony T. Cai's
EstNull.func R function, estimating null distribution from data.
Sun, W., & Cai, T. T. (2007). Oracle and Adaptive Compound Decision
Rules for False Discovery Rate Control.
Journal of the American Statistical Association,
102(479), 901–912.}
\usage{
EstNull(x, gamma = 0.1)
}
\arguments{
\item{x}{Input vector of all coexpression values}

\item{gamma}{Parameter setting the stopping threshold}
}
\value{
List of mean and std
}
\description{
EstNull
This function is a Rcpp version of Wenguang Sun and Tony T. Cai's
EstNull.func R function, estimating null distribution from data.
Sun, W., & Cai, T. T. (2007). Oracle and Adaptive Compound Decision
Rules for False Discovery Rate Control.
Journal of the American Statistical Association,
102(479), 901–912.
}
\author{
Qi Gao
}
