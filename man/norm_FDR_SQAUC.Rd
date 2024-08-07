% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/norm_FDR_SQAUC.R
\name{norm_FDR_SQAUC}
\alias{norm_FDR_SQAUC}
\title{norm_FDR_SQAUC
The function controls the false discovery rate (FDR) of coexpression patterns
using SQAUC method
Jichun Xie and Ruosha Li. "False discovery rate control for high dimensional
networks of quantile associations conditioning on covariates".
In: J R Stat Soc Series B Stat Methodol (2018). doi: 10.1111/rssb.12288.}
\usage{
norm_FDR_SQAUC(value, sam_mean, sam_sd, alpha, n, p)
}
\arguments{
\item{value}{the vector of coexpression patterns}

\item{sam_mean}{the estimated sample mean}

\item{sam_sd}{the estimated sample sd}

\item{alpha}{the type I error rate}

\item{n}{the number of cells}

\item{p}{the number of genes}
}
\value{
lower bound threshold for genes to be significantly coexpressed
}
\description{
norm_FDR_SQAUC
The function controls the false discovery rate (FDR) of coexpression patterns
using SQAUC method
Jichun Xie and Ruosha Li. "False discovery rate control for high dimensional
networks of quantile associations conditioning on covariates".
In: J R Stat Soc Series B Stat Methodol (2018). doi: 10.1111/rssb.12288.
}
