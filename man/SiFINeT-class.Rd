% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_SiFINeT_object.R
\docType{class}
\name{SiFINeT-class}
\alias{SiFINeT-class}
\alias{SiFINeT}
\title{The SiFINeT Class}
\description{
The SiFINeT Class
}
\section{Slots}{

\describe{
\item{\code{data}}{a list of cell (row) by gene (column) count matrix, either regular or sparse matrix}

\item{\code{sparse}}{whether the count matrix should be analyzed as sparse matrix}

\item{\code{meta.data}}{matrix of meta data, the number of rows should equal to the number of cells}

\item{\code{gene.name}}{a vector of names of genes with length equal to the number of genes}

\item{\code{data.name}}{name of the dataset}

\item{\code{n}}{number of cells in the dataset}

\item{\code{p}}{number of genes in the dataset}

\item{\code{data.thres}}{binarized count matrix}

\item{\code{coexp}}{matrix of genes coexpression}

\item{\code{est_ms}}{estimated mean and sd of coexpression values}

\item{\code{thres}}{lower bound of coexpression (or absolute value of coexpression) for network edge assignment}

\item{\code{q5}}{50\% quantile for each gene}

\item{\code{kset}}{index of kept genes after the filtering step}

\item{\code{conn}}{list of connectivities in absolute network}

\item{\code{conn2}}{list of connectivities in positive sub-network}

\item{\code{fg_id}}{index of the candidate feature genes}

\item{\code{uni_fg_id}}{index of the candidate unique feature genes}

\item{\code{uni_cluster}}{cluster result of the candidate unique feature genes}

\item{\code{selected_cluster}}{selected unique feature gene clusters}

\item{\code{featureset}}{detected set of feature genes}
}}

