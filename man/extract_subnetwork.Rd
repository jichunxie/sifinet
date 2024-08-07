% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_subnetwork.R
\name{extract_subnetwork}
\alias{extract_subnetwork}
\title{extract_subnetwork}
\usage{
extract_subnetwork(
  so,
  target_gene_name = NULL,
  target_gene_id = NULL,
  positive = TRUE
)
}
\arguments{
\item{so}{a SiFINeT object}

\item{target_gene_name}{the names of the target genes in the output network}

\item{target_gene_id}{the indices of the target genes in the output network, not used when target_gene_name is not Null}

\item{positive}{whether only positive (default) co-expressions or all co-expressions are considered in assigning edges}
}
\value{
an adjacency matrix of the output subnetwork
}
\description{
The function extract a subnetwork from the co-expression network
}
