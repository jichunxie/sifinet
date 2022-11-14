#' The SiFINeT Class
#'
#' @slot data a list of cell (row) by gene (column) count matrix, either regular or sparse matrix
#' @slot sparse whether the count matrix should be analyzed as sparse matrix
#' @slot meta.data matrix of meta data, the number of rows should equal to the number of cells
#' @slot gene.name a vector of names of genes with length equal to the number of genes 
#' @slot data.name name of the dataset
#' @slot n number of cells in the dataset
#' @slot p number of genes in the dataset
#' @slot data.thres binarized count matrix 
#' @slot coexp matrix of genes coexpression
#' @slot est_ms estimated mean and sd of coexpression values
#' @slot thres lower bound of coexpression (or absolute value of coexpression) for network edge assignment
#' @slot q5 50% quantile for each gene
#' @slot kset index of kept genes after the filtering step
#' @slot conn list of connectivities in absolute network
#' @slot conn2 list of connectivities in positive sub-network
#' @slot fg_id index of the candidate feature genes
#' @slot uni_fg_id index of the candidate unique feature genes
#' @slot uni_cluster cluster result of the candidate unique feature genes
#' @slot selected_cluster selected unique feature gene clusters
#' @slot featureset detected set of feature genes
#' 
#' @name SiFINeT-class
#' @exportClass SiFINeT
#'
SiFINeT <- setClass(
  Class = 'SiFINeT',
  slots = c(
    data = 'list',
    sparse = 'logical',
    meta.data = 'matrix',
    gene.name = 'vector',
    data.name = 'character',
    n = 'numeric',
    p = 'numeric',
    data.thres = 'list',
    coexp = 'matrix',
    est_ms = 'list',
    thres = 'numeric',
    q5 = 'numeric',
    kset = 'integer',
    conn = 'list',
    conn2 = 'list',
    fg_id = 'integer',
    uni_fg_id = 'integer',
    uni_cluster = 'numeric',
    selected_cluster = 'numeric',
    featureset = 'list'
  )
)

#' create_SiFINeT_object
#' 
#' The function classifies count data based on thresholds 
#' defined by quantile regression
#' @param counts count matrix
#' @param gene.name name of the features
#' @param meta.data data.frame of meta data
#' @param data.name name of dataset
#' @param sparse whether the count matrix should be analyzed as sparse matrix
#' @param rowfeature whether the count matrix is feature (row) by cell (column) 
#' @return a SiFINeT object
#' @export
#' 
create_SiFINeT_object <- function(counts, gene.name = NULL, 
                                meta.data = NULL, data.name = NULL, 
                                sparse = FALSE, rowfeature = TRUE){
  if (rowfeature == TRUE){
  	counts <- t(counts)
  }
  if (is.null(gene.name)){
    gene.name <- colnames(counts)
  }
  if (is.null(gene.name)){
    gene.name <- 1:ncol(counts)
  }
  if (is.null(data.name)){
    data.name <- "data1"
  }
  if (is.null(meta.data)){
    meta.data <- matrix(0, nrow(counts), 0)
  }
  data <- list(counts)
  names(data) <- data.name
  object <- new(
    Class = 'SiFINeT',
    data = data,
    sparse = sparse,
    meta.data = meta.data,
    gene.name = gene.name,
    data.name = data.name,
    n = nrow(counts),
    p = ncol(counts),
    q5 = apply(counts, 2, quantile, 0.5),
    kset = 1:ncol(counts),
    featureset = list(unique = list(),
                      shared = list(),
                      enriched = list())
  )
  return(object)
}
