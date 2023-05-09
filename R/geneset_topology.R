#' geneset_topology
#' 
#' The function plots the topology network of the feature gene sets found by SiFINeT.
#' @param so a SiFINeT object
#' @param weightthres edges between nodes (feature gene sets) with weight greater than weightthres would be shown in the plot
#' @param edge_method SiFINeT provides 2 methods of calculating edge weight. 
#' The number of shared feature genes between feature gene sets would be used when edge_method = 1; 
#' while the edge proportion between feature gene sets would be applied if edge_method = 2.
#' @param node_color color of nodes. Should have either length 1 or same length as the number of feature gene sets.
#' @param shiftsize set the distance between center of label and the corresponding feature gene sets node. 
#' @param boundsize set the size of the boundary region. 
#' @param prefix the prefix of the labels
#' @param set_name name of the gene sets
#' @return A ggraph (ggplot) object
#' @details This function visualizes the output feature gene sets of SiFINeT in the form of network.
#' Number of shared feature genes or proportion of edges between feature gene sets could be used to weight the edges.
#' The layout of the nodes is created by create_layout function in ggraph package.
#' @references Thomas Lin Pedersen (2022). ggraph: An Implementation of 
#' Grammar of Graphics for Graphs and Networks. R package 
#' version 2.0.6. https://CRAN.R-project.org/package=ggraph
#' @import ggraph
#' @import ggplot2
#' @import igraph
#' @export
#' 
geneset_topology <- function(so, weightthres = 0.3, edge_method = 2, node_color = "black",
                             shiftsize = 0.05, boundsize = 0.3, prefix = "", set_name = NULL){
  if ((edge_method == 1) & (weightthres < 1)){
    weightthres <- 5
  }
  id_list <- list()
  for (i in 1:length(so@featureset$unique)){
    id_list[[i]] <- match(c(so@featureset$unique[[i]], 
                            so@featureset$shared[[i]], 
                            so@featureset$enriched[[i]]), so@gene.name)
  }
  
  nodes <- data.frame(1:length(id_list), sapply(id_list, length))
  if (is.null(set_name)){
    nodes[,1] <- paste(prefix, nodes[,1], sep = "")
  } else{
    nodes[,1] <- set_name
  }
  edge <- data.frame(t(combn(1:length(id_list), 2)))
  if (edge_method == 1){
    edge$edge_weight <- apply(edge, 1, function(x){length(intersect(id_list[[x[1]]], id_list[[x[2]]]))})
  } else {
    edge_mat <- 1 * (abs(so@coexp - so@est_ms$mean) >= so@thres)
    edge$edge_weight <- apply(edge, 1, function(x){mean(edge_mat[id_list[[x[1]]], id_list[[x[2]]]])})
  }
  if (is.null(set_name)){
    edge[,1] <- paste(prefix, edge[,1], sep = "")
    edge[,2] <- paste(prefix, edge[,2], sep = "")
  } else{
    edge[,1] <- set_name[edge[,1]]
    edge[,2] <- set_name[edge[,2]]
  }
  edge <- edge[edge$edge_weight > weightthres, ]
  me <- max(edge$edge_weight)
  edge$scaled_edge_weight <- edge$edge_weight / me
  colnames(nodes) <- c("node_id", "gene_count")
  
  g <- graph_from_data_frame(edge, directed = FALSE, nodes)
  
  lay_temp <- create_layout(g, layout = "stress")
  shift <- shiftsize * (max(lay_temp$y) - min(lay_temp$y) + max(lay_temp$x) - min(lay_temp$x))
  bound <- boundsize * (max(lay_temp$y) - min(lay_temp$y) + max(lay_temp$x) - min(lay_temp$x))
  gt <- ggraph(lay_temp) + 
    geom_edge_link(aes(width = scaled_edge_weight))  + 
    scale_edge_width(range=c(0.5,2)) +
    geom_node_point(aes(size = 3), shape = 16, color = node_color) +
    geom_label(aes(x = x, y = y + shift, label = name)) + 
    theme_void() + 
    theme(legend.position="none", 
          legend.text=element_blank(), 
          legend.title=element_blank()) + 
    xlim(c(min(lay_temp$x) - bound, 
           max(lay_temp$x) + bound)) + 
    ylim(c(min(lay_temp$y) - bound, 
           max(lay_temp$y) + bound))
  return(gt)
}