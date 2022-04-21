
# Description -------------------------------------------------------------

# Shared internal helper functions used by multiple package functions




# Shared Helpers ----------------------------------------------------------


# Simple wrapper to igraph functions:
# Use adjacency matrix to generate network graph
genNetworkGraph <- function(adjacency_matrix) {

  cat("Generating network from adjacency matrix...\n")

  net <- igraph::graph_from_adjacency_matrix(adjacency_matrix, weighted = TRUE)

  net <- igraph::as.undirected(igraph::simplify(net,
                                                remove.multiple = T,
                                                remove.loops = T))

  return(net)

}
