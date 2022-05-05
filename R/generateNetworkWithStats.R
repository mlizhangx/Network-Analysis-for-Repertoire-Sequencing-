# Main function -----------------------------------------------------------

# Input bulk rep-seq clonotype data; construct repertoire network based on
# desired distance metric, with node and cluster level network characteristics
generateNetworkWithStats <- function(
  data,      # data frame containing req-seq data
  clone_col, # name or number of column of `data` containing clone sequences
  count_col, # name or number of column of `data` containing clone counts
  frac_col,  # name or number of column of `data` containing clone fractions
  dist_type = "hamming", # supports "levenshtein", "hamming", "euclidean_on_atchley"
  edge_dist = 1, # maximum dist threshold for network edges/adjacency
  min_seq_length = 3, # minimum clone sequence length to include
  aggregate_counts = FALSE, # aggregate counts/fracs by unique clone seq?
  group_col = NULL, # optional name or number of column of `data` containing additional grouping variable, if aggregating counts
  output_dir = getwd() # if NULL, output is not saved to file
) {
  # Create output directory if applicable
  if (!is.null(output_dir)) { .createOutputDir(output_dir) }

  # Aggregate counts by unique clone seq if applicable
  if (aggregate_counts) {
    data <- aggregateCountsByAminoAcidSeq(
      data, clone_col, count_col, frac_col, group_col)
    clone_col <- "cloneSeq"
    count_col <- "cloneCount"
  }
  # Remove sequences below  specified
  data <- filterDataBySequenceLength(data, clone_col,
                                     min_length = min_seq_length)
  # stop if not enough unique clonotype sequences
  if (nrow(data) < 2) { stop("Insufficient clone sequences to build network (at least two are needed).") }

  # Generate adjacency matrix for network
  adjacency_matrix <-
    generateNetworkFromClones(data[[clone_col]], dist_type, edge_dist,
                              contig_ids = rownames(data),
                              return_type = "adjacency_matrix")
  if (dist_type != "euclidean_on_atchley") {
    # Subset of data corresponding to the nodes in the network
    data <- data[dimnames(adjacency_matrix)[[1]], ]
  }

  # Generate network from adjacency matrix
  net <- generateNetworkFromAdjacencyMat(adjacency_matrix)

  # Compute node-level network characteristics
  node_data <- computeNodeNetworkStats(net, data)

  # Compute cluster-level network characteristics
  cluster_data <- computeClusterNetworkStats(
    adjacency_matrix, node_data, clone_col, count_col, deg_col = "deg")

  # Plot of Network Graph
  if (!is.null(group_col)) {
    color_nodes_by <- node_data[[group_col]]
    if (is.numeric(group_col)) {
      color_legend_title <- names(node_data)[[group_col]]
    } else { color_legend_title <- group_col }
  } else {
    color_nodes_by <- node_data$deg
    color_legend_title <- "Network degree" }
  graph_plot <- plotNetworkGraph(
    net,
    edge_width = 0.3,
    title = paste0("Network based on distance type: ", dist_type),
    subtitle = paste0("Max edge distance: ", edge_dist),
    color_nodes_by = color_nodes_by,
    size_nodes_by = node_data[[count_col]],
    color_legend_title = color_legend_title,
    size_legend_title = "Clone count") +
    ggplot2::scale_size(range =
                          c(0.1, log(max(node_data[[count_col]])) / 2.5))

  # Write results to disk
  if (!is.null(output_dir)) {
    .saveNetworkResults(node_data, cluster_data, net, adjacency_matrix,
                        graph_plot, output_dir, dist_type, edge_dist,
                        min_seq_length, aggregate_counts, group_col)
  }
  cat("All tasks complete.\n")
  return(list("settings" = list("distance_type" = dist_type,
                                "edge_dist" = edge_dist,
                                "min_seq_length" = min_seq_length,
                                "aggregate_counts" = aggregate_counts,
                                "group_col" = group_col),
              "node_data" = node_data,
              "cluster_data" = cluster_data,
              "adjacency_matrix" = adjacency_matrix,
              "network_graph" = net,
              "graph_plot" = graph_plot))
}



# Helper functions --------------------------------------------------------



.saveNetworkResults <- function(node_data, cluster_data, net, adjacency_matrix,
                                graph_plot, output_dir, dist_type, edge_dist,
                                min_seq_length, aggregate_counts, group_col) {

  cat(paste0("Saving results to ", output_dir, "\n"))

  # Save settings used to generate network
  settings <- data.frame("distance_type" = dist_type,
                         "edge_dist" = edge_dist,
                         "min_seq_length" = min_seq_length,
                         "aggregate_counts" = aggregate_counts)
  if (!is.null(group_col)) settings$group_col <- group_col
  utils::write.csv(settings,
                   file = file.path(output_dir, "settings.csv"),
                   row.names = FALSE)

  # Save Network igraph using edgelist format
  igraph::write_graph(net,
                      file = file.path(output_dir, "network_graph_edgelist.txt"),
                      format = "edgelist")
  cat("Network igraph object saved in edgelist format as 'network_graph_edgelist.txt'\n")

  # Save cell-level info with network characteristics
  utils::write.csv(node_data,
                   file = file.path(output_dir, "node_data.csv"),
                   row.names = FALSE)
  cat("Node-level data and network characteristics saved as 'node_data.csv'\n")

  # Save cluster-level stats
  utils::write.csv(cluster_data,
                   file = file.path(output_dir, "cluster_data.csv"),
                   row.names = FALSE)
  cat("Cluster-level network characteristics saved as 'cluster_data.csv'\n")

  # Save network graph plot
  grDevices::pdf(file.path(output_dir, "network_graph_plot.pdf"),
                 width = 12, height = 8)
  print(graph_plot)
  grDevices::dev.off()
  cat("Plot of network graph saved as 'network_graph_plot.pdf'\n")

  # Save adjacency matrix
  if (dist_type == "euclidean_on_atchley") {
    utils::write.csv(adjacency_matrix,
                     file.path(output_dir, "adjacency_matrix.csv"),
                     row.names = FALSE)
    cat("Adjacency matrix saved as 'adjacency_matrix.csv'\n")
  } else {
    Matrix::writeMM(adjacency_matrix,
                    file.path(output_dir, "adjacency_matrix.mtx"))
    cat("Adjacency matrix saved in sparse format as 'adjacency_matrix.mtx'\n")
  }

}

# # input arguments to buildNetwork other than 'data'
# # check argument types for invalid inputs
# .checkArguments <- function(clonotypes, counts, frequencies, group_labels,
#                             sample_name, distance_type, clone_seq_type) {
#
#   # clonotypes is a character vector
#   if (length(clonotypes) == 0) stop("argument 'clonotypes' has zero length")
#   if (!is.vector(clonotypes) | !is.character(clonotypes)) {
#     stop(paste0("argument 'clonotypes' must be a character vector"))
#   }
#
#   # counts is numeric vector with no NA/Inf values and correct length
#   if (length(counts) == 0) stop("argument 'counts' has zero length")
#   if (!is.vector(counts) | !is.numeric(counts)) {
#     stop(paste0("argument 'counts' must be a numeric vector"))
#   }
#   if (sum(is.na(counts)) != 0) stop("'counts' contains NA/NaN values")
#   if (sum(is.finite(counts)) != length(counts)) {
#     stop("'counts' contains Inf/-Inf values")
#   }
#   if (length(counts) != length(clonotypes)) {
#     stop("length of 'counts' does not match length of 'clonotypes'")
#   }
#
#   # frequency is numeric vector with no NA/Inf values and correct length
#   if (length(frequencies) == 0) stop("argument 'frequencies' has zero length")
#   if (!is.vector(frequencies) | !is.numeric(frequencies)) {
#     stop(paste0("argument 'frequencies' must be a numeric vector"))
#   }
#   if (sum(is.finite(frequencies)) != length(frequencies)) {
#     stop("'frequencies' contains Inf/-Inf values")
#   }
#   if (length(frequencies) != length(clonotypes)) {
#     stop("length of 'frequencies' does not match length of 'clonotypes'")
#   }
#
#   # group_labels is null or vector of correct length, coercible to factor with > 1 level
#   if (!is.null(group_labels)) {
#     if (length(group_labels) == 0) stop("argument 'group_labels' has zero length; if you are trying to omit the group variable, use the default value of `group_labels = NULL`. Otherwise, 'group_labels' must be a vector of the same length as 'clonotypes' containing at least two distinct values")
#     if (!is.vector(group_labels) ) stop("'group_labels' must be a vector or 'NULL'")
#     if (length(group_labels) != length(clonotypes)) stop("length of 'group_labels' does not match length of 'clonotypes'")
#     if (sum(is.na(group_labels)) != 0) stop("'group_labels' contains NA/NaN values")
#     if (length(levels(as.factor(group_labels))) < 2) stop("'group_labels' must contain at least two distinct values")
#   }
#
#   # Valid distance_type arguments
#   programmed_distance_types <- c("levenshtein", "hamming")
#
#   # Check distance_type argument
#   if (length(distance_type) != 1) {
#     stop(paste0("'distance_type' must be one of the following: ",
#                 programmed_distance_types))
#   }
#   if (!distance_type %in% programmed_distance_types) {
#     stop(paste0("'distance_type' must be one of the following: ",
#                 programmed_distance_types))
#   }
#
#   # Check sample_name argument
#   if (length(sample_name) != 1) {
#     stop("argument 'sample_name' must be a character string")
#   }
#   if (!is.character(sample_name)) {
#     stop("argument 'sample_name' must be a character string")
#   }
#
#   # Check clone_type argument
#   if (length(clone_seq_type) != 1) {
#     stop("argument 'clone_seq_type' must be a character string")
#   }
#   if (!is.character(clone_seq_type)) {
#     stop("argument 'clone_seq_type' must be a character string")
#   }
#
#
# }






