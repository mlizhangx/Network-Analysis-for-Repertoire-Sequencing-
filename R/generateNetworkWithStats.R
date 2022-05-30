# Main function -----------------------------------------------------------

# INPUT:
#   bulk rep-seq data
#   network specifications
# DO:
#   aggregate rows by unique clone seq (if desired)
#   filter clone seqs by minimum length
#   build network to spec (dist type, edge dist)
#   compute network stats (if desired)
#   compute clusters (if desired)
#   generate network graph plot

buildRepSeqNetwork <- function(
  data,      # data frame containing req-seq data
  nucleo_col,
  amino_col,
  count_col,
  freq_col,
  vgene_col,
  dgene_col,
  jgene_col,
  cdr3length_col,
  other_cols = NULL,
  clone_seq_type = "amino_acid",
  dist_type = "hamming", # supports "levenshtein", "hamming", "euclidean_on_atchley"
  edge_dist = 1, # maximum dist threshold for network edges/adjacency
  min_seq_length = 3, # minimum clone sequence length to include
  size_nodes_by = count_col,
  color_nodes_by = NULL,
  color_scheme = "default",
  node_stats = FALSE,
  cluster_stats = FALSE,
  node_stat_settings = node_stat_settings(cluster_id = cluster_stats),
  aggregate_reads = FALSE,
  grouping_cols = NULL,
  keep_igraph = TRUE,
  keep_matrix = FALSE,
  output_dir = NULL # if NULL, output is not saved to file
) {
  ### INPUT CHECKS ###
  # Atchley factor embedding only applicable to amino acid sequences
  if (dist_type == "euclidean_on_atchley" & clone_seq_type != "amino_acid") {
    stop("distance type 'euclidean_on_atchley' only applicable to amino acid sequences") }


  ### PREPARE WORKING ENVIRONMENT ###
  # Create output directory if applicable
  if (!is.null(output_dir)) { .createOutputDir(output_dir) }

  # Convert input columns to character if not already
  if (is.numeric(nucleo_col)) { nucleo_col <- names(data)[nucleo_col] }
  if (is.numeric(amino_col)) { amino_col <- names(data)[amino_col] }
  if (is.numeric(count_col)) { count_col <- names(data)[count_col] }
  if (is.numeric(freq_col)) { freq_col <- names(data)[freq_col] }
  if (is.numeric(size_nodes_by)) { size_nodes_by <- names(data)[size_nodes_by] }
  if (!is.null(color_nodes_by)) {
    if (is.numeric(color_nodes_by)) {
      color_nodes_by <- names(data)[color_nodes_by] } }
  if (!aggregate_counts) {
    if (is.numeric(vgene_col)) { vgene_col <- names(data)[vgene_col] }
    if (is.numeric(dgene_col)) { dgene_col <- names(data)[dgene_col] }
    if (is.numeric(jgene_col)) { jgene_col <- names(data)[jgene_col] }
    if (is.numeric(other_cols)) { other_cols <- names(data)[other_cols] }
  } else {
    if (is.numeric(grouping_cols)) {
      grouping_cols <- names(data)[grouping_cols] }
  }

  # Add columns in color_nodes_by to other_cols if not present
  other_cols <- unique(c(other_cols, color_nodes_by))

  # Designate amino acid or nucleotide for clone sequence
  clone_seq_col <- amino_col
  if (clone_seq_type %in% c("nucleo", "nucleotide")) {
    clone_seq_col <- nucleo_col }


  ### FORMAT DATA ###
  if (aggregate_reads) { # Aggregate the counts if specified
    data <- aggregateReads(data, clone_seq_col,
                           count_col, freq_col, grouping_cols)
  } else { # Copy the relevant columns from the input data
    data <-
      data[ , # Keep only the relevant columns:
            c(nucleo_col, amino_col, count_col, freq_col,
              vgene_col, dgene_col, jgene_col, cdr3length_col, other_cols)] }

  # Remove sequences below specified length
  data <- filterDataBySequenceLength(data, clone_seq_col,
                                     min_length = min_seq_length)
  if (nrow(out) < 2) { stop("Insufficient clone sequences to build network (at least two are needed).") }


  ### BUILD NETWORK ###
  # Generate adjacency matrix for network
  adjacency_matrix <-
    generateNetworkFromClones(data[ , clone_seq_col],
                              dist_type, edge_dist,
                              contig_ids = rownames(data),
                              return_type = "adjacency_matrix")

  # Subset data to keep only those clones in the network (nonzero degree)
  if (dist_type != "euclidean_on_atchley") {
    data <- data[dimnames(adjacency_matrix)[[1]], ] }

  # Generate network from adjacency matrix
  net <- generateNetworkFromAdjacencyMat(adjacency_matrix)


  ### NODE/CLUSTER STATS ###
  # Add node-level network characteristics
  if (node_stats) {
    data <- addNodeNetworkStats(data, net, node_stat_settings) }

  # Compute cluster-level network characteristics
  if (cluster_stats) {
    cluster_info <- getClusterStats(
      data, adjacency_matrix, clone_seq_col, count_col,
      cluster_id_col = ifelse("cluster_id" %in% names(data),
                              yes = "cluster_id", no = NULL),
      degree_col = ifelse("degree" %in% names(data),
                          yes = "degree", no = NULL)) }


  ## PLOT(S) OF NETWORK GRAPH ###
  # default variable to color nodes by
  if (is.null(color_nodes_by)) {
    if ("cluster_id" %in% names(data)) {
      color_nodes_by <- "cluster_id"
      color_legend_title <- "cluster"
    } else if ("degree" %in% names(data)) {
      color_nodes_by <- "degree"
      color_legend_title <- "degree"
    } else {
      color_nodes_by <- count_col
      color_legend_title <- "clone count" }
  } else {
    if ("cluster_id" %in% names(data) & !"cluster_id" %in% color_nodes_by) {
      color_nodes_by <- c(color_nodes_by, "cluster_id")
    }
  }

  # If multiple coloring variables, extend color scheme to vector if needed
  if (length(color_nodes_by) > 1 & length(color_scheme) == 1) {
    color_scheme <- rep(color_scheme, length(color_nodes_by)) }

  # If multiple plots (coloring variables), create list to store plots
  if (length(color_nodes_by) > 1) { plots <- list() }

  graph_plot <- plotNetworkGraph(
    net,
    edge_width = 0.3,
    title = paste0("Network based on distance type: ", dist_type),
    subtitle = paste0("Max edge distance: ", edge_dist),
    color_nodes_by = color_nodes_by,
    size_nodes_by = node_data[ , count_col],
    color_legend_title = color_legend_title,
    size_legend_title = "Clone count") +
    ggplot2::scale_size(range =
                          c(0.1, log(max(node_data[ , count_col])) / 2.5))
  print(graph_plot)

  # # # #

  temp_plotlist <- list()
  for (j in 1:length(color_nodes_by)) {
    cat(paste0("Creating cluster graph with nodes colored by ",
               color_nodes_by[[j]], "..."))
    newplot <-
      plotNetworkGraph(
        network, title = plot_title,
        subtitle = paste0("Nodes colored by ", color_nodes_by[[j]]),
        color_nodes_by = data_current_cluster[ , color_nodes_by[[j]]],
        size_nodes_by = data_current_cluster[ , freq_col],
        color_legend_title = color_nodes_by[[j]],
        size_legend_title = "freq in own sample",
        color_scheme = color_scheme[[j]])
    print(newplot) # print to R
    temp_plotlist$newplot <- newplot
    names(temp_plotlist)[[length(names(temp_plotlist))]] <- color_nodes_by[[j]]
    cat("Done.\n")
  }
  if (save_plots & !is.null(output_dir)) {
    grDevices::pdf(file = file.path(output_dir,
                                    paste0("cluster_", i, ".pdf")))
    for (j in 1:length(color_nodes_by)) { print(temp_plotlist[[j]]) }
    grDevices::dev.off()
  }
  if (return_plots) {
    plots$newcluster <- temp_plotlist
    names(plots)[[length(names(plots))]] <- selected_clones[[i]]
  }

  ### SAVE RESULTS ###
  if (!aggregate_reads) {  # Rename data columns
    colnames(data)[1:8] <- c(
      "nucleotideSeq", "aminoAcidSeq", "cloneCount", "cloneFrequency",
      "VGene", "DGene", "JGene", "CDR3Length") }

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






