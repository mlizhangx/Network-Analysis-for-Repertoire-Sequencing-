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

  # Input Data and Columns
  data,      # data frame containing req-seq data
  nucleo_col,
  amino_col,
  count_col,
  freq_col,
  vgene_col, # ignored if aggregate_reads = TRUE
  dgene_col, # ignored if aggregate_reads = TRUE
  jgene_col, # ignored if aggregate_reads = TRUE
  cdr3length_col, # ignored if aggregate_reads = TRUE
  other_cols = NULL, # other cols to keep; ignored if aggregate_reads = TRUE

  # Clone Sequence Settings
  clone_seq_type = "amino acid", # or "nucleotide"
  min_seq_length = 3, # min clone seq length
  drop_chars = NULL, # regular expression, e.g. "[*|_]"
  aggregate_reads = FALSE,
  grouping_cols = NULL,

  # Network Settings
  dist_type = "hamming", # or "levenshtein", "hamming", "euclidean_on_atchley"
  edge_dist = 1, # max dist for edges

  # Network Stats
  node_stats = FALSE,
  cluster_stats = FALSE,
  node_stat_settings = node_stat_settings(cluster_id = cluster_stats),

  # Plot Settings
  plot_title = paste("RepSeq network by", dist_type, "distance on CDR3",
                     clone_seq_type, "sequence"),
  plot_subtitle = ifelse(dist_type == "euclidean_on_atchley",
                         yes = paste("Clone sequences embedded in Euclidean 30-space based on Atchley factor representation using deep learning\nEdges based on a maximum Euclidean distance of", edge_dist, "between embedded values\n"),
                         no = paste("Edges based on a maximum", dist_type, "distance of", edge_dist, "\n")),
  size_nodes_by = count_col, # can use a double, e.g., 1.0, for fixed size
  node_size_limits = NULL, # numeric, length 2
  color_nodes_by = NULL, # use NULL to automatically determine
  color_scheme = "default",

  # Output Settings
  output_dir = getwd(), # if NULL, output is not saved to file
  plot_outfile = "network_graph.pdf",
  data_outfile = "node_data.csv",
  igraph_outfile = NULL, # .txt
  matrix_outfile = NULL, # .mtx (.csv for euclidean on atchley)
  return_all = FALSE # if false, only the node data is returned (unless cluster_stats = TRUE and output_dir = NULL, in which case a list containing the node data and cluster info is returned, with a warning)

) {

  #### INPUT CHECKS ####

  # Atchley factor embedding only applicable to amino acid sequences
  if (dist_type == "euclidean_on_atchley" & clone_seq_type != "amino acid") {
    stop("distance type 'euclidean_on_atchley' only applicable to amino acid sequences") }


  #### PREPARE WORKING ENVIRONMENT ####
  # Create output directory if applicable
  if (!is.null(output_dir)) { .createOutputDir(output_dir) }

  # return type
  return_type <- ifelse(return_all,
                        yes = "all",
                        no = ifelse(cluster_stats & is.null(output_dir),
                                    yes = "node_and_cluster_data",
                                    no = "node_data_only"))

  # Convert input columns to character if not already
  if (is.numeric(nucleo_col)) { nucleo_col <- names(data)[nucleo_col] }
  if (is.numeric(amino_col)) { amino_col <- names(data)[amino_col] }
  if (is.numeric(count_col)) { count_col <- names(data)[count_col] }
  if (is.numeric(freq_col)) { freq_col <- names(data)[freq_col] }
  if (!is.null(color_nodes_by)) {
    if (is.numeric(color_nodes_by)) {
      color_nodes_by <- names(data)[color_nodes_by] } }
  if (is.integer(size_nodes_by)) { size_nodes_by <- names(data)[size_nodes_by] }
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
  if (clone_seq_type == "nucleotide") { clone_seq_col <- nucleo_col }


  #### FORMAT AND FILTER DATA ####
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
  if (nrow(out) < 2) { stop("fewer than two clone sequences meet the specified minimum length") }

  # Drop sequences with specified chars
  if (!is.null(drop_chars)) {
    data <- data[-grep(drop_chars, data[ , clone_seq_col]), ] }


  #### BUILD NETWORK ####
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


  #### NODE/CLUSTER STATS ####
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


  ### PLOT(S) OF NETWORK GRAPH ####
  # if color_nodes_by is NULL, determine default color variable
  if (is.null(color_nodes_by)) {
    if ("cluster_id" %in% names(data)) { # use cluster ID if present
      color_nodes_by <- "cluster_id"
      color_legend_title <- "cluster"
    } else if ("degree" %in% names(data)) { # else try network degree
      color_nodes_by <- "degree"
      color_legend_title <- "degree"
    } else { # if all else fails, just color the nodes by clone count
      color_nodes_by <- count_col
      color_legend_title <- "clone count" }
  } else { # if color variables specified, add cluster ID if applicable
    if ("cluster_id" %in% names(data) & !"cluster_id" %in% color_nodes_by) {
      color_nodes_by <- c(color_nodes_by, "cluster_id") } }

  # If multiple coloring variables, extend color scheme to vector if needed
  if (length(color_nodes_by) > 1 & length(color_scheme) == 1) {
    color_scheme <- rep(color_scheme, length(color_nodes_by)) }

  # Ensure size_nodes_by is a vector or fixed value to use for node sizes
  size_legend_title <- NULL
  if (is.character(size_nodes_by)) {
    size_legend_title <- size_nodes_by
    size_nodes_by <- data[ , size_nodes_by] }

  # Create one plot for each variable used to color the nodes
  temp_plotlist <- list()
  for (j in 1:length(color_nodes_by)) {
    cat(paste0("Generating graph plot with nodes colored by ",
               color_nodes_by[[j]], "..."))
    temp_plotlist$newplot <-
      plotNetworkGraph(
        net, title = plot_title, subtitle = plot_subtitle,
        color_nodes_by = data[ , color_nodes_by[[j]]],
        size_nodes_by = size_nodes_by,
        color_legend_title = color_nodes_by[[j]],
        size_legend_title = size_legend_title,
        color_scheme = color_scheme[[j]])
    if (!is.double(size_nodes_by) & !is.null(node_size_limits)) {
      temp_plotlist$newplot <- temp_plotlist$newplot +
        ggplot2::scale_size(range = node_size_limits) }
    print(temp_plotlist$newplot) # print to R
    names(temp_plotlist)[[length(names(temp_plotlist))]] <- color_nodes_by[[j]]
    cat("Done.\n") }


  #### SAVE RESULTS ####
  if (!aggregate_reads) {  # Rename data columns
    colnames(data)[1:8] <- c(
      "nucleotideSeq", "aminoAcidSeq", "cloneCount", "cloneFrequency",
      "VGene", "DGene", "JGene", "CDR3Length") }

  # Save node [& cluster] data
  if (!is.null(output_dir)) {
    if (!is.null(data_outfile)) {
      utils::write.csv(data, file = file.path(output_dir, data_outfile),
                       row.names = FALSE)
      cat(paste0("Node-level data saved to file:\n", data_outfile, "\n")) }
    if (cluster_stats) {
      utils::write.csv(cluster_info, file = file.path(output_dir,
                                                      "cluster_info.csv"))
      cat(paste0("Cluster-level data saved to file:\n",
                 file.path(output_dir, "cluster_info.csv"), "\n")) } }

  # Save plots to a single pdf
  if (!is.null(output_dir) & !is.null(plot_outfile)) {
    grDevices::pdf(file = file.path(output_dir, plot_outfile))
    for (j in 1:length(color_nodes_by)) { print(temp_plotlist[[j]]) }
    grDevices::dev.off() }

  # Save igraph
  if (!is.null(output_dir) & keep_igraph) {
    igraph::write_graph(net,
                        file = file.path(output_dir, "network_edge_list.txt"),
                        format = "edgelist")
    cat(paste0("Network igraph saved in edgelist format to file:\n",
               file.path(output_dir, "network_edge_list.txt"), "\n")) }

  # Save adjacency matrix
  if (!is.null(output_dir) & keep_matrix) {
    if (dist_type == "euclidean_on_atchley") {
      utils::write.csv(adjacency_matrix,
                       file.path(output_dir, "adjacency_matrix.csv"),
                       row.names = FALSE)
      cat(paste0("Adjacency matrix saved to file:\n",
                 file.path(output_dir, "adjacency_matrix.csv"), "\n"))
    } else {
      Matrix::writeMM(adjacency_matrix,
                      file.path(output_dir, "adjacency_matrix.mtx"))
      cat(paste0("Adjacency matrix saved to file:\n",
                 file.path(output_dir, "adjacency_matrix.mtx"), "\n")) } }


  #### RETURN OUTPUT ####
  if (return_type == "node_data_only") {

    cat("All tasks complete.\n")
    return(data)

  } else {
    out <- list("node_data" = data)
    if (cluster_stats) { out$cluster_info <- cluster_info }
    if (length(temp_plotlist == 1)) { out$plot <- temp_plotlist[[1]]
    } else { out$plot <- temp_plotlist }
    if (keep_igraph) { out$igraph <- net }
    if (keep_matrix) { out$adjacency_matrix <- adjacency_matrix }
    cat("All tasks complete.\n")
    return(out) }

}



# Helper functions --------------------------------------------------------



# .saveNetworkResults <- function(node_data, cluster_data, net, adjacency_matrix,
#                                 graph_plot, output_dir, dist_type, edge_dist,
#                                 min_seq_length, aggregate_counts, group_col) {
#
#   cat(paste0("Saving results to ", output_dir, "\n"))
#
#   # Save settings used to generate network
#   settings <- data.frame("distance_type" = dist_type,
#                          "edge_dist" = edge_dist,
#                          "min_seq_length" = min_seq_length,
#                          "aggregate_counts" = aggregate_counts)
#   if (!is.null(group_col)) settings$group_col <- group_col
#   utils::write.csv(settings,
#                    file = file.path(output_dir, "settings.csv"),
#                    row.names = FALSE)
#
#   # Save Network igraph using edgelist format
#   igraph::write_graph(net,
#                       file = file.path(output_dir, "network_graph_edgelist.txt"),
#                       format = "edgelist")
#   cat("Network igraph object saved in edgelist format as 'network_graph_edgelist.txt'\n")
#
#   # Save cell-level info with network characteristics
#   utils::write.csv(node_data,
#                    file = file.path(output_dir, "node_data.csv"),
#                    row.names = FALSE)
#   cat("Node-level data and network characteristics saved as 'node_data.csv'\n")
#
#   # Save cluster-level stats
#   utils::write.csv(cluster_data,
#                    file = file.path(output_dir, "cluster_data.csv"),
#                    row.names = FALSE)
#   cat("Cluster-level network characteristics saved as 'cluster_data.csv'\n")
#
#   # Save network graph plot
#   grDevices::pdf(file.path(output_dir, "network_graph_plot.pdf"),
#                  width = 12, height = 8)
#   print(graph_plot)
#   grDevices::dev.off()
#   cat("Plot of network graph saved as 'network_graph_plot.pdf'\n")
#
#   # Save adjacency matrix
#   if (dist_type == "euclidean_on_atchley") {
#     utils::write.csv(adjacency_matrix,
#                      file.path(output_dir, "adjacency_matrix.csv"),
#                      row.names = FALSE)
#     cat("Adjacency matrix saved as 'adjacency_matrix.csv'\n")
#   } else {
#     Matrix::writeMM(adjacency_matrix,
#                     file.path(output_dir, "adjacency_matrix.mtx"))
#     cat("Adjacency matrix saved in sparse format as 'adjacency_matrix.mtx'\n")
#   }
#
# }

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






