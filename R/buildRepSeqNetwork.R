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
  nucleo_col = "nucleotideSeq",
  amino_col = "aminoAcidSeq",
  count_col = "cloneCount",
  freq_col = "cloneFreqInSample",
  vgene_col = "VGene",# ignored if aggregate_reads = TRUE
  dgene_col = "DGene",# ignored if aggregate_reads = TRUE
  jgene_col = "JGene",# ignored if aggregate_reads = TRUE
  cdr3length_col = "CDR3Length",# ignored if aggregate_reads = TRUE
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
  node_stats = FALSE,
  stats_to_include = node_stat_settings(), # can also select "all" or "cluster_id_only"
  cluster_stats = FALSE,

  # Plot Settings
  plot_title = paste("Network on CDR3", clone_seq_type, "sequence"),
  plot_subtitle = ifelse(dist_type == "euclidean_on_atchley",
                         yes = paste("Clone sequences embedded in Euclidean 30-space based on Atchley factor representation using deep learning\nEdges based on a maximum Euclidean distance of", edge_dist, "between embedded values\n"),
                         no = paste("Edges based on a maximum", dist_type, "distance of", edge_dist, "\n")),
  edge_width = 0.3,
  size_nodes_by = count_col, # can use a double, e.g., 1.0, for fixed size
  node_size_limits = NULL, # numeric, length 2
  custom_size_legend = NULL, # custom legend title
  color_nodes_by = NULL, # use NULL to automatically determine
  color_scheme = "default",
  custom_color_legend = NULL, # custom title (length must match color_nodes_by)

  # Output Settings
  output_dir = NULL, # if NULL, output is not saved to file
  save_all = FALSE, # by default, only save pdf of plot and csv of node data
  data_outfile = "node_data.csv",
  plot_outfile = "network_graph.pdf",
  plot_width = 12, # passed to pdf()
  plot_height = 10, # passed to pdf()
  cluster_outfile = "cluster_info.csv", # only saved if save_all = TRUE
  igraph_outfile = "network_edgelist.txt", # .txt
  matrix_outfile = ifelse(dist_type == "euclidean_on_atchley",
                          yes = "adjacency_matrix.csv",
                          no = "adjacency_matrix.mtx"), # .mtx (.csv for euclidean on atchley)
  return_all = FALSE # if false, only the node data is returned (unless cluster_stats = TRUE and output_dir = NULL, in which case a list containing the node data and cluster info is returned, with a warning)

) {

  #### INPUT CHECKS ####

  # Atchley factor embedding only applicable to amino acid sequences
  if (dist_type == "euclidean_on_atchley" & clone_seq_type != "amino acid") {
    stop("distance type 'euclidean_on_atchley' only applicable to amino acid sequences") }

  # aggregate_reads only applicable to amino acid sequences

  # each elem of color_nodes_by is a data column or a node stat to be added

  # for each elem of other_cols not present in data, warn that not found


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
  # if (is.integer(size_nodes_by)) { size_nodes_by <- names(data)[size_nodes_by] }
  if (!aggregate_reads) {
    if (is.numeric(vgene_col)) { vgene_col <- names(data)[vgene_col] }
    if (is.numeric(dgene_col)) { dgene_col <- names(data)[dgene_col] }
    if (is.numeric(jgene_col)) { jgene_col <- names(data)[jgene_col] }
    if (is.numeric(cdr3length_col)) { cdr3length_col <- names(data)[cdr3length_col] }
    if (is.numeric(other_cols)) { other_cols <- names(data)[other_cols] }
  } else {
    if (is.numeric(grouping_cols)) {
      grouping_cols <- names(data)[grouping_cols] }
  }

  # Designate amino acid or nucleotide for clone sequence
  clone_seq_col <- amino_col
  if (clone_seq_type == "nucleotide") { clone_seq_col <- nucleo_col }


  #### FORMAT AND FILTER DATA ####
  cat(paste0("Input data contains ", nrow(data), " rows.\n"))

  # Remove sequences below specified length
  if (!is.null(min_seq_length)) {
    cat(paste0("Filtering for minimum clone sequence length (", min_seq_length, ")..."))
    data <- filterDataBySequenceLength(data, clone_seq_col,
                                       min_length = min_seq_length)
    cat(paste0(" Done. ", nrow(data), " rows remaining.\n"))
  }

  # Drop sequences with specified chars
  if (!is.null(drop_chars)) {
    cat(paste0("Filtering out sequences with special characters (", drop_chars, ")..."))
    data <- data[-grep(drop_chars, data[ , clone_seq_col]), ]
    cat(paste0(" Done. ", nrow(data), " rows remaining.\n")) }

  if (aggregate_reads) { # Aggregate the counts if specified
    data <- aggregateReads(data, clone_seq_col,
                           count_col, freq_col, grouping_cols)

    # Update column name references for count and freq
    if (size_nodes_by == count_col) { size_nodes_by <- "aggCloneCount" }
    if (size_nodes_by == freq_col) { size_nodes_by <- "aggCloneFreq" }
    if (count_col %in% color_nodes_by) {
      color_nodes_by[which(color_nodes_by == count_col)] <- "aggCloneCount" }
    if (freq_col %in% color_nodes_by) {
      color_nodes_by[which(color_nodes_by == freq_col)] <- "aggCloneFreq" }
    count_col <- "aggCloneCount"
    freq_col <- "aggCloneFreq"

  } else { # Copy the relevant columns from the input data
    data <-
      data[ , # Keep only the relevant columns:
            intersect(
              unique(
                c(nucleo_col, amino_col, count_col, freq_col,
                  vgene_col, dgene_col, jgene_col, cdr3length_col,
                  other_cols, color_nodes_by)),
              names(data))] }

  if (nrow(data) < 2) { stop("insufficient clone sequences remaining; at least two are needed") }

  #### BUILD NETWORK ####
  # Generate adjacency matrix for network
  adjacency_matrix <-
    generateNetworkFromClones(data[ , clone_seq_col],
                              dist_type, edge_dist,
                              contig_ids = rownames(data),
                              return_type = "adjacency_matrix")

  # Subset data to keep only those clones in the network (nonzero degree)
  if (dist_type != "euclidean_on_atchley") {
    data <- data[as.numeric(dimnames(adjacency_matrix)[[1]]), ] }

  # Generate network from adjacency matrix
  net <- generateNetworkFromAdjacencyMat(adjacency_matrix)


  #### NODE/CLUSTER STATS ####
  # Add node-level network characteristics
  if (node_stats) {
    if (typeof(stats_to_include) != "list") {
      if (stats_to_include == "cluster_id_only") {
        stats_to_include <- node_stat_settings(
          degree = FALSE, cluster_id = TRUE, transitivity = FALSE,
          eigen_centrality = FALSE, centrality_by_eigen = FALSE,
          betweenness = FALSE, centrality_by_betweenness = FALSE,
          authority_score = FALSE, coreness = FALSE, page_rank = FALSE)
      }
    }
    data <- addNodeNetworkStats(data, net, stats_to_include) }

  # Compute cluster-level network characteristics
  if (cluster_stats) {

    if (!"cluster_id" %in% names(data)) {
      data <- addClusterMembership(data, net) }

    degree_col <- NULL
    if ("degree" %in% names(data)) { degree_col <- "degree" }

    cluster_info <- getClusterStats(data, adjacency_matrix, clone_seq_col,
                                    count_col, "cluster_id", degree_col) }


  ### PLOT(S) OF NETWORK GRAPH ####
  # if color_nodes_by is NULL, determine default color variable
  if (is.null(color_nodes_by)) {
    if ("degree" %in% names(data)) { # use network degree if available
      color_nodes_by <- "degree"
    } else { # if degree unavailable, color the nodes by clone count
      color_nodes_by <- count_col } }

  # Add cluster ID to node color variables if applicable
  # if ("cluster_id" %in% names(data) & !"cluster_id" %in% color_nodes_by) {
  #   color_nodes_by <- c(color_nodes_by, "cluster_id") }

  # Vector of legend titles as node variable names
  color_legend_title <- color_nodes_by

  # Force consistent legend labels for count/freq if used for node colors
  # (ensures size and color share same legend if using same variable)
  # if (count_col %in% color_nodes_by) {
  #   if (aggregate_reads) {
  #     color_legend_title[which(color_legend_title == count_col)] <-
  #       "agg clone count"
  #   } else {
  #     color_legend_title[which(color_legend_title == count_col)] <-
  #       "clone count" } }
  # if (freq_col %in% color_nodes_by) {
  #   if (aggregate_reads) {
  #     color_legend_title[which(color_legend_title == count_col)] <-
  #       "agg clone freq" }  }

  # If multiple coloring variables, extend color scheme to vector if needed
  if (length(color_nodes_by) > 1 & length(color_scheme) == 1) {
    color_scheme <- rep(color_scheme, length(color_nodes_by)) }

  # Ensure size_nodes_by is a vector or fixed value to use for node sizes
  # also ensure size legend title matches color legend title for same variable
  size_legend_title <- NULL # default for fixed node size
  if (is.character(size_nodes_by)) {
    # if (aggregate_reads) {
    #   if (size_nodes_by == "aggCloneCount") {
    #     size_legend_title <- "agg clone count"
    #   } else if (size_nodes_by == "aggCloneFreq") {
    #     size_legend_title <- "agg clone freq" }
    # } else if (size_nodes_by == count_col) {
    #   size_legend_title <- "clone count"
    # } else {
    size_legend_title <- size_nodes_by
    # }

    size_nodes_by <- data[ , size_nodes_by]
  }

  # Override legend titles with custom values if provided
  if (!is.null(custom_color_legend)) {
    color_legend_title <- custom_color_legend }
  if (!is.null(custom_size_legend)) {
    size_legend_title <- custom_size_legend }

  # Create one plot for each variable used to color the nodes
  temp_plotlist <- list()
  for (j in 1:length(color_nodes_by)) {
    cat(paste0("Generating graph plot with nodes colored by ",
               color_nodes_by[[j]], "..."))
    temp_plotlist$newplot <-
      plotNetworkGraph(
        net, edge_width, title = plot_title, subtitle = plot_subtitle,
        color_nodes_by = data[ , color_nodes_by[[j]]],
        size_nodes_by = size_nodes_by,
        color_legend_title = color_legend_title[[j]],
        size_legend_title = size_legend_title,
        color_scheme = color_scheme[[j]],
        node_size_limits = node_size_limits)
    print(temp_plotlist$newplot)
    names(temp_plotlist)[[length(names(temp_plotlist))]] <- color_nodes_by[[j]]
    cat(" Done.\n") }


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
      cat(paste0("Node-level data saved to file:\n  ", data_outfile, "\n")) }
    if (cluster_stats) {
      utils::write.csv(cluster_info,
                       file = file.path(output_dir, cluster_outfile),
                       row.names = FALSE)
      cat(paste0("Cluster-level data saved to file:\n  ",
                 file.path(output_dir, cluster_outfile), "\n")) } }

  # Save plots to a single pdf
  if (!is.null(output_dir) & !is.null(plot_outfile)) {
    grDevices::pdf(file = file.path(output_dir, plot_outfile),
                   width = plot_width, height = plot_height)
    for (j in 1:length(color_nodes_by)) { print(temp_plotlist[[j]]) }
    grDevices::dev.off()
    cat(paste0("Network graph plot saved to file:\n  ",
               file.path(output_dir, plot_outfile), "\n")) }

  # Save igraph
  if (!is.null(output_dir) & save_all & !is.null(igraph_outfile)) {
    igraph::write_graph(net,
                        file = file.path(output_dir, igraph_outfile),
                        format = "edgelist")
    cat(paste0("Network igraph saved in edgelist format to file:\n  ",
               file.path(output_dir, igraph_outfile), "\n")) }

  # Save adjacency matrix
  if (!is.null(output_dir) & save_all & !is.null(matrix_outfile)) {
    if (dist_type == "euclidean_on_atchley") {
      utils::write.csv(adjacency_matrix,
                       file.path(output_dir, matrix_outfile),
                       row.names = FALSE)
      cat(paste0("Adjacency matrix saved to file:\n  ",
                 file.path(output_dir, matrix_outfile), "\n"))
    } else {
      Matrix::writeMM(adjacency_matrix,
                      file.path(output_dir, matrix_outfile))
      cat(paste0("Adjacency matrix saved to file:\n  ",
                 file.path(output_dir, matrix_outfile), "\n")) } }


  #### RETURN OUTPUT ####
  cat("Finished building network.\n")
  if (return_type == "node_data_only") {

    return(data)

  } else {
    out <- list("node_data" = data)
    if (cluster_stats) { out$cluster_stats <- cluster_info }
    if (return_type == "all") {
      out$plots <- temp_plotlist
      out$adjacency_matrix <- adjacency_matrix
      out$igraph <- net }
    # cat(paste0("All tasks complete. Returning a list containing the following items:\n  ",
    #            paste(names(out), collapse = ", "), "\n"))

    return(out) }

}



# Helper functions --------------------------------------------------------



# .saveNetworkResults <- function(node_data, cluster_data, net, adjacency_matrix,
#                                 graph_plot, output_dir, dist_type, edge_dist,
#                                 min_seq_length, aggregate_reads, group_col) {
#
#   cat(paste0("Saving results to ", output_dir, "\n"))
#
#   # Save settings used to generate network
#   settings <- data.frame("distance_type" = dist_type,
#                          "edge_dist" = edge_dist,
#                          "min_seq_length" = min_seq_length,
#                          "aggregate_reads" = aggregate_reads)
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






