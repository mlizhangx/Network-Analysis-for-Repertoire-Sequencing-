
buildRepSeqNetwork <- function(

  ## Input ##
  data, seq_col, count_col = NULL,
  subset_cols = NULL,
  min_seq_length = 3, drop_matches = NULL,

  ## Network ##
  dist_type = "hamming", dist_cutoff = 1,
  drop_isolated_nodes = TRUE,
  node_stats = FALSE, stats_to_include = node_stat_settings(),
  cluster_stats = FALSE,

  ## Visualization ##
  plots = TRUE, print_plots = TRUE,
  plot_title = "auto", plot_subtitle = "auto",
  color_nodes_by = "auto", ...,

  ## Output ##
  output_dir = getwd(), output_type = "individual",
  output_name = "MyRepSeqNetwork",
  pdf_width = 12, pdf_height = 10

) {
  .createOutputDir(output_dir)

  # Convert column references to character if not already
  seq_col <- .convertColRef(seq_col, data)
  count_col <- .convertColRef(count_col, data)
  color_nodes_by <- .convertColRef(color_nodes_by, data)
  subset_cols <- .convertColRef(subset_cols, data)

  data <- .filterInputData(data, seq_col, min_seq_length, drop_matches,
                           subset_cols, count_col, color_nodes_by)
  if (nrow(data) < 2) {
    warning("insufficient remaining receptor sequences; at least two needed")
    return(invisible(NULL))
  }
  ### BUILD NETWORK ###
  out <- .generateNetworkObjects(
    data, seq_col, dist_type, dist_cutoff, drop_isolated_nodes)

  ### NODE/CLUSTER STATS ###
  if (node_stats) { out$node_data <- addNodeNetworkStats(
      out$node_data, out$igraph, stats_to_include) }
  if (cluster_stats) { out$cluster_data <-
    .getClusterStats2(out$node_data, out$adjacency_matrix, seq_col, count_col) }

  ### GRAPH PLOTS ###
  if (plots) {
    out$plots <- .generateNetworkGraphPlotsGuarded(
      out$igraph, out$node_data, print_plots,
      .makePlotTitle(plot_title, seq_col = seq_col),
      .makePlotSubtitle(plot_subtitle, seq_col = seq_col,
                        dist_type = dist_type, dist_cutoff = dist_cutoff),
      .passColorNodesBy(color_nodes_by, out$node_data, count_col),
      ...)
    if (is.null(out$plots)) { plots <- FALSE }
  }

  ### SAVE & RETURN RESULTS ###
  .saveNetwork(out, output_dir, output_type, output_name,
               plots, pdf_width, pdf_height, cluster_stats, dist_type)
  return(invisible(out))
}



# Helper functions --------------------------------------------------------

# ## Input checks ##
# # Atchley factor embedding only applicable to amino acid sequences
# if (dist_type == "euclidean_on_atchley" & clone_seq_type != "amino acid") {
#   stop("distance type 'euclidean_on_atchley' only applicable to amino acid sequences") }
# aggregate_identical_clones only applicable to amino acid sequences
# each elem of color_nodes_by is a data column or a node stat to be added
# for each elem of other_cols not present in data, warn that not found
# if color_nodes_by is not a scalar/vector of column names/numbers from data, force to "auto" with warning

# .saveNetworkResults <- function(node_data, cluster_data, net, adjacency_matrix,
#                                 graph_plot, output_dir, dist_type, edge_dist,
#                                 min_seq_length, aggregate_identical_clones, group_col) {
#
#   cat(paste0("Saving results to ", output_dir, "\n"))
#
#   # Save settings used to generate network
#   settings <- data.frame("distance_type" = dist_type,
#                          "edge_dist" = edge_dist,
#                          "min_seq_length" = min_seq_length,
#                          "aggregate_identical_clones" = aggregate_identical_clones)
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







# Previous version -----------------------------------------------------------
# Prior to changes on 9/16/2022:
# specific columns like vgene_col that are not used have been removed
# by default, all input columns are now retained and keep their original names
# (if other_cols is non-NULL, only the other_cols and the sequence column are retained)



# buildRepSeqNetwork <- function(
    #
#   # Input Data and Columns
#   data,      # data frame containing req-seq data
#   nucleo_col,
#   amino_col,
#   count_col,
#   freq_col,
#   vgene_col = NULL,# ignored if aggregate_identical_clones = TRUE
#   dgene_col = NULL,# ignored if aggregate_identical_clones = TRUE
#   jgene_col = NULL,# ignored if aggregate_identical_clones = TRUE
#   cdr3length_col = NULL,# ignored if aggregate_identical_clones = TRUE
#   other_cols = NULL, # other cols to keep; ignored if aggregate_identical_clones = TRUE
#
#   # Clone Sequence Settings
#   clone_seq_type = "amino acid", # or "nucleotide"
#   min_seq_length = 3, # min clone seq length
#   drop_chars = NULL, # regular expression, e.g. "[*|_]"
#   aggregate_identical_clones = FALSE,
#   grouping_cols = NULL,
#
#   # Network Settings
#   dist_type = "hamming", # or "levenshtein", "hamming", "euclidean_on_atchley"
#   edge_dist = 1, # max dist for edges
#   drop_isolated_nodes = TRUE,
#   node_stats = FALSE,
#   stats_to_include = node_stat_settings(), # can also select "all" or "cluster_id_only"
#   cluster_stats = FALSE,
#
#   # Plot Settings
#   plot_title = "auto",
#   plot_subtitle = "auto",
#   color_nodes_by = "auto", # uses degree if available, else clone count
#   color_scheme = "default", #  (accepts vector of same length as color_nodes_by)
#   color_legend = TRUE,
#   color_title = "auto", # custom title (accepts vector of same length as color_nodes_by)
#   edge_width = 0.3,
#   size_nodes_by = count_col, # can use a double, e.g., 1.0, for fixed size
#   node_size_limits = "auto", # numeric, length 2
#   size_title = "auto", # custom legend title
#
#   # Output Settings
#   output_dir = NULL, # if NULL, output is not saved to file
#   save_all = FALSE, # by default, only save pdf of plot and csv of node data
#   data_outfile = "node_data.csv",
#   plot_outfile = "network_graph.pdf",
#   plot_width = 12, # passed to pdf()
#   plot_height = 10, # passed to pdf()
#   cluster_outfile = "cluster_info.csv", # only saved if save_all = TRUE
#   igraph_outfile = "network_edgelist.txt", # .txt
#   matrix_outfile = "auto",
#   return_all = FALSE # if false, only the node data is returned (unless cluster_stats = TRUE and output_dir = NULL, in which case a list containing the node data and cluster info is returned, with a warning)
#
# ) {
#
#   #### INPUT CHECKS ####
#
#   # Atchley factor embedding only applicable to amino acid sequences
#   if (dist_type == "euclidean_on_atchley" & clone_seq_type != "amino acid") {
#     stop("distance type 'euclidean_on_atchley' only applicable to amino acid sequences") }
#
#   # aggregate_identical_clones only applicable to amino acid sequences
#
#   # each elem of color_nodes_by is a data column or a node stat to be added
#
#   # for each elem of other_cols not present in data, warn that not found
#
#   # if color_nodes_by is not a scalar/vector of column names/numbers from data, force to "auto" with warning
#
#   #### PREPARE WORKING ENVIRONMENT ####
#   # Create output directory if applicable
#   if (!is.null(output_dir)) { .createOutputDir(output_dir) }
#
#   # Convert input columns to character if not already
#   if (is.numeric(nucleo_col)) { nucleo_col <- names(data)[nucleo_col] }
#   if (is.numeric(amino_col)) { amino_col <- names(data)[amino_col] }
#   if (is.numeric(count_col)) { count_col <- names(data)[count_col] }
#   if (is.numeric(freq_col)) { freq_col <- names(data)[freq_col] }
#   if (is.numeric(color_nodes_by)) {
#     color_nodes_by <- names(data)[color_nodes_by]
#   }
#   if (!aggregate_identical_clones) {
#     if (is.numeric(vgene_col)) { vgene_col <- names(data)[vgene_col] }
#     if (is.numeric(dgene_col)) { dgene_col <- names(data)[dgene_col] }
#     if (is.numeric(jgene_col)) { jgene_col <- names(data)[jgene_col] }
#     if (is.numeric(cdr3length_col)) { cdr3length_col <- names(data)[cdr3length_col] }
#     if (is.numeric(other_cols)) { other_cols <- names(data)[other_cols] }
#   } else {
#     if (is.numeric(grouping_cols)) {
#       grouping_cols <- names(data)[grouping_cols] }
#   }
#
#   # Designate amino acid or nucleotide for clone sequence
#   clone_seq_col <- amino_col
#   if (clone_seq_type == "nucleotide") { clone_seq_col <- nucleo_col }
#
#
#   #### FORMAT AND FILTER DATA ####
#   cat(paste0("Input data contains ", nrow(data), " rows.\n"))
#
#   # Filter by seq length
#   if (!is.null(min_seq_length)) {
#     cat(paste0("Filtering for minimum clone sequence length (", min_seq_length, ")..."))
#     data <- filterClonesBySequenceLength(data, clone_seq_col,
#                                          min_length = min_seq_length)
#     cat(paste0(" Done. ", nrow(data), " rows remaining.\n"))
#   }
#
#   # Filter seqs with special chars
#   if (!is.null(drop_chars)) {
#     cat(paste0("Filtering out sequences with special characters (", drop_chars, ")..."))
#     drop_matches <- grep(drop_chars, data[ , clone_seq_col])
#     if (length(drop_matches) > 0) { data <- data[-drop_matches, ] }
#     cat(paste0(" Done. ", nrow(data), " rows remaining.\n")) }
#
#   # Format input data
#   if (aggregate_identical_clones) { # Aggregate the counts if specified
#     data <- aggregateIdenticalClones(data, clone_seq_col,
#                                      count_col, freq_col, grouping_cols)
#     # Update column name references for count and freq
#     if (size_nodes_by == count_col) { size_nodes_by <- "AggregatedCloneCount" }
#     if (size_nodes_by == freq_col) { size_nodes_by <- "AggregatedCloneFrequency" }
#     if (count_col %in% color_nodes_by) {
#       color_nodes_by[which(color_nodes_by == count_col)] <- "AggregatedCloneCount" }
#     if (freq_col %in% color_nodes_by) {
#       color_nodes_by[which(color_nodes_by == freq_col)] <- "AggregatedCloneFrequency" }
#     count_col <- "AggregatedCloneCount"
#     freq_col <- "AggregatedCloneFrequency"
#   } else { # Copy the relevant columns from the input data
#     data <-
#       data[ , # Keep only the relevant columns:
#             intersect(
#               unique(
#                 c(nucleo_col, amino_col, count_col, freq_col,
#                   vgene_col, dgene_col, jgene_col, cdr3length_col,
#                   other_cols, color_nodes_by)),
#               names(data))] }
#
#   if (nrow(data) < 2) { stop("insufficient clone sequences remaining; at least two are needed") }
#
#   #### BUILD NETWORK ####
#   # Generate adjacency matrix for network
#   adjacency_matrix <-
#     generateNetworkFromClones(data[ , clone_seq_col],
#                               dist_type, edge_dist,
#                               contig_ids = rownames(data),
#                               return_type = "adjacency_matrix",
#                               drop_isolated_nodes = drop_isolated_nodes)
#
#   # Subset data to keep only those clones in the network (nonzero degree)
#   if (drop_isolated_nodes & dist_type != "euclidean_on_atchley") {
#     data <- data[as.numeric(dimnames(adjacency_matrix)[[1]]), ] }
#
#   # Generate network from adjacency matrix
#   net <- generateNetworkFromAdjacencyMat(adjacency_matrix)
#
#
#   #### NODE/CLUSTER STATS ####
#   # Add node-level network characteristics
#   if (node_stats) {
#     if (typeof(stats_to_include) != "list") {
#       if (stats_to_include == "cluster_id_only") {
#         stats_to_include <- node_stat_settings(
#           degree = FALSE, cluster_id = TRUE, transitivity = FALSE,
#           eigen_centrality = FALSE, centrality_by_eigen = FALSE,
#           betweenness = FALSE, centrality_by_betweenness = FALSE,
#           authority_score = FALSE, coreness = FALSE, page_rank = FALSE)
#       }
#     }
#     data <- addNodeNetworkStats(data, net, stats_to_include)
#   }
#
#   # Compute cluster-level network characteristics
#   if (cluster_stats) {
#     if (!"cluster_id" %in% names(data)) {
#       data <- addClusterMembership(data, net)
#     }
#     degree_col <- NULL
#     if ("degree" %in% names(data)) { degree_col <- "degree" }
#     cluster_info <- getClusterStats(data, adjacency_matrix, clone_seq_col,
#                                     count_col, "cluster_id", degree_col)
#   }
#
#
#   ### PLOT(S) OF NETWORK GRAPH ####
#   # Determine plot title/subtitle
#   if (!is.null(plot_title)) {
#     if (plot_title == "auto") {
#       plot_title <- paste("Network on", clone_seq_type, "sequence")
#     }
#   }
#   if (!is.null(plot_subtitle)) {
#     if (plot_subtitle == "auto") {
#       if (dist_type == "euclidean_on_atchley") {
#         plot_subtitle <- paste(
#           "Clone sequences encoded as numeric vectors using deep learning based on Atchley factor representation\nEdges based on a maximum Euclidean distance of", edge_dist, "between encoded vectors\n")
#       } else {
#         plot_subtitle <- paste("Edges based on a maximum", dist_type, "distance of", edge_dist, "\n")
#       }
#     }
#   }
#
#   # Assign default variable for node colors if applicable
#   if (length(color_nodes_by) == 1) {
#     if (color_nodes_by == "auto") {
#       if ("degree" %in% names(data)) { # use network degree if available
#         color_nodes_by <- "degree"
#       } else { # if degree unavailable, color the nodes by clone count
#         color_nodes_by <- count_col }
#     }
#   }
#
#   # If multiple coloring variables, extend color scheme and legend title to vectors if needed
#   if (length(color_nodes_by) > 1) {
#     if (length(color_scheme) == 1) { # extend to vector
#       color_scheme <- rep(color_scheme, length(color_nodes_by)) }
#     if (!is.null(color_title)) {
#       if (length(color_title) == 1) { # extend to vector
#         color_title <- rep(color_title, length(color_nodes_by)) }
#     } else { # color_title is NULL
#       # vector of empty titles (hack, since can't have NULL vector entries)
#       color_title <- rep("", length(color_nodes_by))
#     }
#   }
#
#   # Set default color legend title if applicable
#   if (!is.null(color_title)) {
#     for (i in 1:length(color_title)) {
#       if (color_title[[i]] == "auto") {
#         color_title[[i]] <- color_nodes_by[[i]]
#       }
#     }
#   }
#   # Set default size legend title if applicable
#   if (!is.null(size_title)) {
#     if (size_title == "auto") {
#       if (is.numeric(size_nodes_by)) {
#         size_title <- NULL
#       }
#       if (is.character(size_nodes_by)) {
#         size_title <- size_nodes_by
#       }
#     }
#   }
#
#   # If using a variable to size nodes, get the vector by column name
#   if (is.character(size_nodes_by)) {
#     size_nodes_by <- data[ , size_nodes_by]
#   }
#
#   # Create one plot for each variable used to color the nodes
#   temp_plotlist <- list()
#   for (j in 1:length(color_nodes_by)) {
#     if (color_nodes_by[[j]] == "auto") { cat("Generating graph plot...")
#     } else {
#       cat(paste0("Generating graph plot with nodes colored by ",
#                  color_nodes_by[[j]], "..."))
#     }
#     temp_plotlist$newplot <-
#       plotNetworkGraph(
#         net, edge_width, title = plot_title, subtitle = plot_subtitle,
#         color_nodes_by = data[ , color_nodes_by[[j]]],
#         color_legend_title = color_title[[j]],
#         color_scheme = color_scheme[[j]],
#         show_color_legend = color_legend,
#         size_nodes_by = size_nodes_by,
#         size_legend_title = size_title,
#         node_size_limits = node_size_limits)
#     print(temp_plotlist$newplot)
#     names(temp_plotlist)[[length(names(temp_plotlist))]] <- color_nodes_by[[j]]
#     cat(" Done.\n") }
#
#
#   #### SAVE RESULTS ####
#   # Rename data columns
#   if (aggregate_identical_clones) {
#     colnames(data)[[1]] <- ifelse(clone_seq_type == "nucleotide",
#                                   yes = "NucleotideSeq",
#                                   no = "AminoAcidSeq")
#   } else {
#     colnames(data)[colnames(data) == nucleo_col] <- "NucleotideSeq"
#     colnames(data)[colnames(data) == amino_col] <- "AminoAcidSeq"
#     colnames(data)[colnames(data) == freq_col] <- "CloneFrequency"
#     colnames(data)[colnames(data) == count_col] <- "CloneCount"
#     if (!is.null(vgene_col)) {
#       colnames(data)[colnames(data) == vgene_col] <- "VGene"
#     }
#     if (!is.null(dgene_col)) {
#       colnames(data)[colnames(data) == dgene_col] <- "DGene"
#     }
#     if (!is.null(jgene_col)) {
#       colnames(data)[colnames(data) == jgene_col] <- "JGene"
#     }
#     if (!is.null(cdr3length_col)) {
#       colnames(data)[colnames(data) == cdr3length_col] <- "CDR3Length"
#     }
#   }
#
#   # Save node [& cluster] data
#   if (!is.null(output_dir)) {
#     if (!is.null(data_outfile)) {
#       utils::write.csv(data, file = file.path(output_dir, data_outfile),
#                        row.names = FALSE)
#       cat(paste0("Node-level data saved to file:\n  ", data_outfile, "\n")) }
#     if (cluster_stats) {
#       utils::write.csv(cluster_info,
#                        file = file.path(output_dir, cluster_outfile),
#                        row.names = FALSE)
#       cat(paste0("Cluster-level data saved to file:\n  ",
#                  file.path(output_dir, cluster_outfile), "\n")) } }
#
#   # Save plots to a single pdf
#   if (!is.null(output_dir) & !is.null(plot_outfile)) {
#     grDevices::pdf(file = file.path(output_dir, plot_outfile),
#                    width = plot_width, height = plot_height)
#     for (j in 1:length(color_nodes_by)) { print(temp_plotlist[[j]]) }
#     grDevices::dev.off()
#     cat(paste0("Network graph plot saved to file:\n  ",
#                file.path(output_dir, plot_outfile), "\n")) }
#
#   # Save igraph
#   if (!is.null(output_dir) & save_all & !is.null(igraph_outfile)) {
#     igraph::write_graph(net,
#                         file = file.path(output_dir, igraph_outfile),
#                         format = "edgelist")
#     cat(paste0("Network igraph saved in edgelist format to file:\n  ",
#                file.path(output_dir, igraph_outfile), "\n")) }
#
#   # Save adjacency matrix
#   if (!is.null(output_dir) & save_all & !is.null(matrix_outfile)) {
#     if (matrix_outfile == "auto") {
#       if (dist_type == "euclidean_on_atchley") {
#         matrix_outfile <- "adjacency_matrix.csv"
#       } else {
#         matrix_outfile <- "adjacency_matrix.mtx"
#       }
#     }
#     if (dist_type == "euclidean_on_atchley") {
#       utils::write.csv(adjacency_matrix,
#                        file.path(output_dir, matrix_outfile),
#                        row.names = FALSE)
#       cat(paste0("Adjacency matrix saved to file:\n  ",
#                  file.path(output_dir, matrix_outfile), "\n"))
#     } else {
#       Matrix::writeMM(adjacency_matrix,
#                       file.path(output_dir, matrix_outfile))
#       cat(paste0("Adjacency matrix saved to file:\n  ",
#                  file.path(output_dir, matrix_outfile), "\n")) } }
#
#
#   #### RETURN OUTPUT ####
#   cat("Finished building network.\n")
#   return_type <- ifelse(return_all,
#                         yes = "all",
#                         no = ifelse(cluster_stats & is.null(output_dir),
#                                     yes = "node_and_cluster_data",
#                                     no = "node_data_only"))
#   if (return_type == "node_data_only") {
#     return(data)
#   } else {
#     out <- list("node_data" = data)
#     if (cluster_stats) { out$cluster_stats <- cluster_info }
#     if (return_type == "all") {
#       out$plots <- temp_plotlist
#       out$adjacency_matrix <- adjacency_matrix
#       out$igraph <- net
#     }
#     # cat(paste0("All tasks complete. Returning a list containing the following items:\n  ",
#     #            paste(names(out), collapse = ", "), "\n"))
#     return(out)
#   }
# }


