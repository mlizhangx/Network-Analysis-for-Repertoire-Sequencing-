# Single-cell TCR network combining similarities in alpha & beta chain sequences

# Assumptions for input data frame:
#   Separate columns for TRA and TRB sequences
#   One cell per row

buildDualChainNetwork <- function(
    data,
    a_col, # alpha chain clone sequences
    b_col,  # beta chain clone sequences
    count_col = NULL, # optional column name or number containing measurements of clonal abundance
    other_cols = NULL, # other cols to keep (if NULL, all are kept); ignored if aggregate_identical_clones = TRUE
    min_seq_length = 3, # min clone seq length (applied to both chains)
    drop_chars = NULL, # regular expression, e.g. "[*|_]"

    # Network Settings
    dist_type = "hamming", # or "levenshtein", "hamming", "euclidean_on_atchley"
    dist_cutoff = 1, # max dist for edges
    drop_isolated_nodes = TRUE,  # from final network
    node_stats = FALSE,
    stats_to_include = node_stat_settings(), # can also select "all" or "cluster_id_only"
    # cluster_stats = FALSE,

    # Plot Settings
    plot_title = "auto",
    plot_subtitle = "auto",
    color_nodes_by = "auto", # uses degree if available, else clone count if available, else nothing
    color_scheme = "default", #  (accepts vector of same length as color_nodes_by)
    color_legend = TRUE,
    color_title = "auto", # custom title (accepts vector of same length as color_nodes_by)
    edge_width = 0.1,
    size_nodes_by = 0.5, # can use a double, e.g., 1.0, for fixed size
    node_size_limits = "auto", # numeric, length 2
    size_title = "auto", # custom legend title

    # Output Settings
    print_plots = TRUE,
    output_dir = NULL, # if NULL, output is not saved to file
    save_all = FALSE, # by default, only save pdf of plot and csv of node data
    data_outfile = "node_data.csv",
    plot_outfile = "network_graph.pdf",
    plot_width = 12, # passed to pdf()
    plot_height = 10, # passed to pdf()
    cluster_outfile = "cluster_info.csv", # only saved if save_all = TRUE
    igraph_outfile = "network_edgelist.txt", # .txt
    matrix_outfile = "auto",
    return_all = FALSE # if false, only the node data is returned (unless cluster_stats = TRUE and output_dir = NULL, in which case a list containing the node data and cluster info is returned, with a warning)


) {

  stopifnot("dist_type must be 'hamming' or 'levenshtein'" =
              dist_type %in% c("hamming", "levenshtein"))

  #### PREPARE WORKING ENVIRONMENT ####
  # Create output directory if applicable
  if (!is.null(output_dir)) { .createOutputDir(output_dir) }

  # Convert column references to character if not already
  if (is.numeric(a_col)) { a_col <- names(data)[a_col] }
  if (is.numeric(b_col)) { b_col <- names(data)[b_col] }
  if (is.numeric(count_col)) { count_col <- names(data)[count_col] }
  if (is.numeric(other_cols)) { other_cols <- names(data)[other_cols] }
  if (is.numeric(color_nodes_by)) {
    color_nodes_by <- names(data)[color_nodes_by]
  }

  #### FORMAT AND FILTER DATA ####
  # Coerce sequence columns to character if needed
  if (!is.character(data[ , a_col])) {
    data[ , a_col] <- as.character(data[ , a_col]) }
  if (!is.character(data[ , b_col])) {
    data[ , b_col] <- as.character(data[ , b_col])

  cat(paste0("Input data contains ", nrow(data), " rows.\n"))

  # Filter by seq length
  if (!is.null(min_seq_length)) {
    cat(paste0("Removing rows containing sequences with length less than ", min_seq_length, "..."))
    data <- filterClonesBySequenceLength(data, a_col, min_seq_length)
    data <- filterClonesBySequenceLength(data, b_col, min_seq_length)
    cat(paste0(" Done. ", nrow(data), " rows remaining.\n"))
  }

  # Filter seqs with special chars
  if (!is.null(drop_chars)) {
    cat(paste0("Removing rows with sequences containing matches to the expression '", drop_chars, "'..."))
    drop_matches <- union(grep(drop_chars, data[ , a_col]),
                          grep(drop_chars, data[ , b_col]))
    if (length(drop_matches) > 0) { data <- data[-drop_matches, ] }
    cat(paste0(" Done. ", nrow(data), " rows remaining.\n")) }

  # Format input data
  if (!is.null(other_cols)) {
    data <-
      data[ , # Keep only the relevant columns:
            intersect(
              unique(c(a_col, b_col, count_col, other_cols, color_nodes_by)),
              names(data))]
  }

  if (nrow(data) < 2) { stop("insufficient data rows remaining; at least two are needed") }


  #### GENERATE DUAL CHAIN NEWTORK ####
  # adjacency matrix for alpha chain
  cat("Computing graph adjacency based on alpha chain sequences:\n")
  adj_mat_a <- sparseAdjacencyMatFromSeqs(
    data[ , a_col], dist_type, dist_cutoff, drop_isolated_nodes = FALSE)
  # adjacency matrix for beta chain
  cat("Computing graph adjacency based on beta chain sequences:\n")
  adj_mat_b <- sparseAdjacencyMatFromSeqs(
    data[ , b_col], dist_type, dist_cutoff, drop_isolated_nodes = FALSE)

  # Combine adjacency matrices for both chains
  # (only edges present for both chains will become edges in the combined graph)
  cat("Intersecting the adjacencies from both chains...")
  adjacency_matrix <- adj_mat_a + adj_mat_b
  adjacency_matrix[adjacency_matrix == 1] <- 0
  adjacency_matrix[adjacency_matrix == 2] <- 1
  cat(" Done.\n")

  # Generate network from combined adjacency matrix
  cat("Building network based on the combined adjacencies... ")
  net <- generateNetworkFromAdjacencyMat(adjacency_matrix); cat(" Done.\n")

  # Drop isolated nodes from final network if specified
  if (drop_isolated_nodes) { cat("Dropping isolated nodes...")
    nodes_to_keep <- igraph::degree(net) > 0
    adjacency_matrix <- adjacency_matrix[nodes_to_keep, nodes_to_keep]
    data <- data[nodes_to_keep, ]
    # regenerate network without isolated nodes
    net <- generateNetworkFromAdjacencyMat(adjacency_matrix)
    cat(paste0(" Done. ", nrow(data), " nodes remaining.\n"))
  }

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
    data <- addNodeNetworkStats(data, net, stats_to_include)
  }

  # # Compute cluster-level network characteristics
  # if (cluster_stats) {
  #   if (!"cluster_id" %in% names(data)) {
  #     data <- addClusterMembership(data, net)
  #   }
  #   degree_col <- NULL
  #   if ("degree" %in% names(data)) { degree_col <- "degree" }
  #   cluster_info <- getClusterStats(data, adjacency_matrix, seq_col,
  #                                   count_col, "cluster_id", degree_col)
  # }


  ### PLOT(S) OF NETWORK GRAPH ####
  # plot title/subtitle
  if (!is.null(plot_title)) { if (plot_title == "auto") {
    plot_title <- paste("Immune Repertoire Network\nby Combined Alpha and Beta Chain Similarity") } }
  if (!is.null(plot_subtitle)) { if (plot_subtitle == "auto") {
    plot_subtitle <- paste("Each node denotes a single TCR/BCR cell\nEdges denote a maximum", dist_type, "distance of", dist_cutoff, "between sequences in corresponding chains\n") } }

  # node color variable
  if (length(color_nodes_by) == 1) { if (color_nodes_by == "auto") {
    if ("degree" %in% names(data)) { color_nodes_by <- "degree"
    } else if (!is.null(count_col)) { color_nodes_by <- count_col
    } else { color_nodes_by <- NULL } } }

  # node color palette and legend title
  if (is.null(color_nodes_by)) { color_title <- NULL
  } else { # color_nodes_by is non NULL
    if (length(color_nodes_by) > 1) { # extend objects to vectors if needed
      if (length(color_scheme) == 1) {
        color_scheme <- rep(color_scheme, length(color_nodes_by)) }
      if (!is.null(color_title)) { if (length(color_title) == 1) {
        color_title <- rep(color_title, length(color_nodes_by)) }
      } else { # color_title is NULL
        color_title <- rep("", length(color_nodes_by)) } } # (hack, since can't have NULL vector entries)
    if (!is.null(color_title)) { # Set default color legend title if applicable
      for (i in 1:length(color_title)) { if (color_title[[i]] == "auto") {
        color_title[[i]] <- color_nodes_by[[i]] } } } }

  # size legend title
  if (!is.null(size_title)) { if (size_title == "auto") { # default size title
    if (is.numeric(size_nodes_by)) { size_title <- NULL } # fixed node sizes
    if (is.character(size_nodes_by)) { size_title <- size_nodes_by } } }

  # node size variable (if applicable)
  if (is.character(size_nodes_by)) { size_nodes_by <- data[ , size_nodes_by] }

  # Create one plot for each variable used to color the nodes
  temp_plotlist <- list()
  if (is.null(color_nodes_by)) {
    cat("Generating graph plot...")
    temp_plotlist$uniform_color <-
      plotNetworkGraph(
        net, edge_width, title = plot_title, subtitle = plot_subtitle,
        color_nodes_by = NULL,
        color_legend_title = NULL,
        color_scheme = color_scheme,
        show_color_legend = FALSE,
        size_nodes_by = size_nodes_by,
        size_legend_title = size_title,
        node_size_limits = node_size_limits)
    if (print_plots) { print(temp_plotlist$uniform_color) }
    cat(" Done.\n")
  } else {
    for (j in 1:length(color_nodes_by)) {
      cat(paste0("Generating graph plot with nodes colored by ",
                 color_nodes_by[[j]], "..."))
      temp_plotlist$newplot <-
        plotNetworkGraph(
          net, edge_width, title = plot_title, subtitle = plot_subtitle,
          color_nodes_by = data[ , color_nodes_by[[j]]],
          color_legend_title = color_title[[j]],
          color_scheme = color_scheme[[j]],
          show_color_legend = color_legend,
          size_nodes_by = size_nodes_by,
          size_legend_title = size_title,
          node_size_limits = node_size_limits)
      if (print_plots) { print(temp_plotlist$newplot) }
      names(temp_plotlist)[[length(names(temp_plotlist))]] <- color_nodes_by[[j]]
      cat(" Done.\n")
    }
  }


  #### SAVE RESULTS ####
  # Save node [& cluster] data
  if (!is.null(output_dir)) {
    if (!is.null(data_outfile)) {
      utils::write.csv(data, file = file.path(output_dir, data_outfile),
                       row.names = FALSE)
      cat(paste0("Node-level data saved to file:\n  ", data_outfile, "\n")) }
    # if (cluster_stats) {
    #   utils::write.csv(cluster_info,
    #                    file = file.path(output_dir, cluster_outfile),
    #                    row.names = FALSE)
    #   cat(paste0("Cluster-level data saved to file:\n  ",
    #              file.path(output_dir, cluster_outfile), "\n")) }
  }

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
    if (matrix_outfile == "auto") { matrix_outfile <- "adjacency_matrix.mtx" }
    Matrix::writeMM(adjacency_matrix, file.path(output_dir, matrix_outfile))
    cat(paste0("Adjacency matrix saved to file:\n  ",
               file.path(output_dir, matrix_outfile), "\n")) }


  #### RETURN OUTPUT ####
  cat("Finished building network.\n")
  if (!return_all) { return(data)
  } else {
    out <- list("node_data" = data)
      out$plots <- temp_plotlist
      out$adjacency_matrix <- adjacency_matrix
      out$igraph <- net
      out$alpha_matrix <- adj_mat_a
      out$beta_matrix <- adj_mat_b
    return(out)
  }
}

