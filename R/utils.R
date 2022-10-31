
# Description -------------------------------------------------------------

# General utility and helper functions, both public and internal


# File Management ---------------------------------------------------------

.createOutputDir <- function(dirname) {
  dir.create(dirname, showWarnings = FALSE, recursive = TRUE)
  # Confirm successful creation of output directory
  if (!dir.exists(dirname)) {
    stop(paste0("Unable to create directory ", dirname,
                ". Check to confirm that a valid directory name was provided."))
  }
}


# Data Manipulation -------------------------------------------------------

# INPUT:
#   rep-seq data w/ specification for count column and one or more grouping cols
#   (clone seq is typically the grouping variable used)
# DO:
#   aggregate the counts by group based on the grouping columns
#   add variable counting the number of reads for each
aggregateIdenticalClones <- function(
    data, # data frame containing columns below
    clone_col,
    count_col, # name or number of column of `data` containing clone counts
    freq_col,
    grouping_cols = NULL # optional integer or character vector specifying additional grouping columns
) {

  # Convert column specifications from numeric to character if not already
  if (is.numeric(clone_col)) { clone_col <- names(data)[clone_col] }
  if (is.numeric(count_col)) { count_col <- names(data)[count_col] }
  if (is.numeric(freq_col)) { freq_col <- names(data)[freq_col] }

  # Define grouping variable(s)
  grouping_variables <- list(data[ , clone_col])
  names(grouping_variables) <- clone_col
  if (!is.null(grouping_cols)) {
    if (is.numeric(grouping_cols)) { grouping_cols <- names(data)[grouping_cols] }
    for (i in 1:length(grouping_cols)) {
      grouping_variables$newvar <- data[ , grouping_cols[[i]]]
      names(grouping_variables)[[length(grouping_variables)]] <-
        grouping_cols[[i]] } }

  # aggregate the reads by group
  cat("Aggregating reads (rows) by unique clone sequence...")
  data_to_aggregate <- list("AggregatedCloneCount" = data[ , c(count_col)],
                            "AggregatedCloneFrequency" = data[ , c(freq_col)])
  agg_counts <- stats::aggregate(data_to_aggregate,
                                 by = grouping_variables, FUN = sum)

  # add variable for num reads (row count)
  groups <- as.data.frame(grouping_variables)
  names(groups)[[1]] <- "temporary_placeholder_name" # for summarize function
  num_reads <- dplyr::summarize(dplyr::group_by_all(groups),
                                UniqueCloneCount = length(temporary_placeholder_name))
  names(num_reads)[[1]] <- clone_col # replace placeholder name with orig

  # Merge aggregate counts with num reads
  out <- merge(agg_counts, num_reads, by = c(clone_col, grouping_cols))
  cat(paste0(" Done. ", nrow(out), " unique clone sequences found.\n"))

  return(out)
}



# FUNCTION: Filter rep-seq data to remove rows for clonotype sequences with
# length below the specified cutoff
filterClonesBySequenceLength <- function(data, seq_col, min_length = 3) {
  if (ncol(data) == 1) {
    out <- as.data.frame(data[nchar(data[ , seq_col]) >= min_length, ])
    colnames(out) <- colnames(data)
    return(out)
  } else {
    return(data[nchar(data[ , seq_col]) >= min_length, ])
  }
}


# FUNCTION: EXTRACT DATA SUBSET FOR ALL SEQUENCES WITHIN SPECIFIED RADIUS OF
# TARGET SEQUENCE BY SPECIFIED DISTANCE TYPE
# If a sample_col is provided, only samples that possess the target sequence
# will be included
getSimilarClones <- function(
    target_seq, # specified candidate sequence for the neighborhood
    data, # data frame containing rep seq data, possibly from multiple samples
    seq_col, # col name/# containing clone sequences
    sample_col = NULL, # optional col name/# containing sample IDs (only samples possessing target seq will be included)
    dist_type = "hamming", # options are "hamming" and "levenshtein"
    max_dist = 1, # Maximum Levenshtein distance allowed for inclusion in neighborhood
    drop_chars = NULL # regular expression for chars to filter sequences by
) {
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }
  # If sample_id is supplied, subset data keeping only samples with target seq
  if (is.null(sample_col)) {
    data_samples_w_targetseq <- data
  } else {
    ### SUBSET DATA: SAMPLES WITH TARGET SEQ ###
    cat("Finding all samples that possess the target sequence...")
    # Get row ids of merged data corresponding to target seq
    rows_for_targetseq <- grep(pattern = paste0("^", target_seq, "$"),
                               x = data[ , seq_col])
    # Extract rows of merged data corresponding to samples with target seq
    data_samples_w_targetseq <-
      data[
        data[ , sample_col] %in% data[rows_for_targetseq, sample_col], ]
    cat(" Done.\n")
  }
  # Remove sequences that match expression in `drop_chars`
  if (!is.null(drop_chars)) {
    drop_matches <- grep(drop_chars, data_samples_w_targetseq[ , seq_col])
    if (length(drop_matches) > 0) {
      data_samples_w_targetseq <- data_samples_w_targetseq[-drop_matches, ]
    }
  }
  ### SUBSET DATA: NEIGHBORHOOD OF TARGET SEQUENCE ###
  # Compute list of bounded distances between target seq and seqs
  # possessed by samples with target seq (values are -1 where bound is exceeded)
  # returned vector will be of type integer; names will be the sequences
  cat("Gathering all cells/clones with receptor sequences similar to the target sequence...")
  if (dist_type == "levenshtein") { dist_fun <- levDistBounded
  } else if (dist_type == "hamming") { dist_fun <- hamDistBounded
  } else { stop("invalid option for `dist_type`") }
  dists_to_targetseq <- sapply(
    X = data_samples_w_targetseq[ , seq_col],
    FUN = dist_fun, b = target_seq, k = max_dist)
  # get data for sequences within the specified radius
  data_targetseq_neighborhood <-
    data_samples_w_targetseq[dists_to_targetseq != -1, ]
  cat(paste0(" Done. ", nrow(data_targetseq_neighborhood), " similar cells/clones found.\n"))
  return(data_targetseq_neighborhood)
}


# Network Building --------------------------------------------------------



# FUNCTION: GENERATE NETWORK FOR A LIST OF RECEPTOR SEQS USING SPECIFIED DISTANCE TYPE AND THRESHOLD
generateNetworkFromSeqs <- function(
    seqs, # character vector of receptor sequences
    dist_type = "hamming", # supports "levenshtein", "hamming", "euclidean_on_atchley"
    dist_cutoff = 1, # max dist threshold for edges
    drop_isolated_nodes = TRUE, # forced to FALSE for dist_type = "euclidean_on_atchley"
    contig_ids = seq_along(seqs), # for dist_type = "euclidean_on_atchley"
    outfile_adjacency_matrix = NULL, # save file for adjacency matrix
    outfile_distance_matrix = NULL, # save file for distance matrix (only for Euclidean on Atchley)
    return_type = "network" # can use "adjacency_matrix" to return the adjacency mat
) {
  ### COMPUTE ADJACENCY MATRIX ###
  if (dist_type %in% c("levenshtein", "hamming")) {
    adjacency_matrix <-
      sparseAdjacencyMatFromSeqs(seqs = seqs,
                                 dist_type = dist_type,
                                 max_dist = dist_cutoff,
                                 drop_isolated_nodes = drop_isolated_nodes)
    if (!is.null(outfile_adjacency_matrix)) {
      Matrix::writeMM(adjacency_matrix, outfile_adjacency_matrix)
    }
  } else if (dist_type == "euclidean_on_atchley") {
    adjacency_matrix <-
      adjacencyMatAtchleyFromSeqs(
        seqs = seqs,
        contig_ids = contig_ids,
        max_dist = dist_cutoff,
        outfile_distance_matrix = outfile_distance_matrix)
  } else { stop("invalid option for argument `dist_type`") }
  if (return_type == "adjacency_matrix") {
    return(adjacency_matrix)
  } else {
    network <- generateNetworkFromAdjacencyMat(adjacency_matrix)
    return(network)
  }
}

# Simple wrapper to igraph functions:
# Use adjacency matrix to generate network graph
generateNetworkFromAdjacencyMat <- function(adjacency_matrix) {
  set.seed(9999)
  net <- igraph::graph_from_adjacency_matrix(adjacency_matrix, weighted = TRUE)
  net <- igraph::as.undirected(
    igraph::simplify(net, remove.multiple = T, remove.loops = T))
  return(net)
}



# Computing Network Statistics --------------------------------------------


# input network and corresponding metadata;
# augment metadata with node-level network info and return
addNodeNetworkStats <- function(
    data, # rep-seq data corresponding to the network
    net, # igraph network object
    stats_to_include = node_stat_settings()
) {

  if (typeof(stats_to_include) != "list")  {
    if (stats_to_include == "all") {
      stats_to_include <- node_stat_settings(all_stats = TRUE) } }
  if (stats_to_include$degree | stats_to_include$all_stats) {
    data$degree <- igraph::degree(net) }

  if (stats_to_include$cluster_id | stats_to_include$all_stats) {
    cat("Computing cluster membership within the network...")
    data$cluster_id <-
      as.factor(as.integer(igraph::cluster_fast_greedy(net)$membership))
    cat(" Done.\n") }

  cat(paste0("Computing node-level network statistics..."))
  if (stats_to_include$transitivity | stats_to_include$all_stats) {
    data$transitivity <- igraph::transitivity(net, type = "local")
  }
  if (stats_to_include$closeness | stats_to_include$all_stats) {
    data$closeness <- igraph::closeness(net, mode = "all", weights = NA)
  }
  if (stats_to_include$centrality_by_closeness | stats_to_include$all_stats) {
    data$centrality_by_closeness <-
      igraph::centr_clo(net, mode = "all", normalized = T)$res
  }
  if (stats_to_include$eigen_centrality | stats_to_include$all_stats) {
    data$eigen_centrality <-
      igraph::eigen_centrality(net, directed = T, weights = NA)$vector
  }
  if (stats_to_include$centrality_by_eigen | stats_to_include$all_stats) {
    data$centrality_by_eigen <-
      igraph::centr_eigen(net, directed = T, normalized = T)$vector
  }
  if (stats_to_include$betweenness | stats_to_include$all_stats) {
    data$betweenness <- igraph::betweenness(net, directed = T, weights = NA) }

  if (stats_to_include$centrality_by_betweenness | stats_to_include$all_stats) {
    data$centrality_by_betweenness <-
      igraph::centr_betw(net, directed = T, normalized = T)$res
  }
  if (stats_to_include$authority_score | stats_to_include$all_stats) {
    data$authority_score <- igraph::authority_score(net, weights = NA)$vector
  }
  if (stats_to_include$coreness | stats_to_include$all_stats) {
    data$coreness <- igraph::coreness(net, mode = "all")
  }
  if (stats_to_include$page_rank | stats_to_include$all_stats) {
    data$page_rank <- igraph::page_rank(net)$vector
  }
  cat(" Done.\n")
  return(data)
}

# Create list of which node-level network statistics to compute; to be used as
# value for argument `stats_to_include` of function `computeNodeNetworkStats()`
node_stat_settings <- function(
    degree = TRUE,
    cluster_id = FALSE,
    transitivity = TRUE,
    closeness = FALSE,
    centrality_by_closeness = FALSE,
    eigen_centrality = TRUE,
    centrality_by_eigen = TRUE,
    betweenness = TRUE,
    centrality_by_betweenness = TRUE,
    authority_score = TRUE,
    coreness = TRUE,
    page_rank = TRUE,
    all_stats = FALSE
) {
  list(degree = degree,
       cluster_id = cluster_id,
       transitivity = transitivity,
       closeness = closeness,
       centrality_by_closeness = centrality_by_closeness,
       eigen_centrality = eigen_centrality,
       centrality_by_eigen = centrality_by_eigen,
       betweenness = betweenness,
       centrality_by_betweenness = centrality_by_betweenness,
       authority_score = authority_score,
       coreness = coreness,
       page_rank = page_rank,
       all_stats = all_stats)
}

# Computing Cluster Statistics --------------------------------------------

# FUNCTION: Compute the clusters for a network and augment the corresponding
# data with a variable containing the cluster membership ID
addClusterMembership <- function(data, net) {
  cat("Computing cluster membership within the network...")
  data$cluster_id <-
    as.factor(as.integer(igraph::cluster_fast_greedy(net)$membership))
  cat(" Done.\n")
  return(data)
}

# FUNCTION: Compute cluster-level network stats
getClusterStats <- function(
    data, # rep-seq data for network, with node-level network stats
    adjacency_matrix, # adjacency matrix for network
    seq_col, # name or number of column of `data` containing the clone sequences
    count_col = NULL, # name or number of column of `data` containing the clone counts
    cluster_id_col = NULL, # optional name or number of column of `data` containing the cluster IDs
    degree_col = NULL, # optional name or number of column of `data` containing the network degree
    seq_length_col = NULL # optional name or number of col containing seq lengths
) {
  if (is.null(cluster_id_col) | is.null(degree_col)) {
    net <- generateNetworkFromAdjacencyMat(adjacency_matrix)
    if (is.null(cluster_id_col)) { # compute cluster id if not provided
      cluster_id_col <- "cluster_id"
      data <- addClusterMembership(data, net)
    }
    if (is.null(degree_col)) { # compute deg if not provided
      degree_col <- "deg"
      data$deg <- igraph::degree(net)
    }
  }
  if (is.null(seq_length_col)) { # compute seq length if not provided
    seq_length_col <- "seq_length"
    data$seq_length <- nchar(data[ , seq_col])
  }

  # Tabulate the number of nodes in each cluster
  out <- as.data.frame(table(data[ , cluster_id_col]))
  colnames(out) <- c("cluster_id", "node_count")
  num_clusters <- nrow(out) # Total number of clusters
  cat(paste0("Computing statistics for the ", num_clusters, " clusters in the network..."))

  ### INITIALIZE VALUES ###
  out$mean_seq_length <- 0
  out$mean_degree <- 0
  out$max_degree <- ""
  out$seq_w_max_degree <- ""
  out$agg_clone_count <- NA
  out$max_clone_count <- NA
  out$seq_w_max_count <- NA
  out$diameter_length <- 0
  out$assortativity <- 0
  out$transitivity <- 0
  out$edge_density <- 0
  out$degree_centrality_index <- 0
  out$closeness_centrality_index <- 0
  out$eigen_centrality_index <- 0
  out$eigen_centrality_eigenvalue <- 0

  ### COMPUTE STATS FOR EACH CLUSTER ###
  for (i in 1:num_clusters) {
    cluster_row <- which(out$cluster_id == i) # current row of cluster data
    node_ids <- data$cluster_id == i  # Rows of node data for current cluster

    # Mean sequence length in cluster
    out$mean_seq_length[[cluster_row]] <-
      round(mean(data[node_ids, seq_length_col]), 2)

    # Mean degree in cluster
    out$mean_degree[[cluster_row]] <-
      round(mean(data[node_ids, degree_col]), 2)

    # Maximum degree (and corresponding seq) within cluster
    max_deg <- max(data[node_ids, degree_col])
    out$max_degree[[cluster_row]] <- max_deg

    node_id_max_deg <- which(node_ids & data[ , degree_col] == max_deg)
    out$seq_w_max_degree[[cluster_row]] <-
      as.character(data[node_id_max_deg, seq_col][[1]])

    if (!is.null(count_col)) {
      # Total aggregate clonotype count in cluster
      out$agg_clone_count[[cluster_row]] <- sum(data[node_ids, count_col])

      # Maximum clonotype count (and corresponding seq) within cluster
      max_count <- max(data[node_ids, count_col])
      out$max_clone_count[[cluster_row]] <- max_count

      node_id_max_count <- which(node_ids & data[ , count_col] == max_count)
      out$seq_w_max_count[[cluster_row]] <-
        as.character(data[node_id_max_count, seq_col][[1]])
    }

    # Build cluster network to get network properties for the cluster
    cluster_adjacency_matrix <- as.matrix(adjacency_matrix[node_ids, node_ids])
    cluster <- generateNetworkFromAdjacencyMat(cluster_adjacency_matrix)

    # Diameter (longest geodesic distance)
    out$diameter_length[[cluster_row]] <-
      length(igraph::get_diameter(cluster, directed = T))

    # Assortativity
    out$assortativity[[cluster_row]] <-
      igraph::assortativity_degree(cluster, directed = F)

    # Transitivity
    out$transitivity[[cluster_row]] <-
      igraph::transitivity(cluster, type = "global")  # cluster is treated as an undirected network

    # Density: The proportion of present edges from all possible ties.
    out$edge_density[[cluster_row]] <-
      igraph::edge_density(cluster, loops = F)

    # Centralization on degree
    out$degree_centrality_index[[cluster_row]] <-
      igraph::centr_degree(cluster, mode = "in", normalized = T)$centralization

    # Centralization on Closeness (centrality based on distance to others in the graph)
    out$closeness_centrality_index[[cluster_row]] <-
      igraph::centr_clo(cluster, mode = "all", normalized = T)$centralization

    # Centralization on Eigenvector (centrality proportional to the sum of connection centralities)
    #  (values of the first eigenvector of the graph adjacency matrix)
    out$eigen_centrality_index[[cluster_row]] <-
      igraph::centr_eigen(cluster, directed = T, normalized = T)$centralization

    out$eigen_centrality_eigenvalue[[cluster_row]] <-
      igraph::eigen_centrality(cluster, directed = T, weights = NA)$value
  }
  cat(" Done.\n")

  return(out)
}




# Plotting Network Graphs -------------------------------------------------


plotNetworkGraph <- function(network, edge_width = 0.3,
                             title = NULL,
                             subtitle = NULL,
                             color_nodes_by = NULL,
                             color_scheme = "default",
                             show_color_legend = TRUE,
                             color_legend_title = "auto",
                             size_nodes_by = NULL,
                             node_size_limits = NULL,
                             size_legend_title = "auto",
                             outfile = NULL
) {
  set.seed(9999)
  layout <- igraph::layout_components(network)

  graph_plot <-
    ggraph::ggraph(network, layout = layout) +
    ggraph::geom_edge_link0(width = edge_width, colour = "grey")

  if (!is.null(color_nodes_by)) {
    # Custom node color scheme
    if (!is.null(size_nodes_by)) {
      # Custom node size scheme
      if (is.numeric(size_nodes_by) & length(size_nodes_by) == 1) {
        # Fixed node sizes
        graph_plot <- graph_plot +
          ggraph::geom_node_point(
            ggplot2::aes(color = color_nodes_by), size = size_nodes_by)
      } else if (length(size_nodes_by) > 1) {
        # size_nodes_by is a numeric vector specifying the size of each node
        graph_plot <- graph_plot +
          ggraph::geom_node_point(
            ggplot2::aes(color = color_nodes_by, size = size_nodes_by))
        if (length(node_size_limits) == 2) {
          # Rescale node sizes if specified
          graph_plot <-
            graph_plot + ggplot2::scale_size(range = node_size_limits)
        }
      }
    } else { # size_nodes_by is null or invalid
      graph_plot <- graph_plot +
        ggraph::geom_node_point(ggplot2::aes(color = color_nodes_by))

    }
  } else { # color_nodes_by is null or invalid
    if (!is.null(size_nodes_by)) {
      if (is.numeric(size_nodes_by) & length(size_nodes_by) == 1) {
        graph_plot <- graph_plot +
          ggraph::geom_node_point(size = size_nodes_by)

      } else if (length(size_nodes_by) > 1) {
        graph_plot <- graph_plot +
          ggraph::geom_node_point(ggplot2::aes(size = size_nodes_by))
        if (!is.null(node_size_limits)) {
          graph_plot <-
            graph_plot + ggplot2::scale_size(range = node_size_limits)
        }
      }
    } else { # size_nodes_by null or invalid
      graph_plot <- graph_plot + ggraph::geom_node_point()
    }
  }

  graph_plot <- graph_plot +
    ggraph::theme_graph(base_family = "sans") +
    ggplot2::labs(title = title, subtitle = subtitle)

  if (show_color_legend) {
    if (is.null(color_legend_title)) {
      graph_plot <- graph_plot +
        ggplot2::guides(color = ggplot2::guide_legend(title = color_legend_title))
    } else if (color_legend_title != "auto") {
      graph_plot <- graph_plot +
        ggplot2::guides(color = ggplot2::guide_legend(title = color_legend_title))
    }
  } else {
    graph_plot <- graph_plot + ggplot2::guides(color = "none")
  }

  if (is.null(size_legend_title)) {
    graph_plot <- graph_plot +
      ggplot2::guides(size = ggplot2::guide_legend(title = size_legend_title))
  } else if (size_legend_title != "auto") {
    graph_plot <- graph_plot +
      ggplot2::guides(size = ggplot2::guide_legend(title = size_legend_title))
  }

  if (is.null(color_nodes_by)) {
    color_type <- "continuous"
  } else {
    color_type <- ggplot2::scale_type(color_nodes_by)[[1]]
  }

  # Convert node-color variable to factor if discrete
  # if (color_type == "discrete") { color_nodes_by <- as.factor(color_nodes_by) }

  if (color_scheme != "default") {

    if (color_type == "continuous") {
      if (color_scheme %in% c("A", "B", "C", "F", "G",
                              "magma", "inferno", "plasma", "rocket", "mako")) {
        graph_plot <- graph_plot +
          ggraph::scale_color_viridis(option = color_scheme,
                                      begin = 0.2, end = 0.8)
      } else if (color_scheme %in% c("D", "viridis")) {
        graph_plot <- graph_plot +
          ggraph::scale_color_viridis(option = color_scheme,
                                      begin = 0, end = 0.9)
      } else if (color_scheme %in% c("E", "H", "cividis",  "turbo")) {
        graph_plot <- graph_plot +
          ggraph::scale_color_viridis(option = color_scheme)
      } else if (color_scheme %in% c("A-1", "B-1", "C-1", "F-1", "G-1",
                                     "magma-1", "inferno-1", "plasma-1", "rocket-1", "mako-1")) {
        graph_plot <- graph_plot +
          ggraph::scale_color_viridis(option = strsplit(color_scheme, "-1")[[1]],
                                      begin = 0.2, end = 0.8, direction = -1)
      } else if (color_scheme %in% c("D-1", "viridis-1")) {
        graph_plot <- graph_plot +
          ggraph::scale_color_viridis(option = strsplit(color_scheme, "-1")[[1]],
                                      begin = 0, end = 0.9, direction = -1)
      } else if (color_scheme %in% c("E-1", "H-1", "cividis-1",  "turbo-1")) {
        graph_plot <- graph_plot +
          ggraph::scale_color_viridis(option = strsplit(color_scheme, "-1")[[1]], direction = -1)
      } else { warning("value for 'color_scheme' is not a valid option for continuous variables; using default color scheme instead") }

    } else { # discrete color scheme
      if (color_scheme %in% c("A", "B", "C", "F", "G",
                              "magma", "inferno", "plasma", "rocket", "mako")) {
        graph_plot <- graph_plot +
          ggraph::scale_color_viridis(option = color_scheme, discrete = TRUE,
                                      begin = 0.2, end = 0.8)
      } else if (color_scheme %in% c("D", "viridis")) {
        graph_plot <- graph_plot +
          ggraph::scale_color_viridis(option = color_scheme, discrete = TRUE,
                                      begin = 0, end = 0.9)
      } else if (color_scheme %in% c("E", "H", "cividis",  "turbo")) {
        graph_plot <- graph_plot +
          ggraph::scale_color_viridis(option = color_scheme, discrete = TRUE,)
      } else if (color_scheme %in% c("A-1", "B-1", "C-1", "F-1", "G-1",
                                     "magma-1", "inferno-1", "plasma-1", "rocket-1", "mako-1")) {
        graph_plot <- graph_plot +
          ggraph::scale_color_viridis(option = strsplit(color_scheme, "-1")[[1]], discrete = TRUE,
                                      begin = 0.2, end = 0.8, direction = -1)
      } else if (color_scheme %in% c("D-1", "viridis-1")) {
        graph_plot <- graph_plot +
          ggraph::scale_color_viridis(option = strsplit(color_scheme, "-1")[[1]], discrete = TRUE,
                                      begin = 0, end = 0.9, direction = -1)
      } else if (color_scheme %in% c("E-1", "H-1", "cividis-1",  "turbo-1")) {
        graph_plot <- graph_plot +
          ggraph::scale_color_viridis(option = strsplit(color_scheme, "-1")[[1]], discrete = TRUE,
                                      direction = -1)
      } else if (color_scheme %in% grDevices::hcl.pals()) {
        graph_plot <- graph_plot +
          ggplot2::scale_color_manual(
            values = grDevices::hcl.colors(n = length(unique(color_nodes_by)),
                                           palette = color_scheme))
      } else { warning("value for 'color_scheme' is not a valid option for discrete variables; using default color scheme instead") } }
  }


  if (!is.null(outfile)) {
    grDevices::pdf(file = outfile, width = 12, height = 8)
    print(graph_plot)
    grDevices::dev.off()
    cat(paste0("Plot of network graph saved to file:\n  ", outfile, "\n"))
  }
  return(graph_plot)
}




# Adjacency and Distance Matrices -----------------------------------------

# FUNCTION: COMPUTE ADJACENCY MATRIX FOR LEVENSHTEIN OR HAMMING DISTANCE
# Intended for use with a large network where the adjacency matrix is sparse
# Returns sparse matrix, includes only nodes with positive network degree
sparseAdjacencyMatFromSeqs <- function(
    seqs, # List of tcr/clonotype sequences
    dist_type = "hamming", # supports "levenshtein" and "hamming"
    max_dist = 1, # Maximum distance threshold for edge/adjacency between two sequences
    drop_isolated_nodes = TRUE # Drop sequences/nodes with zero degree?
) {
  # attempt to coerce seqs to character vector
  if (length(seqs) == 0) stop("'seqs' has zero length")
  seqs <- as.vector(seqs, mode = "character")
  if (!is.character(seqs)) stop("'seqs' must be cocercible to a character vector")
  if (!is.vector(seqs)) stop("'seqs' must be cocercible to a character vector")

  # Compute adjacency matrix
  if (dist_type %in% c("levenshtein", "Levenshtein, lev, Lev, l, L")) {
    cat(paste0("Computing network edges based on a max ", dist_type, " distance of ", max_dist, "..."))
    out <- levAdjacencyMatSparse(seqs, max_dist, drop_isolated_nodes)
  } else if (dist_type %in% c("hamming", "Hamming", "ham", "Ham", "h", "H")) {
    cat(paste0("Computing network edges based on a max ", dist_type, " distance of ", max_dist, "..."))
    out <- hamAdjacencyMatSparse(seqs, max_dist, drop_isolated_nodes)
  } else {
    stop('invalid option for `dist_type`')
  }
  cat(" Done.\n")
  # Number of nodes with positive network degree
  num_nodes <- dim(out)[[1]]
  if (num_nodes == 0) {
    warning(paste0(
      "No edges exist using the specified distance cutoff; try a greater value of `max_dist`"))
  } else {
    if (drop_isolated_nodes) {
      cat(paste0(num_nodes, " nodes in the network (after removing nodes with degree zero).\n"))
      # Import record of selected column IDs and use for matrix row names
      clone_ids <- utils::read.table("col_ids.txt")
      dimnames(out)[[1]] <- clone_ids$V1
      dimnames(out)[[2]] <- seqs[clone_ids$V1]
      # cat(paste0("The row names of the adjacency matrix contain the original index values of the corresponding sequences; the column names contain the sequences themselves. They can be accessed using `dimnames()`\n"))

      # Remove temporary file of column ids
      file.remove("col_ids.txt")
    } else {
      cat(paste0(num_nodes, " nodes in the network.\n"))
    }
  }
  return(out)
}


# Adjacency Matrix: Euclidean Distance on Atchley Factor Embedding
# (Only applicable to TCR CDR3 Amino Acid Sequences)
# This function is intended for building the network for a single cluster, where
# the adjacency matrix is typically dense
adjacencyMatAtchleyFromSeqs <- function(
    seqs, # List of TCR CDR3 amino acid sequences corresponding to the seqs
    contig_ids = seq_along(seqs), # used by BriseisEncoder to perform the Atchley-factor embedding of the TCR sequences
    max_dist, # Maximum Euclidean distance threshold for edge/adjacency between two sequences
    return_type = "adjacency_matrix", # can be set to "distance_matrix" to return the distance matrix instead
    outfile_distance_matrix = NULL # savefile for Euclidean distance matrix
) {
  # Embed amino acid seqs in Euclidean 30-space by Atchley factor representation
  embedded_values <- embedTCRSeqsByAtchleyFactor(seqs, contig_ids)

  # Compute Euclidean distance matrix on embedded sequence values
  cat("Computing Euclidean distances between the embedded values...")
  distance_matrix <- as.matrix(stats::dist(embedded_values[ , -1]))
  cat(" Done.\n")

  if (!is.null(outfile_distance_matrix)) {
    # Save distance matrix to file
    utils::write.csv(distance_matrix, outfile_distance_matrix)
    cat(paste0("Distance matrix saved to file:\n  ", outfile_distance_matrix,
               "\n")) }

  if (return_type == "distance_matrix") {
    return(distance_matrix)
  } else {
    # Convert distance matrix to adjacency matrix using specified bound
    cat(paste0("Generating adjacency matrix based on a maximum distance of ",
               max_dist, "..."))
    adjacency_matrix <-
      matrix(1, nrow = nrow(distance_matrix), ncol = ncol(distance_matrix))
    adjacency_matrix[distance_matrix > max_dist] <- 0
    cat(" Done.\n")
    return(adjacency_matrix)
  }
}



# Atchley Factor Embedding ------------------------------------------------

# Embed TCR CDR3 amino acid sequences in Euclidean 30-space based on the Atchley
# factor representations of their elements, using a trained encoding model
embedTCRSeqsByAtchleyFactor <- function(
    cdr3_AA, # List of TCR CDR3 amino acid sequences
    contig_ids = seq_along(cdr3_AA) # used by BriseisEncoder
) {
  .checkPythonModules()
  if (length(cdr3_AA) != length(contig_ids)) {
    stop("length of `cdr3_AA` and `contig_ids` must match")
  }
  # Write sequences and contig_ids to temporary file
  tempfile_clones <- file.path(getwd(),
                               "Atchley_factor_tcr_only_tmp.csv")
  utils::write.csv(data.frame("contig_id" = contig_ids, "cdr3" = cdr3_AA),
                   file = tempfile_clones, row.names = FALSE)

  # Write files for trained encoder and Atchley factor table to working dir
  file_trained_encoder <-
    system.file(file.path("python", "TrainedEncoder.h5"),
                package = "NAIR")
  file_atchley_table <-
    system.file(file.path("python", "Atchley_factors.csv"),
                package = "NAIR")
  utils::write.csv(data.frame("sysfiles" = c(file_atchley_table,
                                             file_trained_encoder)),
                   file.path(getwd(), "temp_sysfiles.csv"),
                   row.names = FALSE)

  # Run BriseisEncoder python function
  cat("Embedding TCR CDR3 amino acid sequences in 30-dimensional Euclidean space based on Atchley factor representation and a trained encoding model using deep learning routines from the keras & tensorflow python modules. This may produce some warnings...\n")
  reticulate::py_run_file(
    system.file(file.path("python", "BriseisEncoder_modified.py"),
                package = "NAIR"))

  # Import embedded values and remove temp files
  embedded_values <- utils::read.csv("temp_atchley_factors_encoded.csv")
  file.remove("Atchley_factor_tcr_only_tmp.csv", "temp_sysfiles.csv",
              "temp_atchley_factors_encoded.csv")
  cat("Embedding complete.\n")
  warning("the encoder was trained on TCR CDR3 sequences; results not valid for other amino acid sequences")
  return(embedded_values)
}
