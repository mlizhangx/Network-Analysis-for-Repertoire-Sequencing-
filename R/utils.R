
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
aggregateCounts <- function(
  data, # data frame containing columns below
  count_col, # name or number of column of `data` containing clone counts
  group_cols = NULL # optional integer or character vector specifying additional grouping columns
) {

  # Convert column specifications from numeric to character if not already
  if (is.numeric(clone_col)) { clone_col <- names(data)[clone_col] }
  if (is.numeric(count_col)) { count_col <- names(data)[count_col] }

  # Define grouping variable(s)
  grouping_variables <- list(data[ , clone_col])
  names(grouping_variables) <- clone_col
  if (!is.null(group_cols)) {
    if (is.numeric(group_cols)) { group_cols <- names(data)[group_cols] }
    for (i in 1:length(group_cols)) {
      grouping_variables$newvar <- data[ , group_cols[[i]]]
      names(grouping_variables)[[length(grouping_variables)]] <- group_cols[[i]]
    }
  }

  # aggregate the counts by group
  counts <- list(data[ , count_col])
  names(counts) <- count_col
  agg_counts <- stats::aggregate(counts, by = grouping_variables, FUN = sum)

  # add variable for num reads (row count)
  groups <- as.data.frame(data[ , group_cols])
  names(groups)[[1]] <- "temporary_placeholder_name" # for summarize function
  num_reads <- dplyr::summarize(dplyr::group_by_all(groups),
                                numReads = length(temporary_placeholder_name))
  names(num_reads)[[1]] <- group_cols[[1]] # replace placeholder name with orig

  # Merge aggregate counts with num reads
  out <- merge(agg_counts, num_reads, by = group_cols)

  return(out)
}



# FUNCTION: Filter rep-seq data to remove rows for clonotype sequences with
# length below the specified cutoff
filterDataBySequenceLength <- function(data, clone_col, min_length = 3) {
  return(data[nchar(data[ , clone_col]) >= min_length, ])
}


# FUNCTION: EXTRACT DATA SUBSET FOR ALL SEQUENCES WITHIN SPECIFIED RADIUS OF
# TARGET SEQUENCE BY SPECIFIED DISTANCE TYPE
# If a sample_col is provided, only samples that possess the target sequence
# will be included
getSimilarClones <- function(
  target_seq, # specified candidate sequence for the neighborhood
  data, # data frame containing rep seq data, possibly from multiple samples
  clone_col, # col name/# containing clone sequences
  sample_col = NULL, # optional col name/# containing sample IDs (only samples possessing target seq will be included)
  dist_type = "hamming", # options are "hamming" and "levenshtein"
  max_dist = 2 # Maximum Levenshtein distance allowed for inclusion in neighborhood
) {
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }
  # If sample_id is supplied, subset data keeping only samples with target seq
  if (is.null(sample_col)) {
    data_samples_w_targetseq <- data
  } else {
    ### SUBSET DATA: SAMPLES WITH TARGET SEQ ###
    cat("Finding samples that possess the specified sequence...")
    # Get row ids of merged data corresponding to target seq
    rows_for_targetseq <- grep(pattern = paste0("^", target_seq, "$"),
                               x = data[ , clone_col])
    # Extract rows of merged data corresponding to samples with target seq
    data_samples_w_targetseq <-
      data[
        data[ , sample_col] %in% data[rows_for_targetseq, sample_col], ]
    cat("Done.\n")
  }
  # remove seq with * and _
  data_samples_w_targetseq <-
    data_samples_w_targetseq[
      -grep("[*|_]", data_samples_w_targetseq[ , clone_col]), ]
  ### SUBSET DATA: NEIGHBORHOOD OF TARGET SEQUENCE ###
  # Compute list of bounded distances between target seq and seqs
  # possessed by samples with target seq (values are -1 where bound is exceeded)
  # returned vector will be of type integer; names will be the sequences
  cat("Extracting clone sequences similar to the specified sequence...")
  if (dist_type == "levenshtein") { dist_fun <- levDistBounded
  } else if (dist_type == "hamming") { dist_fun <- hamDistBounded
  } else { stop("invalid option for `dist_type`") }
  dists_to_targetseq <- sapply(
    X = data_samples_w_targetseq[ , clone_col],
    FUN = dist_fun, b = target_seq, k = max_dist)
  # get data for sequences within the specified radius
  data_targetseq_neighborhood <-
    data_samples_w_targetseq[dists_to_targetseq != -1, ]
  cat("Done.\n")
  return(data_targetseq_neighborhood)
}


# Network Building --------------------------------------------------------



# FUNCTION: GENERATE NETWORK FOR A LIST OF CLONE SEQS USING SPECIFIED DISTANCE TYPE AND THRESHOLD
generateNetworkFromClones <- function(
  clones, # list of clonotype sequences
  dist_type = "levenshtein", # supports "levenshtein", "hamming", "euclidean_on_atchley"
  edge_dist = 1, # max dist threshold for edges
  contig_ids = seq_along(clones), # for "euclidean_on_atchley" dist_type only
  outfile_adjacency_matrix = NULL, # save file for adjacency matrix
  outfile_distance_matrix = NULL, # save file for distance matrix (only for Euclidean on Atchley)
  return_type = "network" # can use "adjacency_matrix" to return the adjacency mat
) {
  ### COMPUTE ADJACENCY MATRIX ###
  if (dist_type %in% c("levenshtein", "hamming")) {
    adjacency_matrix <-
      sparseAdjacencyMatFromClones(clones = clones,
                                   dist_type = dist_type,
                                   max_dist = edge_dist)
    if (!is.null(outfile_adjacency_matrix)) {
      Matrix::writeMM(adjacency_matrix, outfile_adjacency_matrix)
    }
  } else if (dist_type == "euclidean_on_atchley") {
    adjacency_matrix <-
      adjacencyMatAtchleyFromClones(
        clones = clones,
        contig_ids = contig_ids,
        max_dist = edge_dist,
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
computeNodeNetworkStats <- function(
  net, # igraph network object
  data, # rep-seq data corresponding to the network
  stats_to_include = node_stat_settings()
) {
  cat(paste0("Computing node-level network characteristics...\n"))
  if (stats_to_include$degree) data$degree <- igraph::degree(net)

  if (stats_to_include$cluster_id) data$cluster_id <-
      igraph::cluster_fast_greedy(net)$membership

  if (stats_to_include$transitivity) data$transitivity <-
      igraph::transitivity(net, type = "local")

  if (stats_to_include$closeness) data$closeness <-
      igraph::closeness(net, mode = "all", weights = NA)

  if (stats_to_include$centrality_by_closeness) data$centrality_by_closeness <-
      igraph::centr_clo(net, mode = "all", normalized = T)$res

  if (stats_to_include$eigen_centrality) data$eigen_centrality <-
      igraph::eigen_centrality(net, directed = T, weights = NA)$vector

  if (stats_to_include$centrality_by_eigen) data$centrality_by_eigen <-
      igraph::centr_eigen(net, directed = T, normalized = T)$vector

  if (stats_to_include$betweenness) data$betweenness <-
      igraph::betweenness(net, directed = T, weights = NA)

  if (stats_to_include$centrality_by_betweenness) {
    data$centrality_by_betweenness <-
      igraph::centr_betw(net, directed = T, normalized = T)$res }

  if (stats_to_include$authority_score) data$authority_score <-
      igraph::authority_score(net, weights = NA)$vector

  if (stats_to_include$coreness) data$coreness <-
      igraph::coreness(net, mode = "all")

  if (stats_to_include$page_rank) data$page_rank <-
      igraph::page_rank(net)$vector

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
  page_rank = TRUE
) {
  list(deg = deg,
       cluster_id = cluster_id,
       transitivity = transitivity,
       closeness = closeness,
       centr_clo_res = centr_clo_res,
       eigen_centrality = eigen_centrality,
       centr_eigen = centr_eigen,
       betweenness = betweenness,
       centr_betw = centr_betw,
       authority_score = authority_score,
       coreness = coreness,
       page_rank = page_rank)
}

# Computing Cluster Statistics --------------------------------------------

# FUNCTION: Compute the clusters for a network and augment the corresponding
# data with a variable containing the cluster membership ID
addClusterMembership <- function(net, data) {
  data$cluster_id <- igraph::cluster_fast_greedy(net)$membership
  return(data)
}

# FUNCTION: Compute cluster-level network stats
computeClusterNetworkStats <- function(
  adjacency_matrix, # adjacency matrix for network
  data, # rep-seq data for network, with node-level network stats
  clone_col, # name or number of column of `data` containing the clone sequences
  count_col, # name or number of column of `data` containing the clone counts
  cluster_id_col = NULL, # optional name or number of column of `data` containing the cluster IDs
  deg_col = NULL, # optional name or number of column of `data` containing the network degree
  seq_length_col = NULL # optional name or number of col containing seq lengths
) {
  if (is.null(cluster_id_col) | is.null(deg_col)) {
    net <- generateNetworkFromAdjacencyMat(adjacency_matrix)
    if (is.null(cluster_id_col)) { # compute cluster id if not provided
      cluster_id_col <- "cluster_id"
      data <- addClusterMembership(net, data)
    }
    if (is.null(deg_col)) { # compute deg if not provided
      deg_col <- "deg"
      data$deg <- igraph::degree(net)
    }
  }
  if (is.null(seq_length_col)) { # compute seq length if not provided
    seq_length_col <- "seq_length"
    data$seq_length <- nchar(data[ , clone_col])
  }

  # Tabulate the number of nodes in each cluster
  out <- as.data.frame(table(data[ , cluster_id_col]))
  colnames(out) <- c("cluster_id", "node_count")
  num_clusters <- nrow(out) # Total number of clusters
  cat(paste0(num_clusters, " clusters present in network. Computing cluster characteristics...\n"))

  ### INITIALIZE VALUES ###
  out$mean_seq_length <- 0
  out$deg_mean <- 0
  out$deg_max <- ""
  out$motif_w_max_deg <- ""
  out$agg_clone_count <- 0
  out$max_clone_count <- ""
  out$motif_w_max_count <- ""
  out$diam_length <- 0
  out$assortativity <- 0
  out$transitivity <- 0
  out$edge_density <- 0
  out$centr_degree <- 0
  out$centr_clo <- 0
  out$eigen_centrality <- 0
  out$centr_eigen <- 0

  ### COMPUTE STATS FOR EACH CLUSTER ###
  for(i in 1:num_clusters) {
    cluster_index <- which(out$cluster_id == i) # current row of cluster data
    node_ids <- data$cluster_id == i  # Rows of node data for current cluster


    # Mean sequence length in cluster
    out[cluster_index, ]$mean_seq_length <-
      round(mean(data[node_ids, seq_length_col]), 2)

    # Mean degree in cluster
    out[cluster_index, ]$deg_mean <- round(mean(data[node_ids, ]$deg), 2)

    # Maximum degree (and corresponding seq) within cluster
    max_deg <- max(data[node_ids, ]$deg)
    out[cluster_index, ]$deg_max <- max_deg
    node_id_max_deg <- which(node_ids & data$deg == max_deg)
    out[cluster_index, ]$motif_w_max_deg  <-
      as.character(data[node_id_max_deg, clone_col][[1]])

    # Total aggregate clonotype count in cluster
    out[cluster_index, ]$agg_clone_count <- sum(data[node_ids, count_col])

    # Maximum clonotype count (and corresponding seq) within cluster
    max_count <- max(data[node_ids, count_col])
    out[cluster_index, ]$max_clone_count <- max_count
    node_id_max_count <- which(node_ids & data[ , count_col] == max_count)
    out[cluster_index, ]$motif_w_max_count <-
      as.character(data[node_id_max_count, clone_col][[1]])

    # Build cluster network to get network properties for the cluster
    cluster_adjacency_matrix <- as.matrix(adjacency_matrix[node_ids, node_ids])
    cluster <- generateNetworkFromAdjacencyMat(cluster_adjacency_matrix)
    # Diameter (longest geodesic distance)
    out[cluster_index, ]$diam_length <-
      length(igraph::get_diameter(cluster, directed = T))
    # Assortativity
    out[cluster_index, ]$assortativity <-
      igraph::assortativity_degree(cluster, directed = F)
    # Transitivity
    out[cluster_index, ]$transitivity <-
      igraph::transitivity(cluster, type = "global")  # cluster is treated as an undirected network
    # Density: The proportion of present edges from all possible ties.
    out[cluster_index, ]$edge_density <-
      igraph::edge_density(cluster, loops = F)
    # Centralization on degree
    out[cluster_index, ]$centr_degree <-
      igraph::centr_degree(cluster, mode = "in", normalized = T)$centralization
    # Centralization on Closeness (centrality based on distance to others in the graph)
    out[cluster_index, ]$centr_clo <-
      igraph::centr_clo(cluster, mode = "all", normalized = T)$centralization
    # Centralization on Eigenvector (centrality proportional to the sum of connection centralities)
    #  (values of the first eigenvector of the graph adjacency matrix)
    out[cluster_index, ]$eigen_centrality <-
      igraph::eigen_centrality(cluster, directed = T, weights = NA)$value
    out[cluster_index, ]$centr_eigen <-
      igraph::centr_eigen(cluster, directed = T, normalized = T)$centralization
  }
  # Return cluster-level info
  return(out)
}




# Plotting Network Graphs -------------------------------------------------


plotNetworkGraph <- function(network, edge_width = 0.3,
                             title = NULL,
                             subtitle = NULL,
                             color_nodes_by,
                             size_nodes_by,
                             color_legend_title = NULL,
                             size_legend_title = NULL,
                             color_scheme = "default",
                             outfile = NULL
) {
  set.seed(9999)
  layout <- igraph::layout_components(network)
  graph_plot <-
    ggraph::ggraph(network, layout = layout) +
    ggraph::geom_edge_link0(width = edge_width, colour = "grey") +
    ggraph::geom_node_point(
      ggplot2::aes(color = color_nodes_by, size = size_nodes_by)) +
    ggraph::theme_graph(base_family = "sans") +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::guides(color = ggplot2::guide_legend(title = color_legend_title),
                    size = ggplot2::guide_legend(title = size_legend_title))
  if (color_scheme != "default") {
    if (ggplot2::scale_type(color_nodes_by)[[1]] == "continuous") {
      if (color_scheme %in% c("A", "B", "C", "D", "E", "F", "G", "H",
                               "magma", "inferno", "plasma", "viridis",
                               "cividis", "rocket", "mako", "turbo")) {
        graph_plot <- graph_plot +
          ggraph::scale_color_viridis(option = color_scheme)
      } else { warning("'color_nodes_by' is continuous; 'color_scheme' must be 'default' or a viridis color map option (see `?viridis`); using default color scheme instead") } } else {
      if (color_scheme %in% grDevices::hcl.pals()) {
        graph_plot <- graph_plot +
          ggplot2::scale_color_manual(
            values = grDevices::hcl.colors(n = length(color_nodes_by),
                                           palette = color_scheme))
      } else { warning("'color_scheme' must be 'default' or one of the values contained in `grDevices::hcl.pals()`; using default color scheme instead") } }
  }
  if (!is.null(outfile)) {
    grDevices::pdf(file = outfile, width = 12, height = 8)
    print(graph_plot)
    grDevices::dev.off()
    cat(paste0("Plot of network graph saved as '", outfile, "'\n"))
  }
  return(graph_plot)
}




# Adjacency and Distance Matrices -----------------------------------------

# FUNCTION: COMPUTE ADJACENCY MATRIX FOR LEVENSHTEIN OR HAMMING DISTANCE
# Intended for use with a large network where the adjacency matrix is sparse
# Returns sparse matrix, includes only nodes with positive network degree
sparseAdjacencyMatFromClones <- function(
  clones, # List of tcr/clonotype sequences
  dist_type = "hamming", # supports "levenshtein" and "hamming"
  max_dist = 1 # Maximum distance threshold for edge/adjacency between two sequences
  # drop_isolated_nodes = TRUE # Drop sequences/nodes with zero degree?
) {
  # attempt to coerce clones to character vector
  if (length(clones) == 0) stop("'clones' has zero length")
  clones <- as.vector(clones, mode = "character")
  if (!is.character(clones)) stop("'clones' must be cocercible to a character vector")
  if (!is.vector(clones)) stop("'clones' must be cocercible to a character vector")

  # Compute adjacency matrix
  if (dist_type %in% c("levenshtein", "Levenshtein, lev, Lev, l, L")) {
    cat(paste0("Computing network edges based on a max ", dist_type, " distance of ", max_dist, "..."))
    out <- levAdjacencyMatSparse(clones, max_dist)
  } else if (dist_type %in% c("hamming", "Hamming", "ham", "Ham", "h", "H")) {
    cat(paste0("Computing network edges based on a max ", dist_type, " distance of ", max_dist, "..."))
    out <- hamAdjacencyMatSparse(clones, max_dist)
  } else {
    stop('invalid option for `dist_type`')
  }
  cat("Done.\n")
  # Number of nodes with positive network degree
  num_nodes <- dim(out)[[1]]
  if (num_nodes == 0) {
    warning(paste0(
      "No edges exist using the specified distance cutoff; try a greater value of `max_dist`"))
  } else {
    cat(paste0(num_nodes, " nodes in the network (after removing nodes with degree zero).\n"))
    # Import record of selected column IDs and use for matrix row names
    clone_ids <- utils::read.table("col_ids.txt")
    dimnames(out)[[1]] <- clone_ids$V1
    dimnames(out)[[2]] <- clones[clone_ids$V1]
    # cat(paste0("The row names of the adjacency matrix contain the original index values of the corresponding sequences; the column names contain the sequences themselves. They can be accessed using `dimnames()`\n"))
  }
  # Remove temporary file of column ids
  file.remove("col_ids.txt")
  return(out)
}


# Adjacency Matrix: Euclidean Distance on Atchley Factor Embedding
# (Only applicable to TCR CDR3 Amino Acid Sequences)
# This function is intended for building the network for a single cluster, where
# the adjacency matrix is typically dense
adjacencyMatAtchleyFromClones <- function(
  clones, # List of TCR CDR3 amino acid sequences corresponding to the clones
  contig_ids = seq_along(clones), # used by BriseisEncoder to perform the Atchley-factor embedding of the TCR sequences
  max_dist = 1.5, # Maximum Euclidean distance threshold for edge/adjacency between two sequences
  return_type = "adjacency_matrix", # can be set to "distance_matrix" to return the distance matrix instead
  outfile_distance_matrix = NULL # savefile for Euclidean distance matrix
) {
  # Embed amino acid seqs in Euclidean 30-space by Atchley factor representation
  embedded_values <- embedClonesByAtchleyFactor(clones, contig_ids)

  # Compute Euclidean distance matrix on embedded sequence values
  distance_matrix <- as.matrix(stats::dist(embedded_values[ , -1]))

  if (!is.null(outfile_distance_matrix)) {
    # Save distance matrix to file
    utils::write.csv(distance_matrix, outfile_distance_matrix)
  }
  if (return_type == "distance_matrix") {
    return(distance_matrix)
  } else {
    # Convert distance matrix to adjacency matrix using specified bound
    adjacency_matrix <-
      matrix(1, nrow = nrow(distance_matrix), ncol = ncol(distance_matrix))
    adjacency_matrix[distance_matrix > max_dist] <- 0
    return(adjacency_matrix)
  }
}



# Atchley Factor Embedding ------------------------------------------------

# Embed TCR CDR3 amino acid sequences in Euclidean 30-space based on the Atchley
# factor representations of their elements, using a trained encoding model
embedClonesByAtchleyFactor <- function(
  clones, # List of TCR CDR3 amino acid sequences corresponding to the clones
  contig_ids = seq_along(clones) # used by BriseisEncoder
) {
  .checkPythonModules()
  if (length(clones) != length(contig_ids)) {
    stop("length of `clones` and `contig_ids` must match")
  }
  # Write sequences and contig_ids to temporary file
  tempfile_clones <- file.path(getwd(),
                               "Atchley_factor_tcr_only_tmp.csv")
  utils::write.csv(data.frame("contig_id" = contig_ids, "cdr3" = clones),
                   file = tempfile_clones, row.names = FALSE)

  # Write files for trained encoder and Atchley factor table to working dir
  file_trained_encoder <-
    system.file(file.path("python", "TrainedEncoder.h5"),
                package = "RepSeqNetworkAnalysis")
  file_atchley_table <-
    system.file(file.path("python", "Atchley_factors.csv"),
                package = "RepSeqNetworkAnalysis")
  utils::write.csv(data.frame("sysfiles" = c(file_atchley_table,
                                             file_trained_encoder)),
                   file.path(getwd(), "temp_sysfiles.csv"),
                   row.names = FALSE)

  # Run BriseisEncoder python function
  reticulate::py_run_file(
    system.file(file.path("python", "BriseisEncoder_modified.py"),
                package = "RepSeqNetworkAnalysis"))

  # Import embedded values and remove temp files
  embedded_values <- utils::read.csv("temp_atchley_factors_encoded.csv")
  file.remove("Atchley_factor_tcr_only_tmp.csv", "temp_sysfiles.csv",
              "temp_atchley_factors_encoded.csv")
  return(embedded_values)
}
