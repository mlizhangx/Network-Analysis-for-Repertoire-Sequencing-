
# Adjacency Matrices ------------------------------------------------------
sparseAdjacencyMatFromSeqs <- function(
    seqs,
    dist_type = "hamming",
    max_dist = 1,
    drop_isolated_nodes = TRUE
) {
  .checkargs.sparseAdjacencyMatFromSeqs(
    seqs, dist_type, max_dist, drop_isolated_nodes
  )
  seqs <- as.vector(seqs, mode = "character")
  if (dist_type %in% c("levenshtein", "Levenshtein", "lev", "Lev", "l", "L")) {
    cat(paste0(
      "Computing network edges based on a max ", dist_type, " distance of ",
      max_dist, "..."
    ))
    out <- .levAdjacencyMatSparse(
      seqs, max_dist, drop_isolated_nodes, tempdir()
    )
  } else if (dist_type %in% c("hamming", "Hamming", "ham", "Ham", "h", "H")) {
    cat(paste0(
      "Computing network edges based on a max ", dist_type, " distance of ",
      max_dist, "..."
    ))
    out <- .hamAdjacencyMatSparse(
      seqs, max_dist, drop_isolated_nodes, tempdir()
    )
  } else {
    stop('invalid option for `dist_type`')
  }
  cat(" Done.\n")
  num_nodes <- dim(out)[[1]]  # with positive network degree
  if (num_nodes == 0) {
    warning("No edges exist using the specified distance cutoff")
  } else {
    if (drop_isolated_nodes) {
      cat(paste("Network contains",
                num_nodes, "nodes (after removing isolated nodes).\n"
      ))
      clone_ids <- utils::read.table(file.path(tempdir(), "col_ids.txt"))
      dimnames(out)[[1]] <- clone_ids$V1
      dimnames(out)[[2]] <- seqs[clone_ids$V1]
    } else {
      cat(paste("Network contains", num_nodes, "nodes.\n"))
    }
  }
  if (file.exists(file.path(tempdir(), "col_ids.txt"))) {
    file.remove(file.path(tempdir(), "col_ids.txt"))
  }
  out
}

# Network Building --------------------------------------------------------
generateNetworkFromAdjacencyMat <- function(adjacency_matrix) {
  set.seed(9999)
  net <- igraph::graph_from_adjacency_matrix(adjacency_matrix, weighted = TRUE)
  net <- igraph::as.undirected(
    igraph::simplify(net, remove.multiple = T, remove.loops = T)
  )
  return(net)
}

.generateNetworkFromSeqs <- function(
    seqs,
    dist_type = "hamming",
    dist_cutoff = 1,
    drop_isolated_nodes = TRUE,
    contig_ids = seq_along(seqs),
    outfile_adjacency_matrix = NULL,
    outfile_distance_matrix = NULL,
    return_type = "network"
) {
  dist_keys <- c("levenshtein", "Levenshtein", "lev", "Lev", "l", "L",
                 "hamming", "Hamming", "ham", "Ham", "h", "H"
  )
  if (dist_type %in% dist_keys) {
    adjacency_matrix <-
      sparseAdjacencyMatFromSeqs(seqs = seqs,
                                 dist_type = dist_type,
                                 max_dist = dist_cutoff,
                                 drop_isolated_nodes = drop_isolated_nodes
      )
    if (!is.null(outfile_adjacency_matrix)) {
      Matrix::writeMM(adjacency_matrix, outfile_adjacency_matrix)
    }
    if (sum(dim(adjacency_matrix)) == 0) {
      return(NULL)   # no nodes connected (empty matrix)
    }
  } else { stop("invalid option for argument `dist_type`") }

  if (return_type == "adjacency_matrix") {
    return(adjacency_matrix)
  }
  generateNetworkFromAdjacencyMat(adjacency_matrix)
}

.generateSingleChainNetwork <- function(
    data, seq_col, dist_type, dist_cutoff, drop_isolated_nodes
) {
  adjacency_matrix <- .generateNetworkFromSeqs(
    data[[seq_col]], dist_type, dist_cutoff, contig_ids = rownames(data),
    return_type = "adjacency_matrix", drop_isolated_nodes = drop_isolated_nodes
  )
  if (is.null(adjacency_matrix)) { return(NULL) }
  net <- generateNetworkFromAdjacencyMat(adjacency_matrix)
  if (drop_isolated_nodes) {
    data <- .subsetDataForAdjacencyMatrix(data, adjacency_matrix)
  }
  list("igraph" = net, "adjacency_matrix" = adjacency_matrix,
       "node_data" = as.data.frame(data)
  )
}

.generateDualChainNetwork <- function(
    data, a_col, b_col, dist_type, dist_cutoff, drop_isolated_nodes
) {
  cat("Computing graph adjacency based on sequences in first chain:\n")
  adj_mat_a <- sparseAdjacencyMatFromSeqs(
    data[[a_col]], dist_type, dist_cutoff, drop_isolated_nodes = FALSE
  )
  if (sum(dim(adj_mat_a)) == 0) {
    warning(paste("No edges exist in the network for the first chain",
                  "using the specified distance type and cutoff"
    ))
    return(NULL)
  }
  cat("Computing graph adjacency based on sequences in second chain:\n")
  adj_mat_b <- sparseAdjacencyMatFromSeqs(
    data[[b_col]], dist_type, dist_cutoff, drop_isolated_nodes = FALSE
  )
  if (sum(dim(adj_mat_b)) == 0) {
    warning(paste("No edges exist in the network for the second chain",
                  "using the specified distance type and cutoff"
    ))
    return(NULL)
  }
  cat("Intersecting the adjacencies from both chains...")
  adjacency_matrix <- adj_mat_a + adj_mat_b
  adjacency_matrix[adjacency_matrix == 1] <- 0
  adjacency_matrix[adjacency_matrix == 2] <- 1
  cat(" Done.\n")
  cat("Building network based on the combined adjacencies... ")
  net <- generateNetworkFromAdjacencyMat(adjacency_matrix)
  cat(" Done.\n")
  if (drop_isolated_nodes) {
    cat("Dropping isolated nodes...")
    nodes_to_keep <- igraph::degree(net) > 0
    if (sum(nodes_to_keep) == 0) {
      warning(paste("No edges exist in the combined network for both chains",
                    "using the specified distance type and cutoff"
      ))
      return(NULL)
    }
    adjacency_matrix <- adjacency_matrix[nodes_to_keep, nodes_to_keep]
    data <- data[nodes_to_keep, , drop = FALSE]
    net <- generateNetworkFromAdjacencyMat(adjacency_matrix)
    cat(" Done.\n")
  }
  cat(paste("Network contains", nrow(data), "nodes.\n"))
  list("igraph" = net,
       "adjacency_matrix" = adjacency_matrix,
       "adj_mat_a" = adj_mat_a,
       "adj_mat_b" = adj_mat_b,
       "node_data" = as.data.frame(data)
  )
}

generateNetworkObjects <- function(
    data,
    seq_col,
    dist_type = "hamming",
    dist_cutoff = 1,
    drop_isolated_nodes = TRUE
) {
  data <- as.data.frame(data)
  .checkargs.generateNetworkObjects(
    data, seq_col, dist_type, dist_cutoff, drop_isolated_nodes
  )
  if (length(seq_col) == 1) {
    return(.generateSingleChainNetwork(
      data, seq_col,
      dist_type, dist_cutoff, drop_isolated_nodes
    ))
  } else if (length(seq_col) == 2) {
    return(.generateDualChainNetwork(
      data, seq_col[[1]], seq_col[[2]],
      dist_type, dist_cutoff, drop_isolated_nodes
    ))
  }
}


# Network Properties ------------------------------------------------------

addNodeNetworkStats <- function(
    data,
    net,
    stats_to_include = chooseNodeStats(),
    cluster_fun = cluster_fast_greedy
) {

  .checkargs.addNodeNetworkStats(
    data, net, stats_to_include, cluster_fun
  )
  if (!typeof(stats_to_include) %in% c("list", "logical"))  {
    if (stats_to_include == "all") {
      stats_to_include <- chooseNodeStats(all_stats = TRUE)
    } else if (stats_to_include == "cluster_id_only") {
      stats_to_include <- exclusiveNodeStats(cluster_id = TRUE)
    }
  }
  if (stats_to_include[["degree"]] || stats_to_include[["all_stats"]]) {
    data$degree <- igraph::degree(net)
  }
  if (stats_to_include[["cluster_id"]] || stats_to_include[["all_stats"]]) {
    cat("Computing cluster membership within the network...")
    data$cluster_id <- as.factor(as.integer(cluster_fun(net)$membership))
    cat(" Done.\n")
  }
  cat(paste0("Computing node-level network statistics..."))
  if (stats_to_include[["transitivity"]] || stats_to_include[["all_stats"]]) {
    data$transitivity <- igraph::transitivity(net, type = "local")
  }
  if (stats_to_include[["closeness"]] || stats_to_include[["all_stats"]]) {
    data$closeness <- igraph::closeness(net, mode = "all", weights = NA)
  }
  if (
    stats_to_include[["centrality_by_closeness"]] ||
    stats_to_include[["all_stats"]]
  ) {
    data$centrality_by_closeness <-
      igraph::centr_clo(net, mode = "all", normalized = T)$res
  }
  if (
    stats_to_include[["eigen_centrality"]] || stats_to_include[["all_stats"]]
  ) {
    data$eigen_centrality <-
      igraph::eigen_centrality(net, directed = T, weights = NA)$vector
  }
  if (
    stats_to_include[["centrality_by_eigen"]] || stats_to_include[["all_stats"]]
  ) {
    data$centrality_by_eigen <-
      igraph::centr_eigen(net, directed = T, normalized = T)$vector
  }
  if (stats_to_include[["betweenness"]] || stats_to_include[["all_stats"]]) {
    data$betweenness <- igraph::betweenness(net, directed = T, weights = NA)
  }
  if (
    stats_to_include[["centrality_by_betweenness"]] ||
    stats_to_include[["all_stats"]]
  ) {
    data$centrality_by_betweenness <-
      igraph::centr_betw(net, directed = T, normalized = T)$res
  }
  if (
    stats_to_include[["authority_score"]] || stats_to_include[["all_stats"]]
  ) {
    data$authority_score <- igraph::authority_score(net, weights = NA)$vector
  }
  if (stats_to_include[["coreness"]] || stats_to_include[["all_stats"]]) {
    data$coreness <- igraph::coreness(net, mode = "all")
  }
  if (stats_to_include[["page_rank"]] || stats_to_include[["all_stats"]]) {
    data$page_rank <- igraph::page_rank(net)$vector
  }
  cat(" Done.\n")
  data
}

chooseNodeStats <- function(
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
  .checkargs.chooseNodeStats(
    degree, cluster_id, transitivity, closeness, centrality_by_closeness,
    eigen_centrality, centrality_by_eigen, betweenness,
    centrality_by_betweenness, authority_score, coreness, page_rank, all_stats
  )
  c("degree" = degree,
    "cluster_id" = cluster_id,
    "transitivity" = transitivity,
    "closeness" = closeness,
    "centrality_by_closeness" = centrality_by_closeness,
    "eigen_centrality" = eigen_centrality,
    "centrality_by_eigen" = centrality_by_eigen,
    "betweenness" = betweenness,
    "centrality_by_betweenness" = centrality_by_betweenness,
    "authority_score" = authority_score,
    "coreness" = coreness,
    "page_rank" = page_rank,
    "all_stats" = all_stats
  )
}

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
  lifecycle::deprecate_warn(
    when = "0.0.9035",
    what = "node_stat_settings()",
    with = "chooseNodeStats()",
    details =
      "all node-level network properties are now automatically computed"
  )
  chooseNodeStats(
    degree, cluster_id, transitivity, closeness, centrality_by_closeness,
    eigen_centrality, centrality_by_eigen, betweenness,
    centrality_by_betweenness, authority_score, coreness, page_rank, all_stats
  )
}

exclusiveNodeStats <- function(
    degree = FALSE,
    cluster_id = FALSE,
    transitivity = FALSE,
    closeness = FALSE,
    centrality_by_closeness = FALSE,
    eigen_centrality = FALSE,
    centrality_by_eigen = FALSE,
    betweenness = FALSE,
    centrality_by_betweenness = FALSE,
    authority_score = FALSE,
    coreness = FALSE,
    page_rank = FALSE
) {
  c("degree" = degree,
    "cluster_id" = cluster_id,
    "transitivity" = transitivity,
    "closeness" = closeness,
    "centrality_by_closeness" = centrality_by_closeness,
    "eigen_centrality" = eigen_centrality,
    "centrality_by_eigen" = centrality_by_eigen,
    "betweenness" = betweenness,
    "centrality_by_betweenness" = centrality_by_betweenness,
    "authority_score" = authority_score,
    "coreness" = coreness,
    "page_rank" = page_rank,
    "all_stats" = FALSE
  )
}

addClusterMembership <- function(
    data,
    net,
    fun = cluster_fast_greedy
) {
  .checkargs.addClusterMembership(data, net, fun)
  fun <- match.fun(fun)
  cat("Computing cluster membership within the network...")
  data$cluster_id <- as.factor(as.integer(fun(net)$membership))
  cat(" Done.\n")
  data
}

getClusterStats <- function(
    data,
    adjacency_matrix,
    seq_col = NULL,
    count_col = NULL,
    cluster_id_col = NULL,
    degree_col = NULL,
    cluster_fun = cluster_fast_greedy
) {

  .checkargs.getClusterStats(
    data, adjacency_matrix, seq_col, count_col, cluster_id_col,
    degree_col, cluster_fun
  )
  if (is.null(cluster_id_col) || is.null(degree_col)) {
    net <- generateNetworkFromAdjacencyMat(adjacency_matrix)
    if (is.null(cluster_id_col)) {
      cluster_id_col <- "cluster_id"
      data <- addClusterMembership(data, net, cluster_fun)
    }
    if (is.null(degree_col)) {
      degree_col <- "deg"
      data$deg <- igraph::degree(net)
    }
  }
  if (!is.null(seq_col) && length(seq_col) > 1) {
    return(.computeClusterStatsDualChain(
      data, adjacency_matrix, seq_col, count_col, cluster_id_col, degree_col
    ))
  }
  .computeClusterStats(
    data, adjacency_matrix, seq_col, count_col, cluster_id_col, degree_col
  )
}

.computeClusterStats <- function(
    data, adjacency_matrix,
    seq_col, count_col, cluster_id_col, degree_col
) {
  if (!is.null(seq_col)) {
    seq_lengths <- nchar(data[[seq_col]])
  }
  out <- as.data.frame(table(data[[cluster_id_col]]))
  colnames(out) <- c("cluster_id", "node_count")
  num_clusters <- nrow(out)
  cat(paste0(
    "Computing statistics for the ", num_clusters, " clusters in the network..."
  ))
  out$eigen_centrality_eigenvalue <- out$eigen_centrality_index <-
    out$closeness_centrality_index <- out$degree_centrality_index <-
    out$edge_density <- out$assortativity <- out$global_transitivity <-
    out$diameter_length <- out$seq_w_max_count <- out$max_count <-
    out$agg_count <- out$seq_w_max_degree <- out$max_degree <-
    out$mean_degree <- out$mean_seq_length <- NA
  for (i in 1:num_clusters) {
    cluster_row <- which(out$cluster_id == i)
    node_ids <- data$cluster_id == i
    if (!is.null(seq_col)) {
      out$mean_seq_length[[cluster_row]] <- round(
        mean(seq_lengths[node_ids]), 2
      )
    }
    out$mean_degree[[cluster_row]] <- round(mean(data[node_ids, degree_col]), 2)
    max_deg <- max(data[node_ids, degree_col])
    out$max_degree[[cluster_row]] <- max_deg
    if (!is.null(seq_col)) {
      node_id_max_deg <- which(node_ids & data[[degree_col]] == max_deg)[[1]]
      out$seq_w_max_degree[[cluster_row]] <- as.character(
        data[node_id_max_deg, seq_col]
      )
    }
    if (!is.null(count_col)) {
      out$agg_count[[cluster_row]] <- sum(data[node_ids, count_col])
      max_count <- max(data[node_ids, count_col])
      out$max_count[[cluster_row]] <- max_count
      if (!is.null(seq_col)) {
        node_id_max_count <- which(
          node_ids & data[[count_col]] == max_count
        )[[1]]
        out$seq_w_max_count[[cluster_row]] <- as.character(
          data[node_id_max_count, seq_col]
        )
      }
    }
    cluster_adjacency_matrix <- as.matrix(adjacency_matrix[node_ids, node_ids])
    cluster <- generateNetworkFromAdjacencyMat(cluster_adjacency_matrix)
    out$diameter_length[[cluster_row]] <-
      length(igraph::get_diameter(cluster, directed = T))
    out$assortativity[[cluster_row]] <-
      igraph::assortativity_degree(cluster, directed = F)
    out$global_transitivity[[cluster_row]] <-
      igraph::transitivity(cluster, type = "global")
    out$edge_density[[cluster_row]] <-
      igraph::edge_density(cluster, loops = F)
    out$degree_centrality_index[[cluster_row]] <-
      igraph::centr_degree(cluster, mode = "in", normalized = T)$centralization
    out$closeness_centrality_index[[cluster_row]] <-
      igraph::centr_clo(cluster, mode = "all", normalized = T)$centralization
    out$eigen_centrality_index[[cluster_row]] <-
      igraph::centr_eigen(cluster, directed = T, normalized = T)$centralization
    out$eigen_centrality_eigenvalue[[cluster_row]] <-
      igraph::eigen_centrality(cluster, directed = T, weights = NA)$value
  }
  cat(" Done.\n")
  out
}

.computeClusterStatsDualChain <- function(
    data, adjacency_matrix, seq_col, count_col,
    cluster_id_col, degree_col
) {
  A_seq_lengths <- nchar(data[[seq_col[[1]]]])
  B_seq_lengths <- nchar(data[[seq_col[[2]]]])
  out <- as.data.frame(table(data[[cluster_id_col]]))
  colnames(out) <- c("cluster_id", "node_count")
  num_clusters <- nrow(out)
  cat(paste0(
    "Computing statistics for the ", num_clusters, " clusters in the network..."
  ))
  out$eigen_centrality_eigenvalue <- out$eigen_centrality_index <-
    out$closeness_centrality_index <- out$degree_centrality_index <-
    out$edge_density <- out$assortativity <- out$global_transitivity <-
    out$diameter_length <- out$B_seq_w_max_count <- out$A_seq_w_max_count <-
    out$max_count <- out$agg_count <- out$B_seq_w_max_degree <-
    out$A_seq_w_max_degree <- out$max_degree <- out$mean_degree <-
    out$mean_B_seq_length <- out$mean_A_seq_length <- NA
  for (i in 1:num_clusters) {
    cluster_row <- which(out$cluster_id == i)
    node_ids <- data$cluster_id == i
    out$mean_A_seq_length[[cluster_row]] <- round(
      mean(A_seq_lengths[node_ids]), 2
    )
    out$mean_B_seq_length[[cluster_row]] <- round(
      mean(B_seq_lengths[node_ids]), 2
    )
    out$mean_degree[[cluster_row]] <- round(mean(data[node_ids, degree_col]), 2)
    max_deg <- max(data[node_ids, degree_col])
    out$max_degree[[cluster_row]] <- max_deg
    node_id_max_deg <- which(node_ids & data[[degree_col]] == max_deg)[[1]]
    out$A_seq_w_max_degree[[cluster_row]] <-
      as.character(data[node_id_max_deg, seq_col[[1]]])
    out$B_seq_w_max_degree[[cluster_row]] <-
      as.character(data[node_id_max_deg, seq_col[[2]]])
    if (!is.null(count_col)) {
      out$agg_count[[cluster_row]] <- sum(data[node_ids, count_col])
      max_count <- max(data[node_ids, count_col])
      out$max_count[[cluster_row]] <- max_count
      node_id_max_count <- which(node_ids & data[[count_col]] == max_count)[[1]]
      out$A_seq_w_max_count[[cluster_row]] <-
        as.character(data[node_id_max_count, seq_col[[1]]])
      out$B_seq_w_max_count[[cluster_row]] <-
        as.character(data[node_id_max_count, seq_col[[2]]])
    }
    cluster_adjacency_matrix <- as.matrix(adjacency_matrix[node_ids, node_ids])
    cluster <- generateNetworkFromAdjacencyMat(cluster_adjacency_matrix)
    out$diameter_length[[cluster_row]] <-
      length(igraph::get_diameter(cluster, directed = T))
    out$assortativity[[cluster_row]] <-
      igraph::assortativity_degree(cluster, directed = F)
    out$global_transitivity[[cluster_row]] <-
      igraph::transitivity(cluster, type = "global")
    out$edge_density[[cluster_row]] <-
      igraph::edge_density(cluster, loops = F)
    out$degree_centrality_index[[cluster_row]] <-
      igraph::centr_degree(cluster, mode = "in", normalized = T)$centralization
    out$closeness_centrality_index[[cluster_row]] <-
      igraph::centr_clo(cluster, mode = "all", normalized = T)$centralization
    out$eigen_centrality_index[[cluster_row]] <-
      igraph::centr_eigen(cluster, directed = T, normalized = T)$centralization
    out$eigen_centrality_eigenvalue[[cluster_row]] <-
      igraph::eigen_centrality(cluster, directed = T, weights = NA)$value
  }
  cat(" Done.\n")
  out
}


.getClusterStats2 <- function(
    data, adjacency_matrix, seq_col, count_col,
    cluster_fun = cluster_fast_greedy
) {
  # check to avoid redundant cluster analysis after computing node level stats
  cluster_id_col <- degree_col <- NULL
  if ("cluster_id" %in% names(data)) { cluster_id_col <- "cluster_id" }
  if ("degree" %in% names(data)) { degree_col <- "degree" }
  getClusterStats(
    data, adjacency_matrix, seq_col,
    count_col, cluster_id_col, degree_col, cluster_fun
  )
}
