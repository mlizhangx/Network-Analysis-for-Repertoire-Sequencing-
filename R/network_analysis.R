# NAIR: Network Analysis of Immune Repertoire
# Copyright (C) 2023 Li Zhang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


# Adjacency Matrices ------------------------------------------------------
generateAdjacencyMatrix <- function(
    seqs,
    dist_type = "hamming",
    dist_cutoff = 1,
    drop_isolated_nodes = TRUE,
    method = "default",
    verbose = FALSE
) {
  .MUST.isValidSeqVector(seqs)
  seqs <- as.vector(seqs, mode = "character")
  dist_type <- .checkDistType(dist_type)
  dist_cutoff <- .check(dist_cutoff, .isNonneg, 1)
  drop_isolated_nodes <- .checkTF(drop_isolated_nodes, TRUE)
  method <- .checkMethod(method, dist_cutoff)
  msg <- .makemsg(verbose)
  tmpfile <- tempfile(pattern = "col_ids", fileext = ".txt")
  if (dist_type == "levenshtein") {
    msg("Computing network edges based on a max ", dist_type, " distance of ",
        dist_cutoff, "...", newline = FALSE
    )
    if (method == "pattern") {
      out <- RPatternJoin::similarityJoin(
        seqs, dist_cutoff, "Levenshtein", "partition_pattern",
        drop_deg_one = drop_isolated_nodes,
        special_chars = FALSE,
        output_format = "adj_matrix"
      )
    } else if (method == "sort") {
      out <- .sortAdjacencyMatSparse(seqs, dist_cutoff, "L",
                                     drop_isolated_nodes, tmpfile
      )
    }  else {
      out <- .levAdjacencyMatSparse(seqs, dist_cutoff, drop_isolated_nodes,
                                    tmpfile
      )
    }
  } else if (dist_type == "hamming") {
    dist_type <- "hamming"
    msg("Computing network edges based on a max ", dist_type, " distance of ",
        dist_cutoff, "...", newline = FALSE
    )
    if (method == "pattern") {
      out <- RPatternJoin::similarityJoin(
        seqs, dist_cutoff, "Hamming", "partition_pattern",
        drop_deg_one = drop_isolated_nodes,
        special_chars = FALSE,
        output_format = "adj_matrix"
      )
    } else if (method == "sort") {
      out <- .sortAdjacencyMatSparse(seqs, dist_cutoff, "H",
                                     drop_isolated_nodes, tmpfile
      )
    } else {
      out <- .hamAdjacencyMatSparse(seqs, dist_cutoff, drop_isolated_nodes,
                                    tmpfile
      )
    }
  } else {
    stop("invalid option for ", sQuote("dist_type"))
  }
  msg(" Done.")
  num_nodes <- dim(out)[[1]]  # with positive network degree
  if (num_nodes == 0) {
    warning("No edges exist using the specified distance cutoff")
  } else {
    if (drop_isolated_nodes && method != "pattern") {
      msg("Network contains ", num_nodes,
          " nodes (after removing isolated nodes)."
      )
      clone_ids <- scan(tmpfile, quiet = TRUE)
      dimnames(out)[[1]] <- clone_ids
      dimnames(out)[[2]] <- seqs[clone_ids]
    } else {
      msg("Network contains ", num_nodes, " nodes.")
    }
  }
  unlink(tmpfile)
  out
}

sparseAdjacencyMatFromSeqs <- function(
    seqs,
    dist_type = "hamming",
    dist_cutoff = 1,
    drop_isolated_nodes = TRUE,
    method = "default",
    verbose = FALSE,
    max_dist = deprecated()
) {
  lifecycle::deprecate_warn(
    when = "1.0.1",
    what = "sparseAdjacencyMatFromSeqs()",
    with = "generateAdjacencyMatrix()"
  )
  if (lifecycle::is_present(max_dist)) { dist_cutoff <- max_dist }
  generateAdjacencyMatrix(
    seqs, dist_type, dist_cutoff, drop_isolated_nodes, method, verbose
  )
}
# Network Building --------------------------------------------------------

generateNetworkObjects <- function(
    data,
    seq_col,
    dist_type = "hamming",
    dist_cutoff = 1,
    drop_isolated_nodes = TRUE,
    method = "default",
    verbose = FALSE
) {
  data_name <- deparse(substitute(data))
  data <- as.data.frame(data)
  .MUST.isDataFrame(data, data_name)
  .MUST.hasMultipleRows(data, data_name)
  .MUST.isSeqColrefs(seq_col, data, deparse(substitute(seq_col)), data_name)
  dist_type <- .checkDistType(dist_type, "hamming")
  dist_cutoff <- .check(dist_cutoff, .isNonneg, 1)
  drop_isolated_nodes <- .checkTF(drop_isolated_nodes, TRUE)
  method <- .checkMethod(method, dist_cutoff)
  msg <- .makemsg(verbose)
  if (length(seq_col) == 1) {
    net <- .generateSingleChainNetwork(
      data, seq_col,
      dist_type, dist_cutoff, drop_isolated_nodes, method, verbose
    )
  } else if (length(seq_col) == 2) {
    net <- .generateDualChainNetwork(
      data, seq_col[[1]], seq_col[[2]],
      dist_type, dist_cutoff, drop_isolated_nodes, method, verbose
    )
  }
  msg("Network objects and node metadata stored in a list")
  net
}


generateNetworkGraph <- function(adjacency_matrix) {
  .MUST.isAdjacencyMatrix(adjacency_matrix)
  igraph <- igraph::graph_from_adjacency_matrix(adjacency_matrix,
                                                weighted = TRUE
  )
  igraph <- igraph::as.undirected(
    igraph::simplify(igraph, remove.multiple = T, remove.loops = T)
  )
  igraph
}

generateNetworkFromAdjacencyMat <- function(adjacency_matrix) {
  lifecycle::deprecate_warn(
    when = "1.0.1",
    what = "generateNetworkFromAdjacencyMat()",
    with = "generateNetworkGraph()"
  )
  generateNetworkGraph(adjacency_matrix)
}

.generateNetworkGraphFromSeqs <- function(
    seqs,
    dist_type = "hamming",
    dist_cutoff = 1,
    drop_isolated_nodes = TRUE,
    method = "default",
    contig_ids = seq_along(seqs),
    outfile_adjacency_matrix = NULL,
    outfile_distance_matrix = NULL,
    return_type = "igraph",
    verbose = FALSE
) {
  adjacency_matrix <-
    generateAdjacencyMatrix(seqs = seqs,
                            dist_type = dist_type,
                            dist_cutoff = dist_cutoff,
                            drop_isolated_nodes = drop_isolated_nodes,
                            method = method,
                            verbose = verbose
    )
  if (!is.null(outfile_adjacency_matrix)) {
    Matrix::writeMM(adjacency_matrix, outfile_adjacency_matrix)
  }
  if (sum(dim(adjacency_matrix)) == 0) {
    return(NULL)   # no nodes connected (empty matrix)
  }
  if (return_type == "adjacency_matrix") {
    return(adjacency_matrix)
  }
  generateNetworkGraph(adjacency_matrix)
}


.generateSingleChainNetwork <- function(
    data, seq_col, dist_type, dist_cutoff, drop_isolated_nodes,
    method = "default", verbose = FALSE
) {
  adjacency_matrix <- .generateNetworkGraphFromSeqs(
    data[[seq_col]], dist_type, dist_cutoff, contig_ids = rownames(data),
    return_type = "adjacency_matrix", drop_isolated_nodes = drop_isolated_nodes,
    method = method, verbose = verbose
  )
  if (is.null(adjacency_matrix)) {
    return(NULL)
  }
  igraph <- generateNetworkGraph(adjacency_matrix)
  if (drop_isolated_nodes) {
    data <- .subsetDataForAdjacencyMatrix(data, adjacency_matrix)
  }
  details <- list(
    seq_col = seq_col,
    dist_type = dist_type,
    dist_cutoff = dist_cutoff,
    drop_isolated_nodes = drop_isolated_nodes,
    nodes_in_network = nrow(data)
  )
  list("details" = details,
       "igraph" = igraph,
       "adjacency_matrix" = adjacency_matrix,
       "node_data" = as.data.frame(data)
  )
}

.generateDualChainNetwork <- function(
    data, a_col, b_col, dist_type, dist_cutoff, drop_isolated_nodes,
    method = "default", verbose = FALSE
) {
  msg <- .makemsg(verbose)
  msg("Computing graph adjacency based on sequences in first chain...")
  adj_mat_a <- generateAdjacencyMatrix(
    data[[a_col]], dist_type, dist_cutoff, drop_isolated_nodes = FALSE,
    method = method, verbose = verbose
  )
  if (sum(dim(adj_mat_a)) == 0) {
    warning("No edges exist in the network for the first chain ",
            "using the specified distance type and cutoff"
    )
    return(NULL)
  }
  msg("Computing graph adjacency based on sequences in second chain:")
  adj_mat_b <- generateAdjacencyMatrix(
    data[[b_col]], dist_type, dist_cutoff, drop_isolated_nodes = FALSE,
    method = method, verbose = verbose
  )
  if (sum(dim(adj_mat_b)) == 0) {
    warning("No edges exist in the network for the second chain ",
            "using the specified distance type and cutoff"
    )
    return(NULL)
  }
  msg("Intersecting the adjacencies from both chains...", newline = FALSE)
  adjacency_matrix <- adj_mat_a + adj_mat_b
  adjacency_matrix[adjacency_matrix == 1] <- 0
  adjacency_matrix[adjacency_matrix == 2] <- 1
  msg(" Done.")
  msg("Building network based on the combined adjacencies...", newline = FALSE)
  igraph <- generateNetworkGraph(adjacency_matrix)
  msg(" Done.")
  if (drop_isolated_nodes) {
    msg("Dropping isolated nodes...", newline = FALSE)
    nodes_to_keep <- igraph::degree(igraph) > 0
    if (sum(nodes_to_keep) == 0) {
      warning("No edges exist in the combined network for both chains ",
              "using the specified distance type and cutoff"
      )
      return(NULL)
    }
    adjacency_matrix <- adjacency_matrix[nodes_to_keep, nodes_to_keep]
    data <- data[nodes_to_keep, , drop = FALSE]
    igraph <- generateNetworkGraph(adjacency_matrix)
    msg(" Done.")
  }
  msg("Network contains ", nrow(data), " nodes.")
  details <- list(
    seq_col = c("a_col" = a_col, "b_col" = b_col),
    dist_type = dist_type,
    dist_cutoff = dist_cutoff,
    drop_isolated_nodes = drop_isolated_nodes,
    nodes_in_network = nrow(data)
  )
  list("details" = details,
       "igraph" = igraph,
       "adjacency_matrix" = adjacency_matrix,
       "adj_mat_a" = adj_mat_a,
       "adj_mat_b" = adj_mat_b,
       "node_data" = as.data.frame(data)
  )
}


# Network Properties ------------------------------------------------------

addNodeStats <- function(
    net,
    stats_to_include = chooseNodeStats(),
    cluster_fun = "fast_greedy",
    cluster_id_name = "cluster_id",
    overwrite = FALSE,
    verbose = FALSE,
    ...
) {
  .MUST.isBaseNetworkOutput(net)
  stats_to_include <- .checkStatsToInclude(stats_to_include)
  if (stats_to_include[["degree"]]) {
    net$node_data$degree <- igraph::degree(net$igraph)
  }
  if (stats_to_include[["cluster_id"]]) {
    net <- addClusterMembership(
      net = net,
      cluster_fun = cluster_fun, cluster_id_name = cluster_id_name,
      overwrite = overwrite, verbose = verbose, ... = ...
    )
  }
  cluster_fun <- .check(cluster_fun, .isClusterFun, "fast_greedy")
  cluster_id_name <- make.names(
    .check(cluster_id_name, .isString, "cluster_id")
  )
  overwrite <- .checkTF(overwrite, FALSE)

  msg <- .makemsg(verbose)
  msg("Computing node-level network properties...", newline = FALSE)
  if (stats_to_include[["transitivity"]]) {
    net$node_data$transitivity <-
      igraph::transitivity(net$igraph, type = "local")
  }
  if (stats_to_include[["closeness"]]) {
    net$node_data$closeness <-
      igraph::closeness(net$igraph, mode = "all", weights = NA)
  }
  if (stats_to_include[["centrality_by_closeness"]]) {
    net$node_data$centrality_by_closeness <-
      igraph::centr_clo(net$igraph, mode = "all", normalized = T)$res
  }
  if (stats_to_include[["eigen_centrality"]]) {
    net$node_data$eigen_centrality <-
      igraph::eigen_centrality(net$igraph, directed = T, weights = NA)$vector
  }
  if (stats_to_include[["centrality_by_eigen"]]) {
    net$node_data$centrality_by_eigen <-
      igraph::centr_eigen(net$igraph, directed = T, normalized = T)$vector
  }
  if (stats_to_include[["betweenness"]]) {
    net$node_data$betweenness <-
      igraph::betweenness(net$igraph, directed = T, weights = NA)
  }
  if (stats_to_include[["centrality_by_betweenness"]]) {
    net$node_data$centrality_by_betweenness <-
      igraph::centr_betw(net$igraph, directed = T, normalized = T)$res
  }
  if (stats_to_include[["authority_score"]]) {
    net$node_data$authority_score <-
      igraph::authority_score(net$igraph, weights = NA)$vector
  }
  if (stats_to_include[["coreness"]]) {
    net$node_data$coreness <- igraph::coreness(net$igraph, mode = "all")
  }
  if (stats_to_include[["page_rank"]]) {
    net$node_data$page_rank <- igraph::page_rank(net$igraph)$vector
  }
  msg(" Done.")
  msg("Network properties added to node metadata")
  net
}


addNodeNetworkStats <- function(
    data,
    net,
    stats_to_include = chooseNodeStats(),
    cluster_fun = "fast_greedy",
    cluster_id_name = "cluster_id",
    overwrite = FALSE,
    verbose = FALSE,
    ...
) {
  lifecycle::deprecate_warn(
    when = "1.0.1",
    what = "addNodeNetworkStats()",
    with = "addNodeStats()"
  )
  data_name <- deparse(substitute(data))
  data <- as.data.frame(data)
  .MUST.isDataFrame(data, data_name)
  .MUST.isIgraph(net, deparse(substitute(net)))
  .MUST.doesIgraphMatchData(net, data, deparse(substitute(net)), data_name)
  stats_to_include <- .checkStatsToInclude(stats_to_include)
  if (stats_to_include[["degree"]]) {
    data$degree <- igraph::degree(net)
  }
  if (stats_to_include[["cluster_id"]]) {
    data <- addClusterMembership(
      data = data, net = net,
      cluster_fun = cluster_fun, cluster_id_name = cluster_id_name,
      overwrite = overwrite, verbose = verbose, ... = ...
    )
  }
  cluster_fun <- .check(cluster_fun, .isClusterFun, "fast_greedy")
  cluster_id_name <- make.names(
    .check(cluster_id_name, .isString, "cluster_id")
  )
  overwrite <- .checkTF(overwrite, FALSE)

  msg <- .makemsg(verbose)
  msg("Computing node-level network properties...", newline = FALSE)
  if (stats_to_include[["transitivity"]]) {
    data$transitivity <- igraph::transitivity(net, type = "local")
  }
  if (stats_to_include[["closeness"]]) {
    data$closeness <- igraph::closeness(net, mode = "all", weights = NA)
  }
  if (stats_to_include[["centrality_by_closeness"]]) {
    data$centrality_by_closeness <-
      igraph::centr_clo(net, mode = "all", normalized = T)$res
  }
  if (stats_to_include[["eigen_centrality"]]) {
    data$eigen_centrality <-
      igraph::eigen_centrality(net, directed = T, weights = NA)$vector
  }
  if (stats_to_include[["centrality_by_eigen"]]) {
    data$centrality_by_eigen <-
      igraph::centr_eigen(net, directed = T, normalized = T)$vector
  }
  if (stats_to_include[["betweenness"]]) {
    data$betweenness <- igraph::betweenness(net, directed = T, weights = NA)
  }
  if (stats_to_include[["centrality_by_betweenness"]]) {
    data$centrality_by_betweenness <-
      igraph::centr_betw(net, directed = T, normalized = T)$res
  }
  if (stats_to_include[["authority_score"]]) {
    data$authority_score <- igraph::authority_score(net, weights = NA)$vector
  }
  if (stats_to_include[["coreness"]]) {
    data$coreness <- igraph::coreness(net, mode = "all")
  }
  if (stats_to_include[["page_rank"]]) {
    data$page_rank <- igraph::page_rank(net)$vector
  }
  msg(" Done.")
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
  degree <- .checkTF(degree, TRUE)
  cluster_id <- .checkTF(cluster_id, FALSE)
  transitivity <- .checkTF(transitivity, TRUE)
  closeness <- .checkTF(closeness, FALSE)
  centrality_by_closeness <- .checkTF(centrality_by_closeness, FALSE)
  eigen_centrality <- .checkTF(eigen_centrality, TRUE)
  centrality_by_eigen <- .checkTF(centrality_by_eigen, TRUE)
  betweenness <- .checkTF(betweenness, TRUE)
  centrality_by_betweenness <- .checkTF(centrality_by_betweenness, TRUE)
  authority_score <- .checkTF(authority_score, TRUE)
  coreness <- .checkTF(coreness, TRUE)
  all_stats <- .checkTF(all_stats, FALSE)
  out <-
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
      "page_rank" = page_rank
    )
  if (all_stats) { out[1:length(out)] <- TRUE }
  out
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
  lifecycle::deprecate_soft(
    when = "0.0.9035",
    what = "node_stat_settings()",
    with = "chooseNodeStats()",
    details = "new function name is clearer as to the function's purpose"
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
  degree <- .checkTF(degree, FALSE)
  cluster_id <- .checkTF(cluster_id, FALSE)
  transitivity <- .checkTF(transitivity, FALSE)
  closeness <- .checkTF(closeness, FALSE)
  centrality_by_closeness <- .checkTF(centrality_by_closeness, FALSE)
  eigen_centrality <- .checkTF(eigen_centrality, FALSE)
  centrality_by_eigen <- .checkTF(centrality_by_eigen, FALSE)
  betweenness <- .checkTF(betweenness, FALSE)
  centrality_by_betweenness <- .checkTF(centrality_by_betweenness, FALSE)
  authority_score <- .checkTF(authority_score, FALSE)
  coreness <- .checkTF(coreness, FALSE)
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
    "page_rank" = page_rank
  )
}

addClusterMembership <- function(
    net,
    cluster_fun = "fast_greedy",
    cluster_id_name = "cluster_id",
    overwrite = FALSE,
    verbose = FALSE,
    ...,
    data = deprecated(),
    fun = deprecated()
) {
  if (lifecycle::is_present(fun)) {
    lifecycle::deprecate_soft(
      when = "1.0.1",
      what = "addClusterMembership(fun)",
      with = "addClusterMembership(cluster_fun)",
      details = "Conforms to corresponding argument name in other functions."
    )
    cluster_fun <- fun
  }
  if (lifecycle::is_present(data)) {
    if (isTRUE(inherits(net, "igraph"))) {
      lifecycle::deprecate_soft(
        when = "1.0.1",
        what = "addClusterMembership(data)",
        details = paste(
          "Instead of passing the node metadata to ", sQuote("data"),
          "and the igraph to ", sQuote("net"), ", ",
          "please pass the list containing both objects to ", sQuote("net"),
          "and discontinue use of ", sQuote("data")
        )
      )
      data_name <- deparse(substitute(data))
      data <- as.data.frame(data)
      .MUST.isDataFrame(data, data_name)
      .MUST.hasMultipleRows(data, data_name)
      .MUST.isIgraph(net, deparse(substitute(net)))
      .MUST.doesIgraphMatchData(net, data, deparse(substitute(net)), data_name)
      cluster_fun <- .check(cluster_fun, .isClusterFun, "fast_greedy")
      cluster_id_name <- make.names(
        .check(cluster_id_name, .isString, "cluster_id")
      )
      overwrite <- .checkTF(overwrite, FALSE)
      data <- .addClusterMembershipWithoutList(
        data, deparse(substitute(data)), net, cluster_fun, cluster_id_name,
        overwrite, verbose, ...
      )
      return(data)
    } else {
      lifecycle::deprecate_soft(
        when = "1.0.1",
        what = "addClusterMembership(data)",
        details = paste(
          "Please pass the list of network objects containing the",
          "igraph and node metadata to the ", sQuote("net"), " argument",
          "and discontinue use of the ", sQuote("data"), " argument"
        )
      )
    }
  }
  .MUST.isBaseNetworkOutput(net, "net")
  cluster_fun <- .check(cluster_fun, .isClusterFun, "fast_greedy")
  fun <- .getClusterFun(cluster_fun)
  cluster_id_name <- make.names(
    .check(cluster_id_name, .isString, "cluster_id")
  )
  overwrite <- .checkTF(overwrite, FALSE)
  msg <- .makemsg(verbose)
  if (!overwrite && cluster_id_name %in% names(net$node_data)) {
    message(
      "Node metadata already contains a variable named ",
      sQuote(cluster_id_name), ".\n",
      "To recompute cluster membership, call ",
      sQuote("addClusterMembership()"), " with ",
      sQuote("overwrite = TRUE"), ".\n",
      "To add a new cluster membership variable, ",
      "use a different value for ", sQuote("cluster_id_name"), "."
    )
    return(net)
  }
  overwrite <- overwrite && cluster_id_name %in% names(net$node_data)
  msg("Partitioning the network graph into clusters...", newline = FALSE)
  net$node_data[[cluster_id_name]] <- as.factor(as.integer(
    fun(net$igraph, ...)$membership
  ))
  msg(" Done.")
  msg("Cluster membership variable ", sQuote(cluster_id_name),
      " added to node metadata."
  )
  net <- .updateClusterDetails(net, cluster_id_name, cluster_fun, overwrite)
  net
}

addClusterStats <- function(
    net,
    cluster_id_name = "cluster_id",
    seq_col = NULL,
    count_col = NULL,
    degree_col = "degree",
    cluster_fun = "fast_greedy",
    overwrite = FALSE,
    verbose = FALSE,
    ...
) {
  .MUST.isBaseNetworkOutput(net, "net")
  seq_col <- .check(seq_col, .isSeqColrefs, NULL, data = net$node_data,
                    ornull = TRUE
  )
  if (is.null(seq_col) && .hasDetails(net) &&
      .hasElem(net$details, "seq_col") && is.character(net$details$seq_col) &&
      .isSeqColrefs(net$details$seq_col, net$node_data)
  ) {
    seq_col <- net$details$seq_col
  }
  seq_col <- .convertColRef(seq_col, net$node_data)
  count_col <- .checkCountCol(count_col, net$node_data, NULL)
  count_col <- .convertColRef(count_col, net$node_data)
  degree_col <- .check(degree_col, .isString, "degree", ornull = TRUE)
  if (degree_col %in% names(net$node_data)) {
    .stopifnot(.isNonnegIntegerVector(net$node_data[[degree_col]]),
               "degree_col",
               "specifies an existing column in the node metadata",
               "that does not contain nonnegative integer values"
    )
  }
  cluster_fun <- .check(cluster_fun, .isClusterFun, "fast_greedy")
  cluster_id_name <- make.names(
    .check(cluster_id_name, .isString, "cluster_id")
  )
  if (cluster_id_name %in% names(net$node_data)) {
    .stopifnot(
      .isPosIntegerVector(net$node_data[[cluster_id_name]], factor_ok = TRUE),
      name = NULL,
      "the column", dQuote(cluster_id_name), "is present in",
      sQuote(paste0(deparse(substitute(net)), "$node_data")),
      "but does not contain positive integer values"
    )
  }
  msg <- .makemsg(verbose)
  msg("Obtaining cluster properties...")
  if (!overwrite && "cluster_data" %in% names(net)) {
    message(sQuote(paste0(deparse(substitute(net)), "$cluster_data")),
            " already exists.\n",
            "To overwrite, call ", sQuote("addClusterStats()"),
            " with ", sQuote("overwrite = TRUE")
    )
    return(net)
  }
  if (!overwrite && isTRUE(cluster_id_name %in% names(net$node_data))) {
    msg("Using cluster ID variable ", sQuote(cluster_id_name),
        " present in node metadata")
  } else {
    net <- addClusterMembership(
      net = net, cluster_fun = cluster_fun, cluster_id_name = cluster_id_name,
      overwrite = overwrite, verbose = verbose, ... = ...
    )
  }
  if (is.null(degree_col)) {
    net$node_data$degree <- igraph::degree(net$igraph)
    degree_col <- "degree"
    msg("Variable ", sQuote("degree"), " added to node metadata")
  } else if (!isTRUE(degree_col %in% names(net$node_data))) {
    msg("Variable ", sQuote(degree_col), " specified by ", sQuote("degree_col"),
        " is not present in node metadata"
    )
    net$node_data[[degree_col]] <- igraph::degree(net$igraph)
    msg("Variable ", sQuote(degree_col),
        " for network degree added to node metadata"
    )
  }

  net$cluster_data <- .computeClusterStats(
    net$node_data,
    net$adjacency_matrix,
    cluster_id_name,
    degree_col,
    seq_col = seq_col,
    count_col = count_col,
    verbose = FALSE
  )
  msg("Data frame ", sQuote("cluster_data"),
      " added to list of network objects."
  )
  if (isTRUE("details" %in% names(net))) {
    net$details$cluster_data_goes_with <- cluster_id_name
    net$details$count_col_for_cluster_data <-
      ifelse(is.null(count_col), yes = as.character(NA), no = count_col)
  }
  net
}

getClusterStats <- function(
    data,
    adjacency_matrix,
    seq_col = NULL,
    count_col = NULL,
    cluster_id_col = "cluster_id",
    degree_col = NULL,
    cluster_fun = deprecated(),
    verbose = FALSE
) {
  if (lifecycle::is_present(cluster_fun)) {
    lifecycle::deprecate_soft(
      when = "1.0.1",
      what = "getClusterStats(cluster_fun)",
      details =
        "Cluster membership must be computed prior to calling getClusterStats()"
    )
  }
  data_name <- deparse(substitute(data))
  mat_name <- deparse(substitute(adjacency_matrix))
  data <- as.data.frame(data)
  .MUST.isDataFrame(data, data_name)
  .MUST.isAdjacencyMatrix(adjacency_matrix, mat_name)
  .MUST.doesDataMatchMatrix(data, adjacency_matrix, data_name, mat_name)
  .MUST.isDataColref(cluster_id_col, data)
  cluster_id_col <- .convertColRef(cluster_id_col, data)
  .MUST.isPosIntegerVector(data[[cluster_id_col]], factor_ok = TRUE)
  seq_col <- .check(seq_col, .isSeqColrefs, NULL, data = data, ornull = TRUE)
  seq_col <- .convertColRef(seq_col, data)
  count_col <- .checkCountCol(count_col, data, NULL)
  count_col <- .convertColRef(count_col, data)
  degree_col <- .check(degree_col, .isStringOrPosInt, "degree", ornull = TRUE)
  if (is.numeric(degree_col)) {
    degree_col <- .check(degree_col, .isDataColref, "degree", data = data)
    degree_col <- .convertColRef(degree_col, data)
  }
  if (isTRUE(degree_col %in% names(data))) {
    .stopifnot(.isNonnegIntegerVector(data[[degree_col]]),
               "degree_col",
               "specifies an existing column in the node metadata",
               "that does not contain nonnegative integer values"
    )
  }
  if (is.null(degree_col) || !isTRUE(degree_col %in% names(data))) {
    degree_col <- "degree"
    data$degree <- igraph::degree(generateNetworkGraph(adjacency_matrix))
  }
  .computeClusterStats(
    data,
    adjacency_matrix,
    cluster_id_col,
    degree_col,
    seq_col = seq_col,
    count_col = count_col,
    verbose = FALSE
  )
}

.getClusterFun <- function(name) {
  if (.isString(name)) {
    if (pmatch(name,
               c("cluster_edge_betweenness", "edge_betweenness", "betweenness"),
               nomatch = 0
    )) {
      return(igraph::cluster_edge_betweenness)
    }
    if (pmatch(name, c("cluster_fast_greedy", "fast_greedy", "greedy"), 0)) {
      return(igraph::cluster_fast_greedy)
    }
    if (pmatch(name, c("cluster_infomap", "infomap"), 0)) {
      return(igraph::cluster_infomap)
    }
    if (pmatch(name, c("cluster_label_prop", "label_prop", "prop"), 0)) {
      return(igraph::cluster_label_prop)
    }
    if (pmatch(name, c("cluster_leading_eigen", "leading_eigen", "eigen"), 0)) {
      return(igraph::cluster_leading_eigen)
    }
    if (pmatch(name, c("cluster_leiden", "leiden"), 0)) {
      return(igraph::cluster_leiden)
    }
    if (pmatch(name, c("cluster_louvain", "louvain"), 0)) {
      return(igraph::cluster_louvain)
    }
    if (pmatch(name, c("cluster_optimal", "optimal"), 0)) {
      return(igraph::cluster_optimal)
    }
    if (pmatch(name, c("cluster_spinglass", "spinglass"), 0)) {
      return(igraph::cluster_spinglass)
    }
    if (pmatch(name, c("cluster_walktrap", "walktrap"), 0)) {
      return(igraph::cluster_walktrap)
    }
  }
  warning("invalid specification for clustering algorithm, defaulting to",
          dQuote("fast_greedy")
  )
  igraph::cluster_fast_greedy
}

.addClusterMembershipWithoutList <- function(
    data,
    data_name,
    net,
    cluster_fun = "fast_greedy",
    cluster_id_name = "cluster_id",
    overwrite = FALSE,
    verbose = FALSE,
    ...
) {
  fun <- .getClusterFun(cluster_fun)
  msg <- .makemsg(verbose)
  if (!overwrite && cluster_id_name %in% names(data)) {
    message(
      sQuote(data_name), " already contains a variable named ",
      sQuote(cluster_id_name), ".\n",
      "To recompute cluster membership, call ",
      sQuote("addClusterMembership()"),
      "with ", sQuote("overwrite = TRUE"), ".\n",
      "To add a new cluster membership variable, ",
      "use a different value for ", sQuote("cluster_id_name")
    )
    return(data)
  }
  msg("Partitioning the network graph into clusters...", newline = FALSE)
  data[[cluster_id_name]] <- as.factor(as.integer(fun(net, ...)$membership))
  msg(" Done.\nVariable ", sQuote(cluster_id_name),
      " added to ", sQuote(data_name), "."
  )
  data
}

.calc_clustersnodes_components <- function(
  data, out, adjacency_matrix, cluster_id_col
) {
  cluster_id2node_ids <- split(
    seq_along(data[[cluster_id_col]]),
    data[[cluster_id_col]]
  )
  cluster_ids <- as.character(out$cluster_id)
  clusters_nodes <- cluster_id2node_ids[cluster_ids]
  clusters_adj <- lapply(
    clusters_nodes,
    function(node_ids) adjacency_matrix[node_ids, node_ids]
  )
  components <- lapply(
    clusters_adj,
    igraph::graph_from_adjacency_matrix,
    mode = "undirected",
    diag = FALSE
  )
  return(list(clusters_nodes = clusters_nodes, components = components))
}

.calc_seq_w_max_count <- function(
  clusters_nodes, out, clusters_counts, data, seq_col
) {
  vapply(
    seq_along(clusters_nodes),
    function(i) {
      max_count <- out$max_count[[i]]
      id_incluster <- match(max_count, clusters_counts[[i]])
      id_global <- clusters_nodes[[i]][id_incluster]
      data[[seq_col]][id_global]
    },
    character(1)
  )
}

.calc_igraph_stats <- function(components, out) {
  calculate_igraph_stat <- function(
    stat_function, output_value = NULL, ...
  ) {
    vapply(
      components,
      function(component) {
        result <- stat_function(component, ...)
        if (!is.null(output_value)) {
          result <- result[[output_value]]
        }
        result
      },
      numeric(1)
    )
  }
  out$eigen_centrality_eigenvalue <- calculate_igraph_stat(
    igraph::eigen_centrality,
    output_value = "value",
    directed = TRUE, # why directed?
    weights = NULL
  )
  out$eigen_centrality_index <- calculate_igraph_stat(
    igraph::centr_eigen,
    output_value = "centralization",
    directed = TRUE, # why directed?
    normalized = TRUE
  )
  out$closeness_centrality_index <- calculate_igraph_stat(
    igraph::centr_clo,
    output_value = "centralization",
    mode = "all",
    normalized = TRUE
  )
  out$degree_centrality_index <- calculate_igraph_stat(
    igraph::centr_degree,
    output_value = "centralization",
    mode = "in",
    normalized = TRUE
  )
  out$edge_density <- calculate_igraph_stat(
    igraph::edge_density
  )
  out$global_transitivity <- calculate_igraph_stat(
    igraph::transitivity,
    type = "global"
  )
  out$assortativity <- calculate_igraph_stat(
    igraph::assortativity_degree
  )
  out$diameter_length <-  calculate_igraph_stat(
    igraph::diameter,
    directed = TRUE # why directed?
  ) + 1
  out
}

.calc_seq_stats <- function(
  data, seq_col, clusters_nodes, clusters_degrees
) {
  seq_lens <- nchar(data[[seq_col]])
  clusters_seq_lens <- lapply(clusters_nodes, function(nodes) seq_lens[nodes])
  mean_seq_length <- round(vapply(clusters_seq_lens, mean, numeric(1)), 2)
  seq_w_max_degree <- vapply(
    seq_along(clusters_nodes),
    function(i) {
      max_degree <- max(clusters_degrees[[i]])
      id_incluster <- match(max_degree, clusters_degrees[[i]])
      id_global <- clusters_nodes[[i]][id_incluster]
      data[[seq_col]][id_global]
    },
    character(1)
  )

  list(mean_seq_length = mean_seq_length, seq_w_max_degree = seq_w_max_degree)
}

.computeClusterStats <- function(
  data,
  adjacency_matrix,
  cluster_id_col,
  degree_col,
  seq_col = NULL,
  count_col = NULL,
  verbose = FALSE
) {
  msg <- .makemsg(verbose)
  out <- as.data.frame(table(data[[cluster_id_col]]))
  colnames(out) <- c("cluster_id", "node_count")
  num_clusters <- nrow(out)
  msg("Computing statistics for the ", num_clusters,
      " clusters in the network...",
      newline = FALSE
  )

  clust_comp <- .calc_clustersnodes_components(
    data, out, adjacency_matrix, cluster_id_col
  )
  clusters_nodes <- clust_comp$clusters_nodes
  components <- clust_comp$components

  out <- .calc_igraph_stats(components, out)

  clusters_degrees <- lapply(
    clusters_nodes, function(nodes) data[[degree_col]][nodes]
  )
  out$max_degree <- vapply(clusters_degrees, max, numeric(1))
  out$mean_degree <- round(vapply(clusters_degrees, mean, numeric(1)), 2)

  if (length(seq_col) == 1) {
    seq_stat <- .calc_seq_stats(data, seq_col, clusters_nodes, clusters_degrees)
    out$mean_seq_length <- seq_stat$mean_seq_length
    out$seq_w_max_degree <- seq_stat$seq_w_max_degree
  } else if (length(seq_col) == 2) {
    seq_a_stat <- .calc_seq_stats(
      data, seq_col[[1]], clusters_nodes, clusters_degrees
    )
    seq_b_stat <- .calc_seq_stats(
      data, seq_col[[2]], clusters_nodes, clusters_degrees
    )
    out$mean_A_seq_length <- seq_a_stat$mean_seq_length
    out$mean_B_seq_length <- seq_b_stat$mean_seq_length
    out$A_seq_w_max_degree <- seq_a_stat$seq_w_max_degree
    out$B_seq_w_max_degree <- seq_b_stat$seq_w_max_degree
  }

  if (!is.null(count_col)) {
    clusters_counts <- lapply(
      clusters_nodes, function(idxs) data[[count_col]][idxs]
    )
    out$max_count <- vapply(clusters_counts, max, numeric(1))
    out$agg_count <- vapply(clusters_counts, sum, numeric(1))
    if (length(seq_col) == 1) {
      out$seq_w_max_count <- .calc_seq_w_max_count(
        clusters_nodes, out, clusters_counts, data, seq_col
      )
    } else if (length(seq_col) == 2) {
      out$A_seq_w_max_count <- .calc_seq_w_max_count(
        clusters_nodes, out, clusters_counts, data, seq_col[[1]]
      )
      out$B_seq_w_max_count <- .calc_seq_w_max_count(
        clusters_nodes, out, clusters_counts, data, seq_col[[2]]
      )
    }
  }

  msg(" Done.")
  out
}

.updateClusterDetails <- function(
    net, cluster_id_name, cluster_fun, overwrite
) {
  if (!isTRUE("details" %in% names(net))) {
    return(net)
  }
  if (is.list(net$details) && "clusters_in_network" %in% names(net$details)) {
    last_entry <- length(net$details$clusters_in_network)
    if (overwrite) {
      net$details$clusters_in_network[[last_entry]] <-
        length(unique(net$node_data[[cluster_id_name]]))
      net$details$cluster_id_variable[[last_entry]] <- cluster_id_name
      names(net$details$clusters_in_network)[[last_entry]] <- cluster_fun
      names(net$details$cluster_id_variable)[[last_entry]] <- cluster_fun
    } else {
      net$details$clusters_in_network <- append(
        net$details$clusters_in_network,
        length(unique(net$node_data[[cluster_id_name]]))
      )
      net$details$cluster_id_variable <- append(net$details$cluster_id_variable,
                                                cluster_id_name
      )
      names(net$details$clusters_in_network)[[last_entry + 1]] <- cluster_fun
      names(net$details$cluster_id_variable)[[last_entry + 1]] <- cluster_fun
    }
  } else {
    net$details$clusters_in_network <-
      length(unique(net$node_data[[cluster_id_name]]))
    net$details$cluster_id_variable <- cluster_id_name
    names(net$details$clusters_in_network) <- cluster_fun
    names(net$details$cluster_id_variable) <- cluster_fun
  }
  net
}
