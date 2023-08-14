
# Top-level functions -----------------------------------------------------

findPublicClusters <- function(

  file_list,
  input_type,
  data_symbols = NULL,
  header = TRUE, sep = "",
  sample_ids = paste0("Sample", 1:length(file_list)),
  seq_col,
  count_col = NULL,
  min_seq_length = 3,
  drop_matches = "[*|_]",
  top_n_clusters = 20,
  min_node_count = 10,
  min_clone_count = 100,
  plots = FALSE,
  print_plots = FALSE,
  plot_title = "auto",
  color_nodes_by = "cluster_id",
  output_dir =
    file.path(getwd(), "public_clusters"),
  output_type = "rds",
  output_dir_unfiltered = NULL,
  output_type_unfiltered = "rds",
  ...

) {

  .checkargs.findPublicClusters(
    file_list, input_type, data_symbols, header, sep, sample_ids, seq_col,
    count_col, min_seq_length, drop_matches, top_n_clusters, min_node_count,
    min_clone_count, plots, print_plots, plot_title, color_nodes_by, output_dir,
    output_type, output_dir_unfiltered, output_type_unfiltered
  )
  color_nodes_by <- .colorNodesBy.findPublicClusters(color_nodes_by)
  sample_ids <- as.character(sample_ids)
  if (all(plots, !print_plots, is.null(output_dir_unfiltered))) {
    plots <- FALSE
  }
  .ensureOutputDir(output_dir)
  cat(paste0("<<< Beginning search for public clusters >>>\n"))
  for (i in 1:length(file_list)) {
    cat(paste0("Processing sample ", i, " of ", length(file_list), ": ",
               sample_ids[[i]], "\n"
    ))
    .findPublicClustersOneSample(
      top_n_clusters, min_node_count, min_clone_count,
      file_list[[i]], sample_ids[[i]], input_type, data_symbols, header, sep,
      seq_col, count_col, min_seq_length, drop_matches, plots, print_plots,
      plot_title, color_nodes_by, output_dir, output_type,
      output_dir_unfiltered, output_type_unfiltered, ...
    )
    cat(
      "----------------------------------------------------------------------\n"
    )
  }
  cat(paste0("All samples complete. ",
             "Filtered data is located in the following directory:",
             "\n  ", output_dir, "\n"
  ))
  if (!is.null(output_dir_unfiltered)) {
    cat(paste0("Unfiltered data for full sample networks ",
               "is located in the following directory:",
               "\n  ", output_dir_unfiltered, "\n"
    ))
  }
  invisible(TRUE)

}


buildPublicClusterNetwork <- function(
    file_list =
      list.files(file.path(getwd(), "public_clusters", "node_meta_data")),
    input_type = "rds",
    data_symbols = "ndat",
    header = TRUE, sep = "",
    seq_col,
    drop_isolated_nodes = FALSE,
    node_stats = deprecated(),
    stats_to_include = deprecated(),
    cluster_stats = deprecated(),
    color_nodes_by = "SampleID",
    color_scheme = "turbo",
    plot_title = "Global Network of Public Clusters",
    output_dir = file.path(getwd(), "public_clusters"),
    output_name = "PublicClusterNetwork",
    ...
) {
  .checkDeprecated.buildPublicClusterNetwork(
    node_stats, stats_to_include, cluster_stats
  )
  .checkargs.buildPublicClusterNetwork(
    file_list, input_type, data_symbols, header, sep, seq_col,
    drop_isolated_nodes,
    color_nodes_by, color_scheme, plot_title, output_dir, output_name
  )
  .createOutputDir(output_dir)
  color_nodes_by <- .colorNodesBy.buildPublicClusterNetwork(color_nodes_by)
  data <- loadDataFromFileList(
    file_list, input_type, data_symbols, header, sep
  )
  cat("Building network of public clusters:\n")
  net <- buildRepSeqNetwork(
    data = data,
    seq_col = seq_col,
    drop_isolated_nodes = drop_isolated_nodes,
    node_stats = TRUE,
    stats_to_include = "all",
    cluster_stats = TRUE,
    color_nodes_by = color_nodes_by,
    color_scheme = color_scheme,
    output_dir = output_dir,
    output_name = output_name,
    ...
  )
  if (is.null(net)) {
    return(NULL)
  }
  names(net$node_data) <- .renamePublicNodeStats(names(net$node_data))
  invisible(net)
}


buildPublicClusterNetworkByRepresentative <- function(
    file_list,
    input_type = "rds",
    data_symbols = "cdat",
    header = TRUE, sep = "",
    seq_col = "seq_w_max_count",
    count_col = "agg_count",
    dist_type = "hamming",
    dist_cutoff = 1,
    cluster_fun = cluster_fast_greedy,
    plots = TRUE,
    print_plots = TRUE,
    plot_title = "auto",
    plot_subtitle = "auto",
    color_nodes_by = "SampleID",
    color_scheme = "turbo",
    ...,
    output_dir = file.path(getwd(), "public_clusters"),
    output_type = "rda",
    output_name = "PubClustByRepresentative",
    pdf_width = 12, pdf_height = 10
) {
  .checkargs.buildPublicClusterNetworkByRep(
    file_list, input_type, data_symbols, header, sep, seq_col, count_col,
    dist_type, dist_cutoff, cluster_fun, plots, print_plots, plot_title,
    plot_subtitle, color_nodes_by, color_scheme, output_dir, output_type,
    output_name, pdf_width, pdf_height
  )
  if (!is.null(color_nodes_by)) {
    color_nodes_by[color_nodes_by == "cluster_id"] <- "ClusterIDPublic"
  }
  .createOutputDir(output_dir)
  data <- loadDataFromFileList(file_list, input_type, data_symbols, header, sep)
  cat(paste("Building network of public clusters using a representative",
            "sequence from each cluster:\n"
  ))
  net <- buildRepSeqNetwork(
    data = data,
    seq_col = seq_col,
    count_col = count_col,
    subset_cols = NULL,
    min_seq_length = NULL,
    drop_matches = NULL,
    dist_type = dist_type,
    dist_cutoff = dist_cutoff,
    drop_isolated_nodes = FALSE,
    node_stats = TRUE,
    stats_to_include = "all",
    cluster_fun = cluster_fun,
    cluster_stats = FALSE,
    plots = FALSE,
    output_dir = NULL
  )
  if (is.null(net)) {
    return(NULL)
  }
  ndat <- net$node_data
  ndat$RepresentativeSeq <- ndat[[seq_col]]
  names(ndat)[names(ndat) == "cluster_id"] <- "ClusterIDPublic"
  net$node_data <- ndat
  cdat <- as.data.frame(table(ndat[["ClusterIDPublic"]]))
  colnames(cdat) <- c("ClusterIDPublic", "NodeCount")
  cat(paste0(
    "Performing clustering on the nodes in the new network resulted in ",
    nrow(cdat), " clusters.",
    "\nComputing network properties of the new clusters..."
  ))
  net$cluster_data <- .addClusterOnClusterStats(
    ndat, cdat, net$adjacency_matrix
  )
  cat(" Done.\n")
  if (plots) {
    net$plots <- .generateNetworkGraphPlotsGuarded(
      net$igraph,
      net$node_data,
      print_plots,
      .makePlotTitle(plot_title, "pub_clust_rep"),
      .makePlotSubtitle(plot_subtitle, "pub_clust_rep", seq_col),
      color_nodes_by = color_nodes_by,
      color_scheme = color_scheme,
      ...
    )
  }
  saveNetwork(net, output_dir, output_type, output_name, pdf_width, pdf_height)
  cat("All tasks complete.\n")
  invisible(net)
}



# Helper functions --------------------------------------------------------
.findPublicClustersOneSample <- function(
    top_n_clusters, min_node_count, min_clone_count,
    input_file, sample_id, input_type, data_symbols, header, sep,
    seq_col, count_col, min_seq_length, drop_matches, plots, print_plots,
    plot_title, color_nodes_by, output_dir, output_type, output_dir_unfiltered,
    output_type_unfiltered, ...
) {
  data <- .loadDataFromFile(input_file, input_type, data_symbols, header, sep)
  if (!is.null(plot_title) && plot_title == "auto") {
    plot_title <- paste0("Sample: ", sample_id)
  }
  net <- buildRepSeqNetwork(
    data = data, seq_col = seq_col, count_col = count_col,
    min_seq_length = min_seq_length, drop_matches = drop_matches,
    node_stats = TRUE, stats_to_include = "all", cluster_stats = TRUE,
    plots = plots, print_plots = print_plots, plot_title = plot_title,
    color_nodes_by = color_nodes_by, output_dir = output_dir_unfiltered,
    output_type = output_type_unfiltered, output_name = sample_id, ...
  )
  if (is.null(net)) {
    warning("couldn't build network for current sample. Proceeding to next")
    return(invisible(NULL))
  }
  cat(">>> Filtering clusters in the current sample...")
  ndat <- net$node_data
  cdat <- net$cluster_data
  if (is.null(count_col)) { min_clone_count <- NULL }
  cdat <- .filterClustersBySize(cdat, top_n_clusters,
                                min_node_count, min_clone_count
  )
  ndat <- ndat[ndat$cluster_id %in% cdat$cluster_id , ]
  cat(paste0(" Done.
             \n* ", nrow(cdat), " clusters (", nrow(ndat), " nodes) remain. "
  ))
  names(ndat) <- .renameSampleLevelNodeStats(names(ndat))
  names(cdat)[names(cdat) == "cluster_id"] <- "ClusterIDInSample"
  names(cdat)[names(cdat) == "transitivity"] <- "SampleLevelTransitivity"
  ndat$SampleID <- as.character(sample_id)
  cdat$SampleID <- as.character(sample_id)
  .savePublicClustersOneSample(ndat, cdat, output_type, output_dir, sample_id)
}

.filterClustersBySize <- function(
    cdat, top_n_clusters, min_node_count, min_clone_count
) {
  if (nrow(cdat) < top_n_clusters) {
    return(cdat)
  }
  ids_top_n_clusters <- order(-cdat$node_count)[1:top_n_clusters]
  ids_by_node_count <- which(cdat$node_count >= min_node_count)
  filtered_ids <- union(ids_top_n_clusters, ids_by_node_count)
  if (!is.null(min_clone_count)) {
    ids_by_clone_count <- which(cdat$agg_count >= min_clone_count)
    filtered_ids <- union(filtered_ids, ids_by_clone_count)
  }
  cdat[filtered_ids, ]
}

.renameSampleLevelNodeStats <- function(column_names) {
  column_names[column_names == "cluster_id"] <-
    "ClusterIDInSample"
  column_names[column_names == "degree"] <-
    "SampleLevelNetworkDegree"
  column_names[column_names == "transitivity"] <-
    "SampleLevelTransitivity"
  column_names[column_names == "closeness"] <-
    "SampleLevelCloseness"
  column_names[column_names == "centrality_by_closeness"] <-
    "SampleLevelCentralityByCloseness"
  column_names[column_names == "eigen_centrality"] <-
    "SampleLevelEigenCentrality"
  column_names[column_names == "centrality_by_eigen"] <-
    "SampleLevelCentralityByEigen"
  column_names[column_names == "betweenness"] <-
    "SampleLevelBetweenness"
  column_names[column_names == "centrality_by_betweenness"] <-
    "SampleLevelCentralityByBetweenness"
  column_names[column_names == "authority_score"] <-
    "SampleLevelAuthorityScore"
  column_names[column_names == "coreness"] <-
    "SampleLevelCoreness"
  column_names[column_names == "page_rank"] <-
    "SampleLevelPageRank"
  column_names
}

.renamePublicNodeStats <- function(column_names) {
  column_names[column_names == "cluster_id"] <-
    "ClusterIDPublic"
  column_names[column_names == "degree"] <-
    "PublicNetworkDegree"
  column_names[column_names == "transitivity"] <-
    "PublicTransitivity"
  column_names[column_names == "closeness"] <-
    "PublicCloseness"
  column_names[column_names == "centrality_by_closeness"] <-
    "PublicCentralityByCloseness"
  column_names[column_names == "eigen_centrality"] <-
    "PublicEigenCentrality"
  column_names[column_names == "centrality_by_eigen"] <-
    "PublicCentralityByEigen"
  column_names[column_names == "betweenness"] <-
    "PublicBetweenness"
  column_names[column_names == "centrality_by_betweenness"] <-
    "PublicCentralityByBetweenness"
  column_names[column_names == "authority_score"] <-
    "PublicAuthorityScore"
  column_names[column_names == "coreness"] <-
    "PublicCoreness"
  column_names[column_names == "page_rank"] <-
    "PublicPageRank"
  column_names
}

.savePublicClustersOneSample <- function(
    ndat, cdat, output_type, output_dir, sample_id
) {
  cat("Saving results...")
  node_dir <- file.path(output_dir, "node_meta_data")
  cluster_dir <- file.path(output_dir, "cluster_meta_data")
  .createOutputDir(node_dir)
  .createOutputDir(cluster_dir)
  node_file <- file.path(output_dir, "node_meta_data",
                         paste0(sample_id, ".", output_type)
  )
  cluster_file <- file.path(output_dir, "cluster_meta_data",
                            paste0(sample_id, ".", output_type)
  )
  if (output_type == "rda") {
    save(ndat, file = node_file)
    save(cdat, file = cluster_file)
  } else if (output_type == "csv") {
    utils::write.csv(ndat, file = node_file, row.names = FALSE)
    utils::write.csv(cdat, file = cluster_file, row.names = FALSE)
  } else {
    saveRDS(ndat, file = node_file)
    saveRDS(cdat, file = cluster_file)
  }
  cat(" Done.\n")
}

.addClusterOnClusterStats <- function(ndat, cdat, adjacency_matrix) {

  cdat$TotalSampleLevelNodes <- 0
  cdat$TotalCloneCount <- 0
  cdat$MeanOfMeanSeqLength <- 0
  cdat$MeanDegreeInPublicNet <- 0
  cdat$MaxDegreeInPublicNet <- 0
  cdat$SeqWithMaxDegree <- ""
  cdat$MaxCloneCount <- 0
  cdat$SampleWithMaxCloneCount <- ""
  cdat$SeqWithMaxCloneCount <- ""
  cdat$MaxAggCloneCount <- 0
  cdat$SampleWithMaxAggCloneCount <- ""
  cdat$SeqWithMaxAggCloneCount <- ""
  cdat$DiameterLength <- 0
  cdat$Assortativity <- 0
  cdat$GlobalTransitivity <- 0
  cdat$EdgeDensity <- 0
  cdat$DegreeCentralityIndex <- 0
  cdat$ClosenessCentralityIndex <- 0
  cdat$EigenCentralityIndex <- 0
  cdat$EigenCentralityEigenvalue <- 0
  for (i in 1:nrow(cdat)) {
    cluster_id <- which(cdat$ClusterIDPublic == i)
    node_ids <- ndat$ClusterIDPublic == i
    cdat$TotalSampleLevelNodes[[cluster_id]] <-
      sum(ndat[node_ids, "node_count"])
    cdat$TotalCloneCount[[cluster_id]] <- sum(ndat[node_ids, "agg_count"])
    cdat$MeanOfMeanSeqLength[[cluster_id]] <- round(
      mean(ndat[node_ids, "mean_seq_length"]), 2
    )
    cdat$MeanDegreeInPublicNet[[cluster_id]] <- round(
      mean(ndat[node_ids, "degree"]), 2
    )
    max_deg <- max(ndat[node_ids, "degree"])
    cdat$MaxDegreeInPublicNet[[cluster_id]] <- max_deg
    node_id_max_deg <- which(
      node_ids & ndat[["degree"]] == max_deg
    )[[1]]
    cdat$SeqWithMaxDegree[[cluster_id]] <- as.character(
      ndat[[node_id_max_deg, "RepresentativeSeq"]]
    )
    max_count <- max(ndat[node_ids, "max_count"])
    cdat$MaxCloneCount[[cluster_id]] <- max_count
    node_id_max_count <- which(
      node_ids & ndat[["max_count"]] == max_count
    )[[1]]
    cdat$SampleWithMaxCloneCount[[cluster_id]] <-
      ndat[[node_id_max_count, "SampleID"]]
    cdat$SeqWithMaxCloneCount[[cluster_id]] <-
      ndat[[node_id_max_count, "RepresentativeSeq"]]
    max_agg_count <- max(ndat[node_ids, "agg_count"])
    cdat$MaxAggCloneCount[[cluster_id]] <- max_agg_count
    node_id_max_agg_count <- which(
      node_ids & ndat[["agg_count"]] == max_agg_count
    )[[1]]
    cdat$SampleWithMaxAggCloneCount[[cluster_id]] <-
      ndat[node_id_max_agg_count, "SampleID"]
    cdat$SeqWithMaxAggCloneCount[[cluster_id]] <-
      ndat[node_id_max_agg_count, "RepresentativeSeq"]

    cluster <- generateNetworkFromAdjacencyMat(
      as.matrix(adjacency_matrix[node_ids, node_ids]))
    cdat$DiameterLength[[cluster_id]] <- length(
      igraph::get_diameter(cluster, directed = T)
    )
    cdat$Assortativity[[cluster_id]] <-
      igraph::assortativity_degree(cluster, directed = F)
    cdat$GlobalTransitivity[[cluster_id]] <-
      igraph::transitivity(cluster, type = "global")
    cdat$EdgeDensity[[cluster_id]] <-
      igraph::edge_density(cluster, loops = F)
    cdat$DegreeCentralityIndex[[cluster_id]] <-
      igraph::centr_degree(cluster, mode = "in", normalized = T)$centralization
    cdat$ClosenessCentralityIndex[[cluster_id]] <-
      igraph::centr_clo(cluster, mode = "all", normalized = T)$centralization
    cdat$EigenCentralityIndex[[cluster_id]] <-
      igraph::centr_eigen(cluster, directed = T, normalized = T)$centralization
    cdat$EigenCentralityEigenvalue[[cluster_id]] <-
      igraph::eigen_centrality(cluster, directed = T, weights = NA)$value
  }
  cdat
}


.colorNodesBy.findPublicClusters <- function(arg) {

  if (!is.null(arg)) {
    arg[arg == "ClusterIDInSample"] <-
      "cluster_id"
    arg[arg == "SampleLevelNetworkDegree"] <-
      "degree"
    arg[arg == "SampleLevelTransitivity"] <-
      "transitivity"
    arg[arg == "SampleLevelCloseness"] <-
      "closeness"
    arg[arg == "SampleLevelCentralityByCloseness"] <-
      "centrality_by_closeness"
    arg[arg == "SampleLevelEigenCentrality"] <-
      "eigen_centrality"
    arg[arg == "SampleLevelCentralityByEigen"] <-
      "centrality_by_eigen"
    arg[arg == "SampleLevelBetweenness"] <-
      "betweenness"
    arg[arg == "SampleLevelCentralityByBetweenness"] <-
      "centrality_by_betweenness"
    arg[arg == "SampleLevelAuthorityScore"] <-
      "authority_score"
    arg[arg == "SampleLevelCoreness"] <-
      "coreness"
    arg[arg == "SampleLevelPageRank"] <-
      "page_rank"
  }
  arg
}


.colorNodesBy.buildPublicClusterNetwork <- function(arg) {

  if (!is.null(arg)) {
    arg[arg == "ClusterIDPublic"] <-
      "cluster_id"
    arg[arg == "PublicNetworkDegree"] <-
      "degree"
    arg[arg == "PublicTransitivity"] <-
      "transitivity"
    arg[arg == "PublicCloseness"] <-
      "closeness"
    arg[arg == "PublicCentralityByCloseness"] <-
      "centrality_by_closeness"
    arg[arg == "PublicEigenCentrality"] <-
      "eigen_centrality"
    arg[arg == "PublicCentralityByEigen"] <-
      "centrality_by_eigen"
    arg[arg == "PublicBetweenness"] <-
      "betweenness"
    arg[arg == "PublicCentralityByBetweenness"] <-
      "centrality_by_betweenness"
    arg[arg == "PublicAuthorityScore"] <-
      "authority_score"
    arg[arg == "PublicCoreness"] <-
      "coreness"
    arg[arg == "PublicPageRank"] <-
      "page_rank"
  }
  arg
}
