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

# Top-level functions -----------------------------------------------------

findPublicClusters <- function(

  file_list,
  input_type,
  data_symbols = NULL,
  header, sep, read.args,
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
  output_dir,
  output_type = "rds",
  output_dir_unfiltered = NULL,
  output_type_unfiltered = "rds",
  verbose = FALSE,
  ...

) {

  # Process arguments
  if (missing(header)) {
    header <-  switch(input_type, "txt" = FALSE, "table" = FALSE, TRUE)
  }
  if (missing(sep)) {
    sep <- switch(input_type, "csv" = ",", "csv2" = ";", "tsv" = "\t", "")
  }
  if (missing(read.args)) { read.args <- NULL }
  .checkargs.InputFiles(file_list, input_type, data_symbols,
                        header, sep, read.args
  )
  if (input_type %in% c("csv", "csv2", "tsv", "txt", "table")) {
    read.args <- .checkReadArgs(read.args, header, sep)
  }
  .MUST.isCharOrIntegerVector(seq_col, "seq_col")
  .MUST.hasLength(seq_col, c(1, 2))
  .MUST.isString(output_dir, "output_dir")
  .createOutputDir(output_dir)
  .requireOutputDir(output_dir)
  output_type <- .checkOutputType(output_type, "findPublicClusters")
  sample_ids <- .checkIDs(sample_ids, length(file_list),
                          default = paste0("Sample", 1:length(file_list))
  )
  sample_ids <- as.vector(sample_ids, mode = "character")
  count_col <- .check(count_col, .isStringOrPosInt, NULL, ornull = TRUE)
  min_seq_length <- .check(min_seq_length, .isNonneg, NULL, ornull = TRUE)
  drop_matches <- .check(drop_matches, .isString, NULL, ornull = TRUE)
  top_n_clusters <- .check(top_n_clusters, .isPosInt, 20)
  min_node_count <- .check(min_node_count, .isPosInt, 10)
  min_clone_count <- .check(min_clone_count, .isPos, NULL, ornull = TRUE)
  plots <- .checkTF(plots, FALSE)
  if (isTRUE(plots)) {
    print_plots <- .checkTF(print_plots, FALSE)
    plot_title <- .check(plot_title, .isString, "auto", ornull = TRUE)
    color_nodes_by <- .check(color_nodes_by, .isCharVector, "cluster_id",
                             ornull = TRUE
    )
    color_nodes_by <- .colorNodesBy.findPublicClusters(color_nodes_by)
  }
  output_dir_unfiltered <- .check(output_dir_unfiltered, .isString, NULL,
                                  ornull = TRUE
  )
  .createOutputDir(output_dir_unfiltered)
  if (!is.null(output_dir_unfiltered) &&
      !isTRUE(dir.exists(output_dir_unfiltered))
  ) {
    warning("directory ", dQuote(output_dir_unfiltered),
            " specified for ", sQuote("output_dir_unfiltered"),
            " does not exist and could not be created. ",
            "Unfiltered data for sample networks will not be saved"
    )
    output_dir_unfiltered <- NULL
  }
  if (!is.null(output_dir_unfiltered)) {
    output_type_unfiltered <- .check(output_type_unfiltered, .isOutputType,
                                     "rds"
    )
  }
  if (all(plots, !print_plots, is.null(output_dir_unfiltered))) {
    warning(
      "ignoring `plots = TRUE` ",
      "since ", sQuote("print_plots"), " is ", dQuote("FALSE"), " and ",
      sQuote("output_dir_unfiltered"), " is ", dQuote("NULL")
    )
    plots <- FALSE
  }
  msg <- .makemsg(verbose)

  # Execute Tasks
  msg("Beginning search for public clusters...")
  for (i in 1:length(file_list)) {
    msg("Processing sample ", i, " of ", length(file_list),
        " (", sample_ids[[i]], ")..."
    )
    .findPublicClustersOneSample(top_n_clusters,
                                 min_node_count, min_clone_count,
                                 i, file_list[[i]], sample_ids[[i]],
                                 input_type, data_symbols,
                                 header, sep, read.args,
                                 seq_col, count_col,
                                 min_seq_length, drop_matches,
                                 plots, print_plots,
                                 plot_title, color_nodes_by,
                                 output_dir, output_type,
                                 output_dir_unfiltered, output_type_unfiltered,
                                 verbose, msg, ...
    )
  }
  msg("All samples complete. Filtered data is located in the following ",
      "directory:\n  ", output_dir
  )
  if (!is.null(output_dir_unfiltered)) {
    msg("Unfiltered data for full sample networks is located in the following ",
        "directory:\n  ", output_dir_unfiltered
    )
  }
  invisible(TRUE)

}


buildPublicClusterNetwork <- function(
    file_list,
    input_type = "rds",
    data_symbols = "ndat",
    header = TRUE, sep,
    read.args = list(row.names = 1),
    seq_col,
    drop_isolated_nodes = FALSE,
    node_stats = deprecated(),
    stats_to_include = deprecated(),
    cluster_stats = deprecated(),
    color_nodes_by = "SampleID",
    color_scheme = "turbo",
    plot_title = "Global Network of Public Clusters",
    output_dir = NULL,
    output_name = "PublicClusterNetwork",
    verbose = FALSE,
    ...
) {
  .checkDeprecated.buildPublicClusterNetwork(
    node_stats, stats_to_include, cluster_stats
  )
  data <- loadDataFromFileList(file_list, input_type, data_symbols,
                               header, sep, read.args
  )
  tmp_rownames <- rownames(data)
  tmp_rownames <- sapply(
    tmp_rownames, # Strip file*. prefix
    function(x) { substr(x, start = 1 + regexpr("\\.", x), stop = nchar(x)) },
    USE.NAMES = FALSE
  )
  rownames(data) <- paste0(data$SampleID, ".", tmp_rownames)
  color_nodes_by <- .check(color_nodes_by, .isCharVector, "SampleID",
                           ornull = TRUE
  )
  color_nodes_by <- .colorNodesBy.buildPublicClusterNetwork(color_nodes_by)
  color_scheme <- .checkColorScheme(color_scheme, color_nodes_by, "Viridis")
  plot_title <- .check(plot_title, .isString,
                       "Global Network of Public Clusters", ornull = TRUE)
  if (!is.null(output_dir)) {
    output_name <- .checkOutputName(output_name, "PublicClusterNetwork")
  }
  msg <- .makemsg(verbose)
  msg("Building global network of public clusters...")
  net <- buildRepSeqNetwork(
    data = data,
    seq_col = seq_col,
    drop_isolated_nodes = drop_isolated_nodes,
    node_stats = TRUE,
    stats_to_include = "all",
    cluster_stats = TRUE,
    cluster_id_name = "ClusterIDPublic",
    color_nodes_by = color_nodes_by,
    color_scheme = color_scheme,
    plot_title = plot_title,
    output_dir = output_dir,
    output_name = output_name,
    verbose = verbose,
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
    header = TRUE, sep, read.args,
    seq_col = "seq_w_max_count",
    count_col = "agg_count",
    dist_type = "hamming",
    dist_cutoff = 1,
    cluster_fun = "fast_greedy",
    plots = TRUE,
    print_plots = FALSE,
    plot_title = "auto",
    plot_subtitle = "auto",
    color_nodes_by = "SampleID",
    color_scheme = "turbo",
    ...,
    output_dir = NULL,
    output_type = "rds",
    output_name = "PubClustByRepresentative",
    pdf_width = 12, pdf_height = 10,
    verbose = FALSE
) {
  data <- loadDataFromFileList(file_list, input_type, data_symbols,
                               header, sep, read.args
  )
  output_dir <- .check(output_dir, .isString, NULL, ornull = TRUE)
  .createOutputDir(output_dir)
  output_dir <- .checkOutputDir(output_dir)
  plots <- .checkTF(plots, TRUE)
  if (isTRUE(plots)) {
    print_plots <- .checkTF(print_plots, TRUE)
    plot_title <- .check(plot_title, .isString, "auto", ornull = TRUE)
    plot_subtitle <- .check(plot_subtitle, .isString, "auto", ornull = TRUE)
    if (!is.null(color_nodes_by)) {
      color_nodes_by[color_nodes_by == "cluster_id"] <- "ClusterIDPublic"
    }
    color_nodes_by <- .checkColorNodesBy(
      color_nodes_by, data, node_stats = TRUE, cluster_stats = FALSE,
      plots = plots, cluster_id_name = "ClusterIDPublic",
      stats_to_include = "all", default = "SampleID"
    )
    color_nodes_by <- .check(color_nodes_by, .isCharVector, "SampleID",
                             ornull = TRUE
    )
    if (!is.null(output_dir)) {
      pdf_width <- .check(pdf_width, .isPos, 12)
      pdf_height <- .check(pdf_width, .isPos, 10)
    }
  }
  if (isTRUE(plots) || !is.null(output_dir)) {
    output_name <- .checkOutputName(output_name, "PubClustByRepresentative")
  }
  if (!is.null(output_dir)) {
    output_type <- .check(output_type, .isOutputType, "rds")
  }
  msg <- .makemsg(verbose)
  msg("Building network of public clusters using a representative",
      "sequence from each cluster:"
  )
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
    cluster_id_name = "ClusterIDPublic",
    cluster_stats = FALSE,
    plots = FALSE,
    output_dir = NULL,
    verbose = verbose
  )
  if (is.null(net)) {
    return(NULL)
  }
  net$node_data$RepresentativeSeq <- net$node_data[[seq_col]]
  net$cluster_data <- as.data.frame(table(net$node_data$ClusterIDPublic))
  colnames(net$cluster_data) <- c("cluster_id", "node_count")
  msg(
    "Performing clustering on the nodes in the new network resulted in ",
    nrow(net$cluster_data), " clusters.",
    "\nComputing network properties of the new clusters...",
    newline = FALSE
  )
  net <- .addClusterOnClusterStats(net)
  msg(" Done.")
  if (plots) {
    net <- addPlots(
      net,
      print_plots,
      .makePlotTitle(plot_title, "pub_clust_rep"),
      .makePlotSubtitle(plot_subtitle, "pub_clust_rep", seq_col,
                        dist_type, dist_cutoff
      ),
      color_nodes_by = color_nodes_by,
      color_scheme = color_scheme, verbose = verbose,
      ...
    )
  }
  saveNetwork(net, output_dir, output_type, output_name, pdf_width, pdf_height,
              verbose
  )
  msg("All tasks complete.")
  invisible(net)
}



# Helper functions --------------------------------------------------------
.findPublicClustersOneSample <- function(
    top_n_clusters, min_node_count, min_clone_count,
    sample_index, input_file, sample_id, input_type, data_symbols, header, sep,
    read.args,
    seq_col, count_col, min_seq_length, drop_matches, plots, print_plots,
    plot_title, color_nodes_by, output_dir, output_type, output_dir_unfiltered,
    output_type_unfiltered, verbose, msg, ...
) {
  data <- .loadDataFromFile(input_file, input_type, data_symbols, header, sep,
                            read.args
  )
  if (!is.null(plot_title) && plot_title == "auto") {
    plot_title <- paste0("Sample: ", sample_id)
  }
  current_name <- paste0("sample_", sample_index)
  if (nchar(.sanitizeFilenamePart(sample_id)) > 0) {
    current_name <- paste0(current_name, "_", .sanitizeFilenamePart(sample_id))
  }
  net <- buildRepSeqNetwork(
    data = data, seq_col = seq_col, count_col = count_col,
    min_seq_length = min_seq_length, drop_matches = drop_matches,
    node_stats = TRUE, stats_to_include = "all", cluster_stats = TRUE,
    cluster_id_name = "ClusterIDInSample",
    plots = plots, print_plots = print_plots, plot_title = plot_title,
    color_nodes_by = color_nodes_by, output_dir = output_dir_unfiltered,
    output_type = output_type_unfiltered,
    output_name = current_name,
    verbose = verbose,
    ...
  )
  if (is.null(net)) {
    warning("unable to build network for current sample. Proceeding to next")
    return(invisible(NULL))
  }
  msg("Filtering clusters in the current sample...", newline = FALSE)
  ndat <- net$node_data
  cdat <- net$cluster_data
  if (is.null(count_col)) { min_clone_count <- NULL }
  cdat <- .filterClustersBySize(cdat, top_n_clusters,
                                min_node_count, min_clone_count
  )
  ndat <- ndat[ndat$ClusterIDInSample %in% cdat$cluster_id, , drop = FALSE]
  msg(" Done.\n", nrow(cdat), " clusters (", nrow(ndat), " nodes) remain.",
      newline = FALSE
  )
  names(ndat) <- .renameSampleLevelNodeStats(names(ndat))
  names(cdat)[names(cdat) == "cluster_id"] <- "ClusterIDInSample"
  names(cdat)[names(cdat) == "transitivity"] <- "SampleLevelTransitivity"
  ndat$SampleID <- as.character(sample_id)
  cdat$SampleID <- as.character(sample_id)
  .savePublicClustersOneSample(ndat, cdat, output_type, output_dir, sample_id,
                               sample_index, msg
  )
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
    ndat, cdat, output_type, output_dir, sample_id, sample_index, msg
) {
  msg(" Saving results...", newline = FALSE)
  node_dir <- file.path(output_dir, "node_meta_data")
  cluster_dir <- file.path(output_dir, "cluster_meta_data")
  .createOutputDir(node_dir)
  .createOutputDir(cluster_dir)
  file_stem <- sample_index
  if (nchar(.sanitizeFilenamePart(sample_id)) > 0) {
    file_stem <- paste0(file_stem, "_", .sanitizeFilenamePart(sample_id))
  } else {
    file_stem <- paste0("sample_", file_stem)
  }
  node_file <- file.path(output_dir, "node_meta_data",
                         paste0(file_stem, ".", output_type)
  )
  cluster_file <- file.path(output_dir, "cluster_meta_data",
                            paste0(file_stem, ".", output_type)
  )
  if (output_type == "rda") {
    save(ndat, file = node_file)
    save(cdat, file = cluster_file)
  } else if (output_type == "csv") {
    utils::write.csv(ndat, file = node_file)
    utils::write.csv(cdat, file = cluster_file, row.names = FALSE)
  } else {
    saveRDS(ndat, file = node_file)
    saveRDS(cdat, file = cluster_file)
  }
  msg(" Done.")
}

.addClusterOnClusterStats <- function(net) {
  ndat <- net$node_data
  cdat <- net$cluster_data
  adjacency_matrix <- net$adjacency_matrix

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
    cluster_id <- which(cdat$cluster_id == i)
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

    cluster <- generateNetworkGraph(
      as.matrix(adjacency_matrix[node_ids, node_ids])
    )
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
  net$cluster_data <- cdat
  net
}


.colorNodesBy.findPublicClusters <- function(arg) {

  if (!is.null(arg)) {
    arg[arg == "cluster_id"] <-
      "ClusterIDInSample"
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
    arg[arg == "cluster_id"] <-
      "ClusterIDPublic"
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
