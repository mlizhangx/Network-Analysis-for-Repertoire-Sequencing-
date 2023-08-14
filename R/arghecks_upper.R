
# buildRepSeqNetwork ------------------------------------------------------


.checkArgs.buildRepSeqNetwork <- function(
    data, seq_col, count_col, subset_cols, min_seq_length, drop_matches,
    dist_type, dist_cutoff, drop_isolated_nodes, node_stats, stats_to_include,
    cluster_stats, cluster_fun, plots, print_plots, plot_title,
    plot_subtitle, color_nodes_by, output_dir, output_type, output_name,
    pdf_width, pdf_height
) {

  .isDataFrame(data, "data")
  .hasAtLeastTwoRows(data)
  .isSeqCol(data, seq_col)
  .isDataColOrNull(data, count_col, "count_col")
  .isDataColsOrNull(data, subset_cols, "subset_cols")
  .isTF(drop_isolated_nodes, "drop_isolated_nodes")
  .isTF(node_stats, "node_stats")
  .isTF(cluster_stats, "cluster_stats")
  .isTF(plots, "plots")
  .isTF(print_plots, "print_plots")
  .isString(output_name, "output_name")
  .orNull(.isString, plot_title, "plot_title")
  .orNull(.isString, plot_subtitle, "plot_subtitle")
  .orNull(.isString, drop_matches, "drop_matches")
  .orNull(.isPosInt, min_seq_length, "min_seq_length")
  .isNonneg(dist_cutoff, "dist_cutoff")
  .isPos(pdf_width, "pdf_width")
  .isPos(pdf_height, "pdf_height")
  .isDistType(dist_type)
  .orNull(.isString, output_dir, "output_dir")
  if (!is.null(output_dir)) {
    .isOutputType(output_type)
  }
  # color_nodes_by must be NULL, "auto", or valid column reference
  .checkColorNodesBy(color_nodes_by, data, node_stats, plots)
  .checkStatsToInclude(stats_to_include)
  .checkClusterFun(cluster_fun)

}



# Public Clusters ---------------------------------------------------------


.checkargs.findPublicClusters <- function(
    file_list, input_type, data_symbols, header, sep, sample_ids, seq_col,
    count_col, min_seq_length, drop_matches, top_n_clusters, min_node_count,
    min_clone_count, plots, print_plots, plot_title, color_nodes_by, output_dir,
    output_type, output_dir_unfiltered, output_type_unfiltered
) {

  .checkargs.InputFiles(
    file_list, input_type, data_symbols, header, sep
  )
  .isCharOrNumericVector(sample_ids, "sample_ids")
  .isCharOrNumericVector(seq_col, "seq_col")
  .orNull(.isCharOrNumericScalar, count_col, "count_col")
  .orNull(.isNonneg, min_seq_length, "min_seq_length")
  .orNull(.isString, drop_matches, "drop_matches")
  .isPosInt(top_n_clusters, "top_n_clusters")
  .isPosInt(min_node_count, "min_node_count")
  .orNull(.isPos, min_clone_count, "min_clone_count")
  .isTF(plots, "plots")
  .isTF(print_plots, "print_plots")
  .orNull(.isString, plot_title, "plot_title")
  .orNull(.isCharOrNumericVector, color_nodes_by, "color_nodes_by")
  .orNull(.isString, output_dir, "output_dir")
  .isOutputType(output_type, type = "findPublicClusters")
  .orNull(.isString, output_dir_unfiltered, "output_dir")
  if (!is.null(output_dir_unfiltered)) {
    .isOutputType(output_type_unfiltered)
  }
  if (all(plots, !print_plots, is.null(output_dir_unfiltered))) {
    warning("ignoring `plots = TRUE` since `print_plots = FALSE` and `output_dir_unfiltered` is NULL")
  }
  stopifnot("lengths of file_list and sample_ids must match" =
              length(file_list) == length(sample_ids))

}

.checkargs.buildPublicClusterNetwork <- function(
    file_list, input_type, data_symbols, header, sep, seq_col,
    drop_isolated_nodes,
    color_nodes_by, color_scheme, plot_title, output_dir, output_name
) {

  .checkargs.InputFiles(
    file_list, input_type, data_symbols, header, sep
  )
  .isCharOrNumericVector(seq_col, "seq_col")
  .isTF(drop_isolated_nodes, "drop_isolated_nodes")
  .isCharOrNumericVector(color_nodes_by, "color_nodes_by")
  .isCharVector(color_scheme, "color_scheme")
  .orNull(.isString, plot_title, "plot_title")
  .orNull(.isString, output_dir, "output_dir")
  .isString(output_name, "output_name")

}

.checkargs.buildPublicClusterNetworkByRep <- function(
    file_list, input_type, data_symbols, header, sep, seq_col, count_col,
    dist_type, dist_cutoff, cluster_fun, plots, print_plots, plot_title,
    plot_subtitle, color_nodes_by, color_scheme, output_dir, output_type,
    output_name, pdf_width, pdf_height
) {

  .checkargs.InputFiles(
    file_list, input_type, data_symbols, header, sep
  )
  .isCharOrNumericVector(seq_col, "seq_col")
  .orNull(.isCharOrNumericScalar, count_col, "count_col")
  .isDistType(dist_type)
  .isNonneg(dist_cutoff, "dist_cutoff")
  .checkClusterFun(cluster_fun)
  .isTF(plots, "plots")
  .isTF(print_plots, "print_plots")
  .orNull(.isString, plot_title, "plot_title")
  .orNull(.isString, plot_subtitle, "plot_subtitle")
  .isCharOrNumericVector(color_nodes_by, "color_nodes_by")
  .checkColorScheme(color_scheme, color_nodes_by, plots)
  .orNull(.isString, output_dir, "output_dir")
  if (!is.null(output_dir)) {
    .isOutputType(output_type)
  }
  .isString(output_name, "output_name")
  .isPos(pdf_width, "pdf_width")
  .isPos(pdf_height, "pdf_height")

}


# Associated Clusters -----------------------------------------------------


.checkargs.findAssociatedSeqs <- function(
    file_list, input_type, data_symbols, header, sep,
    subject_ids, group_ids, seq_col, freq_col,
    min_seq_length, drop_matches, min_sample_membership, pval_cutoff, outfile
) {
  .checkargs.InputFiles(
    file_list, input_type, data_symbols, header, sep
  )
  .orNull(.isCharOrNumericVector, subject_ids, "subject_ids")
  .isCharOrNumericVector(group_ids, "group_ids")
  .isCharOrNumericVector(seq_col, "seq_col")
  .orNull(.isCharOrNumericScalar, freq_col, "freq_col")
  .orNull(.isNonneg, min_seq_length, "min_seq_length")
  .orNull(.isString, drop_matches, "drop_matches")
  .orNull(.isNonneg, min_sample_membership, "min_sample_membership")
  .isPos(pval_cutoff, "pval_cutoff")
  .orNull(.isString, outfile, "outfile")
  if (!is.null(subject_ids)) {
    stopifnot(
      "file_list and subject_ids have non-matching lengths" =
        length(file_list) == length(subject_ids)
    )
  }
  stopifnot(
    "lengths of file_list and group_ids must match" =
      length(file_list) == length(group_ids)
  )
  stopifnot(
    "file_list contains duplicate values" =
      length(file_list) == length(unique(file_list))
  )
}

.checkargs.findAssociatedSeqs2 <- function(
    data, seq_col, sample_col, subject_col, group_col, groups, freq_col,
    min_seq_length, drop_matches, min_sample_membership, pval_cutoff, outfile
) {

  if (!is.null(groups)) {
    warning(paste(
      "groups argument is deprecated; group labels are now determined from",
      "the unique values of the group_ids argument. Avoid using the groups",
      "argument to avoid errors in future versions of NAIR"
    ))
  }
  .hasAtLeastTwoRows(data)
  .isSeqCol(data, seq_col)
  .isDataCol(data, group_col, "group_col")
  .isDataCol(data, sample_col, "sample_col")
  .isDataCol(data, subject_col, "subject_col")
  .isDataColOrNull(data, freq_col, "freq_col")
  .orNull(.isNonneg, min_seq_length, "min_seq_length")
  .orNull(.isString, drop_matches, "drop_matches")
  .orNull(.isNonneg, min_sample_membership, "min_sample_membership")
  .isPos(pval_cutoff, "pval_cutoff")
  .orNull(.isString, outfile, "outfile")
  if (!is.null(groups)) {
    warning("`groups` argument is deprecated; group labels are now determined from the unique values of the `group_ids` argument. Avoid using the `groups` argument to avoid errors in future versions of NAIR")
  }
}

.checkargs.findAssociatedClones <- function(
    file_list, input_type, data_symbols, header, sep,
    sample_ids, subject_ids, group_ids, seq_col,
    assoc_seqs, nbd_radius, dist_type, min_seq_length, drop_matches,
    subset_cols, output_dir, output_type, verbose
) {

  .checkargs.InputFiles(
    file_list, input_type, data_symbols, header, sep
  )
  .orNull(.isCharOrNumericVector, subject_ids, "subject_ids")
  .isCharOrNumericVector(group_ids, "group_ids")
  .isCharOrNumericVector(sample_ids, "sample_ids")
  .isCharOrNumericVector(seq_col, "seq_col")
  .isCharOrNumericVector(assoc_seqs, "assoc_seqs")
  .isNonneg(nbd_radius, "nbd_radius")
  .isDistType(dist_type)
  .orNull(.isNonneg, min_seq_length, "min_seq_length")
  .orNull(.isString, drop_matches, "drop_matches")
  .orNull(.isCharOrNumericVector, subset_cols, "subset_cols")
  .isString(output_dir, "output_dir")
  .isOutputType(output_type, "findAssociatedClones")
  .isTF(verbose, "verbose")

}

.checkargs.buildAssociatedClusterNetwork <- function(
    file_list, input_type, data_symbols, header, sep,
    seq_col, min_seq_length, drop_matches, drop_isolated_nodes,
    node_stats, stats_to_include, cluster_stats, color_nodes_by, output_name
) {

  .checkargs.InputFiles(
    file_list, input_type, data_symbols, header, sep
  )
  .isCharOrNumericVector(seq_col, "seq_col")
  .orNull(.isNonneg, min_seq_length, "min_seq_length")
  .orNull(.isString, drop_matches, "drop_matches")
  .isTF(drop_isolated_nodes, "drop_isolated_nodes")
  .isTF(node_stats, "node_stats")
  .isTF(cluster_stats, "cluster_stats")
  .checkStatsToInclude(stats_to_include)
  .orNull(.isCharOrNumericVector, color_nodes_by, "color_nodes_by")
  .isString(output_name, "output_name")

}


# Plotting Functions ------------------------------------------------------

.checkargs.generateNetworkGraphPlots <- function(
    igraph, data, print_plots, plot_title, plot_subtitle, color_nodes_by,
    color_scheme, color_legend, color_title, edge_width, size_nodes_by,
    node_size_limits, size_title
) {

  .isIgraph(igraph, "igraph")
  .isDataFrame(data, "data")
  .checkIgraphAgainstData(igraph, data)
  .isTF(print_plots, "print_plots")
  .orNull(.isString, plot_title, "plot_title")
  .orNull(.isString, plot_subtitle, "plot_subtitle")
  .checkColorNodesBy(color_nodes_by, data)
  .checkColorScheme(color_scheme, color_nodes_by)
  .isTFOrAuto(color_legend, "color_legend")
  .checkColorTitle(color_title, color_nodes_by)
  .isPos(edge_width, "edge_width")
  .checkSizeNodesBy(size_nodes_by, data)
  .checkNodeSizeLimits(node_size_limits)
  .orNull(.isString, size_title, "size_title")

}

# Mid-level Functions -----------------------------------------------------

.checkargs.generateNetworkObjects <- function(
    data, seq_col, dist_type, dist_cutoff, drop_isolated_nodes
) {

  .isDataFrame(data, "data")
  .hasAtLeastTwoRows(data)
  .isSeqCol(data, seq_col)
  .isDistType(dist_type)
  .isNonneg(dist_cutoff, "dist_cutoff")
  .isTF(drop_isolated_nodes, "drop_isolated_nodes")

}

.checkargs.InputFiles <- function(
    file_list, input_type, data_symbols, header, sep
) {
  .isCharVector(file_list, "file_list")
  .isInputType(input_type)
  if (input_type == "rda") {
    .isString(data_symbols, "data_symbols")
  }
  if (input_type %in% c("csv", "table", "tsv", "txt")) {
    .isTF(header, "header")
    .isString(sep, "sep")
  }
}

.checkargs.CombineSamples <- function(
    file_list, input_type, data_symbols, header, sep,
    seq_col, min_seq_length, drop_matches, subset_cols,
    sample_ids, subject_ids, group_ids
) {
  .checkargs.InputFiles(
    file_list, input_type, data_symbols, header, sep
  )
  .isCharOrNumericVector(seq_col, "seq_col")
  .orNull(.isNonneg, min_seq_length, "min_seq_length")
  .orNull(.isString, drop_matches, "drop_matches")
  .orNull(.isCharOrNumericVector, subset_cols, "subset_cols")
  .orNull(.isCharOrNumericVector, sample_ids, "sample_ids")
  .orNull(.isCharOrNumericVector, subject_ids, "subject_ids")
  .orNull(.isCharOrNumericVector, group_ids, "group_ids")
}

.checkargs.filterInputData <- function(
    data, seq_col, min_seq_length, drop_matches, subset_cols, count_col
) {

  .isDataFrame(data, "data")
  .isSeqCol(data, seq_col)
  .orNull(.isNonneg, min_seq_length, "min_seq_length")
  .orNull(.isString, drop_matches, "drop_matches")
  .isDataColsOrNull(data, subset_cols, "subset_cols")
  .isDataColOrNull(data, count_col, "count_col")

}

.checkargs.addClusterLabels <- function(
    plot, net, top_n_clusters, cluster_id_col, criterion,
    size, color, greatest_values
) {

  .isGgraph(plot, "plot")
  .hasNodeAndClusterData(net, "net")
  .isPosInt(top_n_clusters, "top_n_clusters")
  .isDataCol(net$node_data, cluster_id_col, "cluster_id_col")
  .isDataCol(net$cluster_data, criterion, "criterion")
  .isPos(size, "size")
  .isString(color, "color")
  .isTF(greatest_values, "greatest_values")

}

.checkargs.saveNetwork <- function(
    net, output_dir, output_type, output_filename, pdf_width, pdf_height
) {

  .isBaseNetworkOutput(net, "net")
  .orNull(.isString, output_dir, "output_dir")
  .isOutputType(output_type)
  .isString(output_filename, "output_filename")
  .isPos(pdf_width, "pdf_width")
  .isPos(pdf_height, "pdf_height")

}

.checkargs.saveNetworkPlots <- function(
    plotlist, outfile, pdf_width, pdf_height
) {
  .isPlotlist(plotlist, "plotlist")
  .isString(outfile, "outfile")
  .isPos(pdf_width, "pdf_width")
  .isPos(pdf_height, "pdf_height")
}


.checkargs.chooseNodeStats <- function(
    degree, cluster_id, transitivity, closeness, centrality_by_closeness,
    eigen_centrality, centrality_by_eigen, betweenness,
    centrality_by_betweenness, authority_score, coreness, page_rank, all_stats
) {
  .isTF(degree, "degree")
  .isTF(cluster_id, "cluster_id")
  .isTF(transitivity, "transitivity")
  .isTF(closeness, "closeness")
  .isTF(centrality_by_closeness, "centrality_by_closeness")
  .isTF(eigen_centrality, "eigen_centrality")
  .isTF(centrality_by_eigen, "centrality_by_eigen")
  .isTF(betweenness, "betweenness")
  .isTF(centrality_by_betweenness, "centrality_by_betweenness")
  .isTF(authority_score, "authority_score")
  .isTF(coreness, "coreness")
  .isTF(page_rank, "page_rank")
  .isTF(all_stats, "all_stats")
}

.checkargs.addNodeNetworkStats <- function(
  data, net, stats_to_include, cluster_fun
) {
  .isDataFrame(data, "data")
  .isIgraph(net, "net")
  .checkIgraphAgainstData(net, data)
  .checkStatsToInclude(stats_to_include)
  .checkClusterFun(cluster_fun)
}
.checkargs.addClusterMembership <- function(data, net, fun) {
  .isDataFrame(data, "data")
  .isIgraph(net, "net")
  .checkIgraphAgainstData(net, data)
  .checkClusterFun(fun, "fun")
}

.checkargs.getClusterStats <- function(
    data, adjacency_matrix, seq_col, count_col, cluster_id_col,
    degree_col, cluster_fun
) {
  .isDataFrame(data, "data")
  .isAdjacencyMatrix(adjacency_matrix, "adjacency_matrix")
  .checkDataAgainstMatrix(data, adjacency_matrix)
  .isSeqCol(data, seq_col)
  .isDataColOrNull(data, count_col, "count_col")
  .isDataColOrNull(data, cluster_id_col, "cluster_id_col")
  .isDataColOrNull(data, degree_col, "degree_col")
  .checkClusterFun(cluster_fun)
}

.checkargs.sparseAdjacencyMatFromSeqs <- function(
    seqs, dist_type, max_dist, drop_isolated_nodes
) {
  .isValidSeqVector(seqs)
  .isDistType(dist_type)
  .isNonneg(max_dist, "max_dist")
  .isTF(drop_isolated_nodes, "drop_isolated_nodes")
}


# Deprecated Arguments ----------------------------------------------------
.checkDeprecated.buildPublicClusterNetwork <- function(
    node_stats, stats_to_include, cluster_stats,
    env = caller_env(),
    user_env = caller_env(2)
) {
  if (lifecycle::is_present(node_stats)) {
    lifecycle::deprecate_warn(
      when = "0.0.9038",
      what = "buildPublicClusterNetwork(node_stats)",
      details =
        "all node-level network properties are now automatically computed",
      env = env,
      user_env = user_env
    )
  }
  if (lifecycle::is_present(stats_to_include)) {
    lifecycle::deprecate_warn(
      when = "0.0.9038",
      what = "buildPublicClusterNetwork(stats_to_include)",
      details =
        "all node-level network properties are now automatically computed",
      env = env,
      user_env = user_env
    )
  }
  if (lifecycle::is_present(cluster_stats)) {
    lifecycle::deprecate_warn(
      when = "0.0.9038",
      what = "buildPublicClusterNetwork(cluster_stats)",
      details =
        "all cluster-level network properties are now automatically computed",
      env = env,
      user_env = user_env
    )
  }
}

.checkDeprecated.findAssociatedSeqs <- function(
    sample_ids, groups,
    env = caller_env(),
    user_env = caller_env(2)
) {
  if (lifecycle::is_present(groups)) {
    lifecycle::deprecate_warn(
      when = "0.0.9038",
      what = "findAssociatedSeqs(groups)",
      details =
        "group labels are now determined from the unique values of group_ids",
      env = env,
      user_env = user_env
    )
  }
  if (lifecycle::is_present(sample_ids)) {
    lifecycle::deprecate_warn(
      when = "0.0.9038",
      what = "findAssociatedSeqs(sample_ids)",
      details =
        "custom sample IDs are not relevant to this function",
      env = env,
      user_env = user_env
    )
  }
}