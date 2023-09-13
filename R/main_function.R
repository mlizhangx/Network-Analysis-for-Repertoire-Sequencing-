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

buildRepSeqNetwork <- function(
    data,
    seq_col,
    count_col = NULL,
    subset_cols = NULL,
    min_seq_length = 3,
    drop_matches = NULL,
    dist_type = "hamming",
    dist_cutoff = 1,
    drop_isolated_nodes = TRUE,
    node_stats = FALSE,
    stats_to_include = chooseNodeStats(),
    cluster_stats = FALSE,
    cluster_fun = "fast_greedy",
    cluster_id_name = "cluster_id",
    plots = TRUE,
    print_plots = FALSE,
    plot_title = "auto",
    plot_subtitle = "auto",
    color_nodes_by = "auto",
    ...,
    output_dir = NULL,
    output_type = "rds",
    output_name = "MyRepSeqNetwork",
    pdf_width = 12,
    pdf_height = 10,
    verbose = FALSE
) {

  # Process arguments
  data_name <- deparse(substitute(data))
  data <- as.data.frame(data)
  .MUST.isDataFrame(data, data_name)
  .MUST.hasMultipleRows(data, data_name)
  .MUST.isSeqColrefs(seq_col, data, deparse(substitute(seq_col)), data_name)
  count_col <- .checkCountCol(count_col, data, NULL)
  subset_cols <- .checkDataColrefs(subset_cols, data, NULL)
  min_seq_length <- .check(min_seq_length, .isNonneg, NULL, ornull = TRUE)
  drop_matches <- .check(drop_matches, .isString, NULL, ornull = TRUE)
  dist_type <- .checkDistType(dist_type, "hamming")
  dist_cutoff <- .check(dist_cutoff, .isNonneg, 1)
  drop_isolated_nodes <- .checkTF(drop_isolated_nodes, TRUE)
  node_stats <- .checkTF(node_stats, FALSE)
  if (isTRUE(node_stats)) {
    stats_to_include <- .checkStatsToInclude(stats_to_include)
  }
  cluster_stats <- .checkTF(cluster_stats, FALSE)
  cluster_fun <- .check(cluster_fun, .isClusterFun, "fast_greedy")
  cluster_id_name <- make.names(
    .check(cluster_id_name, .isString, "cluster_id")
  )
  output_dir <- .check(output_dir, .isString, NULL, ornull = TRUE)
  .createOutputDir(output_dir)
  output_dir <- .checkOutputDir(output_dir)
  plots <- .checkTF(plots, TRUE)
  if (isTRUE(plots)) {
    print_plots <- .checkTF(print_plots, TRUE)
    plot_title <- .check(plot_title, .isString, "auto", ornull = TRUE)
    plot_subtitle <- .check(plot_subtitle, .isString, "auto", ornull = TRUE)
    color_nodes_by <- .checkColorNodesBy(
      color_nodes_by, data, node_stats, cluster_stats, plots, cluster_id_name,
      stats_to_include, default = "auto", auto_ok = TRUE
    )
    if (!is.null(output_dir)) {
      pdf_width <- .check(pdf_width, .isPos, 12)
      pdf_height <- .check(pdf_width, .isPos, 10)
    }
  }
  if (isTRUE(plots) || !is.null(output_dir)) {
    output_name <- .checkOutputName(output_name, "MyRepSeqNetwork")
  }
  if (!is.null(output_dir)) {
    output_type <- .check(output_type, .isOutputType, "rds")
  }
  seq_col <- .convertColRef(seq_col, data)
  count_col <- .convertColRef(count_col, data)
  color_nodes_by <- .convertColRef(color_nodes_by, data)
  subset_cols <- .convertColRef(subset_cols, data)
  subset_cols <- .processSubsetCols(subset_cols, c(count_col, color_nodes_by))
  data <- filterInputData(data, seq_col, min_seq_length, drop_matches,
                          subset_cols, verbose = verbose
  )
  if (nrow(data) < 2) {
    warning(
      "Returning NULL since fewer than two observations remain after ",
      "filtering the data"
    )
    return(NULL)
  }

  # Build network
  net <- generateNetworkObjects(data, seq_col, dist_type, dist_cutoff,
                                drop_isolated_nodes, verbose
  )
  if (is.null(net)) {
    warning("Graph contains no nodes; returning NULL. ",
            "Consider increasing ", sQuote("dist_cutoff")
    )
    return(NULL)
  }
  if (node_stats) {
    net <- addNodeStats(net, stats_to_include, cluster_fun, cluster_id_name,
                        verbose
    )
  }
  if (cluster_stats) {
    net <- addClusterStats(net, cluster_id_name = cluster_id_name,
                           seq_col = seq_col, count_col = count_col,
                           cluster_fun = cluster_fun, verbose = verbose
    )
  }
  if (plots) {
    net <- addPlots(
      net, print_plots,
      .makePlotTitle(plot_title, network_name = output_name),
      .makePlotSubtitle(plot_subtitle, seq_col = seq_col,
                        dist_type = dist_type, dist_cutoff = dist_cutoff
      ),
      .passColorNodesBy(color_nodes_by, net$node_data, count_col),
      verbose = verbose,
      ...
    )
  }
  net$details$min_seq_length <- ifelse(
    is.null(min_seq_length), yes = "NULL", no = min_seq_length
  )
  net$details$drop_matches <- ifelse(
    is.null(drop_matches), yes = "NULL", no = drop_matches
  )
  saveNetwork(net, output_dir, output_type, output_name, pdf_width, pdf_height,
              verbose
  )
  invisible(net)
}

buildNet <- buildRepSeqNetwork