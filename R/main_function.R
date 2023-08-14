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
    cluster_fun = cluster_fast_greedy,
    plots = TRUE,
    print_plots = TRUE,
    plot_title = "auto",
    plot_subtitle = "auto",
    color_nodes_by = "auto",
    ...,
    output_dir = getwd(),
    output_type = "individual",
    output_name = "MyRepSeqNetwork",
    pdf_width = 12,
    pdf_height = 10
) {

  data <- as.data.frame(data)
  .checkArgs.buildRepSeqNetwork(
    data, seq_col, count_col, subset_cols, min_seq_length, drop_matches,
    dist_type, dist_cutoff, drop_isolated_nodes, node_stats, stats_to_include,
    cluster_stats, cluster_fun, plots, print_plots, plot_title, plot_subtitle,
    color_nodes_by, output_dir, output_type, output_name, pdf_width, pdf_height
  )
  seq_col <- .convertColRef(seq_col, data)
  count_col <- .convertColRef(count_col, data)
  color_nodes_by <- .convertColRef(color_nodes_by, data)
  subset_cols <- .convertColRef(subset_cols, data)
  subset_cols <- .processSubsetCols(subset_cols, c(count_col, color_nodes_by))
  data <- filterInputData(data, seq_col, min_seq_length, drop_matches,
                          subset_cols, count_col
  )
  if (nrow(data) < 2) {
    warning(paste(
      "Returning NULL, since fewer than two observations remain after",
      "filtering the data. Check the values used for the min_seq_length",
      "and drop_matches arguments"
    ))
    return(NULL)
  }
  .createOutputDir(output_dir)
  out <- generateNetworkObjects(data, seq_col, dist_type, dist_cutoff,
                                drop_isolated_nodes
  )
  if (is.null(out)) {
    warning(paste(
      "Graph contains no nodes; returning NULL. You may wish to consider",
      "another distance type or a greater distance cutoff, or setting",
      "drop_isolated_nodes = FALSE"
    ))
    return(NULL)
  }
  if (node_stats) {
    out$node_data <- addNodeNetworkStats(out$node_data, out$igraph,
                                         stats_to_include, cluster_fun
    )
  }
  if (cluster_stats) {
    out$node_data <- .addClusterAndDegree(out, cluster_fun)
    out$cluster_data <- .getClusterStats2(out$node_data, out$adjacency_matrix,
                                          seq_col, count_col, cluster_fun
    )
  }
  if (plots) {
    out$plots <- .generateNetworkGraphPlotsGuarded(
      out$igraph, out$node_data, print_plots,
      .makePlotTitle(plot_title, network_name = output_name),
      .makePlotSubtitle(plot_subtitle, seq_col = seq_col,
                        dist_type = dist_type, dist_cutoff = dist_cutoff
      ),
      .passColorNodesBy(color_nodes_by, out$node_data, count_col),
      ...
    )
  }
  saveNetwork(out, output_dir, output_type, output_name, pdf_width, pdf_height)
  invisible(out)
}

.addClusterAndDegree <- function(network, cluster_fun) {
  if (!"cluster_id" %in% names(network$node_data)) {
    network$node_data <-
      addClusterMembership(network$node_data, network$igraph, cluster_fun)
  }
  if (!"degree" %in% names(network$node_data)) {
    network$node_data$deg <- igraph::degree(network$igraph)
  }
  return(network$node_data)
}