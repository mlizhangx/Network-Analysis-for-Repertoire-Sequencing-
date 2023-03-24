
# Top-Level Checking Functions --------------------------------------------

.checkArgs.buildRepSeqNetwork <- function(
  data, seq_col, count_col, subset_cols, min_seq_length, drop_matches,
  dist_type, dist_cutoff, drop_isolated_nodes, node_stats, stats_to_include,
  cluster_stats, plots, print_plots, plot_title,
  plot_subtitle, color_nodes_by, output_type, output_name,
  pdf_width, pdf_height) {

  .isRepSeqData(data)
  .isSeqCol(data, seq_col)
  .isDataColOrNull(data, count_col, "count_col")
  .isDataColsOrNull(data, subset_cols, "subset_cols")

  .isTF(drop_isolated_nodes, "drop_isolated_nodes")
  .isTF(node_stats, "node_stats")
  .isTF(cluster_stats, "cluster_stats")
  .isTF(plots, "plots")
  .isTF(print_plots, "print_plots")

  .isString(output_name, "output_name")
  .isStringOrNull(plot_title, "plot_title")
  .isStringOrNull(plot_subtitle, "plot_subtitle")
  .isStringExprOrNull(drop_matches, "drop_matches")

  .isPosIntOrNull(min_seq_length, "min_seq_length")

  .isPos(dist_cutoff, "dist_cutoff")
  .isPos(pdf_width, "pdf_width")
  .isPos(pdf_height, "pdf_height")

  .isDistType(dist_type)

  .isOutputType(output_type)

  # color_nodes_by must be NULL, "auto", or valid column reference
  .checkColorNodesBy(color_nodes_by, data, node_stats, cluster_stats)

  .checkStatsToInclude(stats_to_include)

}