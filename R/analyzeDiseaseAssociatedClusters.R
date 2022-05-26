# Inputs:
#       Combined RepSeq Data From Top Disease-Associated Clusters
# Do:
#       Build Network
#       Compute Node-level Network Statistics
#       Plot Network Graphs
#       Perform Atchley Factor Embedding of Clone Seqs and Plot Corr Heatmap

analyzeDiseaseAssociatedClusters <- function(
  data, # combined rep-seq data from disease-associated clusters to be included
  clone_col, # column of data containing clone sequences
  size_col, # column of data to use for sizing the nodes (e.g., clone count, clone fraction)
  color_cols, # vector specifying columns of data to use for node coloring (one plot per variable)
  dist_type = "hamming",
  edge_dist = 1) {

  # Build network for combined clusters
  net <- generateNetworkFromClones(data[ , clone_col],
                                   dist_type = dist_type, edge_dist = edge_dist)
  # Compute node-level network characteristics and cluster membership
  data <- computeNodeNetworkStats(net, data,
                                  node_stat_settings(cluster_id = TRUE))
  ## Create plots ##
  plots <- list()
  color_indices <- c("D", "C", "E", "H", "B", "G", "F" ,"A") # viridis colors
  if (is.numeric(size_col)) {
    size_label <- names(data)[[size_label]]
  } else { size_label <- size_col }
  if (is.numeric(color_cols)) {
    color_labels <- names(data)[color_cols]
  } else { color_labels <- color_cols }
  for (i in 1:length(color_cols)) {
    newplot <-
      plotNetworkGraph(
        net,
        title = "Clusters For Top Disease-Associated Clones",
        subtitle = paste0("Nodes colored by ", color_labels[[i]]),
        color_nodes_by = as.factor(data[ , color_cols[[i]]]),
        size_nodes_by = as.numeric(data[ , size_col]),
        color_legend_title = color_labels[[i]],
        size_legend_title = size_label)
    color_code_type <-
      ggplot2::scale_type(as.factor(data[ , color_cols[[i]]]))[[1]]
    if (color_code_type == "continuous") {
      newplot <- newplot +
        ggraph::scale_color_viridis(
          begin = 0.1, end = 0.85, direction = -1,
          option = color_indices[[i]]) }
    if (color_code_type == "discrete") {
      newplot <- newplot +
        ggraph::scale_color_viridis(
          begin = 0.1, end = 0.85, direction = -1,
          option = color_indices[[i]], discrete = TRUE) }
    newplot # print to R
    plots <- c(plots, newplot)
    names(plots)[[length(names(plots))]] <- color_labels[[i]]
  }

  # Perform Atchley Factor Encoding

  # Correlation Plot on Embedded Values

  # Save results

  # Return results
  # return(list("node_data" = data,
  #             "igraph" = net,
  #             "plots" = plots,
  #             "AF_embed" = AF_embedded_values,
  #             "AF_corr" = AF_corr_heatmap))

}