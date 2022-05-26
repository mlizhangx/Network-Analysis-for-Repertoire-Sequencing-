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
  size_col, # column of data to use for sizing the nodes (e.g. clone fraction)
  color_cols = NULL, # character vec specifying columns of data to use for node coloring (one plot per variable)
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
  # Plot network graph with nodes colored by disease status
  plots$disease <- plotNetworkGraph(
    net,
    title = "Clusters For Top Disease-Associated Clones",
    subtitle = paste0("Nodes colored by disease status"),
    color_nodes_by = as.factor(data[ , disease_col]),
    size_nodes_by = as.numeric(data[ , clone_frac_col]),
    color_legend_title = "Disease Status",
    size_legend_title = "Clone Fraction")
  plots$disease
  # Plot network graph with nodes colored by disease-associated cluster
  plots$cluster <- plotNetworkGraph(
    net,
    title = "Clusters For Top Disease-Associated Clones",
    subtitle = paste0("Nodes colored by disease-associated cluster"),
    color_nodes_by = as.factor(data[ , disease_col]),
    size_nodes_by = as.numeric(data[ , clone_frac_col]),
    color_legend_title = "Cluster",
    size_legend_title = "Clone Fraction") +
    ggraph::scale_color_viridis(
      begin = 0, end = 1, direction = -1, option = "turbo", discrete = TRUE)
  plots$cluster
  if (!is.null(color_cols)) {
    color_indices <- c("D", "C", "E", "H", "B", "G", "F" ,"A")
    for (i in 1:length(color_cols)) {
      newplot <-
        plotNetworkGraph(
          net,
          title = "Clusters For Top Disease-Associated Clones",
          subtitle = paste0("Nodes colored by ", color_cols[[i]]),
          color_nodes_by = as.factor(data[ , color_cols[[i]]]),
          size_nodes_by = as.numeric(data[ , clone_frac_col]),
          color_legend_title = color_cols[[i]],
          size_legend_title = "Clone Fraction")
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
      newplot
      plots <- c(plots, newplot)
      names(plots)[[length(names(plots))]] <- color_cols[[i]]
    }
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