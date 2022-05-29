# Inputs:
#       Combined RepSeq Data From Top Disease-Associated Clusters
# Do:
#       Build Network
#       Compute Node-level Network Statistics
#       Plot Network Graphs
#       Perform Atchley Factor Embedding of Clone Seqs and Plot Corr Heatmap

combineClustersIntoNetwork <- function(

  data, # combined rep-seq data from disease-associated clusters to be included
  clone_col, # column of data containing clone sequences
  size_col, # column of data to use for sizing the nodes (e.g., clone count, clone fraction)
  color_cols, # vector specifying columns of data to use for node coloring (one plot per variable)
  dist_type = "hamming",
  edge_dist = 1,
  output_dir = NULL

) {
  # Check that sample_id_col is valid (or that include_atchley = FALSE)
  # TO DO: Add check

  # Build network for combined clusters
  net <- generateNetworkFromClones(data[ , clone_col],
                                   dist_type = dist_type, edge_dist = edge_dist)
  # Compute node-level network characteristics and cluster membership
  data <- computeNodeNetworkStats(net, data,
                                  node_stat_settings(cluster_id = TRUE))
  ## Create plots ##
  plots <- list()
  # color_indices <- c("D", "C", "E", "H", "B", "G", "F" ,"A") # viridis colors
  if (is.numeric(size_col)) {
    size_label <- names(data)[[size_label]]
  } else { size_label <- size_col }
  if (is.numeric(color_cols)) {
    color_labels <- names(data)[color_cols]
  } else { color_labels <- color_cols }
  for (i in 1:length(color_cols)) {
    discrete_colors <- FALSE
    color_code_type <-
      ggplot2::scale_type(as.factor(data[ , color_cols[[i]]]))[[1]]
    if (color_code_type == "discrete") { discrete_colors <- TRUE }
    newplot <-
      plotNetworkGraph(
        net, edge_width = 0.1,
        title = "Clusters For Top Disease-Associated Clones",
        subtitle = paste0("Nodes colored by ", color_labels[[i]]),
        color_nodes_by = as.factor(data[ , color_cols[[i]]]),
        size_nodes_by = as.numeric(data[ , size_col]),
        color_legend_title = color_labels[[i]],
        size_legend_title = size_label) +
      ggraph::scale_color_viridis(
        begin = 0.1, end = 0.85, direction = -1,
        option = "plasma", discrete = discrete_colors)
    print(newplot) # print to R
    plots$newplot <- newplot
    names(plots)[[length(names(plots))]] <- color_labels[[i]]
  }

  # Save results
  if (!is.null(output_dir)) {
    # Save Network igraph using edgelist format
    igraph::write_graph(net,
                        file = file.path(output_dir, "network_graph_edgelist.txt"),
                        format = "edgelist")
    cat("Network igraph object saved in edgelist format as 'network_graph_edgelist.txt'\n")

    # Save node-level data with network characteristics
    utils::write.csv(data,
                     file = file.path(output_dir, "node_data.csv"),
                     row.names = FALSE)
    cat("Node-level data and network characteristics saved as 'node_data.csv'\n")

    # Save network graph plots
    grDevices::pdf(file.path(output_dir, "network_graph_plots.pdf"),
                   width = 12, height = 8)
    for (plot in plots) { print(plot) }
    grDevices::dev.off()
    cat("Plots of network graph saved as 'network_graph_plots.pdf'\n")
  }

  # Return results
  return(list("node_data" = data,
              "igraph" = net,
              "plots" = plots))

}



