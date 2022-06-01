# Inputs:
#       Combined TCR rep-seq data from multiple samples
#       List of selected clone sequences (CDR3 nucleotide or amino acid)
# Do:
#       Build a cluster around each of the selected clones
#       Print the plot for each cluster and return the combined cluster data

buildClustersAroundSelectedClones <- function(

  # Input Data and Columns
  data, # merged bulk rep-seq data from all patients/samples
  nucleo_col,
  amino_col,
  count_col,
  freq_col,
  vgene_col,
  dgene_col,
  jgene_col,
  cdr3length_col,
  sample_col,
  other_cols = NULL,

  selected_clones,
  selected_clone_labels = NULL,

  clone_seq_type = "amino acid",
  drop_chars = "[*|_]", # passed to getSimilarClones() for getting cluster data

  # Network Settings
  dist_type = "hamming", # options are "hamming", "levenshtein", "euclidean_on_atchley"
  cluster_radius = 1,
  edge_dist = 1,

  # Plot Settings
  # main_title = paste("Global cluster network"),
  # main_subtitle = ifelse(dist_type == "euclidean_on_atchley",
  #                        yes = paste("Clone sequences embedded in Euclidean 30-space based on Atchley factor representation using deep learning\nEdges based on a maximum Euclidean distance of", edge_dist, "between embedded values\n"),
  #                        no = paste("Edges based on a maximum", dist_type, "distance of", edge_dist, "\n")),
  size_nodes_by = count_col, # can use a double, e.g., 1.0, for fixed size
  node_size_limits = NULL, # numeric, length 2
  custom_size_legend = NULL, # custom legend title
  edge_width = 0.3,
  color_nodes_by = sample_col, # accepts multiple values (one plot per value)
  color_scheme = "default", # passed to plotNetworkGraph(); accepts multiple values (one per value of color_nodes_by)
  custom_color_legend = NULL, # custom title (length must match color_nodes_by)
  long_captions = TRUE, # should plot subtitles include details on cluster settings like dist_type and edge_dist?

  # Output Settings
  output_dir = getwd(),
  data_outfile = "data_all_clusters.csv",
  save_plots = FALSE,
  single_plot = TRUE, # if false, when saving plots, save one file per cluster
  plot_width = 12, # passed to pdf()
  plot_height = 10, # passed to pdf()
  return_plots = FALSE # should function return a list of dataframe + ggplots, or just print/write plots and return the dataframe?

) {

  ### INPUT CHECKS ###
  # Atchley factor embedding only applicable to amino acid sequences
  if (dist_type == "euclidean_on_atchley" & clone_seq_type != "amino acid") {
    stop("distance type 'euclidean_on_atchley' only applicable to amino acid sequences") }


  #### PREPARE WORKING ENVIRONMENT ####
  # Create output directory if applicable
  if (!is.null(output_dir)) {
    .createOutputDir(output_dir)
    if (save_plots & !single_plot) {
      .createOutputDir(file.path(output_dir, "individual_cluster_plots"))
    }}

  # Convert input columns to character if not already
  if (is.numeric(color_nodes_by)) { color_nodes_by <- names(data)[color_nodes_by] }
  if (is.numeric(nucleo_col)) { nucleo_col <- names(data)[nucleo_col] }
  if (is.numeric(amino_col)) { amino_col <- names(data)[amino_col] }
  if (is.numeric(count_col)) { count_col <- names(data)[count_col] }
  if (is.numeric(freq_col)) { freq_col <- names(data)[freq_col] }
  if (is.numeric(vgene_col)) { vgene_col <- names(data)[vgene_col] }
  if (is.numeric(dgene_col)) { dgene_col <- names(data)[dgene_col] }
  if (is.numeric(jgene_col)) { jgene_col <- names(data)[jgene_col] }
  if (is.numeric(cdr3length_col)) { cdr3length_col <- names(data)[cdr3length_col] }
  if (is.numeric(sample_col)) { sample_col <- names(data)[sample_col] }
  if (is.numeric(other_cols)) { other_cols <- names(data)[other_cols] }
  if (is.integer(size_nodes_by)) { size_nodes_by <- names(data)[size_nodes_by] }

  # Designate amino acid or nucleotide for clone sequence
  clone_seq_col <- amino_col
  if (clone_seq_type == "nucleotide") { clone_seq_col <- nucleo_col }

  # Determine distance type for cluster radius
  cluster_radius_dist_type <- "levenshtein" # applies to lev and atchley
  if (dist_type == "hamming") { cluster_radius_dist_type <- "hamming" }

  # Format the input data
  data <-
    data[ , # Keep only the relevant columns:
          intersect(
            unique(
              c(nucleo_col, amino_col, count_col, freq_col, vgene_col,
                dgene_col, jgene_col, cdr3length_col, sample_col,
                other_cols, color_nodes_by)),
            names(data))]

  # Initialize output directory and objects
  data_all_clusters <- as.data.frame(matrix(nrow = 0, ncol = ncol(data) + 3))
  names(data_all_clusters) <- c(colnames(data), "assocClustID",
                                "assocClustSeq", "degreeInAssocClust")

  if (return_plots | (save_plots & single_plot & !is.null(output_dir))) {
    plots <- list() }

  ### PLOT SETTINGS ###
  # Fixed subtitle content across selected clone sequences
  if (long_captions) {
    subtitle_extra <- paste0("Cluster includes clone sequences with a maximum ", cluster_radius_dist_type, " distance of ", cluster_radius, " from the central sequence\n")
    if (dist_type == "euclidean_on_atchley") {
      subtitle_extra <- paste0(subtitle_extra, "Clone sequences embedded in Euclidean 30-space based on Atchley factor representation using deep learning\nEdges based on a maximum Euclidean distance of ", edge_dist, " between embedded values\n")
    } else {
      subtitle_extra <- paste0(subtitle_extra, "Edges based on a maximum ", dist_type, " distance of ", edge_dist, "\n") }
  }

  # If multiple coloring variables, extend color scheme to vector if needed
  if (length(color_nodes_by) > 1 & length(color_scheme) == 1) {
    color_scheme <- rep(color_scheme, length(color_nodes_by)) }

  # Legend titles for size and color
  size_legend_title <- NULL # default for fixed node size
  if (is.character(size_nodes_by)) {
    if (!is.null(custom_size_legend)) {
      size_legend_title <- custom_size_legend
    } else { size_legend_title <- size_nodes_by }
  }

  if (!is.null(custom_color_legend)) {
    color_legend_title <- custom_color_legend
  } else { color_legend_title <- color_nodes_by }

  #### BUILD CLUSTERS FOR SELECTED CLONES ####
  # Iterate over the selected clones
  for (i in 1:length(selected_clones)) {

    ### BUILD CLUSTER NETWORK ###
    # Get cluster data
    cat(paste0("Gathering cluster data for clone sequence ", i,
               " (", selected_clones[[i]], "):\n"))
    data_current_cluster <-
      getSimilarClones(selected_clones[[i]], data, clone_seq_col, sample_col,
                       cluster_radius_dist_type, max_dist = cluster_radius,
                       drop_chars = drop_chars)

    # Build cluster network
    # Generate adjacency matrix for network
    adjacency_matrix <-
      generateNetworkFromClones(data_current_cluster[ , clone_seq_col],
                                dist_type, edge_dist,
                                contig_ids = rownames(data_current_cluster),
                                return_type = "adjacency_matrix")

    # Subset data to keep only those clones in the network (nonzero degree)
    if (dist_type != "euclidean_on_atchley") {
      data_current_cluster <-
        data_current_cluster[as.numeric(dimnames(adjacency_matrix)[[1]]), ] }

    # Generate network from adjacency matrix
    network <- generateNetworkFromAdjacencyMat(adjacency_matrix)


    # Add variables for cluster ID, central sequence, and degree in cluster
    data_current_cluster$assocClusterID <- i
    data_current_cluster$assocClusterSeq <- selected_clones[[i]]
    data_current_cluster$degreeInAssocCluster <- igraph::degree(network)

    # Add data for current cluster to combined cluster data
    data_all_clusters <- rbind(data_all_clusters, data_current_cluster)


    ### PLOT(S) OF NETWORK GRAPH ###
    # Create labels for plots
    plot_title <- paste0("Cluster ", i, " (", selected_clones[[i]], ")")
    plot_subtitle <- NULL
    if (!is.null(selected_clone_labels)) {
      plot_subtitle <- paste0(selected_clone_labels[[i]], "\n") }
    if (long_captions) {
      plot_subtitle <- paste0(plot_subtitle, subtitle_extra) }

    # Ensure size_nodes_by is a vector or fixed value to use for node sizes
    if (is.character(size_nodes_by)) {
      size_code_vector <- data_current_cluster[ , size_nodes_by] }

    # Create one plot for each variable in color_nodes_by
    temp_plotlist <- list()
    for (j in 1:length(color_nodes_by)) {
      cat(paste0("Creating cluster graph with nodes colored by ",
                 color_nodes_by[[j]], "..."))
      temp_plotlist$newplot <-
        plotNetworkGraph(
          network, edge_width, title = plot_title,
          subtitle = plot_subtitle,
          color_nodes_by = data_current_cluster[ , color_nodes_by[[j]]],
          size_nodes_by = size_code_vector,
          color_legend_title = color_legend_title[[j]],
          size_legend_title = size_legend_title,
          color_scheme = color_scheme[[j]],
          node_size_limits = node_size_limits)
      print(temp_plotlist$newplot) # print to R
      names(temp_plotlist)[[length(names(temp_plotlist))]] <-
        color_nodes_by[[j]]
      cat(" Done.\n")
    }

    # Save plots for current cluster to a single pdf if applicable
    if (save_plots & !single_plot & !is.null(output_dir)) {
      plot_outfile <- file.path(output_dir, "individual_cluster_plots",
                                paste0("cluster_", i, ".pdf"))
      grDevices::pdf(file = plot_outfile,
                     width = plot_width, height = plot_height)
      for (j in 1:length(color_nodes_by)) { print(temp_plotlist[[j]]) }
      grDevices::dev.off()
      cat(paste0("Cluster graph plot saved to file:\n  ", plot_outfile, "\n"))
    }

    # Add plots to output list if applicable
    if (return_plots | (save_plots & single_plot & !is.null(output_dir))) {
      plots$newcluster <- temp_plotlist
      names(plots)[[length(names(plots))]] <- selected_clones[[i]]
    }

  } # done looping over selected clones
  cat("All clusters complete.\n")

  # Rename some columns of combined cluster data
  colnames(data_all_clusters)[1:9] <- c(
    "nucleotideSeq", "aminoAcidSeq", "cloneCount", "cloneFreqInSample",
    "VGene", "DGene", "JGene", "CDR3Length", "sampleID")

  # Save combined cluster data if applicable
  if (!is.null(output_dir) & !is.null(data_outfile)) {
    utils::write.csv(data_all_clusters, file.path(output_dir, data_outfile),
                     row.names = FALSE)
    cat(paste0(
      "Cluster data saved to file:\n  ",
      file.path(output_dir, data_outfile), "\n")) }

  # Save plots for all clusters to a single pdf if applicable
  if (save_plots & single_plot & !is.null(output_dir)) {
    plot_outfile <- file.path(output_dir, "individual_cluster_plots.pdf")
    grDevices::pdf(file = plot_outfile,
                   width = plot_width, height = plot_height)
    for (i in 1:length(plots)) {
      for (j in 1:length(plots[[i]])) { print(plots[[i]][[j]]) } }
    grDevices::dev.off()
    cat(paste0("Individual cluster graph plots saved to file:\n  ",
               plot_outfile, "\n"))
  }

  # Return output
  if (return_plots) {
    cat("All tasks complete. Returning a list containing the plots and combined cluster data.\n")
    return(list("cluster_data" = data_all_clusters,
                "plots" = plots))
  } else {
    cat("All tasks complete. Returning combined cluster data.\n")
    return(data_all_clusters)
  }
}


# Helpers -----------------------------------------------------------------







# .saveResultsForCandidateSeqNetwork <- function(
#   network, data_current_cluster, adjacency_matrix, netplot_disease,
#   netplot_sampleid, netplot_patid, netplot_deg, keep_adjacency_matrix,
#   dist_type, output_dir, outfilestem) {
#
#   # pdf: network graphs
#   plotfile <- paste0(outfilestem, "_network_plots.pdf")
#   grDevices::pdf(file.path(output_dir, plotfile), width = 12, height = 8)
#   if (!is.null(netplot_deg)) { print(netplot_deg) }
#   if (!is.null(netplot_disease)) { print(netplot_disease) }
#   if (!is.null(netplot_sampleid)) { print(netplot_sampleid) }
#   if (!is.null(netplot_patid)) { print(netplot_patid) }
#   grDevices::dev.off()
#
#   # Save metadata for candidiate sequence network
#   utils::write.csv(
#     data_current_cluster,
#     file.path(output_dir, paste0(outfilestem, "_network_metadata.csv")))
#
#   # Save Network igraphs using edgelist format
#   igraph::write_graph(
#     network,
#     file = file.path(
#       output_dir, paste0(outfilestem, "_network_graph_edgelist.txt")),
#     format = "edgelist")
#
#   # Save adjacency matrices
#   if (keep_adjacency_matrix) {
#     matfile <- file.path(output_dir, paste0(outfilestem, "_adjacency_matrix"))
#     if (dist_type == "euclidean_on_atchley") {
#       utils::write.csv(adjacency_matrix, paste0(matfile, ".csv"))
#     } else { #hamming/levenshtein returns sparse matrix
#       Matrix::writeMM(adjacency_matrix, paste0(matfile, ".mtx"))
#     }
#   }
# }