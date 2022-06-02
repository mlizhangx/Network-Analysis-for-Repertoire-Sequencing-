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

  # Clone sequence settings
  clone_seq_type = "amino acid",
  selected_clones,
  selected_clone_labels = NULL, # added to single-cluster plot subtitles
  drop_chars = "[*|_]", # passed to getSimilarClones() for getting cluster data

  # Network Settings
  dist_type = "hamming", # options are "hamming", "levenshtein", "euclidean_on_atchley"
  edge_dist = 1,
  cluster_radius = 1,
  node_stats = TRUE, # cluster_id will always be computed even if FALSE
  stats_to_include = node_stat_settings(), # cluster_id forced to TRUE
  cluster_stats = FALSE,

  # Plot Settings (global cluster network)
  custom_title = NULL, # for plot of global cluster network
  custom_subtitle = NULL,
  edge_width = 0.1,
  size_nodes_by = count_col, # can use a column name/# of data
  node_size_limits = NULL, # numeric length 2
  custom_size_legend = NULL, # custom legend title
  color_nodes_by = sample_col, # accepts multiple values (one plot per value)
  color_scheme = "default", # passed to plotNetworkGraph(); accepts multiple values (one per value of color_nodes_by)
  custom_color_legend = NULL, # custom title (length must match color_nodes_by)

  # Plot Settings (individual cluster plots)
  single_cluster_plots = TRUE,
  sc_edge_width = 0.3,
  sc_size_nodes_by = count_col,
  sc_node_size_limits = NULL,
  sc_custom_size_legend = NULL, # custom legend title
  sc_color_nodes_by = sample_col, # accepts multiple values (one plot per value)
  sc_color_scheme = "default", # passed to plotNetworkGraph(); accepts multiple values (one per value of color_nodes_by)
  sc_custom_color_legend = NULL, # custom title (length must match color_nodes_by)

  # Output Settings
  output_dir = getwd(),
  data_outfile = "data_global_cluster_network.csv",
  global_plot_outfile = "global_cluster_network_graph.pdf",
  cluster_plots_outfile = "individual_cluster_graphs.pdf",
  igraph_outfile = NULL, # network igraph of global cluster network
  matrix_outfile = NULL, # adjacency matrix for global cluster network
  plot_width = 12, # passed to pdf()
  plot_height = 10, # passed to pdf()
  return_all = FALSE # should function return a list, or just the data?

) {


  ### INPUT CHECKS ###
  # Atchley factor embedding only applicable to amino acid sequences
  if (dist_type == "euclidean_on_atchley" & clone_seq_type != "amino acid") {
    stop("distance type 'euclidean_on_atchley' only applicable to amino acid sequences") }

  # Check that each input column is a distinct col of data and meets specs

  # each variable for size/color exists or will be added to data

  #### PREPARE WORKING ENVIRONMENT ####
  # Create output directory if applicable
  if (!is.null(output_dir)) { .createOutputDir(output_dir) }

  # return type
  return_type <- ifelse(
    return_all,
    yes = "all",
    no = ifelse(cluster_stats &
                  (is.null(output_dir) | is.null(cluster_plots_outfile)),
                yes = "node_and_cluster_data",
                no = "node_data_only"))

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
  # if (is.integer(size_nodes_by)) { size_nodes_by <- names(data)[size_nodes_by] }
  if (is.numeric(sc_color_nodes_by)) { sc_color_nodes_by <- names(data)[sc_color_nodes_by] }
  # if (is.integer(sc_size_nodes_by)) { sc_size_nodes_by <- names(data)[sc_size_nodes_by] }

  # Designate amino acid or nucleotide for clone sequence
  clone_seq_col <- amino_col
  if (clone_seq_type == "nucleotide") { clone_seq_col <- nucleo_col }

  # Determine distance type for cluster radius
  cluster_radius_dist_type <- "levenshtein" # applies to lev and atchley
  if (dist_type == "hamming") { cluster_radius_dist_type <- "hamming" }

  # Format the input data
  extra_cols <- intersect(
    unique(c(other_cols, color_nodes_by, sc_color_nodes_by)), names(data))
  data <-
    data[ , # Keep only the relevant columns, in specified order:
          unique(c(nucleo_col, amino_col, count_col, freq_col,
                   vgene_col, dgene_col, jgene_col, cdr3length_col, sample_col,
                   extra_cols))]

  # Initialize output directory and objects
  data_all_clusters <- as.data.frame(matrix(nrow = 0, ncol = ncol(data) + 3))
  names(data_all_clusters) <- c(colnames(data), "assocClustID",
                                "assocClustSeq", "degreeInAssocClust")

  ### GLOBAL NETWORK PLOT SETTINGS ###
  # Title for global cluster network plot
  if (is.null(custom_title)) {
    main_title <- paste("Global network of clusters around",
                        length(selected_clones), "selected clone sequences")
  } else { main_title <- custom_title }

  # Default fixed component of subtitle across all plots
  subtitle_part <- ifelse(
    dist_type == "euclidean_on_atchley",
    yes = paste(
      "\nClone sequences embedded in Euclidean 30-space",
      "based on Atchley factor representation using deep learning",
      "\nEdges based on a maximum Euclidean distance of", edge_dist,
      "between embedded values\n"),
    no = paste(
      "\nEdges based on a maximum", dist_type, "distance of", edge_dist, "\n"))

  # Subtitle for global cluster network plot
  if (is.null(custom_subtitle)) {
    main_subtitle <- paste(
      "Network includes clone sequences with a maximum",
      cluster_radius_dist_type, "distance of", cluster_radius,
      "from one of the selected sequences",
      subtitle_part)
  } else { main_subtitle <- custom_subtitle }

  # # If multiple coloring variables, extend color scheme to vector if needed
  # if (length(color_nodes_by) > 1 & length(color_scheme) == 1) {
  #   color_scheme <- rep(color_scheme, length(color_nodes_by)) }

  # # Legend titles for size and color (global plot)
  # size_legend_title <- NULL # default for fixed node size
  # if (is.character(size_nodes_by)) {
  #   if (!is.null(custom_size_legend)) {
  #     size_legend_title <- custom_size_legend
  #   } else { size_legend_title <- size_nodes_by }
  # }
  # if (!is.null(custom_color_legend)) {
  #   color_legend_title <- custom_color_legend
  # } else { color_legend_title <- color_nodes_by }


  ### SINGLE-CLUSTER PLOT SETTINGS ###
  if (single_cluster_plots) {
    # Single-cluster plot subtitle (fixed part of subtitle across clusters)
    sc_subtitle <- paste(
      "Cluster includes clone sequences with a maximum",
      cluster_radius_dist_type, "distance of", cluster_radius,
      "from the central sequence",
      subtitle_part)

    # If multiple coloring variables, extend color scheme to vector if needed
    if (length(sc_color_nodes_by) > 1 & length(sc_color_scheme) == 1) {
      sc_color_scheme <- rep(sc_color_scheme, length(sc_color_nodes_by)) }

    # Legend titles for size and color (single-cluster plots)
    sc_size_legend_title <- NULL # default for fixed node size
    if (is.character(sc_size_nodes_by)) {
      if (!is.null(sc_custom_size_legend)) {
        sc_size_legend_title <- sc_custom_size_legend
      } else { sc_size_legend_title <- sc_size_nodes_by }
    }
    if (!is.null(sc_custom_color_legend)) {
      sc_color_legend_title <- sc_custom_color_legend
    } else { sc_color_legend_title <- sc_color_nodes_by }

    if (single_cluster_plots &
        (return_all | (!is.null(output_dir) & !is.null(cluster_plots_outfile)))) {
      sc_plots <- list()
    }
  }


  #### BUILD CLUSTERS FOR SELECTED CLONES ####
  # Iterate over the selected clones
  for (i in 1:length(selected_clones)) {

    ### BUILD SINGLE-CLUSTER NETWORK ###
    # Get cluster data
    cat(paste0("Gathering cluster data for clone sequence ", i,
               " (", selected_clones[[i]], "):\n"))
    data_current_cluster <-
      getSimilarClones(selected_clones[[i]], data, clone_seq_col, sample_col,
                       cluster_radius_dist_type, max_dist = cluster_radius,
                       drop_chars = drop_chars)
    if (nrow(data_current_cluster) < 2) {
      warning("not enough clones to build network; proceeding to next clone sequence")
      next
    }

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


    ### PLOT(S) OF SINGLE-CLUSTER NETWORK GRAPH ###
    if (single_cluster_plots) {
      # Create labels for plots
      plot_title <- paste0("Cluster ", i, " (", selected_clones[[i]], ")")
      sc_subtitle_prefix <- NULL
      if (!is.null(selected_clone_labels)) {
        sc_subtitle_prefix <- paste0(selected_clone_labels[[i]], "\n") }
      plot_subtitle <- paste0(sc_subtitle_prefix, sc_subtitle)

      # Ensure size_nodes_by is a vector or fixed value to use for node sizes
      if (is.character(sc_size_nodes_by)) {
        size_code <- data_current_cluster[ , sc_size_nodes_by]
      } else { size_code <- sc_size_nodes_by }

      # Create one plot for each variable in color_nodes_by
      temp_plotlist <- list()
      for (j in 1:length(sc_color_nodes_by)) {
        cat(paste0("Creating single-cluster graph with nodes colored by ",
                   sc_color_nodes_by[[j]], "..."))
        temp_plotlist$newplot <-
          plotNetworkGraph(
            network, sc_edge_width, title = plot_title,
            subtitle = plot_subtitle,
            color_nodes_by = data_current_cluster[ , sc_color_nodes_by[[j]]],
            size_nodes_by = size_code,
            color_legend_title = sc_color_legend_title[[j]],
            size_legend_title = sc_size_legend_title,
            color_scheme = sc_color_scheme[[j]],
            node_size_limits = sc_node_size_limits)
        print(temp_plotlist$newplot)
        names(temp_plotlist)[[length(names(temp_plotlist))]] <-
          sc_color_nodes_by[[j]]
        cat(" Done.\n")
      }

      # # Save plots for current cluster to a single pdf if applicable
      # if (save_plots & !single_plot & !is.null(output_dir)) {
      #   plot_outfile <- file.path(output_dir, "individual_cluster_plots",
      #                             paste0("cluster_", i, ".pdf"))
      #   grDevices::pdf(file = plot_outfile,
      #                  width = plot_width, height = plot_height)
      #   for (j in 1:length(color_nodes_by)) { print(temp_plotlist[[j]]) }
      #   grDevices::dev.off()
      #   cat(paste0("Cluster graph plot saved to file:\n  ", plot_outfile, "\n"))
      # }

      # Add plots to output list if applicable
      if (single_cluster_plots &
          (return_all | (!is.null(output_dir) & !is.null(cluster_plots_outfile)))) {
        sc_plots$newcluster <- temp_plotlist
        names(sc_plots)[[length(names(sc_plots))]] <- selected_clones[[i]]
      }
    }

  } # done looping over selected clones
  cat("All clusters complete.\n")

  # Format additional variables in data
  data_all_clusters$assocClusterID <- as.factor(data_all_clusters$assocClusterID)

  #### BUILD GLOBAL CLUSTER NETWORK ####
  # Ensure cluster ID is computed
  if (!node_stats) {
    stats_to_include <- "cluster_id_only"
  } else if (!stats_to_include$cluster_id) {
    stats_to_include$cluster_id <- TRUE }

  if ("globalClusterID" %in% color_nodes_by) {
    color_nodes_by[which(color_nodes_by == "globalClusterID")] <- "cluster_id"
  }

  cat("Building global cluster network using combined cluster data:\n")
  global_net <- buildRepSeqNetwork(
    data_all_clusters,
    nucleo_col, amino_col, count_col, freq_col, vgene_col, dgene_col, jgene_col,
    cdr3length_col,
    other_cols = c(sample_col, extra_cols, "assocClusterID", "assocClusterSeq",
                   "degreeInAssocCluster"),
    clone_seq_type,
    min_seq_length = NULL, dist_type = dist_type, edge_dist = edge_dist,
    node_stats = TRUE, stats_to_include = stats_to_include,
    plot_title = main_title, plot_subtitle = main_subtitle,
    edge_width = edge_width, size_nodes_by = size_nodes_by,
    node_size_limits = node_size_limits,
    custom_size_legend = custom_size_legend,
    color_nodes_by = color_nodes_by, color_scheme = color_scheme,
    custom_color_legend = custom_color_legend,
    return_all = TRUE)
  data_all_clusters <- global_net$node_data

  # Rename some columns of combined cluster data
  names(data_all_clusters)[
    which(names(data_all_clusters) == sample_col)] <- "sampleID"
  names(data_all_clusters)[
    which(names(data_all_clusters) == "cluster_id")] <- "globalClusterID"
  if ("degree" %in% names(data_all_clusters)) {
    names(data_all_clusters)[which(names(data_all_clusters) == "degree")] <-
      "globalDegree" }
  # colnames(data_all_clusters)[1:9] <- c(
  #   "nucleotideSeq", "aminoAcidSeq", "cloneCount", "cloneFreqInSample",
  #   "VGene", "DGene", "JGene", "CDR3Length", "sampleID")


  #### SAVE RESULTS ####
  if (!is.null(output_dir)) {
    # Save data for global cluster network if applicable
    if (!is.null(data_outfile)) {
      utils::write.csv(data_all_clusters, file.path(output_dir, data_outfile),
                       row.names = FALSE)
      cat(paste0(
        "Global cluster network data and rep-seq data saved to file:\n  ",
        file.path(output_dir, data_outfile), "\n"))
    }

    # Save global network cluster plots to a single pdf if applicable
    if (!is.null(global_plot_outfile)) {
      grDevices::pdf(file = file.path(output_dir, global_plot_outfile),
                     width = plot_width, height = plot_height)
      for (i in 1:length(global_net$plots)) { print(global_net$plots[[i]]) }
      grDevices::dev.off()
      cat(paste0("Plot of global cluster network graph saved to file:\n  ",
                 file.path(output_dir, global_plot_outfile), "\n"))
    }

    # Save all single-cluster plots to a single pdf if applicable
    if (!is.null(cluster_plots_outfile) & single_cluster_plots) {
      grDevices::pdf(file = file.path(output_dir, cluster_plots_outfile),
                     width = plot_width, height = plot_height)
      for (i in 1:length(sc_plots)) {
        for (j in 1:length(sc_plots[[i]])) { print(sc_plots[[i]][[j]]) } }
      grDevices::dev.off()
      cat(paste0("Individual single-cluster graph plots saved to file:\n  ",
                 file.path(output_dir, cluster_plots_outfile), "\n"))
    }

    # Save igraph
    if (!is.null(igraph_outfile)) {
      igraph::write_graph(global_net$igraph,
                          file = file.path(output_dir, igraph_outfile),
                          format = "edgelist")
      cat(paste0("Global cluster network igraph saved in edgelist format to file:\n  ",
                 file.path(output_dir, igraph_outfile), "\n")) }

    # Save adjacency matrix
    if (!is.null(matrix_outfile)) {
      if (dist_type == "euclidean_on_atchley") {
        utils::write.csv(adjacency_matrix,
                         file.path(output_dir, matrix_outfile),
                         row.names = FALSE)
        cat(paste0("Global cluster network adjacency matrix saved to file:\n  ",
                   file.path(output_dir, matrix_outfile), "\n"))
      } else {
        Matrix::writeMM(adjacency_matrix,
                        file.path(output_dir, matrix_outfile))
        cat(paste0("Global cluster network adjacency matrix saved to file:\n  ",
                   file.path(output_dir, matrix_outfile), "\n")) } }
  }


  #### RETURN OUTPUT ####
  if (return_type == "node_data_only") {

    cat("All tasks complete.\n")
    return(data_all_clusters)

  } else {
    out <- list("data" = data_all_clusters)
    if (cluster_stats) { out$cluster_stats <- global_net$cluster_stats }
    if (return_type == "all") {
      out$global_plots <- global_net$plots
      out$cluster_plots <- sc_plots
      out$adjacency_matrix <- adjacency_matrix
      out$igraph <- global_net$igraph }
    cat(paste0("All tasks complete. Returning a list containing the following items:\n  ",
               paste(names(out), collapse = ", "), "\n"))

    return(out)

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