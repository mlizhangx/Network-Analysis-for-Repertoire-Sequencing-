# Inputs:
#       Combined TCR rep-seq data from multiple samples
#       List of selected clone sequences (CDR3 nucleotide or amino acid)
# Do:
#       Build a cluster around each of the selected clones
#       Print the plot for each cluster and return the combined cluster data

getAssociatedClusters <- function(

  # Input Data/Settings
  data, # merged bulk rep-seq data from all patients/samples
  seq_col,
  # vgene_col = NULL,
  # dgene_col = NULL,
  # jgene_col = NULL,
  # cdr3length_col = NULL,
  sample_col,
  target_seqs,
  # count_col = NULL,
  # freq_col = NULL,
  other_cols = NULL,
  target_seq_labels = NULL, # added to single-cluster plot subtitles
  drop_chars = "[*|_]", # passed to getSimilarClones() for getting cluster data

  # Network Settings
  dist_type = "hamming", # options are "hamming", "levenshtein", "euclidean_on_atchley"
  dist_cutoff = 1,
  neighborhood_radius = 1,
  node_stats = TRUE, # cluster_id will always be computed even if FALSE
  stats_to_include = node_stat_settings(), # cluster_id forced to TRUE
  cluster_stats = FALSE,

  # Plot Settings (global network)
  plot_title = "auto", # for plot of global cluster network
  plot_subtitle = "auto",
  color_nodes_by = sample_col, # accepts multiple values (one plot per value)
  color_scheme = "auto", # passed to plotNetworkGraph(); accepts multiple values (one per value of color_nodes_by)
  color_legend = TRUE,
  color_title = "auto", # custom title (length must match color_nodes_by)
  edge_width = 0.1,
  size_nodes_by = 0.5, # can use a column name of data (a numeric value yields fixed node sizes)
  node_size_limits = "auto", # numeric length 2
  size_title = "auto", # custom legend title

  # Plot Settings (individual neighborhoods)
  neighborhood_plots = TRUE,
  nbd_color_nodes_by = sample_col, # accepts multiple values (one plot per value)
  nbd_color_scheme = "default", # passed to plotNetworkGraph(); accepts multiple values (one per value of color_nodes_by)
  nbd_color_legend = TRUE,
  nbd_color_title = "auto", # custom title (length must match color_nodes_by)
  nbd_edge_width = 0.3,
  nbd_size_nodes_by = 1,
  nbd_node_size_limits = "auto",
  nbd_size_title = "auto", # custom legend title

  # Output Settings
  print_plots = TRUE,
  output_dir = file.path(getwd(), "associated_clusters"),
  data_outfile = "global_network_node_data.csv",
  global_plot_outfile = "global_network_graph.pdf",
  nbd_plots_outfile = "neighborhood_graphs.pdf",
  cluster_stats_outfile = "global_network_cluster_data.csv",
  igraph_outfile = NULL, # network igraph of global cluster network
  matrix_outfile = NULL, # adjacency matrix for global cluster network
  plot_width = 12, # passed to pdf()
  plot_height = 10, # passed to pdf()
  return_all = FALSE # should function return a list, or just the data?

) {


  ### INPUT CHECKS ###
  # # Atchley factor embedding only applicable to amino acid sequences
  # if (dist_type == "euclidean_on_atchley" & clone_seq_type != "amino acid") {
  #   stop("distance type 'euclidean_on_atchley' only applicable to amino acid sequences") }

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
                  (is.null(output_dir) | is.null(nbd_plots_outfile)),
                yes = "node_and_cluster_data",
                no = "node_data_only"))

  # Convert column references to character if not already
  if (is.numeric(seq_col)) { seq_col <- names(data)[seq_col] }
  # if (is.numeric(nucleo_col)) { nucleo_col <- names(data)[nucleo_col] }
  # if (is.numeric(amino_col)) { amino_col <- names(data)[amino_col] }
  # if (is.numeric(count_col)) { count_col <- names(data)[count_col] }
  # if (is.numeric(freq_col)) { freq_col <- names(data)[freq_col] }
  # if (is.numeric(vgene_col)) { vgene_col <- names(data)[vgene_col] }
  # if (is.numeric(dgene_col)) { dgene_col <- names(data)[dgene_col] }
  # if (is.numeric(jgene_col)) { jgene_col <- names(data)[jgene_col] }
  # if (is.numeric(cdr3length_col)) { cdr3length_col <- names(data)[cdr3length_col] }
  if (is.numeric(sample_col)) { sample_col <- names(data)[sample_col] }
  if (is.numeric(other_cols)) { other_cols <- names(data)[other_cols] }
  if (is.numeric(color_nodes_by)) { color_nodes_by <- names(data)[color_nodes_by] }
  if (is.numeric(nbd_color_nodes_by)) { nbd_color_nodes_by <- names(data)[nbd_color_nodes_by] }

  # # Designate amino acid or nucleotide for clone sequence
  # clone_seq_col <- amino_col
  # if (clone_seq_type == "nucleotide") { clone_seq_col <- nucleo_col }

  # Determine distance type for cluster radius
  neighborhood_radius_dist_type <- "levenshtein" # applies to lev and atchley
  if (dist_type == "hamming") { neighborhood_radius_dist_type <- "hamming" }

  # Format the input data
  # Coerce sequence column to character if needed
  if (!is.character(data[ , seq_col])) {
    data[ , seq_col] <- as.character(data[ , seq_col]) }

  if (!is.null(other_cols)) {
    extra_cols <- intersect(
      unique(c(other_cols, color_nodes_by, nbd_color_nodes_by)), names(data))
    data <-
      data[ , # Keep only the relevant columns, in specified order:
            unique(c(seq_col,
                     # count_col,
                     # freq_col,
                     # vgene_col, dgene_col, jgene_col, cdr3length_col,
                     sample_col,
                     extra_cols))]
  }

  # Rename frequency column if provided
  # if (!is.null(freq_col)) {
  #   new_freq_colname <- "FreqInSample"
  #   if (size_nodes_by == freq_col) { size_nodes_by <- new_freq_colname }
  #   if (nbd_size_nodes_by == freq_col) { nbd_size_nodes_by <- new_freq_colname }
  #   if (freq_col %in% color_nodes_by) {
  #     color_nodes_by[color_nodes_by == freq_col] <- new_freq_colname
  #   }
  #   if (freq_col %in% nbd_color_nodes_by) {
  #     nbd_color_nodes_by[nbd_color_nodes_by == freq_col] <- new_freq_colname
  #   }
  #   names(data)[names(data) == freq_col] <- new_freq_colname
  #   freq_col <- new_freq_colname
  # }

  # Initialize output directory and objects
  data_all_nbds <- as.data.frame(matrix(nrow = 0, ncol = ncol(data) + 3))
  names(data_all_nbds) <- c(colnames(data), "NeighborhoodID",
                            "NeighborhoodSeq", "DegreeInNeighborhood")

  ### GLOBAL NETWORK PLOT SETTINGS ###
  # Title for global cluster network plot
  if (plot_title == "auto") {
    main_title <- paste("Global network of associated clusters")
  } else { main_title <- plot_title }

  # Default fixed component of subtitle across all plots
  subtitle_part <- ifelse(
    dist_type == "euclidean_on_atchley",
    yes = paste(
      "\nAmino acid sequences embedded in Euclidean 30-space",
      "based on Atchley factor representation using deep learning",
      "\nEdges denote a maximum Euclidean distance of", dist_cutoff,
      "between embedded values\n"),
    no = paste(
      "\nEdges denote a maximum", dist_type, "distance of", dist_cutoff, "between receptor sequences\n"))

  # Subtitle for global cluster network plot
  if (plot_subtitle == "auto") {
    main_subtitle <- paste(
      "Each node represents a TCR/BCR cell or clone\nEach associated cluster is a neighborhood of cells/clones whose receptor sequences have a max",
      neighborhood_radius_dist_type, "distance of", neighborhood_radius,
      "from a target sequence",
      subtitle_part)
  } else { main_subtitle <- plot_subtitle }


  ### SINGLE-CLUSTER PLOT SETTINGS ###
  if (neighborhood_plots) {
    # Single-cluster plot subtitle (fixed part of subtitle across clusters)
    nbd_subtitle <- paste(
      "Each node represents a TCR/BCR cell or clone\nNeighborhood includes cells/clones whose receptor sequences have a max",
      neighborhood_radius_dist_type, "distance of", neighborhood_radius,
      "from the target sequence",
      subtitle_part)

    # If multiple coloring variables, extend color scheme and legend title to vectors if needed
    if (length(nbd_color_nodes_by) > 1) {
      if (length(nbd_color_scheme) == 1) { # extend to vector
        nbd_color_scheme <- rep(nbd_color_scheme, length(nbd_color_nodes_by)) }
      if (!is.null(nbd_color_title)) {
        if (length(nbd_color_title) == 1) { # extend to vector
          nbd_color_title <- rep(nbd_color_title, length(nbd_color_nodes_by)) }
      } else { # nbd_color_title is NULL
        # vector of empty titles (hack, since can't have NULL vector entries)
        nbd_color_title <- rep("", length(nbd_color_nodes_by))
      }
    }

    # Set default color legend title if applicable
    if (!is.null(nbd_color_title)) {
      for (i in 1:length(nbd_color_title)) {
        if (nbd_color_title[[i]] == "auto") {
          nbd_color_title[[i]] <- nbd_color_nodes_by[[i]]
        }
      }
    }
    # Set default size legend title if applicable
    if (!is.null(nbd_size_title)) {
      if (nbd_size_title == "auto") {
        if (is.numeric(nbd_size_nodes_by)) {
          nbd_size_title <- NULL
        }
        if (is.character(nbd_size_nodes_by)) {
          nbd_size_title <- nbd_size_nodes_by
        }
      }
    }

    # Initialize storage for plots if applicable
    if (neighborhood_plots &
        (return_all | (!is.null(output_dir) & !is.null(nbd_plots_outfile)))) {
      nbd_plots <- list()
    }
  }


  #### BUILD NEIGHBORHOODS AROUND TARGET CLONES ####
  # Iterate over the selected clones
  for (i in 1:length(target_seqs)) {

    ### BUILD NEIGHBORHOOD ###
    # Get neighborhood data
    cat(paste0("Gathering neighborhood data for target sequence ", i,
               " (", target_seqs[[i]], "):\n"))
    data_current_nbd <-
      getSimilarClones(target_seqs[[i]], data, seq_col, sample_col,
                       neighborhood_radius_dist_type, max_dist = neighborhood_radius,
                       drop_chars = drop_chars)
    if (nrow(data_current_nbd) < 2) {
      warning("not enough cells/clones to build network; proceeding to next target sequence")
      next
    }

    # Generate adjacency matrix for network
    adjacency_matrix <-
      generateNetworkFromSeqs(data_current_nbd[ , seq_col],
                              dist_type, dist_cutoff,
                              contig_ids = rownames(data_current_nbd),
                              return_type = "adjacency_matrix")

    # Subset data to keep only those clones in the network (nonzero degree)
    if (dist_type != "euclidean_on_atchley") {
      data_current_nbd <-
        data_current_nbd[as.numeric(dimnames(adjacency_matrix)[[1]]), ] }

    # Generate network from adjacency matrix
    network <- generateNetworkFromAdjacencyMat(adjacency_matrix)


    # Add variables for cluster ID, central sequence, and degree in cluster
    data_current_nbd$NeighborhoodID <- i
    data_current_nbd$NeighborhoodSeq <- target_seqs[[i]]
    data_current_nbd$DegreeInNeighborhood <- igraph::degree(network)

    # Add data for current cluster to combined cluster data
    data_all_nbds <- rbind(data_all_nbds, data_current_nbd)


    ### PLOT(S) OF SINGLE-CLUSTER NETWORK GRAPH ###
    if (neighborhood_plots) {
      # Create labels for plots
      current_title <- paste0("Neighborhood ", i, " (", target_seqs[[i]], ")")
      nbd_subtitle_prefix <- NULL
      if (!is.null(target_seq_labels)) {
        nbd_subtitle_prefix <- paste0(target_seq_labels[[i]], "\n") }
      current_subtitle <- paste0(nbd_subtitle_prefix, nbd_subtitle)

      # If using column to size nodes, get column vector from column name
      if (is.character(nbd_size_nodes_by)) {
        size_code <- data_current_nbd[ , nbd_size_nodes_by]
      } else {
        size_code <- nbd_size_nodes_by # numeric value = fixed node size
      }

      # Create one plot for each variable in color_nodes_by
      temp_plotlist <- list()
      if (is.null(nbd_color_nodes_by)) {
        cat("Creating neighborhood network graph plot...")
        temp_plotlist$uniform_color <-
          plotNetworkGraph(
            network, nbd_edge_width,
            title = current_title,
            subtitle = current_subtitle,
            color_nodes_by = NULL,
            color_scheme = nbd_color_scheme[[j]],
            color_legend_title = nbd_color_title[[j]],
            show_color_legend = nbd_color_legend,
            size_nodes_by = size_code,
            size_legend_title = nbd_size_title,
            node_size_limits = nbd_node_size_limits)
        if (print_plots) { print(temp_plotlist$uniform_color) }
        cat(" Done.\n")
      } else {
        for (j in 1:length(nbd_color_nodes_by)) {
          cat(paste0(
            "Creating neighborhood network graph with nodes colored by ",
            nbd_color_nodes_by[[j]], "..."))
          temp_plotlist$newplot <-
            plotNetworkGraph(
              network, nbd_edge_width,
              title = current_title,
              subtitle = current_subtitle,
              color_nodes_by = data_current_nbd[ , nbd_color_nodes_by[[j]]],
              color_scheme = nbd_color_scheme[[j]],
              color_legend_title = nbd_color_title[[j]],
              show_color_legend = nbd_color_legend,
              size_nodes_by = size_code,
              size_legend_title = nbd_size_title,
              node_size_limits = nbd_node_size_limits)
          if (print_plots) { print(temp_plotlist$newplot) }
          names(temp_plotlist)[[length(names(temp_plotlist))]] <-
            nbd_color_nodes_by[[j]]
          cat(" Done.\n")
        }
      }

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
    if (neighborhood_plots &
        (return_all | (!is.null(output_dir) & !is.null(nbd_plots_outfile)))) {
      nbd_plots$newnbd <- temp_plotlist
      names(nbd_plots)[[length(names(nbd_plots))]] <- target_seqs[[i]]
    }
  }

  # done looping over selected sequences
  cat("All neighborhoods complete.\n")

  # Format additional variables in data
  data_all_nbds$NeighborhoodID <- as.factor(data_all_nbds$NeighborhoodID)


  #### BUILD GLOBAL CLUSTER NETWORK ####
  # Ensure cluster ID is computed
  if (!node_stats) {
    stats_to_include <- "cluster_id_only"
  } else if (!stats_to_include$cluster_id) {
    stats_to_include$cluster_id <- TRUE }

  if ("GlobalClusterID" %in% color_nodes_by) {
    color_nodes_by[which(color_nodes_by == "GlobalClusterID")] <- "cluster_id"
  }

  if (!is.null(other_cols)) {
    global_other_cols <- c(sample_col, extra_cols,
                           "NeighborhoodID", "NeighborhoodSeq",
                           "DegreeInNeighborhood")
  } else {
    global_other_cols <- NULL
  }

  cat("Building global network using data from all neighborhoods:\n")
  global_net <- buildRepSeqNetwork(
    data_all_nbds, seq_col,
    other_cols = global_other_cols,
    min_seq_length = NULL, dist_type = dist_type, dist_cutoff = dist_cutoff,
    node_stats = TRUE, stats_to_include = stats_to_include,
    cluster_stats = cluster_stats,
    plot_title = main_title, plot_subtitle = main_subtitle,
    edge_width = edge_width, size_nodes_by = size_nodes_by,
    node_size_limits = node_size_limits,
    size_title = size_title,
    color_nodes_by = color_nodes_by, color_scheme = color_scheme,
    color_title = color_title, color_legend = color_legend,
    print_plots = print_plots, return_all = TRUE)

  # Rename some columns of combined cluster data
  # if (!is.null(freq_col)) {
  #   names(global_net$node_data)[
  #     which(names(global_net$node_data) == "CloneFrequency")] <- new_freq_colname
  # }
  # names(global_net$node_data)[
  #   which(names(global_net$node_data) == sample_col)] <- "SampleID"
  names(global_net$node_data)[
    which(names(global_net$node_data) == "cluster_id")] <- "GlobalClusterID"
  if ("degree" %in% names(global_net$node_data)) {
    names(global_net$node_data)[which(names(global_net$node_data) == "degree")] <-
      "GlobalDegree" }


  #### SAVE RESULTS ####
  if (!is.null(output_dir)) {
    # Save node meta data for global network if applicable
    if (!is.null(data_outfile)) {
      utils::write.csv(global_net$node_data, file.path(output_dir, data_outfile),
                       row.names = FALSE)
      cat(paste0(
        "Node-level meta-data for global network saved to file:\n  ",
        file.path(output_dir, data_outfile), "\n"))
    }

    # Save global network plots to a single pdf if applicable
    if (!is.null(global_plot_outfile)) {
      grDevices::pdf(file = file.path(output_dir, global_plot_outfile),
                     width = plot_width, height = plot_height)
      for (i in 1:length(global_net$plots)) { print(global_net$plots[[i]]) }
      grDevices::dev.off()
      cat(paste0("Plot(s) for global network saved to file:\n  ",
                 file.path(output_dir, global_plot_outfile), "\n"))
    }

    # Save all neighborhood plots to a single pdf if applicable
    if (!is.null(nbd_plots_outfile) & neighborhood_plots) {
      grDevices::pdf(file = file.path(output_dir, nbd_plots_outfile),
                     width = plot_width, height = plot_height)
      for (i in 1:length(nbd_plots)) {
        for (j in 1:length(nbd_plots[[i]])) { print(nbd_plots[[i]][[j]]) } }
      grDevices::dev.off()
      cat(paste0("Plots for individual neighborhood networks saved to file:\n  ",
                 file.path(output_dir, nbd_plots_outfile), "\n"))
    }
    # Save cluster data for global network if applicable
    if (cluster_stats & !is.null(cluster_stats_outfile)) {
      utils::write.csv(global_net$cluster_stats,
                       file.path(output_dir, cluster_stats_outfile),
                       row.names = FALSE)
      cat(paste0(
        "Cluster-level meta-data for global network saved to file:\n  ",
        file.path(output_dir, data_outfile), "\n"))
    }

    # Save igraph
    if (!is.null(igraph_outfile)) {
      igraph::write_graph(global_net$igraph,
                          file = file.path(output_dir, igraph_outfile),
                          format = "edgelist")
      cat(paste0("Global network igraph saved in edgelist format to file:\n  ",
                 file.path(output_dir, igraph_outfile), "\n")) }

    # Save adjacency matrix
    if (!is.null(matrix_outfile)) {
      if (dist_type == "euclidean_on_atchley") {
        utils::write.csv(adjacency_matrix,
                         file.path(output_dir, matrix_outfile),
                         row.names = FALSE)
        cat(paste0("Global network adjacency matrix saved to file:\n  ",
                   file.path(output_dir, matrix_outfile), "\n"))
      } else {
        Matrix::writeMM(adjacency_matrix,
                        file.path(output_dir, matrix_outfile))
        cat(paste0("Global network adjacency matrix saved to file:\n  ",
                   file.path(output_dir, matrix_outfile), "\n")) } }
  }


  #### RETURN OUTPUT ####
  if (return_type == "node_data_only") {

    cat("All tasks complete.\n")
    return(global_net$node_data)

  } else {
    out <- list("data" = global_net$node_data)
    if (cluster_stats) { out$cluster_stats <- global_net$cluster_stats }
    if (return_type == "all") {
      out$global_plots <- global_net$plots
      if (neighborhood_plots) { out$cluster_plots <- nbd_plots }
      out$adjacency_matrix <- adjacency_matrix
      out$igraph <- global_net$igraph
    }
    cat(paste0("All tasks complete. Returning a list containing the following items:\n  ",
               paste(names(out), collapse = ", "), "\n"))

    return(out)

  }

}

