# Inputs:
#       Merged RepSeq Data From Multiple Samples
#       Candidate Disease-Related TCR Sequences
# Do:
#       Build Cluster Networks Around First N Candidates
getClustersForSelectedClones <- function(
  data, # merged bulk rep-seq data from all patients/samples
  nucleo_col,
  amino_col,
  count_col,
  freq_col,
  vgene_col,
  dgene_col,
  jgene_col,
  sample_col,
  other_cols = NULL,
  selected_clones,
  selected_clone_labels = NULL,
  clone_seq_type = "amino_acid",
  dist_type = "hamming", # options are "hamming", "levenshtein", "euclidean_on_atchley"
  cluster_radius = 1,
  edge_dist = 1,
  node_colors = sample_col, # accepts multiple values (one plot per value)
  output_dir = NULL,
  return_plots = FALSE # should function return a list of dataframe + ggplots, or just print/write plots and return the dataframe?
) {
  if (dist_type == "euclidean_on_atchley" & clone_seq_type != "amino_acid") {
    stop("distance type 'euclidean_on_atchley' only applicable to amino acid sequences") }

  # Format the input data
  data <-
    data[ , # Keep only the relevant columns:
          c(nucleo_col, amino_col, count_col, freq_col,
            vgene_col, dgene_col, jgene_col, sample_col, other_cols)]
  colnames(data)[1:8] <- c(
    "nucleotideSeq", "aminoAcidSeq", "cloneCount", "cloneFrequency",
    "VGene", "DGene", "JGene", "sampleID")
  # Designate amino acid or nucleotide for clone sequence
  clone_col <- "aminoAcidSeq"
  if (cloneSeq %in% c("nucleo", "nucleotide")) { clone_col <- "nucleotideSeq" }

  # Determine distance type for cluster radius
  cluster_radius_dist_type <- "levenshtein"
  if (dist_type == "hamming") { cluster_radius_dist_type <- "hamming" }

  # Initialize output directory and objects
  if (!is.null(output_dir)) { .createOutputDir(output_dir) }
  data_all_clusters <- as.data.frame(matrix(nrow = 0, ncol = ncol(data) + 3))
  names(data_all_clusters) <- c(colnames(data), "cluster_seq")
  if (return_plots) { plots <- list() }

  # Iterate over selected clones
  for (i in 1:length(selected_clones)) {

    # Get cluster data
    cat(paste0("Gathering cluster data for selected clone sequence ", i, ": ", selected_clones[[i]], "...\n"))
    data_current_cluster <-
      getSimilarClones(selected_clones[[i]], data, clone_col, sample_col,
                       cluster_radius_dist_type, max_dist = cluster_radius)

    # Build cluster network
    cat("Computing edges between sequences in the cluster...\n")
    # if (keep_adjacency_matrix) {
    #   adjacency_matrix <-
    #     generateNetworkFromClones(data_current_cluster[ , clone_col],
    #                               dist_type = dist_type,
    #                               edge_dist = edge_dist,
    #                               contig_ids = rownames(data_current_cluster),
    #                               return_type = "adjacency_matrix")
    #   network <-
    #     generateNetworkFromAdjacencyMat(adjacency_matrix)
    # } else {
    network <-
      generateNetworkFromClones(data_current_cluster[ , clone_col],
                                dist_type = dist_type,
                                edge_dist = edge_dist,
                                contig_ids = rownames(data_current_cluster))
    # }

    # Compute network degree within current cluster
    data_current_cluster$degree_in_assoc_cluster <- igraph::degree(network)

    # Create labels for plots
    plot_title <- paste0("Cluster ", i, " (", selected_clones[[i]], ")")
    if (dist_type == "euclidean_on_atchley") {
      dist_label <- "sequence embeddings in Euclidean 30-space based on Atchley factor representation "
    } else { dist_label <- paste0(dist_type, " distance") }
    plot_subtitle <- NULL
    if (!is.null(selected_clone_labels)) {
      plot_subtitle <- selected_clone_labels[[i]] }
    plot_subtitle <- paste0(
      plot_subtitle,
      "\nNetwork includes clones with ", cluster_radius_dist_type,
      " distance at most ", cluster_radius, " from disease-associated sequence",
      "\nEdges based on ", dist_label, " (max edge dist = ", edge_dist, ")")

    # Generate plots of network graph
    cat("Generating plot(s)...\n")
    if (!is.null(disease_col)) {
      netplot_disease <- plotNetworkGraph(
        network, title = plot_title,
        subtitle = paste0(plot_subtitle, "\nNodes colored by disease status"),
        color_nodes_by = as.factor(data_current_cluster[ , disease_col]),
        size_nodes_by = as.numeric(data_current_cluster[ , clone_frac_col]),
        color_legend_title = "Disease Status",
        size_legend_title = "Clone Fraction")
      if (display_plots) { print(netplot_disease) }
    }
    if (!is.null(sample_col)) {
      netplot_sampleid <- plotNetworkGraph(
        network, title = plot_title,
        subtitle = paste0(plot_subtitle, "\nNodes colored by Sample ID"),
        color_nodes_by = as.factor(data_current_cluster[ , sample_col]),
        size_nodes_by = as.numeric(data_current_cluster[ , clone_frac_col]),
        color_legend_title = "Sample ID",
        size_legend_title = "Clone Fraction") +
        ggraph::scale_color_viridis(
          begin = 0, end = 1, direction = -1, option = "turbo", discrete = TRUE)
      if (display_plots) { print(netplot_sampleid) }
    }
    if (!is.null(pat_id_col)) {
      netplot_patid <- plotNetworkGraph(
        network, title = plot_title,
        subtitle = paste0(plot_subtitle, "\nNodes colored by Patient ID"),
        color_nodes_by = as.factor(data_current_cluster[ , pat_id_col]),
        size_nodes_by = as.numeric(data_current_cluster[ , clone_frac_col]),
        color_legend_title = "Patient ID",
        size_legend_title = "Clone Fraction") +
        ggraph::scale_color_viridis(
          begin = 0, end = 1, direction = -1, option = "turbo", discrete = TRUE)
      if (display_plots) { print(netplot_patid) }
    }
    if (is.null(netplot_disease) & is.null(netplot_sampleid) &
        is.null(netplot_patid)) {
      netplot_deg <- plotNetworkGraph(
        network, title = plot_title,
        subtitle = paste0(plot_subtitle, "\nNodes colored by Patient ID"),
        color_nodes_by = as.factor(data_current_cluster$deg),
        size_nodes_by = as.numeric(data_current_cluster[ , clone_frac_col]),
        color_legend_title = "Network Degree",
        size_legend_title = "Clone Fraction")
      if (display_plots) { print(netplot_deg) }
    }

    # Save results if applicable
    if (!is.null(output_dir)) {
      .saveResultsForCandidateSeqNetwork(
        network, data_current_cluster, adjacency_matrix, netplot_disease,
        netplot_sampleid, netplot_patid, netplot_deg, keep_adjacency_matrix,
        dist_type, output_dir, outfilestem = paste(i, selected_clones[[i]], sep = "_"))
    }
    # Add data for current cluster to combined cluster data
    data_all_clusters <-
      rbind(
        data_all_clusters,
        cbind(data_current_cluster,
              "assoc_clust_seq" = rep(selected_clones[[i]], nrow(data_current_cluster)),
              "assoc_clust_id" = rep(i, nrow(data_current_cluster))))
    # Add results for current candidate to output list if applicable
    if (return_output) {
      new_sublist <- list(igraph = network,
                          repseq_data = data_current_cluster)
      if (keep_adjacency_matrix) {
        new_sublist$adjacency_matrix <- adjacency_matrix }
      if (!is.null(netplot_disease)) {
        new_sublist$plot_disease <- netplot_disease }
      if (!is.null(netplot_sampleid)) {
        new_sublist$plot_sampleid <- netplot_sampleid }
      if (!is.null(netplot_patid)) {
        new_sublist$plot_patid <- netplot_patid }
      if (!is.null(netplot_deg)) {
        new_sublist$plot <- netplot_deg }
      outlist$new_sublist <- new_sublist
      names(outlist)[[length(names(outlist))]] <- selected_clones[[i]]
    }
  } # end looping over top candidates

  cat("Done building individual clusters for specified disease-associated clone sequences.\n")
  # Save combined cluster data if applicable
  if (!is.null(output_dir)) {
    utils::write.csv(
      data_all_clusters,
      file.path(output_dir, "data_all_clusters.csv"))
    cat(paste0(
      "All results saved to output directory '", output_dir, "'.\n"))
  }
  # Return final output list if applicable
  if (return_output) {
    outlist$combined_cluster_data <- data_all_clusters
    outlist$settings <- list(dist_type = dist_type,
                             cluster_radius = cluster_radius,
                             edge_dist = edge_dist)
    cat("Results combined into a list and returned as output.\n")
    return(outlist)
  }
}


# Helpers -----------------------------------------------------------------


# FUNCTION: EXTRACT DATA SUBSET FOR ALL SEQUENCES WITHIN SPECIFIED RADIUS OF
# TARGET SEQUENCE BY SPECIFIED DISTANCE TYPE
# If a sample_col is provided, only samples that possess the target sequence
# will be included
getSimilarClones <- function(
  target_seq, # specified candidate sequence for the neighborhood
  data, # data frame containing rep seq data, possibly from multiple samples
  clone_col, # col name/# containing clone sequences
  sample_col = NULL, # optional col name/# containing sample IDs (only samples possessing target seq will be included)
  dist_type = "hamming", # options are "hamming" and "levenshtein"
  max_dist = 2 # Maximum Levenshtein distance allowed for inclusion in neighborhood
) {
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }
  # If sample_id is supplied, subset data keeping only samples with target seq
  if (is.null(sample_col)) {
    data_samples_w_targetseq <- data
  } else {
    ### SUBSET DATA: SAMPLES WITH TARGET SEQ ###
    cat("Finding samples that possess the selected sequence...\n")
    # Get row ids of merged data corresponding to target seq
    rows_for_targetseq <- grep(pattern = paste0("^", target_seq, "$"),
                               x = data[ , clone_col])
    # Extract rows of merged data corresponding to samples with target seq
    data_samples_w_targetseq <-
      data[
        data[ , sample_col] %in% data[rows_for_targetseq, sample_col], ]
  }
  # remove seq with * and _
  data_samples_w_targetseq <-
    data_samples_w_targetseq[
      -grep("[*|_]", data_samples_w_targetseq[ , clone_col]), ]
  ### SUBSET DATA: NEIGHBORHOOD OF TARGET SEQUENCE ###
  # Compute list of bounded distances between target seq and seqs
  # possessed by samples with target seq (values are -1 where bound is exceeded)
  # returned vector will be of type integer; names will be the sequences
  cat("Extracting clone sequences similar to the selected sequence...\n")
  if (dist_type == "levenshtein") { dist_fun <- levDistBounded
  } else if (dist_type == "hamming") { dist_fun <- hamDistBounded
  } else { stop("invalid option for `dist_type`") }
  dists_to_targetseq <- sapply(
    X = data_samples_w_targetseq[ , clone_col],
    FUN = dist_fun, b = target_seq, k = max_dist)
  # get data for sequences within the specified radius
  data_targetseq_neighborhood <-
    data_samples_w_targetseq[dists_to_targetseq != -1, ]
  return(data_targetseq_neighborhood)
}




.saveResultsForCandidateSeqNetwork <- function(
  network, data_current_cluster, adjacency_matrix, netplot_disease,
  netplot_sampleid, netplot_patid, netplot_deg, keep_adjacency_matrix,
  dist_type, output_dir, outfilestem) {

  # pdf: network graphs
  plotfile <- paste0(outfilestem, "_network_plots.pdf")
  grDevices::pdf(file.path(output_dir, plotfile), width = 12, height = 8)
  if (!is.null(netplot_deg)) { print(netplot_deg) }
  if (!is.null(netplot_disease)) { print(netplot_disease) }
  if (!is.null(netplot_sampleid)) { print(netplot_sampleid) }
  if (!is.null(netplot_patid)) { print(netplot_patid) }
  grDevices::dev.off()

  # Save metadata for candidiate sequence network
  utils::write.csv(
    data_current_cluster,
    file.path(output_dir, paste0(outfilestem, "_network_metadata.csv")))

  # Save Network igraphs using edgelist format
  igraph::write_graph(
    network,
    file = file.path(
      output_dir, paste0(outfilestem, "_network_graph_edgelist.txt")),
    format = "edgelist")

  # Save adjacency matrices
  if (keep_adjacency_matrix) {
    matfile <- file.path(output_dir, paste0(outfilestem, "_adjacency_matrix"))
    if (dist_type == "euclidean_on_atchley") {
      utils::write.csv(adjacency_matrix, paste0(matfile, ".csv"))
    } else { #hamming/levenshtein returns sparse matrix
      Matrix::writeMM(adjacency_matrix, paste0(matfile, ".mtx"))
    }
  }
}