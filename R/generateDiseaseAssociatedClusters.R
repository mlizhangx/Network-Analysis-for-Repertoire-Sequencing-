# Inputs:
#       Merged RepSeq Data From Multiple Samples
#       Candidate Disease-Related TCR Sequences
# Do:
#       Build Cluster Networks Around First N Candidates
generateDiseaseAssociatedClusters <- function(
  data, # merged bulk rep-seq data from all patients/samples
  clone_col, # column of data containing clone sequences
  clone_frac_col, # column of data containing clone counts/fractions
  pat_id_col = NULL, # optional column of data containing patient ID (to color nodes)
  sample_id_col = NULL, # optional column of data containing sample ID (to color nodes)
  disease_col = NULL, # optional column of data containing disease status (to color nodes)
  assoc_clones_meta, # meta for disease-assoc seqs
  meta_clone_col = 1, # col of assoc_clones_meta with clone seqs
  meta_label_col = ncol(assoc_clones_meta), # optional col of assoc_clones_meta with label for plot subtitle (use NULL to exclude)
  # num_assoc_clones = min(20, nrow(assoc_clones_meta)),
  dist_type = "hamming", # options are "hamming", "levenshtein", "euclidean_on_atchley"
  cluster_radius = 1,
  edge_dist = 1,
  output_dir =  # if NULL, output not saved
    file.path(getwd(), paste(dist_type, "edgedist", edge_dist,
                             "radius", cluster_radius, sep = "_")),
  return_output = TRUE, # return results for all top N candidates when done?
  keep_adjacency_matrix = FALSE, # save/return adjacency matrix for each cluster
  display_plots = TRUE # print plots?
) {
  if (is.null(output_dir) & !return_output) stop("Function requires an output directory when called with 'return_output = FALSE'")

  # Initialize output directory and objects
  data_all_clusters_combined <- matrix(nrow = 0, ncol = ncol(data) + 1)
  colnames(data_all_clusters_combined) <- c(colnames(data), "cluster_seq")
  data_all_clusters_combined <- as.data.frame(data_all_clusters_combined)
  if (!is.null(output_dir)) { .createOutputDir(output_dir) }
  netplot_disease <- netplot_sampleid <- netplot_patid <- netplot_deg <-
    adjacency_matrix <- NULL
  if (return_output) { outlist <- list() }
  cat(paste0("Building cluster networks around each of the ",
             nrow(assoc_clones_meta), " disease-associated clone sequences contained in 'assoc_clones_meta'...\n"))

  # Iterate over specified number of candidates
  for (i in 1:nrow(assoc_clones_meta)) {
    # Get candidate sequence
    assoc_clone <- assoc_clones_meta[i, meta_clone_col]
    cat(paste0("Building cluster network for disease-associated clone sequence ", i, ": ", assoc_clone, "\n"))

    # Subset Data For Sequences Nearby Candidate Seq
    # (only includes samples possessing the candidate seq)
    cluster_radius_dist_type <- "levenshtein"
    if (dist_type == "hamming") { cluster_radius_dist_type <- "hamming" }
    cat(paste0(
      "Collecting samples that share the disease-associated clone sequence. \n",
      "Extracting sequences whose ", cluster_radius_dist_type,
      " distance from the disease-associated clone sequence is at most ", cluster_radius, "...\n"))
    data_current_cluster <-
      subsetDataNearTargetSeq(
        data, clone_col, sample_id_col, assoc_clone,
        cluster_radius_dist_type, max_dist = cluster_radius)

    # Compute Cluster Network Around Candidate Seq
    cat(paste0(
      "Computing network edges based on a maximum ",
      dist_type, " distance of ", edge_dist, "...\n"))
    if (keep_adjacency_matrix) {
      adjacency_matrix <-
        generateNetworkFromClones(data_current_cluster[ , clone_col],
                                  dist_type = dist_type,
                                  edge_dist = edge_dist,
                                  contig_ids = rownames(data_current_cluster),
                                  return_type = "adjacency_matrix")
      network <-
        generateNetworkFromAdjacencyMat(adjacency_matrix)
    } else {
      network <-
        generateNetworkFromClones(data_current_cluster[ , clone_col],
                                  dist_type = dist_type,
                                  edge_dist = edge_dist,
                                  contig_ids = rownames(data_current_cluster))
    }
    # Add degree to data for candidate seq network
    data_current_cluster$deg <- igraph::degree(network)

    # Create labels for plots
    if (dist_type == "euclidean_on_atchley") {
      dist_label <- "Euclidean distance on Atchley factor encoding"
    } else { dist_label <- paste0(dist_type, " distance") }
    plot_title <- paste0(
      "Cluster for Disease-Associated Clone Seq #", i, " (", assoc_clone, ")")
    plot_subtitle <- NULL
    if (!is.null(meta_label_col)) {
      plot_subtitle <- assoc_clones_meta[i, meta_label_col] }
    plot_subtitle <- paste0(
      plot_subtitle,
      "\nNetwork includes clones with ", cluster_radius_dist_type,
      " distance at most ", cluster_radius, " from disease-associated sequence",
      "\nEdges between clones based on ", dist_label, " (max edge dist = ", edge_dist, ")")

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
    if (!is.null(sample_id_col)) {
      netplot_sampleid <- plotNetworkGraph(
        network, title = plot_title,
        subtitle = paste0(plot_subtitle, "\nNodes colored by Sample ID"),
        color_nodes_by = as.factor(data_current_cluster[ , sample_id_col]),
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
        dist_type, output_dir, outfilestem = paste(i, assoc_clone, sep = "_"))
    }
    # Add data for current cluster to combined cluster data
    data_all_clusters_combined <-
      rbind(
        data_all_clusters_combined,
        cbind(data_current_cluster,
              "assoc_clust_seq" = rep(assoc_clone, nrow(data_current_cluster)),
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
      names(outlist)[[length(names(outlist))]] <- assoc_clone
    }
  } # end looping over top candidates

  cat("Done building individual clusters for specified disease-associated clone sequences.\n")
  # Save combined cluster data if applicable
  if (!is.null(output_dir)) {
    utils::write.csv(
      data_all_clusters_combined,
      file.path(output_dir, "data_all_clusters_combined.csv"))
    cat(paste0(
      "All results saved to output directory '", output_dir, "'.\n"))
  }
  # Return final output list if applicable
  if (return_output) {
    outlist$combined_cluster_data <- data_all_clusters_combined
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
# If a sample_id_col is provided, only samples that possess the target sequence
# will be included
subsetDataNearTargetSeq <- function(
  data, # data frame containing rep seq data, possibly from multiple samples
  clone_col, # col name/# containing clone sequences
  sample_id_col = NULL, # optional col name/# containing sample IDs (only samples possessing target seq will be included)
  target_seq, # specified candidate sequence for the neighborhood
  dist_type = "hamming", # options are "hamming" and "levenshtein"
  max_dist = 2 # Maximum Levenshtein distance allowed for inclusion in neighborhood
) {
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }
  # If sample_id is supplied, subset data keeping only samples with target seq
  if (is.null(sample_id_col)) {
    data_samples_w_targetseq <- data
  } else {
    ### SUBSET DATA: SAMPLES WITH TARGET SEQ ###
    # Get row ids of merged data corresponding to target seq
    rows_for_targetseq <- grep(pattern = paste0("^", target_seq, "$"),
                               x = data[ , clone_col])
    # Extract rows of merged data corresponding to samples with target seq
    data_samples_w_targetseq <-
      data[
        data[ , sample_id_col] %in% data[rows_for_targetseq, sample_id_col], ]
  }
  # remove seq with * and _
  data_samples_w_targetseq <-
    data_samples_w_targetseq[
      -grep("[*|_]", data_samples_w_targetseq[ , clone_col]), ]
  ### SUBSET DATA: NEIGHBORHOOD OF TARGET SEQUENCE ###
  # Compute list of bounded distances between target seq and seqs
  # possessed by samples with target seq (values are -1 where bound is exceeded)
  # returned vector will be of type integer; names will be the sequences
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