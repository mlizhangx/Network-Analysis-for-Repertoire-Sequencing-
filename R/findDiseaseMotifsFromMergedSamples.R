# Given Merged RepSeq Data From Multiple Samples, Build Cluster Networks Around
# Top N Candidate TCR Sequences for Disease Relation Ranked By Fisher P-Value
findDiseaseMotifsFromMergedSamples <- function(
  data,
  pat_id_col = 1, sample_id_col = 2, disease_col = 3,
  clone_col = 4, clone_frac_col = 5,
  disease_encoding, # unique values in disease_col corresponding to disease
  candidate_seqs = data[[clone_col]],
  pgen_values = NULL, # optional, length must match candidate_seqs
  top_n_candidates = min(10, length(candidate_seqs)), #to be ranked by fisher pv
  dist_type = "hamming", # options are "hamming", "levenshtein", "euclidean_on_atchley"
  cluster_radius = 2,
  edge_dist = 2,
  color_nodes_by = c("disease", "sample_id", "pat_id"), # character vector containing nonempty subset of these 3 options; 1 plot per value
  output_dir =  # if NULL, output not saved
    file.path(getwd(), paste("top", top_n_candidates, dist_type, "edgedist",
                             edge_dist, "radius", cluster_radius)),
  return_output = TRUE, # return results for all top N candidates when done?
  keep_adjacency_matrix = FALSE, # save/return adjacency matrix for each cluster
  display_plots = FALSE # print plots?
) {
  # Initialize output directory and objects
  if (!is.null(output_dir)) { .createOutputDir(output_dir) }
  netplot_disease <- netplot_sampleid <- netplot_patid <- adjacency_matrix <-
    NULL
  if (return_output) { outlist <- list() }
  # Compute Metadata & Fisher P-value for Candidate Seqs
  candidates_meta <- # rows will be sorted by fisher pv
    computeMetaForCandidateSeqs(
      data, pat_id_col, sample_id_col, disease_col, clone_col, clone_frac_col,
      disease_encoding, candidate_seqs, pgen_values)
  cat(paste0("Building short-range networks/clusters around each of the top ",
             top_n_candidates, " candidate disease motifs ranked by P-value on Fisher's exact test...\n"))
  # Iterate over specified number of top candidates
  for (i in 1:top_n_candidates) {
    # Get candidate sequence
    candseq <- candidates_meta[i, "clone_seq"]
    cat(paste0("Building network for candidate ", i,
               ": ", candseq, "\n"))

    # Subset Data For Sequences Nearby Candidate Seq
    # (only includes samples possessing the candidate seq)
    cluster_radius_dist_type <- "levenshtein"
    if (dist_type == "hamming") { cluster_radius_dist_type <- "hamming" }
    cat(paste0(
      "Collecting all samples possessing the candidate disease motif; extracting from these samples all sequences with ",
      cluster_radius_dist_type, " distance at most ", cluster_radius,
      " from the candidate disease motif...\n"))
    data_candseq_network <-
      subsetDataNearTargetMotif(
        data, clone_col, sample_id_col, candseq,
        cluster_radius_dist_type, max_dist = cluster_radius)

    # Compute Network Around Candidate Motif
    cat(paste0(
      "Computing network edges among the selected sequences based on an (inclusive) ",
      dist_type, " distance threshold of ", edge_dist, "...\n"))
    if (keep_adjacency_matrix) {
      adjacency_matrix <-
        generateNetworkFromClones(data_candseq_network[[clone_col]],
                                  dist_type = dist_type,
                                  edge_dist = edge_dist,
                                  contig_ids = rownames(data_candseq_network),
                                  return_type = "adjacency_matrix")
      set.seed(9999)
      network <-
        generateNetworkFromAdjacencyMat(adjacency_matrix)
    } else {
      network <-
        generateNetworkFromClones(data_candseq_network[[clone_col]],
                                  dist_type = dist_type,
                                  edge_dist = edge_dist,
                                  contig_ids = rownames(data_candseq_network))
    }
    # Add degree and cluster ID to data for candidate seq network
    data_candseq_network$deg <- igraph::degree(network)
    data_candseq_network$cluster_id <-
      igraph::cluster_fast_greedy(network)$membership

    # Create labels for plots
    if (is.null(pgen_values)) { pGen_label <- NULL
    } else { pGen_label <- paste0("pGen = ",
                                  signif(pgen_values[[i]], digits = 3), ", ") }
    fisher_pval <- signif(candidates_meta[i, "pv_fisher"], digits = 3)
    max_clonefrac <- signif(candidates_meta[i, "max_clonefrac"], digits = 3)
    count_shared_samples <- candidates_meta[i, "shared_by_n_samples"]
    count_shared_disease_pt <- candidates_meta[i, "shared_by_n_pt_disease"]
    count_shared_total_pt <-
      count_shared_disease_pt + candidates_meta[i, "absent_from_n_pt_disease"]
    if (dist_type == "euclidean_on_atchley") {
      dist_label <- "Euclidean distance on Atchley factor encoding"
    } else { dist_label <- paste0(dist_type, " distance") }
    plot_title <- paste0(
      "Network for Candidiate Disease Motif #", i, " (", candseq, ")\n",
      "by ", dist_label, " (max edge dist ", edge_dist, ")")
    plot_subtitle <- paste0(
      "Motif present in ", count_shared_samples, " samples and ", count_shared_total_pt, " patients (", count_shared_disease_pt, " of which have positive disease status)",
      "\n", pGen_label, "Fisher P-value = ", fisher_pval, ", Max Clone Fraction in Network = ", max_clonefrac,
      "\nNetwork includes sequences with ", cluster_radius_dist_type, " distance at most ", cluster_radius, " from candidate motif")

    # Generate plots of network graph
    cat("Creating plot(s) of the short-range network graph for the candidate disease motif...\n")
    if ("disease" %in% color_nodes_by) {
      netplot_disease <- plotNetworkGraph(
        network, title = plot_title,
        subtitle = paste0(plot_subtitle, "\nNodes colored by disease status"),
        color_nodes_by = as.factor(data_candseq_network[[disease_col]]),
        size_nodes_by = as.numeric(data_candseq_network[[clone_frac_col]]),
        color_legend_title = "Disease Status", size_legend_title = "Clone Fraction")
      if (display_plots) { print(netplot_disease) }
    }
    if ("sample_id" %in% color_nodes_by) {
      netplot_sampleid <- plotNetworkGraph(
        network, title = plot_title,
        subtitle = paste0(plot_subtitle, "\nNodes colored by Sample ID"),
        color_nodes_by = as.factor(data_candseq_network[[sample_id_col]]),
        size_nodes_by = as.numeric(data_candseq_network[[clone_frac_col]]),
        color_legend_title = "Sample ID",
        size_legend_title = "Clone Fraction") +
        ggraph::scale_color_viridis(
          begin = 0, end = 1, direction = -1, option = "turbo", discrete = TRUE)
      if (display_plots) { print(netplot_sampleid) }
    }
    if ("pat_id" %in% color_nodes_by) {
      netplot_patid <- plotNetworkGraph(
        network, title = plot_title,
        subtitle = paste0(plot_subtitle, "\nNodes colored by Patient ID"),
        color_nodes_by = as.factor(data_candseq_network[[pat_id_col]]),
        size_nodes_by = as.numeric(data_candseq_network[[clone_frac_col]]),
        color_legend_title = "Patient ID",
        size_legend_title = "Clone Fraction") +
        ggraph::scale_color_viridis(
          begin = 0, end = 1, direction = -1, option = "turbo", discrete = TRUE)
      if (display_plots) { print(netplot_patid) }
    }

    # Save results if applicable
    if (!is.null(output_dir)) {
      .saveResultsForCandidateSeqNetwork(
        candseq, network, data_candseq_network, adjacency_matrix, netplot_disease,
        netplot_sampleid, netplot_patid, color_nodes_by, keep_adjacency_matrix,
        dist_type, output_dir, outfilestem = paste(i, candseq, sep = "_"))
    }
    # Add results for current candidate to output list if applicable
    if (return_output) {
      new_sublist <- list(igraph = network,
                          network_meta = data_candseq_network)
      if (keep_adjacency_matrix) {
        new_sublist$adjacency_matrix <- adjacency_matrix }
      if ("disease" %in% color_nodes_by) {
        new_sublist$plot_disease <- netplot_disease }
      if ("sample_id" %in% color_nodes_by) {
        new_sublist$plot_sampleid <- netplot_sampleid }
      if ("pat_id" %in% color_nodes_by) {
        new_sublist$plot_patid <- netplot_patid }
      outlist$new_sublist <- new_sublist
      names(outlist)[[length(names(outlist))]] <- candseq
    }
  } # end looping over top candidates
  cat("All tasks complete.\n")
  # Return final output list if applicable
  if (return_output) {
    outlist$candidates_metadata <- candidates_meta
    outlist$settings <- list(top_n_candidates = top_n_candidates,
                             dist_type = dist_type,
                             cluster_radius = cluster_radius,
                             edge_dist = edge_dist,
                             color_nodes_by = color_nodes_by)
    return(outlist)
  }
}


# Helpers -----------------------------------------------------------------

# FUNCTION: COMPUTE METADATA FOR CANDIDATE SEQUENCES USING MERGED SAMPLE DATA
computeMetaForCandidateSeqs <- function(
  merged_samples,
  pat_id_col = 1, sample_id_col = 2, disease_col = 3,
  clone_col = 4, clone_frac_col = 5,
  disease_encoding,
  candidate_seqs = merged_samples[[clone_col]],
  pgen_values = NULL
) {
  ### INITIALIZE OUTPUT DATA FRAME ###
  candidate_seqs <- unique(candidate_seqs)
  out <- data.frame("clone_seq" = as.vector(as.character(candidate_seqs)))
  out$seq_length <- nchar(out$clone_seq)
  if (!is.null(pgen_values)) {
    if (length(pgen_values) != nrow(out)) stop(
      "length of 'pgen_values' does not match length of 'out'")
    out$pGen <- as.vector(as.numeric(pgen_values))
  }
  # Keep only unique list of candidate sequences
  out <- out[!duplicated(out$clone_seq), ]
  # Add variable for "shared by n samples" to candidate sequence data
  out$shared_by_n_samples <- rep(NA, nrow(out))
  # Add variable for "shared by n healthy patients" to candidate sequence data
  out$shared_by_n_pt_healthy <- rep(NA, nrow(out))
  # Add variable for "shared by n disease patients" to candidate sequence data
  out$shared_by_n_pt_disease <- rep(NA, nrow(out))
  # Add variable for "absent from n healthy patients" to candidate sequence data
  out$absent_from_n_pt_healthy <- rep(NA, nrow(out))
  # Add variable for "absent from n disease patients" to candidate sequence data
  out$absent_from_n_pt_disease <- rep(NA, nrow(out))
  # Add variable for Fisher P-value to candidate sequence data
  out$pv_fisher <- rep(NA, nrow(out))
  # Add variable for max clone fraction to candidate sequence data
  out$max_clonefrac <- rep(NA, nrow(out))

  # Get counts of number of total/disease/healthy patients
  n_pt_total <- length(unique(merged_samples[[pat_id_col]]))
  rowids_disease_pt <- merged_samples[[disease_col]] %in% disease_encoding
  rowids_healthy_pt <- !rowids_disease_pt
  n_pt_disease_total <-
    length(unique(merged_samples[rowids_disease_pt, pat_id_col]))
  n_pt_healthy_total <- n_pt_total - n_pt_disease_total
  cat(paste0("Data contains a total of ", n_pt_total, " patients: ",
             n_pt_disease_total, " disease patients and ",
             n_pt_healthy_total, " healthy patients.\n"))
  cat("Computing metadata for candidate sequences and performing Fisher's exact tests...\n")

  ### ITERATE OVER CANDIDATE SEQUENCES ###
  for (i in 1:nrow(out)) {
    ### GET CANDIDATE SEQUENCE ###
    candseq <- out$clone_seq[[i]]
    ### COMPUTE METADATA FOR CURRENT CANDIDATE SEQUENCE ###
    # logical vector for rows of merged data corresponding to candidate seq
    rowids_candseq <- grepl(pattern = paste0("^", candseq, "$"),
                            x = merged_samples[[clone_col]])
    # Number of samples sharing candidate sequence
    out$shared_by_n_samples[[i]] <-
      length(unique(merged_samples[rowids_candseq, sample_id_col]))
    # Number of disease patients with candidate sequence
    out$shared_by_n_pt_disease[[i]] <-
      length(unique(
        merged_samples[rowids_candseq & rowids_disease_pt, pat_id_col]))
    # Number of disease patients without candidate sequence
    out$absent_from_n_pt_disease[[i]] <-
      n_pt_disease_total - out$shared_by_n_pt_disease[[i]]
    # Number of healthy patients with candidate sequence
    out$shared_by_n_pt_healthy[[i]] <-
      length(unique(
        merged_samples[rowids_candseq & rowids_healthy_pt, pat_id_col]))
    # Number of healthy patients without candidate sequence
    out$absent_from_n_pt_healthy[[i]] <-
      n_pt_healthy_total - out$shared_by_n_pt_healthy[[i]]
    # Compute fisher p-value for candidate sequence
    out$pv_fisher[[i]] <-
      stats::fisher.test(data.frame(
        "healthy" = c(out$shared_by_n_pt_healthy[[i]],
                      out$absent_from_n_pt_healthy[[i]]),
        "disease" = c(out$shared_by_n_pt_disease[[i]],
                      out$absent_from_n_pt_disease[[i]]))
      )$p.value
    # Compute max clone fraction of candidate seq among samples containing it
    out$max_clonefrac[[i]] <- max(merged_samples[rowids_candseq, clone_frac_col])
  }
  # Sort candidate sequence metadata by fisher P-value and return
  return(out[order(out$pv_fisher), ])
}


# FUNCTION: EXTRACT DATA SUBSET FOR ALL SEQUENCES WITHIN SPECIFIED RADIUS OF
# TARGET MOTIF BY SPECIFIED DISTANCE TYPE
# If a sample_id_col is provided, only samples that possess the target sequence
# will be included
subsetDataNearTargetMotif <- function(
  data, # data frame containing rep seq data, possibly from multiple samples
  clone_col, # col name/# containing clone sequences
  sample_id_col = NULL, # optional col name/# containing sample IDs (only samples possessing target seq will be included)
  target_motif, # specified candidate sequence for the neighborhood
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
    rows_for_targetseq <- grep(pattern = paste0("^", target_motif, "$"),
                               x = data[[clone_col]])
    # Extract rows of merged data corresponding to samples with target seq
    data_samples_w_targetseq <-
      data[data[[sample_id_col]] %in% data[rows_for_targetseq, sample_id_col], ]
  }
  # remove seq with * and _
  data_samples_w_targetseq <-
    data_samples_w_targetseq[
      -grep("[*|_]", data_samples_w_targetseq[[clone_col]]), ]
  ### SUBSET DATA: NEIGHBORHOOD OF TARGET SEQUENCE ###
  # Compute list of bounded distances between target seq and seqs
  # possessed by samples with target seq (values are -1 where bound is exceeded)
  # returned vector will be of type integer; names will be the sequences
  if (dist_type == "levenshtein") { dist_fun <- levDistBounded
  } else if (dist_type == "hamming") { dist_fun <- hamDistBounded
  } else { stop("invalid option for `dist_type`") }
  dists_to_targetseq <- sapply(
    X = data_samples_w_targetseq[[clone_col]],
    FUN = dist_fun, b = target_motif, k = max_dist)
  # get data for sequences within the specified radius
  data_targetseq_neighborhood <-
    data_samples_w_targetseq[dists_to_targetseq != -1, ]
  return(data_targetseq_neighborhood)
}




.saveResultsForCandidateSeqNetwork <- function(
  candseq, network, data_candseq_network, adjacency_matrix, netplot_disease,
  netplot_sampleid, netplot_patid, color_nodes_by, keep_adjacency_matrix, dist_type,
  output_dir, outfilestem) {

  cat(paste0(
    "Saving results for short-range network around candidate disease motif ",
    candseq, "...\n"))

  # pdf: network graphs
  plotfile <- paste0(outfilestem, "_network_plots.pdf")
  grDevices::pdf(file.path(output_dir, plotfile), width = 12, height = 8)
  if ("disease" %in% color_nodes_by) { print(netplot_disease) }
  if ("sample_id" %in% color_nodes_by) { print(netplot_sampleid) }
  if ("pat_id" %in% color_nodes_by) { print(netplot_patid) }
  grDevices::dev.off()

  # Save metadata for candidiate sequence network
  utils::write.csv(
    data_candseq_network,
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