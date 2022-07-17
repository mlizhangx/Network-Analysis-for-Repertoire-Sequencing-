# Inputs:
#       rep-seq data from multiple samples (separate files, same directory)
# Do:
#       Identify public clones and build public clone network

findPublicClusters <- function(

  # Input Data and Columns
  input_dir = getwd(),
  file_list,
  sample_id_list,
  group_list = NULL, #exactly one of group_list or group_cols must be non-NULL
  csv_files = FALSE, # use read.csv instead of read.table?
  header = FALSE,
  sep = "",
  nucleo_col,
  amino_col,
  count_col,
  freq_col,
  vgene_col,
  dgene_col,
  jgene_col,
  cdr3length_col,
  group_col = NULL, #exactly one of group_list or group_cols must be non-NULL
  other_cols = NULL,

  # Clone sequence settings
  clone_seq_type = "amino acid",
  min_seq_length = 3, # min clone seq length
  drop_chars = "[*|_]",
  aggregate_identical_clones = FALSE,

  # Network Settings
  dist_type = "hamming", # options are "hamming", "levenshtein", "euclidean_on_atchley"
  edge_dist = 1,

  # Filter pass settings for sample-level clusters
  top_n_clusters = 20,
  min_node_count = 10,
  min_clone_count = 100,

  # Downstream filter settings for public clones
  bayes_factor_col = NULL, # column containing adjusted bayes factor pvalues
  bayes_factor_cutoff = 0.05,
  diff_test_col = NULL, # column containing pvalues from differential testing
  diff_test_cutoff = 0.05,

  # Settings for K-means clustering on Atchley factor
  kmeans_atchley = ifelse(clone_seq_type == "nucleotide", FALSE, TRUE), # include kmeans clustering? only applicable to TCRB CDR3 amino acid seqs
  k = 100, # number of clusters
  k_plot_width = 15,
  k_plot_height = 15,
  k_plot_margin = 15,
  k_plot_viridis = FALSE, # use viridis color palettes for color-blindness robustness


  # Plot Settings (public cluster network)
  custom_title = NULL, #
  custom_subtitle = NULL,
  edge_width = 0.1,
  size_nodes_by = 0.5, # can use a column name of data (a numeric value yields fixed node sizes)
  node_size_limits = NULL, # numeric length 2
  custom_size_legend = NULL, # custom legend title
  color_nodes_by = c(group_col, "cluster_id", "SampleID"), # accepts multiple values (one plot per value)
  color_scheme = c("viridis", "turbo", "turbo"), # passed to plotNetworkGraph(); accepts multiple values (one per value of color_nodes_by)
  custom_color_legend = NULL, # custom title (length must match color_nodes_by)

  # Plot Settings (public cluster-core network)
  cores_edge_width = 0.3,
  cores_size_nodes_by = "node_count",
  cores_node_size_limits = NULL,
  cores_custom_size_legend = NULL, # custom legend title
  cores_color_nodes_by = "SampleID", # accepts multiple values (one plot per value)
  cores_color_scheme = "default", # passed to plotNetworkGraph(); accepts multiple values (one per value of color_nodes_by)
  cores_custom_color_legend = NULL, # custom title (length must match color_nodes_by)

  # Plot Settings (sample-level networks)
  sample_edge_width = 0.3,
  sample_size_nodes_by = count_col,
  sample_node_size_limits = NULL,
  sample_custom_size_legend = NULL, # custom legend title
  sample_color_nodes_by = "cluster_id", # accepts multiple values (one plot per value)
  sample_color_scheme = "default", # passed to plotNetworkGraph(); accepts multiple values (one per value of color_nodes_by)
  sample_custom_color_legend = NULL, # custom title (length must match color_nodes_by)

  # Output Settings
  output_dir = file.path(getwd(), "public_clusters"),
  plot_width = 12, # passed to pdf()
  plot_height = 10 # passed to pdf()

) {


  ### INPUT CHECKS ###
  # Atchley factor embedding only applicable to amino acid sequences
  if (dist_type == "euclidean_on_atchley" & clone_seq_type != "amino acid") {
    stop("distance type 'euclidean_on_atchley' only applicable to amino acid sequences") }
  # warn about computational limits for atchley compared to ham/lev?

  # each variable for size/color exists or will be added to data
  # need to account for aggregate_identical_clones if TRUE


  #### PREPARE WORKING ENVIRONMENT ####
  # Create output directory if applicable
  if (!is.null(output_dir)) { .createOutputDir(output_dir) }

  # Check that #exactly one of group_list or group_cols is non-NULL


  # Convert input columns to character if not already
  if (is.numeric(nucleo_col)) { nucleo_col <- names(data)[nucleo_col] }
  if (is.numeric(amino_col)) { amino_col <- names(data)[amino_col] }
  if (is.numeric(count_col)) { count_col <- names(data)[count_col] }
  if (is.numeric(freq_col)) { freq_col <- names(data)[freq_col] }
  if (is.numeric(vgene_col)) { vgene_col <- names(data)[vgene_col] }
  if (is.numeric(dgene_col)) { dgene_col <- names(data)[dgene_col] }
  if (is.numeric(jgene_col)) { jgene_col <- names(data)[jgene_col] }
  if (is.numeric(cdr3length_col)) { cdr3length_col <- names(data)[cdr3length_col] }
  if (is.numeric(group_col)) { group_col <- names(data)[group_col] }
  if (is.numeric(other_cols)) { other_cols <- names(data)[other_cols] }
  if (is.numeric(color_nodes_by)) { color_nodes_by <- names(data)[color_nodes_by] }
  if (is.numeric(cores_color_nodes_by)) { cores_color_nodes_by <- names(data)[cores_color_nodes_by] }
  if (is.numeric(sample_color_nodes_by)) { sample_color_nodes_by <- names(data)[sample_color_nodes_by] }
  if ("SampleLevelClusterID" %in% sample_color_nodes_by) {
    sample_color_nodes_by[
      which(sample_color_nodes_by == "SampleLevelClusterID")] <- "cluster_id"
  }
  if ("PublicClusterID" %in% color_nodes_by) {
    color_nodes_by[
      which(color_nodes_by == "PublicClusterID")] <- "cluster_id"
  }

  # Designate amino acid or nucleotide for clone sequence
  clone_seq_col <- amino_col
  if (clone_seq_type == "nucleotide") { clone_seq_col <- nucleo_col }

  # Initialize output directory and objects
  data_public_clusters <- data_cores_network <- sample_net <- extra_cols <-
    keep_cols <- NULL #init

  # New name for frequency column
  old_freq_colname <- freq_col
  new_freq_colname <- ifelse(aggregate_identical_clones,
                             yes = "AggCloneFreqInSample",
                             no = "CloneFreqInSample")


  # Apply new name for freq column to color/size variables for plots
  old_sample_color_nodes_by <- sample_color_nodes_by
  old_sample_size_nodes_by <- sample_size_nodes_by
  old_color_nodes_by <- color_nodes_by
  old_size_nodes_by <- size_nodes_by
  if (size_nodes_by == old_freq_colname) { size_nodes_by <- new_freq_colname }
  if (cores_size_nodes_by == old_freq_colname) {
    warning("can't size nodes for cluster-level network by clone frequency since this variable isn't present at the cluster level. Defaulting to sizing nodes by the total clone count in the cluster.")
    cores_size_nodes_by <- "agg_clone_count" }
  if (sample_size_nodes_by == old_freq_colname) { sample_size_nodes_by <- new_freq_colname }
  if (old_freq_colname %in% color_nodes_by) {
    color_nodes_by[color_nodes_by == old_freq_colname] <- new_freq_colname
  }
  if (old_freq_colname %in% cores_color_nodes_by) {

    warning("can't color nodes for cluster-level network by clone frequency since this variable isn't present at the cluster level. Defaulting to coloring nodes by the total clone count in the cluster.")
    cores_color_nodes_by[cores_color_nodes_by == old_freq_colname] <- "agg_clone_count"
  }
  if (old_freq_colname %in% sample_color_nodes_by) {
    sample_color_nodes_by[sample_color_nodes_by == old_freq_colname] <- new_freq_colname
  }

  #### BUILD SAMPLE-LEVEL NETWORKS ####
  cat("Building individual networks for each sample and searching for public clusters...\n")
  for (i in 1:length(file_list)) {
    cat(paste0("Processing data for sample ", i, ": ", sample_id_list[[i]], "\n"))

    cat("Loading and formatting data...")
    # Load data
    data <-
      ifelse(csv_files,
             read.csv(file.path(input_dir, file_list[[i]]), header, sep),
             read.table(file.path(input_dir, file_list[[i]]), header, sep))

    # Check that each input column is a distinct col of data and meets specs
    # Only need to perform this check for i = 1
    #

    # Format the input data
    if (i == 1) {
      extra_cols <-
        intersect(
          unique(c(other_cols, bayes_factor_col, diff_test_col,
                   old_color_nodes_by, old_size_nodes_by,
                   old_sample_color_nodes_by, old_sample_size_nodes_by)),
          names(data))
      keep_cols <-
        unique(c(nucleo_col, amino_col, count_col, old_freq_colname,
                 vgene_col, dgene_col, jgene_col, cdr3length_col, group_col,
                 extra_cols))
    }
    data <- data[ , keep_cols]
    names(data)[names(data) == old_freq_colname] <- new_freq_colname
    data$SampleID <- sample_id_list[[i]]
    cat(" Done.\n")


    ### BUILD SAMPLE NETWORK ###
    cat("Building network for current sample...\n")
    sample_net <- buildRepSeqNetwork(
      data, nucleo_col, amino_col, count_col, new_freq_colname, vgene_col,
      dgene_col, jgene_col, cdr3length_col, c(extra_cols, "SampleID"),
      clone_seq_type, min_seq_length, drop_chars, aggregate_identical_clones,
      grouping_cols = group_col, #only used if aggregate_identical_clones = TRUE
      dist_type = dist_type, edge_dist = edge_dist,
      node_stats = TRUE, stats_to_include = "all",
      cluster_stats = TRUE,
      plot_title = paste("Network for Sample:", sample_id_list[[i]]),
      plot_subtitle = NULL,
      edge_width = sample_edge_width, size_nodes_by = sample_size_nodes_by,
      node_size_limits = sample_node_size_limits,
      custom_size_legend = sample_custom_size_legend,
      color_nodes_by = sample_color_nodes_by,
      color_scheme = sample_color_scheme,
      custom_color_legend = sample_custom_color_legend,
      return_all = TRUE)

    # Rename sample-level node stats
    if (aggregate_identical_clones) {
      names(sample_net$node_data)[
        names(sample_net$node_data) == "AggregatedCloneFrequency"] <-
        new_freq_colname
    }
    names(sample_net$node_data)[names(sample_net$node_data) == "degree"] <-
      "SampleLevelNetworkDegree"
    names(sample_net$node_data)[names(sample_net$node_data) == "cluster_id"] <-
      "SampleLevelClusterID"
    names(sample_net$node_data)[names(sample_net$node_data) == "transitivity"] <-
      "SampleLevelTransitivity"
    names(sample_net$node_data)[names(sample_net$node_data) == "Closeness"] <-
      "SampleLevelCloseness"
    names(sample_net$node_data)[names(sample_net$node_data) == "centrality_by_closeness"] <-
      "SampleLevelCentralityByCloseness"
    names(sample_net$node_data)[names(sample_net$node_data) == "eigen_centrality"] <-
      "SampleLevelEigenCentrality"
    names(sample_net$node_data)[names(sample_net$node_data) == "centrality_by_eigen"] <-
      "SampleLevelCentralityByEigen"
    names(sample_net$node_data)[names(sample_net$node_data) == "betweenness"] <-
      "SampleLevelBetweenness"
    names(sample_net$node_data)[names(sample_net$node_data) == "centrality_by_betweenness"] <-
      "SampleLevelCentralityByBetweenness"
    names(sample_net$node_data)[names(sample_net$node_data) == "authority_score"] <-
      "SampleLevelAuthorityScore"
    names(sample_net$node_data)[names(sample_net$node_data) == "coreness"] <-
      "SampleLevelCoreness"
    names(sample_net$node_data)[names(sample_net$node_data) == "page_rank"] <-
      "SampleLevelPageRank"

    # Rename cluster stats
    names(sample_net$cluster_stats)[names(sample_net$cluster_stats) == "cluster_id"] <-
      "SampleLevelClusterID"
    names(sample_net$cluster_stats)[names(sample_net$cluster_stats) == "transitivity"] <-
      "cluster_transitivity"

    # Add sample ID to cluster stats
    sample_net$cluster_stats$SampleID <- sample_id_list[[i]]

    ## Save node & cluster data & plots ##
    sample_output_dir <- file.path(output_dir, "sample-level networks",
                                   sample_id_list[[i]])
    cat(paste0("Saving output for current sample to:\n  ", sample_output_dir,
               "\n"))
    utils::write.csv(
      sample_net$node_data,
      file = file.path(sample_output_dir, "node_level_meta_data.csv"),
      row.names = FALSE)
    cat("Node-level meta data saved as 'node_level_meta_data.csv'\n")
    utils::write.csv(sample_net$cluster_stats,
                     file = file.path(sample_output_dir, "cluster_info.csv"),
                     row.names = FALSE)
    cat("Cluster meta data saved as 'cluster_info.csv'\n")
    grDevices::pdf(
      file = file.path(sample_output_dir, "network_graph_plot.pdf"),
      width = plot_width, height = plot_height)
    for (j in 1:length(sample_net$plots)) { print(sample_net$plots[[j]]) }
    grDevices::dev.off()
    cat("Network graph plot saved as 'network_graph_plot.pdf'\n")


    ### FIND PUBLIC CLUSTERS AND ADD TO PUBLIC CLUSTER DATA ###
    cat("Finding public clusters in the current sample...")
    # row ids of top n clusters by node count
    ids_top_n_clusters <-
      which(
        sample_net$cluster_stats$SampleLevelClusterID %in%
          sample_net$cluster_stats[
            order(
              sample_net$cluster_stats$node_count, decreasing = TRUE
            )[1:top_n_clusters]
            ,
            SampleLevelClusterID
          ]
      )
    # row ids of clusters with at least min_node_count nodes
    ids_node_count_gt_n <-
      which(sample_net$cluster_stats$node_count >= min_node_count)
    # row ids of clusters with at least min_clone_count total clone count
    ids_clone_count_gt_n <-
      which(sample_net$cluster_stats$agg_clone_count >= min_clone_count)
    # Take the union of the filtered row ids
    filtered_ids <-
      unique(c(ids_top_n_clusters, ids_node_count_gt_n, ids_clone_count_gt_n))
    sample_net$cluster_stats <- sample_net$cluster_stats[filtered_ids, ]
    cat(" Done.\n")

    # (GET NODE-LEVEL DATA FOR FILTERED CLUSTERS AND ADD TO PUBLIC CLONE DATA)
    sample_net$node_data <- # subset rows:
      sample_net$node_data[ # keep if cluster ID appears in filtered cluster data
        sample_net$node_data$SampleLevelClusterID %in%
          sample_net$cluster_stats$SampleLevelClusterID
        , ]
    # Filter using Bayes FDR if supplied
    if (!is.null(bayes_factor_col)) {
      cat("Filtering nodes in public clusters for current sample by FDR-adjusted Bayes Factor P-value...")
      sample_net$node_data <- # subset rows:
        sample_net$node_data[ # keep if Bayes FDR is below cutoff
          sample_net$node_data[ , bayes_factor_col] < bayes_factor_cutoff
          , ]
      cat(" Done.\n")
    }
    # Filter using differential testing P-value if supplied
    if (!is.null(diff_test_col)) {
      cat("Filtering nodes in public clusters for current sample by P-value from differential testing...")
      sample_net$node_data <- # subset rows:
        sample_net$node_data[ # keep if Bayes FDR is below cutoff
          sample_net$node_data[ , diff_test_col] < diff_test_cutoff
          , ]
      cat(" Done.\n")
    }
    cat("Adding public clusters for current sample to existing public cluster data...")
    data_cores_network <- # Add cluster-level data to public cluster-cores network
      ifelse(i == 1,
             yes = sample_net$cluster_stats,
             no = rbind(data_cores_network, sample_net$cluster_stats))
    data_public_clusters <- # Add node-level data to public cluster network
      ifelse(i == 1,
             yes = sample_net$node_data,
             no = rbind(data_public_clusters, sample_net$node_data))
    cat(" Done.\n")


  } # done looping over selected clones
  cat("All samples have been processed.\n")



  #### BUILD NETWORK OF PUBLIC CLUSTERS BY REPRESENTATIVE SEQUENCE ####
  cat("Building public network using a single representative (core) clone from each cluster...\n")
  data_cores_network$empty <- NA # placeholder column
  cores_nucleo_col <- cores_amino_col <- "empty"
  if (clone_seq_type == "nucleotide") {
    cores_nucleo_col <- "seq_w_max_count"
  } else { cores_amino_col <- "seq_w_max_count" }
  cores_extra_cols <-
    c("SampleID", "SampleLevelClusterID", "node_count", "mean_seq_length",
      "mean_degree", "max_degree", "seq_w_max_degree", "max_clone_count",
      "diameter_length", "assortativity", 'cluster_transitivity',
      "edge_density", "degree_centrality_index", "closeness_centrality_index",
      "eigen_centrality_index", "eigen_centrality_eigenvalue")
  # Build network
  cores_network <- buildRepSeqNetwork(
    data_cores_network, cores_nucleo_col, cores_amino_col,
    "agg_clone_count", "empty", "empty", "empty", "empty", "empty",
    cores_extra_cols, clone_seq_type,
    dist_type = dist_type, edge_dist = edge_dist,
    drop_isolated_nodes = FALSE, # (keep zero-degree nodes)
    node_stats = TRUE, stats_to_include = "all",
    cluster_stats = FALSE,
    plot_title =
      "Network Using One Representative Clone From Each Public Cluster",
    plot_subtitle =
      paste0(
        "Includes top ", top_n_clusters,
        " clusters from each sample by node count; clusters with node count > ",
        min_node_count - 1, "; and clusters with total clone count > ",
        min_clone_count - 1),
    edge_width = cores_edge_width, size_nodes_by = cores_size_nodes_by,
    node_size_limits = cores_node_size_limits,
    custom_size_legend = cores_custom_size_legend,
    color_nodes_by = cores_color_nodes_by,
    color_scheme = cores_color_scheme,
    custom_color_legend = cores_custom_color_legend,
    return_all = TRUE)

  # Drop/rename variables for node-level data
  if (clone_seq_type == "nucleotide") {
    drop_cols <- c(2, 4:8)
  } else { drop_cols <- c(1, 4:8) }
  cores_network$node_data <- cores_network$node_data[ , -drop_cols]
  names(cores_network$node_data)[1:2] <-
    c("RepresentativeCloneSeq", "AggCloneCount")
  names(cores_network$node_data)[
    names(cores_network$node_data) == "node_count"] <- "NodeCountInSampleNetwork"
  names(cores_network$node_data)[
    names(cores_network$node_data) == "cluster_id"] <- "ClusterIDInCoresNetwork"
  names(cores_network$node_data)[
    names(cores_network$node_data) == "degree"] <- "DegreeInCoresNetwork"

  # Compute cluster level cluster stats
  # Tabulate the number of nodes in each cluster
  cores_network_cluster_info <-
    as.data.frame(table(cores_network$node_data[ , "ClusterIDInCoresNetwork"]))
  colnames(cores_network_cluster_info) <-
    c("ClusterIDInCoresNetwork", "NodeCountInCoresNetwork")
  num_clusters <- nrow(cores_network_cluster_info) # Total number of clusters
  cat(paste0("Computing statistics for the ", num_clusters, " cluster-level clusters..."))

  ### INITIALIZE VALUES ###
  cores_network_cluster_info$AggNodeCountInSampleNetwork <- 0
  cores_network_cluster_info$AggCloneCountInSampleNetwork <- 0
  cores_network_cluster_info$MeanValueOfMeanSeqLength <- 0
  cores_network_cluster_info$MeanDegreeInCoresNetwork <- 0
  cores_network_cluster_info$MaxDegreeInCoresNetwork <- 0
  cores_network_cluster_info$SeqWithMaxDegreeInCoresNetwork <- ""
  cores_network_cluster_info$MaxCloneCount <- 0
  cores_network_cluster_info$SampleWithMaxCloneCount <- ""
  cores_network_cluster_info$SeqWithMaxCloneCount <- ""
  cores_network_cluster_info$MaxAggCloneCount <- 0
  cores_network_cluster_info$SampleWithMaxAggCloneCount <- ""
  cores_network_cluster_info$SeqWithMaxAggCloneCount <- ""
  cores_network_cluster_info$DiameterLength <- 0
  cores_network_cluster_info$Assortativity <- 0
  cores_network_cluster_info$Transitivity <- 0
  cores_network_cluster_info$EdgeDensity <- 0
  cores_network_cluster_info$DegreeCentralityIndex <- 0
  cores_network_cluster_info$ClosenessCentralityIndex <- 0
  cores_network_cluster_info$EigenCentralityIndex <- 0
  cores_network_cluster_info$EigenCentralityEigenvalue <- 0

  ### COMPUTE STATS FOR EACH CLUSTER ###
  for (i in 1:num_clusters) {

    cores_cluster_id <- which(cores_network_cluster_info$ClusterIDInCoresNetwork == i) # current row of cluster cores_network$node_data
    node_ids <- cores_network$node_data$ClusterIDInCoresNetwork == i  # Rows of node cores_network$node_data for current cluster

    # Aggregate Sample-level Node Count
    cores_network_cluster_info$AggNodeCountInSampleNetwork[[cores_cluster_id]] <-
      sum(cores_network$node_data[node_ids, "NodeCountInSampleNetwork"])

    # Aggregate Sample-level Clone Count
    cores_network_cluster_info$AggCloneCountInSampleNetwork[[cores_cluster_id]] <-
      sum(cores_network$node_data[node_ids, "AggCloneCount"])

    # Mean value of mean sequence length
    cores_network_cluster_info$MeanValueOfMeanSeqLength[[cores_cluster_id]] <-
      round(mean(cores_network$node_data[node_ids, "mean_seq_length"]), 2)

    # Mean degree in cluster
    cores_network_cluster_info$MeanDegreeInCoresNetwork[[cores_cluster_id]] <-
      round(mean(cores_network$node_data[node_ids, "DegreeInCoresNetwork"]), 2)

    # Maximum degree (and corresponding seq) within cluster
    max_deg <- max(cores_network$node_data[node_ids, "DegreeInCoresNetwork"])
    cores_network_cluster_info$MaxDegreeInCoresNetwork[[cores_cluster_id]] <- max_deg
    node_id_max_deg <-
      which(node_ids &
              cores_network$node_data[ , "DegreeInCoresNetwork"] == max_deg)[[1]]
    cores_network_cluster_info$SeqWithMaxDegreeInCoresNetwork[[cores_cluster_id]] <-
      as.character(
        cores_network$node_data[node_id_max_deg, "RepresentativeCloneSeq"])

    # max value of max clone count (and corresponding sample & seq)
    max_count <- max(cores_network$node_data[node_ids, "max_clone_count"])
    cores_network_cluster_info$MaxCloneCount[[cores_cluster_id]] <- max_count
    node_id_max_count <-
      which(node_ids &
              cores_network$node_data[ , "max_clone_count"] == max_count)[[1]]
    cores_network_cluster_info$SampleWithMaxCloneCount[[cores_cluster_id]] <-
      cores_network$node_data[node_id_max_count, "SampleID"]
    cores_network_cluster_info$SeqWithMaxCloneCount[[cores_cluster_id]] <-
      cores_network$node_data[node_id_max_count, "RepresentativeCloneSeq"]

    # max value of agg clone count (and corresponding sample & seq)
    max_agg_count <- max(cores_network$node_data[node_ids, "AggCloneCount"])
    cores_network_cluster_info$MaxAggCloneCount[[cores_cluster_id]] <- max_agg_count
    node_id_max_agg_count <-
      which(node_ids &
              cores_network$node_data[ , "AggCloneCount"] == max_agg_count)[[1]]
    cores_network_cluster_info$SampleWithMaxAggCloneCount[[cores_cluster_id]] <-
      cores_network$node_data[node_id_max_agg_count, "SampleID"]
    cores_network_cluster_info$SeqWithMaxAggCloneCount[[cores_cluster_id]] <-
      cores_network$node_data[node_id_max_agg_count, "RepresentativeCloneSeq"]


    # Build cluster network to get network properties for the cluster
    cluster <- generateNetworkFromAdjacencyMat(
      as.matrix(cores_network$adjacency_matrix[node_ids, node_ids]))

    # Diameter (longest geodesic distance)
    cores_network_cluster_info$DiameterLength[[cores_cluster_id]] <-
      length(igraph::get_diameter(cluster, directed = T))

    # Assortativity
    cores_network_cluster_info$Assortativity[[cores_cluster_id]] <-
      igraph::assortativity_degree(cluster, directed = F)

    # Transitivity
    cores_network_cluster_info$Transitivity[[cores_cluster_id]] <-
      igraph::transitivity(cluster, type = "global")  # cluster is treated as an undirected network

    # Density: The proportion of present edges from all possible ties.
    cores_network_cluster_info$EdgeDensity[[cores_cluster_id]] <-
      igraph::edge_density(cluster, loops = F)

    # Centralization on degree
    cores_network_cluster_info$DegreeCentralityIndex[[cores_cluster_id]] <-
      igraph::centr_degree(cluster, mode = "in", normalized = T)$centralization

    # Centralization on Closeness (centrality based on distance to others in the graph)
    cores_network_cluster_info$ClosenessCentralityIndex[[cores_cluster_id]] <-
      igraph::centr_clo(cluster, mode = "all", normalized = T)$centralization

    # Centralization on Eigenvector (centrality proportional to the sum of connection centralities)
    #  (values of the first eigenvector of the graph adjacency matrix)
    cores_network_cluster_info$EigenCentralityIndex[[cores_cluster_id]] <-
      igraph::centr_eigen(cluster, directed = T, normalized = T)$centralization

    cores_network_cluster_info$EigenCentralityEigenvalue[[cores_cluster_id]] <-
      igraph::eigen_centrality(cluster, directed = T, weights = NA)$value
  }
  cat(" Done.\n")

  ## Save node & cluster data & plots ##
  cores_output_dir <- file.path(output_dir, "network_on_public_cluster_cores")
  cat(paste0("Saving output for network of public cluster cores to:\n  ",
             cores_output_dir, "\n"))
  utils::write.csv(
    cores_network$node_data,
    file = file.path(cores_output_dir, "node_level_meta_data.csv"),
    row.names = FALSE)
  cat("Node-level meta data saved as 'node_level_meta_data.csv'\n")
  utils::write.csv(cores_network_cluster_info,
                   file = file.path(cores_output_dir, "cluster_info.csv"),
                   row.names = FALSE)
  cat("Cluster meta data saved as 'cluster_info.csv'\n")
  grDevices::pdf(
    file = file.path(cores_output_dir, "network_graph_plot.pdf"),
    width = plot_width, height = plot_height)
  for (j in 1:length(cores_network$plots)) { print(cores_network$plots[[j]]) }
  grDevices::dev.off()
  cat("Network graph plot saved as 'network_graph_plot.pdf'\n")



  #### BUILD FULL PUBLIC CLUSTER NETWORK ####
  cat("Building network using all clones in the public clusters...\n")

  pub_cols <- c("SampleID", other_cols,
                bayes_factor_col, diff_test_col,
                color_nodes_by, size_nodes_by,
                "SampleLevelClusterID", "SampleLevelTransitivity",
                "SampleLevelCloseness", "SampleLevelCentralityByCloseness",
                "SampleLevelEigenCentrality", "SampleLevelCentralityByEigen",
                "SampleLevelBetweenness", "SampleLevelCentralityByBetweenness",
                "SampleLevelAuthorityScore", "SampleLevelCoreness",
                "SampleLevelPageRank")
  if (aggregate_identical_clones) { pub_cols <- c("UniqueCloneCount", pub_cols) }
  pub_extra_cols <- intersect(unique(pub_cols), names(data_public_clusters))
  pub_nucleo_col <- "NucleotideSeq"
  pub_amino_col  <-  "AminoAcidSeq"
  if (aggregate_identical_clones) {
    data_public_clusters$empty <- NA # placeholder variable
    if (clone_seq_type == "nucleotide") {
      pub_amino_col <- "empty"
    } else {
      pub_nucleo_col <- "empty"
    }
    pub_count_col <- "AggregatedCloneCount"
    pub_freq_col <- "AggCloneFreqInSample"
    pub_vgene_col <- "empty"
    pub_dgene_col <- "empty"
    pub_jgene_col <- "empty"
    pub_cdr3length_col <- "empty"
  } else {
    pub_count_col <- "CloneCount"
    pub_freq_col <- "CloneFreqInSample"
    pub_vgene_col <- "VGene"
    pub_dgene_col <- "DGene"
    pub_jgene_col <- "JGene"
    pub_cdr3length_col <- "CDR3Length"
  }

  pub_clusters <- buildRepSeqNetwork(
    data = data_public_clusters,
    nucleo_col = pub_nucleo_col, amino_col = pub_amino_col,
    count_col = pub_count_col, freq_col = pub_freq_col,
    vgene_col = pub_vgene_col, dgene_col = pub_dgene_col, jgene_col = pub_jgene_col,
    cdr3length_col = pub_cdr3length_col, other_cols = pub_extra_cols,
    clone_seq_type = clone_seq_type, min_seq_length = NULL, drop_chars = NULL,
    aggregate_identical_clones = FALSE, grouping_cols = NULL,
    dist_type = dist_type, edge_dist = edge_dist,
    drop_isolated_nodes = FALSE,
    node_stats = TRUE, stats_to_include = "all", cluster_stats = TRUE,
    plot_title = "Public Clone Network", plot_subtitle = NULL,
    edge_width = edge_width, size_nodes_by = size_nodes_by,
    node_size_limits = node_size_limits,
    custom_size_legend = custom_size_legend,
    color_nodes_by = color_nodes_by, color_scheme = color_scheme,
    custom_color_legend = custom_color_legend,
    return_all = TRUE)

  # Rename public node stats
  names(pub_clusters$node_data)[names(pub_clusters$node_data) == "degree"] <-
    "PublicNetworkDegree"
  names(pub_clusters$node_data)[names(pub_clusters$node_data) == "cluster_id"] <-
    "PublicClusterID"
  names(pub_clusters$node_data)[names(pub_clusters$node_data) == "transitivity"] <-
    "PublicTransitivity"
  names(pub_clusters$node_data)[names(pub_clusters$node_data) == "Closeness"] <-
    "PublicCloseness"
  names(pub_clusters$node_data)[names(pub_clusters$node_data) == "centrality_by_closeness"] <-
    "PublicCentralityByCloseness"
  names(pub_clusters$node_data)[names(pub_clusters$node_data) == "eigen_centrality"] <-
    "PublicEigenCentrality"
  names(pub_clusters$node_data)[names(pub_clusters$node_data) == "centrality_by_eigen"] <-
    "PublicCentralityByEigen"
  names(pub_clusters$node_data)[names(pub_clusters$node_data) == "betweenness"] <-
    "PublicBetweenness"
  names(pub_clusters$node_data)[names(pub_clusters$node_data) == "centrality_by_betweenness"] <-
    "PublicCentralityByBetweenness"
  names(pub_clusters$node_data)[names(pub_clusters$node_data) == "authority_score"] <-
    "PublicAuthorityScore"
  names(pub_clusters$node_data)[names(pub_clusters$node_data) == "coreness"] <-
    "PublicCoreness"
  names(pub_clusters$node_data)[names(pub_clusters$node_data) == "page_rank"] <-
    "PublicPageRank"

  ## Save node & cluster data & plots ##
  cat(paste0("Saving output for full network of public clusters to:\n  ",
             output_dir, "\n"))
  utils::write.csv(
    pub_clusters$node_data,
    file = file.path(output_dir, "public_cluster_network_node_level_meta_data.csv"),
    row.names = FALSE)
  cat("Node-level meta data saved as 'public_cluster_network_node_level_meta_data.csv'\n")
  utils::write.csv(pub_clusters$cluster_stats,
                   file = file.path(output_dir, "public_cluster_network_cluster_info.csv"),
                   row.names = FALSE)
  cat("Cluster meta data saved as 'public_cluster_network_cluster_info.csv'\n")
  grDevices::pdf(
    file = file.path(output_dir, "public_cluster_network_graph_plot.pdf"),
    width = plot_width, height = plot_height)
  for (j in 1:length(pub_clusters$plots)) { print(pub_clusters$plots[[j]]) }
  grDevices::dev.off()
  cat("Network graph plot saved as 'public_cluster_network_graph_plot.pdf'\n")



  #### ENCODE CLONES BY ATCHLEY FACTOR AND PERFORM K-MEANS CLUSTERING ####
  cat("Embedding the TCR sequences from the public clusters in Euclidean space to perform K-means clustering...\n")
  kmeansAtchley(
    pub_clusters$node_data,
    amino_col = "AminoAcidSeq",
    sample_col = "SampleID",
    group_col = group_col,
    k = 100,
    plot_width = k_plot_width,
    plot_height = k_plot_height,
    margin_size = k_plot_margin,
    use_viridis = k_plot_viridis,
    output_dir = output_dir,
    outfile_heatmap =
      "public_cluster_network_atchley_kmeans_relative_clust_sizes.pdf",
    outfile_corr_heatmap =
      "public_cluster_network_atchley_kmeans_corr_in_relative_clust_sizes.pdf",
    return_output = FALSE)

  cat("All tasks complete.\n")

}
