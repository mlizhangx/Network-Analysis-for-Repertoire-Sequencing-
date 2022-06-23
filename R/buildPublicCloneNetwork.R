# Inputs:
#       rep-seq data from multiple samples (separate files, same directory)
# Do:
#       Identify public clones and build public clone network

buildPublicCloneNetwork <- function(

  # Input Data and Columns
  input_dir = getwd(),
  file_list,
  sample_id_list = seq_along(file_list),
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
  # sample_col,
  other_cols = NULL,

  # Clone sequence settings
  clone_seq_type = "amino acid",
  min_seq_length = 3, # min clone seq length
  drop_chars = "[*|_]",
  aggregate_identical_clones = FALSE,

  # Network Settings
  dist_type = "hamming", # options are "hamming", "levenshtein", "euclidean_on_atchley"
  edge_dist = 1,
  # node_stats = TRUE,
  # stats_to_include = node_stat_settings(), # for final public clone network
  # cluster_stats = FALSE,
  include_atchley_embedding = FALSE, # only applicable to TCRB CDR3 amino acid seqs

  # Filter pass settings for sample-level clusters
  top_n_clusters = 20,
  min_node_count = 10,
  min_clone_count = 100,

  # Downstream filter settings for public clones
  bayes_factor_col = NULL, # column containing adjusted bayes factor pvalues
  bayes_factor_cutoff = 0.05,
  diff_test_col = NULL, # column containing pvalues from differential testing
  diff_test_cutoff = 0.05,

  # Plot Settings (public clone network)
  custom_title = NULL, #
  custom_subtitle = NULL,
  edge_width = 0.1,
  size_nodes_by = 0.5, # can use a column name of data (a numeric value yields fixed node sizes)
  node_size_limits = NULL, # numeric length 2
  custom_size_legend = NULL, # custom legend title
  color_nodes_by = "GlobalClusterID", # accepts multiple values (one plot per value)
  color_scheme = "default", # passed to plotNetworkGraph(); accepts multiple values (one per value of color_nodes_by)
  custom_color_legend = NULL, # custom title (length must match color_nodes_by)

  # Plot Settings (public cluster-level network)
  # single_cluster_plots = TRUE,
  cluster_edge_width = 0.3,
  cluster_size_nodes_by = count_col,
  cluster_node_size_limits = NULL,
  cluster_custom_size_legend = NULL, # custom legend title
  cluster_color_nodes_by = "ClusterNetworkClusterID", # accepts multiple values (one plot per value)
  cluster_color_scheme = "default", # passed to plotNetworkGraph(); accepts multiple values (one per value of color_nodes_by)
  cluster_custom_color_legend = NULL, # custom title (length must match color_nodes_by)

  # Plot Settings (sample-level networks)
  # single_cluster_plots = TRUE,
  sample_edge_width = 0.3,
  sample_size_nodes_by = count_col,
  sample_node_size_limits = NULL,
  sample_custom_size_legend = NULL, # custom legend title
  sample_color_nodes_by = "cluster_id", # accepts multiple values (one plot per value)
  sample_color_scheme = "default", # passed to plotNetworkGraph(); accepts multiple values (one per value of color_nodes_by)
  sample_custom_color_legend = "sample-level cluster ID", # custom title (length must match color_nodes_by)

  # Output Settings
  output_dir = file.path(getwd(), "public_clones_output"),
  plot_width = 12, # passed to pdf()
  plot_height = 10 # passed to pdf()
  # return_all = FALSE # should function return a list, or just the data?

) {


  ### INPUT CHECKS ###
  # Atchley factor embedding only applicable to amino acid sequences
  if (dist_type == "euclidean_on_atchley" & clone_seq_type != "amino acid") {
    stop("distance type 'euclidean_on_atchley' only applicable to amino acid sequences") }

  # each variable for size/color exists or will be added to data
  # need to account for aggregate_identical_clones if TRUE

  # Load data for first sample to validate column names
  data <-
    ifelse(csv_files,
           read.csv(file.path(input_dir, file_list[[1]]), header, sep),
           read.table(file.path(input_dir, file_list[[1]]), header, sep))

  # Check that each input column is a distinct col of data and meets specs
  #
  #

  # Format the input data
  extra_cols <- intersect(
    unique(c(other_cols, bayes_factor_col, diff_test_col,
             color_nodes_by, cluster_color_nodes_by, sample_color_nodes_by)),
    names(data))

  data <-
    data[ , # Keep only the relevant columns, in specified order:
          unique(c(nucleo_col, amino_col, count_col, freq_col,
                   vgene_col, dgene_col, jgene_col, cdr3length_col, extra_cols))]
  data$SampleID <- sample_id_list[[i]]

  #### PREPARE WORKING ENVIRONMENT ####
  # Create output directory if applicable
  if (!is.null(output_dir)) { .createOutputDir(output_dir) }


  # Convert input columns to character if not already
  if (is.numeric(nucleo_col)) { nucleo_col <- names(data)[nucleo_col] }
  if (is.numeric(amino_col)) { amino_col <- names(data)[amino_col] }
  if (is.numeric(count_col)) { count_col <- names(data)[count_col] }
  if (is.numeric(freq_col)) { freq_col <- names(data)[freq_col] }
  if (is.numeric(vgene_col)) { vgene_col <- names(data)[vgene_col] }
  if (is.numeric(dgene_col)) { dgene_col <- names(data)[dgene_col] }
  if (is.numeric(jgene_col)) { jgene_col <- names(data)[jgene_col] }
  if (is.numeric(cdr3length_col)) { cdr3length_col <- names(data)[cdr3length_col] }
  if (is.numeric(other_cols)) { other_cols <- names(data)[other_cols] }
  if (is.numeric(color_nodes_by)) { color_nodes_by <- names(data)[color_nodes_by] }
  if (is.numeric(cluster_color_nodes_by)) { cluster_color_nodes_by <- names(data)[cluster_color_nodes_by] }
  if (is.numeric(sample_color_nodes_by)) { sample_color_nodes_by <- names(data)[sample_color_nodes_by] }

  # Designate amino acid or nucleotide for clone sequence
  clone_seq_col <- amino_col
  if (clone_seq_type == "nucleotide") { clone_seq_col <- nucleo_col }

  # Initialize output directory and objects
  data_cluster_network <- sample_net <- NULL #init


  #### BUILD SAMPLE-LEVEL NETWORKS ####
  # Iterate over the selected clones
  for (i in 1:length(file_list)) {

    cat(paste0("Loading data for sample ", i, ": ", sample_id_list[[i]], "\n"))
    # Load data
    data <-
      ifelse(csv_files,
             read.csv(file.path(input_dir, file_list[[i]]), header, sep),
             read.table(file.path(input_dir, file_list[[i]]), header, sep))

    # Check that each input column is a distinct col of data and meets specs
    # Only need to perform this check for i = 1
    #

    # Format the input data
    extra_cols <- intersect(
      unique(c(other_cols, bayes_factor_col, diff_test_col,
               color_nodes_by, cluster_color_nodes_by, sample_color_nodes_by)),
      names(data))
    data <-
      data[ , # Keep only the relevant columns, in specified order:
            unique(c(nucleo_col, amino_col, count_col, freq_col,
                     vgene_col, dgene_col, jgene_col, cdr3length_col, extra_cols))]
    data$SampleID <- sample_id_list[[i]]

    # Rename frequency column
    old_freq_colname <- freq_col
    old_sample_color_nodes_by <- sample_color_nodes_by
    old_sample_size_nodes_by <- sample_size_nodes_by
    new_freq_colname <- "CloneFreqInSample"
    # if (size_nodes_by == freq_col) { size_nodes_by <- new_freq_colname }
    # if (cluster_size_nodes_by == freq_col) { cluster_size_nodes_by <- new_freq_colname }
    if (sample_size_nodes_by == freq_col) { sample_size_nodes_by <- new_freq_colname }
    # if (freq_col %in% color_nodes_by) {
    #   color_nodes_by[color_nodes_by == freq_col] <- new_freq_colname
    # }
    # if (freq_col %in% cluster_color_nodes_by) {
    #   cluster_color_nodes_by[cluster_color_nodes_by == freq_col] <- new_freq_colname
    # }
    if (freq_col %in% sample_color_nodes_by) {
      sample_color_nodes_by[sample_color_nodes_by == freq_col] <- new_freq_colname
    }
    names(data)[names(data) == freq_col] <- new_freq_colname
    freq_col <- new_freq_colname


    ### BUILD SAMPLE NETWORK ###
    cat(paste0("Building network for sample ", i, ": ", sample_id_list[[i]], "\n"))
    sample_net <- buildRepSeqNetwork(
      data, nucleo_col, amino_col, count_col, freq_col, vgene_col, dgene_col,
      jgene_col, cdr3length_col, c(extra_cols, "SampleID"),
      clone_seq_type, min_seq_length, drop_chars, aggregate_identical_clones,
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


    ### FILTER CLUSTER STATS ###
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
    ids_node_count_gt_n <-
      which(sample_net$cluster_stats$node_count >= min_node_count)
    ids_clone_count_gt_n <-
      which(sample_net$cluster_stats$agg_clone_count >= min_clone_count)
    filtered_ids <-
      unique(c(ids_top_n_clusters, ids_node_count_gt_n, ids_clone_count_gt_n))
    sample_net$cluster_stats <- sample_net$cluster_stats[filtered_ids, ]

    # Add cluster stats to combined data for public cluster network
    data_cluster_network <-
      ifelse(i == 1,
             yes = sample_net$cluster_stats,
             no = rbind(data_cluster_network, sample_net$cluster_stats))

    # Reset freq_col
    freq_col <- old_freq_colname
    sample_color_nodes_by <- old_sample_color_nodes_by
    sample_size_nodes_by <- old_sample_size_nodes_by

  }

} # done looping over selected clones
cat("Done building sample-level networks.\n")







# Format additional variables in data
data_public_clones$AssocClusterID <- as.factor(data_public_clones$AssocClusterID)

#### BUILD GLOBAL CLUSTER NETWORK ####
# Ensure cluster ID is computed
if (!node_stats) {
  stats_to_include <- "cluster_id_only"
} else if (!stats_to_include$cluster_id) {
  stats_to_include$cluster_id <- TRUE }

if ("GlobalClusterID" %in% color_nodes_by) {
  color_nodes_by[which(color_nodes_by == "GlobalClusterID")] <- "cluster_id"
}

cat("Building global cluster network using combined cluster data:\n")
global_net <- buildRepSeqNetwork(
  data_public_clones,
  nucleo_col, amino_col, count_col, freq_col, vgene_col, dgene_col, jgene_col,
  cdr3length_col,
  other_cols = c(sample_col, extra_cols, "AssocClusterID", "AssocClusterSeq",
                 "DegreeInAssocCluster"),
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

# Rename some columns of combined cluster data
names(global_net$node_data)[
  which(names(global_net$node_data) == "CloneFrequency")] <- new_freq_colname
names(global_net$node_data)[
  which(names(global_net$node_data) == sample_col)] <- "SampleID"
names(global_net$node_data)[
  which(names(global_net$node_data) == "cluster_id")] <- "GlobalClusterID"
if ("degree" %in% names(global_net$node_data)) {
  names(global_net$node_data)[which(names(global_net$node_data) == "degree")] <-
    "globalDegree" }
# colnames(global_net$node_data)[1:9] <- c(
#   "NucleotideSeq", "AminoAcidSeq", "CloneCount", "CloneFreqInSample",
#   "VGene", "DGene", "JGene", "CDR3Length", "SampleID")


#### SAVE RESULTS ####
if (!is.null(output_dir)) {
  # Save data for global cluster network if applicable
  if (!is.null(data_outfile)) {
    utils::write.csv(global_net$node_data, file.path(output_dir, data_outfile),
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
  return(global_net$node_data)

} else {
  out <- list("data" = global_net$node_data)
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