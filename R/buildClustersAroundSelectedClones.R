# Inputs:
#       Combined TCR rep-seq data from multiple samples
#       List of selected clone sequences (CDR3 nucleotide or amino acid)
# Do:
#       Build a cluster around each of the selected clones
#       Print the plot for each cluster and return the combined cluster data

buildClustersAroundSelectedClones <- function(

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
  clone_seq_type = "amino_acid",
  dist_type = "hamming", # options are "hamming", "levenshtein", "euclidean_on_atchley"
  cluster_radius = 1,
  edge_dist = 1,
  color_nodes_by = sample_col, # accepts multiple values (one plot per value)
  color_scheme = "default", # passed to plotNetworkGraph(); accepts multiple values (one per value of color_nodes_by)
  long_captions = FALSE, # should plot subtitles include details on cluster settings like dist_type and edge_dist?
  output_dir = NULL,
  save_plots = FALSE,
  return_plots = FALSE # should function return a list of dataframe + ggplots, or just print/write plots and return the dataframe?

) {
  # Create output directory if applicable
  if (!is.null(output_dir)) { .createOutputDir(output_dir) }

  ### INPUT CHECKS ###
  # Atchley factor embedding only applicable to amino acid sequences
  if (dist_type == "euclidean_on_atchley" & clone_seq_type != "amino_acid") {
    stop("distance type 'euclidean_on_atchley' only applicable to amino acid sequences") }

  ### PREPARE WORKING ENVIRONMENT ###
  # Convert input columns to character if not already
  if (is.numeric(color_nodes_by)) { color_nodes_by <- names(data)[color_nodes_by] }
  if (is.numeric(nucleo_col)) { nucleo_col <- names(data)[nucleo_col] }
  if (is.numeric(amino_col)) { amino_col <- names(data)[amino_col] }
  if (is.numeric(count_col)) { count_col <- names(data)[count_col] }
  if (is.numeric(freq_col)) { freq_col <- names(data)[freq_col] }
  if (is.numeric(vgene_col)) { vgene_col <- names(data)[vgene_col] }
  if (is.numeric(dgene_col)) { dgene_col <- names(data)[dgene_col] }
  if (is.numeric(jgene_col)) { jgene_col <- names(data)[jgene_col] }
  if (is.numeric(sample_col)) { sample_col <- names(data)[sample_col] }
  if (is.numeric(other_cols)) { other_cols <- names(data)[other_cols] }

  # Add columns in color_nodes_by to other_cols if not present
  other_cols <- unique(c(other_cols, color_nodes_by))

  # Format the input data
  data <-
    data[ , # Keep only the relevant columns:
          c(nucleo_col, amino_col, count_col, freq_col, vgene_col, dgene_col,
            jgene_col, cdr3length_col, sample_col, other_cols)]

  # Designate amino acid or nucleotide for clone sequence
  clone_col <- amino_col
  if (clone_seq_type %in% c("nucleo", "nucleotide")) { clone_col <- nucleo_col }

  # Determine distance type for cluster radius
  cluster_radius_dist_type <- "levenshtein" # applies to lev and atchley
  if (dist_type == "hamming") { cluster_radius_dist_type <- "hamming" }

  # Initialize output directory and objects
  if (!is.null(output_dir)) { .createOutputDir(output_dir) }
  data_all_clusters <- as.data.frame(matrix(nrow = 0, ncol = ncol(data) + 3))
  names(data_all_clusters) <- c(colnames(data), "assocClustID",
                                "assocClustSeq", "degreeInAssocClust")
  if (return_plots) { plots <- list() }

  ### ITERATE OVER SELECTED CLONES ###
  for (i in 1:length(selected_clones)) {

    # Get cluster data
    cat(paste0("Gathering cluster data for clone sequence ", i, " (", selected_clones[[i]], "):\n"))
    data_current_cluster <-
      getSimilarClones(selected_clones[[i]], data, clone_col, sample_col,
                       cluster_radius_dist_type, max_dist = cluster_radius)

    # Build cluster network
    network <-
      generateNetworkFromClones(data_current_cluster[ , clone_col],
                                dist_type = dist_type,
                                edge_dist = edge_dist,
                                contig_ids = rownames(data_current_cluster))

    # Add variables for cluster ID, central sequence, and degree in cluster
    data_current_cluster$assocClustID <- i
    data_current_cluster$assocClustSeq <- selected_clones[[i]]
    data_current_cluster$degreeInAssocClust <- igraph::degree(network)

    # Add data for current cluster to combined cluster data
    data_all_clusters <- rbind(data_all_clusters, data_current_cluster)

    # Create labels for plots
    plot_title <- paste0("Cluster ", i, " (", selected_clones[[i]], ")")
    plot_subtitle <- NULL
    if (!is.null(selected_clone_labels)) {
      plot_subtitle <- paste0(selected_clone_labels[[i]], "\n") }
    if (long_captions) {
      plot_subtitle <- paste0(plot_subtitle, "\nCluster includes clone sequences with a maximum ", cluster_radius_dist_type, " distance of ", cluster_radius, " from the central sequence")
      if (dist_type == "euclidean_on_atchley") {
        plot_subtitle <- paste0(plot_subtitle, "\nClone sequences embedded in Euclidean 30-space based on Atchley factor representation using deep learning\nEdges based on a maximum Euclidean distance of ", edge_dist, " between embedded values")
      } else {
        plot_subtitle <- paste0(plot_subtitle, "\nEdges based on a maximum ", dist_type, " distance of ", edge_dist, "") }
    }

    # Generate plots of network graph
    if (length(color_nodes_by) > 1 & length(color_scheme) == 1) {
      color_scheme <- rep(color_scheme, length(color_nodes_by)) }
    temp_plotlist <- list()
    for (j in 1:length(color_nodes_by)) {
      cat(paste0("Creating cluster graph with nodes colored by ",
                 color_nodes_by[[j]], "..."))
      newplot <-
        plotNetworkGraph(
          network, title = plot_title,
          subtitle = paste0("Nodes colored by ", color_nodes_by[[j]]),
          color_nodes_by = data_current_cluster[ , color_nodes_by[[j]]],
          size_nodes_by = data_current_cluster[ , freq_col],
          color_legend_title = color_nodes_by[[j]],
          size_legend_title = "freq in own sample",
          color_scheme = color_scheme[[j]])
      print(newplot) # print to R
      temp_plotlist$newplot <- newplot
      names(temp_plotlist)[[length(names(temp_plotlist))]] <- color_nodes_by[[j]]
      cat("Done.\n")
    }
    if (save_plots & !is.null(output_dir)) {
      grDevices::pdf(file = file.path(output_dir,
                                      paste0("cluster_", i, ".pdf")))
      for (j in 1:length(color_nodes_by)) { print(temp_plotlist[[j]]) }
      grDevices::dev.off()
    }
    if (return_plots) {
      plots$newcluster <- temp_plotlist
      names(plots)[[length(names(plots))]] <- selected_clones[[i]]
    }

  } # done looping over selected clones
  cat("All clusters complete.\n")

  # Rename some columns of combined cluster data
  colnames(data_all_clusters)[1:9] <- c(
    "nucleotideSeq", "aminoAcidSeq", "cloneCount", "cloneFrequency",
    "VGene", "DGene", "JGene", "CDR3Length", "sampleID")

  # Save combined cluster data if applicable
  if (!is.null(output_dir)) {
    utils::write.csv(data_all_clusters,
                     file.path(output_dir, "data_all_clusters.csv"))
    cat(paste0(
      "Cluster data saved to file:\n",
      file.path(output_dir, "data_all_clusters.csv"), "\n")) }

  # Return output
  if (return_plots) {
    outlist$cluster_data <- data_all_clusters
    cat("Returning plots and combined cluster data.\n")
    return(outlist)
  } else {
    cat("Returning combined cluster data.\n")
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