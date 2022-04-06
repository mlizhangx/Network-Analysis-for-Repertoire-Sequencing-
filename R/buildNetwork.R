# Main function -----------------------------------------------------------

# Input bulk rep-seq clonotype data; construct repertoire network based on
# desired distance metric, with node and cluster level network characteristics
buildNetwork <- function(clonotypes, counts, frequencies,
                         group_labels = NULL,
                         distance_type = "levenshtein",
                         clone_seq_type = "unspecified clonotype",
                         sample_name = "unknown_sample",
                         group_type = "no group labels provided",
                         save_output = FALSE,
                         report_time = FALSE) {

  ### INITIALIZE ###

  # check for invalid arguments
  .checkArguments(clonotypes, counts, frequencies, group_labels,
                  sample_name, distance_type, clone_seq_type)

  # create output directory
  if (save_output) {
    dir_output <- .createOutputDir(distance_type, clone_seq_type, sample_name)
  }

  # record start time
  if (report_time) start_time <- Sys.time()



  ### PROCESS CLONOTYPE DATA ###
  meta_data <- .processClonotypeData(clonotypes, counts, frequencies,
                                     group_labels)


  ### GENERATE ADJACENCY MATRIX ###
  # call python: generate adjacency matrix for clonotypes with deg > 0
  adjacency_matrix <- adjacencyMatrix(meta_data$cloneSeq, distance_type,
                                      max_dist = 1, sparse = TRUE)

  # clonotypes with deg > 0 (i.e., those corresponding to adjacency matrix)
  clone_ids <- utils::read.table("col_ids.txt")
  meta_data <- meta_data[clone_ids$V1, ]



  ### GENERATE NETWORK AND COMPUTE CHARACTERISTICS ###
  # Set random seed to fix graph output
  set.seed(9999)

  # Generate network from adjacency matrix
  net <- .genNetworkGraph(adjacency_matrix)

  # Compute node-level network characteristics
  cell_level_info <- .genNetworkStats(net, meta_data)

  # Compute cluster-level network characteristics
  cluster_level_info <- .genClusterStats(net, cell_level_info, adjacency_matrix)



  ### VISUAL PLOT OF NETWORK GRAPH ###
  # Create ggraph of network
  if (!is.null(group_labels)) { # color nodes by group label if present
    network_plot <-
      .genGraphPlot(net, cell_level_info$cloneCount,
                    sample_name, distance_type, clone_seq_type,
                    groups = as.character(cell_level_info$groupLabel),
                    group_legend_title = group_type)
  } else { # default to coloring nodes by network degree
    network_plot <-
      .genGraphPlot(net, cell_level_info$cloneCount,
                    sample_name, distance_type, clone_seq_type,
                    groups = as.character(cell_level_info$deg),
                    group_legend_title = "Network Degree")
  }

  # Display graph
  print(network_plot)



  ### SAVE AND RETURN RESULTS ###
  # Write results to disk
  if (save_output) {
    .saveResults(cell_level_info, cluster_level_info, net, network_plot,
                 dir_output)
  }

  # Remove temporary files
  .removeTempFiles()

  # Report completion and elapsed time
  cat("All processes complete.")
  if (report_time) {
    end_time <- Sys.time()
    cat(paste0(" Total elapsed time: ",
               round(difftime(end_time, start_time, units = "secs"), 2),
               " seconds."))
  }
  cat("\n")

  # Return results as list
  return(list("settings" = list("sample_name" = sample_name,
                                "clone_seq_type" = clone_seq_type,
                                "distance_type" = distance_type,
                                "group_labels" = group_type),
              "igraph" = net,
              "cell_info" = cell_level_info,
              "cluster_info" = cluster_level_info,
              "ggraph" = network_plot))

}






# Helper functions for buildNetwork ---------------------------------------





# create output directory named according to sample ID and distance type
# directory created as subfolder of current working directory
.createOutputDir <- function(distance_type, clone_seq_type, sample_name) {

  dir_output <- paste0(getwd(), '/', sample_name, '/', clone_seq_type, '/',
                       distance_type)
  dir.create(dir_output, showWarnings = FALSE, recursive = TRUE)

  # Confirm successful creation of output directory
  if (!dir.exists(dir_output)) {
    stop(paste0("Unable to create directory ", dir_output,
                ". Check that 'sample_name' does not contain characters which are invalid for file and directory names"))
  }

  return(dir_output)

}




# input arguments to buildNetwork other than 'data'
# check argument types for invalid inputs
.checkArguments <- function(clonotypes, counts, frequencies, group_labels,
                            sample_name, distance_type, clone_seq_type) {

  # clonotypes is a character vector
  if (length(clonotypes) == 0) stop("argument 'clonotypes' has zero length")
  if (!is.vector(clonotypes) | !is.character(clonotypes)) {
    stop(paste0("argument 'clonotypes' must be a character vector"))
  }

  # counts is numeric character vector with no NA/Inf values and correct length
  if (length(counts) == 0) stop("argument 'counts' has zero length")
  if (!is.vector(counts) | !is.numeric(counts)) {
    stop(paste0("argument 'counts' must be a numeric vector"))
  }
  if (sum(is.na(counts)) != 0) stop("'counts' contains NA/NaN values")
  if (sum(is.finite(counts)) != length(counts)) {
    stop("'counts' contains Inf/-Inf values")
  }
  if (length(counts) != length(clonotypes)) {
    stop("length of 'counts' does not match length of 'clonotypes'")
  }

  # frequency is numeric character vector with no NA/Inf values and correct length
  if (length(frequencies) == 0) stop("argument 'frequencies' has zero length")
  if (!is.vector(frequencies) | !is.numeric(frequencies)) {
    stop(paste0("argument 'frequencies' must be a numeric vector"))
  }
  if (sum(is.finite(frequencies)) != length(frequencies)) {
    stop("'frequencies' contains Inf/-Inf values")
  }
  if (length(frequencies) != length(clonotypes)) {
    stop("length of 'frequencies' does not match length of 'clonotypes'")
  }

  # group_labels is null or vector of correct length, coercible to factor with > 1 level
  if (!is.null(group_labels)) {
    if (length(group_labels) == 0) stop("argument 'group_labels' has zero length; if you are trying to omit the group variable, use the default value of `group_labels = NULL`. Otherwise, 'group_labels' must be a vector of the same length as 'clonotypes' containing at least two distinct values")
    if (!is.vector(group_labels) ) stop("'group_labels' must be a vector or 'NULL'")
    if (length(group_labels) != length(clonotypes)) stop("length of 'group_labels' does not match length of 'clonotypes'")
    if (sum(is.na(group_labels)) != 0) stop("'group_labels' contains NA/NaN values")
    if (length(levels(as.factor(group_labels))) < 2) stop("'group_labels' must contain at least two distinct values")
  }

  # Valid distance_type arguments
  programmed_distance_types <- c("levenshtein", "hamming")

  # Check distance_type argument
  if (length(distance_type) != 1) {
    stop(paste0("'distance_type' must be one of the following: ",
                programmed_distance_types))
  }
  if (!distance_type %in% programmed_distance_types) {
    stop(paste0("'distance_type' must be one of the following: ",
                programmed_distance_types))
  }

  # Check sample_name argument
  if (length(sample_name) != 1) {
    stop("argument 'sample_name' must be a character string")
  }
  if (!is.character(sample_name)) {
    stop("argument 'sample_name' must be a character string")
  }

  # Check clone_type argument
  if (length(clone_seq_type) != 1) {
    stop("argument 'clone_seq_type' must be a character string")
  }
  if (!is.character(clone_seq_type)) {
    stop("argument 'clone_seq_type' must be a character string")
  }


}




# Input bulk rep-seq data; return cleaned and formatted data for unique
# nontrivial clonotypes; write list of unique nontrivial clonotypes to disk
.processClonotypeData <- function(clonotypes, counts, frequencies,
                                  group_labels) {

  cat("Processing repertoire-sequence data...\n")

  # Does the data include group labels?
  groups <- !is.null(group_labels)

  # aggregate the data for clonotype sequences with duplicate entries
  clone_seq_agg <- stats::aggregate(data.frame(counts, frequencies),
                                    by = list(clonotypes),
                                    FUN = sum)
  colnames(clone_seq_agg) <- c("cloneSeq", "cloneCount", "cloneFraction")

  if (groups) { # Assign group with greatest count to each clone seq

    # Initialize variable 'groupLabel' for aggregate data
    # (ensures type is same as that of 'group_labels';
    # initial values are unimportant)
    clone_seq_agg$groupLabel <- group_labels[1:nrow(clone_seq_agg)]

    # iterate over unique clonotype sequences
    for (i in 1:nrow(clone_seq_agg)) {
      # which rows of original data correspond to current sequence?
      current_seq_indices <- which(clonotypes == clone_seq_agg$cloneSeq[[i]])

      # get group labels and counts from rows in previous step
      labels_sub <- group_labels[current_seq_indices]
      counts_sub <- counts[current_seq_indices]

      # aggregate counts for common group label values
      agg_counts_by_group <- stats::aggregate(data.frame(counts_sub),
                                              by = list(labels_sub),
                                              FUN = sum)
      colnames(agg_counts_by_group) <- c("label", "count")

      # are there nonempty group labels for current sequence?
      nonempty_groups <- sum(agg_counts_by_group$label != "") > 0

      if (nonempty_groups) {
        # maximum count among nonempty group labels for current sequence
        max_count <-
          max(agg_counts_by_group$count[agg_counts_by_group$label != ""])

        # index for corresponding group label
        ind_maxcount <-
          which(agg_counts_by_group$count == max_count &
                  agg_counts_by_group$label != "" )[[1]]

      } else {
        ind_maxcount <- which.max(agg_counts_by_group$count)
      }

      # assign group label with greatest count to current sequence
      clone_seq_agg$groupLabel[[i]] <- agg_counts_by_group$label[[ind_maxcount]]
    }

  }

  # for each sequence, count # of rows for that sequence in input data
  unique_clone_seq <- as.data.frame(table(clonotypes))
  colnames(unique_clone_seq) <- c("cloneSeq", "uniqueCount")

  # merge aggregate sequence data with unique counts
  meta_data <- merge(clone_seq_agg, unique_clone_seq, by = "cloneSeq")

  # remove any remaining duplicate sequences (original comment: count if not AA count but one nucleotide count)
  meta_data <- meta_data[!duplicated(meta_data$cloneSeq), ]

  cat(paste0(nrow(meta_data), " unique clonotype sequences found.\n"))

  # remove sequences with length < 3 (empty "" values are included)
  meta_data$seq_length <- nchar(meta_data$cloneSeq)
  nontrivial_clone_ids <- meta_data$seq_length > 2
  cat(paste0("Removing ", sum(!nontrivial_clone_ids),
             " clonotype sequences containing fewer than 3 characters...\n"))
  meta_data <- meta_data[nontrivial_clone_ids, ]
  cat(paste0(nrow(meta_data), " unique clonotype sequences remaining.\n"))

  # stop if not enough unique clonotype sequences
  if (nrow(meta_data) < 2) {
    stop("Insufficient number of unique non-trivial clonotype sequences to construct network.")
  }

  # Write clonotype sequences to temporary file in working directory
  utils::write.table(meta_data$cloneSeq,
                     file = "temp_clone_seq_list.txt",
                     quote = FALSE,  col.names = FALSE,
                     row.names = FALSE, sep = "\t")

  # Return clonotype data for unique nontrivial clonotype sequences
  return(meta_data)

}




# Use adjacency matrix and corresponding clonotype data to generate network
.genNetworkGraph <- function(adjacency_matrix) {

  cat("Generating network from adjacency matrix...\n")

  net <- igraph::graph_from_adjacency_matrix(adjacency_matrix, weighted = TRUE)

  net <- igraph::as.undirected(igraph::simplify(net,
                                                remove.multiple = T,
                                                remove.loops = T))

  return(net)

}




# input network and corresponding metadata;
# augment metadata with node-level network info and return
.genNetworkStats <- function(net, meta_data) {

  cat(paste0("Computing node-level network characteristics for the ",
             nrow(meta_data), " clonotypes (nodes) in the network...\n"))

  # add degree to data
  meta_data$deg <- igraph::degree(net)

  # add cluster_id to data
  # cfg <- igraph::cluster_fast_greedy(igraph::as.undirected(net)) # as.undirected redundant?
  cfg <- igraph::cluster_fast_greedy(net)
  meta_data$cluster_id <- cfg$membership

  meta_data$transitivity <- igraph::transitivity(net, type = "local")

  # meta_data$closeness <- closeness(net, mode="all", weights=NA)
  # meta_data$centr_clo_res <- igraph::centr_clo(net, mode="all", normalized=T)$res

  meta_data$eigen_centrality <-
    igraph::eigen_centrality(net, directed = T, weights = NA)$vector

  meta_data$centr_eigen <-
    igraph::centr_eigen(net, directed = T, normalized = T)$vector

  meta_data$betweenness <-
    igraph::betweenness(net, directed = T, weights = NA)

  meta_data$centr_betw <-
    igraph::centr_betw(net, directed = T, normalized = T)$res

  meta_data$authority_score <-
    igraph::authority_score(net, weights = NA)$vector

  meta_data$coreness <-
    igraph::coreness(net, mode = "all")

  meta_data$page_rank <-
    igraph::page_rank(net)$vector

  return(meta_data)

}




# input network and corresponding metadata with node-level stats;
# compute cluster-level network info and return
.genClusterStats <- function(net, cell_level_info, adjacency_matrix) {


  # Tabulate the number of cells in each cluster
  out <- as.data.frame(table(cell_level_info$cluster_id))
  colnames(out) <- c("cluster_id", "node_count")

  num_clusters <- nrow(out) # Total number of clusters


  cat(paste0(num_clusters, " clusters in network. Computing cluster-level network characteristics...\n"))

  ### INITIALIZE VALUES ###
  # out$motif_top_deg_alpha <- ""
  out$motif_top_deg_beta <- ""
  out$max_deg_within_cluster <- ""
  out$deg_mean <- 0
  out$motif_top_count_beta <- ""
  out$max_count_within_cluster <- ""
  out$cluster_total_count <- 0
  # out$motif_freq_50_alpha <- ""
  # out$motif_freq_50_beta <- ""

  out$betaCDR3_length <- 0
  # out$alphaCDR3_length <- 0

  # out$Count_PRE_INFUSION <- 0
  # out$Count_DOSE_2 <- 0

  # social network properites
  # out$deg_avg <- 0
  out$diam_length <- 0
  out$assortativity <- 0
  out$transitivity <- 0
  out$edge_density <- 0
  out$centr_degree <- 0
  out$centr_clo <- 0
  out$eigen_centrality <- 0
  out$centr_eigen <- 0

  ### COMPUTE STATS FOR EACH CLUSTER ###
  for(membership_id in 1:num_clusters) {

    # Row of cluster_level_info corresponding to current cluster
    cluster_index <- which(out$cluster_id == membership_id)

    # Rows of cell_level_info corresponding to cells in current cluster
    cell_ids <- cell_level_info$cluster_id == membership_id

    # Mean degree in cluster
    out[cluster_index, ]$deg_mean <-
      round(mean(cell_level_info[cell_ids, ]$deg), 2)

    # Mean sequence length in cluster
    out[cluster_index, ]$betaCDR3_length  <-
      round(mean(cell_level_info[cell_ids, ]$seq_length),2)
    # out[cluster_index,]$alphaCDR3_length  <- round(mean(cell_level_info[cell_ids,]$alphaCDR3_length),2)

    # Total aggregate clonotype count in cluster
    out[cluster_index, ]$cluster_total_count  <-
      sum(cell_level_info[cell_ids, ]$cloneCount)
    # out[cluster_index,]$Count_PRE_INFUSION  <- sum(cell_level_info[cell_ids,]$Count_PRE_INFUSION)
    # out[cluster_index,]$Count_DOSE_2  <- sum(cell_level_info[cell_ids,]$Count_DOSE_2)

    # Maximum degree within cluster
    max_deg <- max(cell_level_info[cell_ids, ]$deg)
    out[cluster_index, ]$max_deg_within_cluster <- max_deg

    # Clonotype sequence corresponding to maximum degree within cluster
    cell_id_top_deg_beta <- which(cell_ids & cell_level_info$deg == max_deg)
    out[cluster_index, ]$motif_top_deg_beta  <-
      as.character(cell_level_info[cell_id_top_deg_beta, 1][[1]])
    # out[cluster_index,]$motif_top_deg_alpha  <- as.character(cell_level_info[cell_ids & cell_level_info$deg_mean == max_deg,]$alphaCDR3[1])

    # Maximum clonotype count within cluster
    max_count <- max(cell_level_info[cell_ids, ]$cloneCount)
    out[cluster_index, ]$max_count_within_cluster <- max_count

    # Clonotype sequence with maximum count within cluster
    cell_id_top_count_beta <-
      which(cell_ids & cell_level_info$cloneCount == max_count)
    out[cluster_index, ]$motif_top_count_beta  <-
      as.character(cell_level_info[cell_id_top_count_beta, 1][[1]])

    ### get the representative motif ###
    # take long time to run, comment out
    # betaCDR3_length_round <- round(out[cluster_index,]$betaCDR3_length)
    # for(i in 1:betaCDR3_length_round) {
    #   # cat(i)
    #   string_list <- cell_level_info[cell_ids & cell_level_info$betaCDR3_length == betaCDR3_length_round,]$aminoAcid
    #   freq_table <- as.data.frame(table(substring(string_list, i,i)))
    #   select_letter <- ifelse(dim(freq_table[freq_table$Freq > length(string_list)/2,])[1] == 1, as.character(freq_table[freq_table$Freq > length(string_list)/2,]$Var1), "*")
    #   # cat(select_letter)
    #   select_letter_index <- paste0("char_", i)
    #   assign(select_letter_index, select_letter)
    # }
    # out[cluster_index,]$motif_freq_50_beta  <- paste(noquote(mget(mixedsort(ls(pattern= "char_")))), collapse = "")
    # rm(list = ls(pattern= "char_"))


    # Build cluster network to get network properties for the cluster
    adjacency_matrix_cluster <- as.matrix(adjacency_matrix[cell_ids, cell_ids])
    # # adjacency_matrix[adjacency_matrix > 1] <- 0
    # matrix_2_net <- as.matrix(adjacency_matrix)
    cluster <- igraph::graph_from_adjacency_matrix(adjacency_matrix_cluster)
    cluster <- igraph::as.undirected(igraph::simplify(cluster,
                                                      remove.multiple = T,
                                                      remove.loops = T))
    deg <- igraph::degree(cluster, mode = "all")
    # table(deg) # there should be no deg== 0.
    # out[cluster_index,]$deg_avg <- round(mean(deg),2)


    # Diameter (longest geodesic distance)
    out[cluster_index, ]$diam_length <-
      length(igraph::get_diameter(cluster, directed = T))

    # Assortativity
    out[cluster_index, ]$assortativity <-
      igraph::assortativity_degree(cluster, directed = F)

    # Transitivity
    out[cluster_index, ]$transitivity <-
      igraph::transitivity(cluster, type = "global")  # cluster is treated as an undirected network

    # Density
    # The proportion of present edges from all possible ties.
    out[cluster_index, ]$edge_density <-
      igraph::edge_density(cluster, loops = F)

    # centralization on degree
    out[cluster_index, ]$centr_degree <-
      igraph::centr_degree(cluster, mode = "in", normalized = T)$centralization

    # centralization on Closeness (centrality based on distance to others in the
    # graph)
    out[cluster_index, ]$centr_clo <-
      igraph::centr_clo(cluster, mode = "all", normalized = T)$centralization

    # centralization on Eigenvector (centrality proportional to the sum of
    # connection centralities) Values of the first eigenvector of the graph
    # adjacency matrix
    out[cluster_index, ]$eigen_centrality <-
      igraph::eigen_centrality(cluster, directed = T, weights = NA)$value

    out[cluster_index, ]$centr_eigen <-
      igraph::centr_eigen(cluster, directed = T, normalized = T)$centralization

  }

  # Return cluster-level info
  return(out)
}




# input network and cell-level data; create ggraph for network
.genGraphPlot <- function(net, clone_counts, sample_name, distance_type,
                          clone_seq_type, groups, group_legend_title) {

  cat("Creating visual plot of network graph...\n")

  l <- igraph::layout_components(net)

  # add cluster id in plot
  # membership_label <- cell_level_info$cluster_id
  # membership_label[duplicated(membership_label)] <- NA
  # membership_label <- as.factor(membership_label)

  network_plot <-
    ggraph::ggraph(net, layout = l) +
    ggraph::geom_edge_link0(colour = "grey") +
    ggraph::geom_node_point(
      ggplot2::aes(color = groups,
                   size = clone_counts)) +
    ggplot2::scale_size(range = c(0.1, log(max(clone_counts)) / 2.5)) +
    # theme(legend.position = "none") +
    # geom_node_text(aes(label = membership_label), size = 1) +
    ggraph::theme_graph(base_family = "sans") +
    ggplot2::labs(color = group_legend_title,
                  size = "Clone Count",
                  title = sample_name,
                  subtitle = paste0(distance_type, " distance on ",
                                    clone_seq_type, " sequence"))

  return(network_plot)
}




# save output to disk and remove temporary files
.saveResults <- function(cell_level_info, cluster_level_info, net, network_plot,
                         dir_output) {

  cat(paste0("Saving results to ", dir_output, "\n"))

  # Save Network igraph using edgelist format
  igraph::write_graph(net,
                      file = paste0(dir_output, "/network_graph_edgelist.txt"),
                      format = "edgelist")
  cat("Network igraph object saved in edgelist format as 'network_graph_edgelist.txt'\n")

  # Save cell-level info with network characteristics
  utils::write.csv(cell_level_info,
                   file = paste0(dir_output, "/cell_level_info.csv"),
                   row.names = FALSE)
  cat("Cell-level data and network characteristics saved as 'cell_level_info.csv'\n")

  # Save cluster-level stats
  utils::write.csv(cluster_level_info,
                   file = paste0(dir_output, "/cluster_level_info.csv"),
                   row.names = FALSE)
  cat("Cluster-level network characteristics saved as 'cluster_level_info.csv'\n")

  # Save ggraph
  grDevices::pdf(file = paste0(dir_output, "/graph_plot.pdf"))
  print(network_plot)
  grDevices::dev.off()
  cat("Plot of network graph saved as 'graph_plot.pdf'\n")

}




.removeTempFiles <- function() {

  # Remove temporary files
  file.remove("temp_clone_seq_list.txt")
  file.remove("temp_adjacency_matrix.mtx")
  file.remove("select_col_id.csv")

}
