

# File Loading ------------------------------------------------------------


# load data from file
.loadDataFromFile <- function(
    input_file, input_type = "rds", data_symbol = NULL, header = TRUE, sep = "")
{
  if (input_type == "csv") {
    data <- utils::read.csv(input_file, header, sep)
  } else if (input_type %in% c("table", "tsv", "txt")) {
    data <- utils::read.table(input_file, header, sep)
  } else if (input_type == "rds") {
    data <- readRDS(input_file)
  } else if (input_type %in% c("rda", "Rda", "Rdata", "rdata")) {
    stopifnot("must specify data symbol via argument `data_symbol`" =
                !is.null(data_symbol))
    load(input_file); assign(x = "data", value = get(data_symbol))
  } else {
    data <- utils::read.table(input_file, header, sep)
  }
  return(data)
}

# load and merge multiple samples from file list
loadDataFromFileList <- function(
    file_list, input_type, data_symbols = NULL, header = TRUE, sep = "")
{
  if (input_type == "csv") {
    data <- plyr::ldply(file_list, utils::read.csv, header = header, sep = sep)
  } else if (input_type %in% c("table", "tsv", "txt")) {
    data <- plyr::ldply(file_list, utils::read.table, header = header, sep = sep)
  } else if (input_type == "rds") {
    data <- plyr::ldply(file_list, readRDS)
  } else if (input_type %in% c("rda", "Rda", "Rdata", "rdata")) {
    stopifnot("must specify data symbols via argument `data_symbols`" =
                !is.null(data_symbols))
    if (length(file_list) > length(data_symbols)) {
      data_symbols <- rep(data_symbols, length.out = length(file_list))
    }
    temp_symbols <- paste0("data_for_sample", 1:length(file_list))
    for (i in 1:length(file_list)) {
      load(file_list[[i]])
      assign(temp_symbols[[i]], get(data_symbols[[i]])); rm(data_symbols[[i]])
    }
    data_list <- mget(temp_symbols); rm(list = temp_symbols)
    data <- do.call(rbind, data_list)
  }
  return(data)
}

# load and compile data from multiple samples, adding variables for
#  sample ID, subject ID, group ID
combineSamples <- function(
    file_list, input_type,
    data_symbols = NULL, header = TRUE, sep = "",
    seq_col,
    min_seq_length = NULL, drop_matches = NULL, subset_cols = NULL,
    sample_ids = NULL, subject_ids = NULL, group_ids = NULL
) {

  cat(">>> Loading and compiling data from all samples:\n")

  if (input_type %in% c("rda", "Rda", "Rdata", "rdata")) {
    stopifnot("must specify data symbols via argument `data_symbols`" =
                !is.null(data_symbols))
    if (length(file_list) > length(data_symbols)) {
      data_symbols <- rep(data_symbols, length.out = length(file_list))
    }
  }

  temp_symbols <- paste0("data_for_sample", 1:length(file_list))

  for (i in 1:length(file_list)) {
    cat(paste0("Loading sample ", i, ": "))
    tmp <- .loadDataFromFile(
      file_list[[i]], input_type, data_symbols[[i]], header, sep)
    if (i == 1) {
      seq_col <- .convertColRef(seq_col, tmp)
      subset_cols <- .convertColRef(subset_cols, tmp)
    }
    tmp <- filterInputData(tmp, seq_col, min_seq_length, drop_matches,
                           subset_cols)
    if (!is.null(sample_ids)) { tmp$SampleID <- sample_ids[[i]] }
    if (!is.null(subject_ids)) { tmp$SubjectID <- subject_ids[[i]] }
    if (!is.null(group_ids)) { tmp$GroupID <- group_ids[[i]] }
    assign(temp_symbols[[i]], tmp); rm(tmp)

  }

  data_list <- mget(temp_symbols); rm(list = temp_symbols)
  data <- do.call(rbind, data_list)
  cat("All samples loaded.\n")
  return(data)

}

# File Saving -------------------------------------------------------------


saveNetwork <- function(
    net, output_dir = getwd(), output_type = "individual",
    output_filename = "MyRepSeqNetwork", pdf_width = 12, pdf_height = 10)
{
  if (is.null(output_dir)) { return(invisible(NULL)) } # exit if no output dir
  .createOutputDir(output_dir)
  if (output_type == "rds") {
    .saveNetworkRDS(net, output_dir, output_filename)
    if ("plots" %in% names(net)) {
      saveNetworkPlots(
        net$plots, file.path(output_dir, paste0(output_filename, ".pdf")),
        pdf_width, pdf_height) }
  } else if (output_type == "individual") {
    .saveNetworkObjects(net, output_dir, output_filename, pdf_width, pdf_height)
  } else {
    .saveNetworkRDA(net, output_dir, output_filename)
    if ("plots" %in% names(net)) {
      saveNetworkPlots(
        net$plots, file.path(output_dir, paste0(output_filename, ".pdf")),
        pdf_width, pdf_height) }
  }
}


saveNetworkPlots <- function(plotlist, outfile = "MyRepSeqNetwork.pdf",
                             pdf_width = 12, pdf_height = 10) {
  grDevices::pdf(file = outfile, width = pdf_width, height = pdf_height)
  for (j in 1:length(plotlist)) { print(plotlist[[j]]) }
  grDevices::dev.off()
  cat(paste0("Network graph plots saved to file:\n  ", outfile, "\n"))
}

.ensureOutputDir <- function(output_dir) {
  stopifnot("output_dir is required" = !is.null(output_dir))
  if (!dir.exists(output_dir)) { .createOutputDir(output_dir) }
  stopifnot("could not create output_dir" = dir.exists(output_dir))
}

.createOutputDir <- function(dirname) {
  if (!is.null(dirname)) {
    if (!dir.exists(dirname)) {
      dir.create(dirname, showWarnings = FALSE, recursive = TRUE)
    }
    if (!dir.exists(dirname)) {
      stop(paste0("Unable to create directory ", dirname,
                  ". Check to confirm that a valid directory name was provided."))
    }
  }
}

.createDirectories <- function(dirs) {
  for (i in 1:length(dirs)) {
    dir.create(dirs[[i]], showWarnings = FALSE, recursive = TRUE)
  }
}


.saveDataGeneric <- function(data, output_dir, output_name, output_type)
{
  if (output_type == "rds") {
    saveRDS(object = data,
            file = file.path(output_dir, paste0(output_name, ".rds")))
  } else if (output_type == "csv") {
    data <- apply(data, MARGIN = 2, FUN = as.character)
    utils::write.csv(data, row.names = FALSE,
                     file = file.path(output_dir, paste0(output_name, ".csv")))
  } else if (output_type %in% c("tsv", "table")) {
    data <- apply(data, MARGIN = 2, FUN = as.character)
    utils::write.table(data, row.names = FALSE,
                       file = file.path(output_dir, paste0(output_name, ".tsv")))
  } else {
    save(data, file = file.path(output_dir, paste0(output_name, ".rda")))
  }
}


.saveNetworkRDA <- function(net, output_dir, output_filename) {
  outfile <- file.path(output_dir, paste0(output_filename, ".rda"))
  save(net, file = outfile)
  cat(paste0("Output saved to file:\n  ", outfile, "\n"))
}


.saveNetworkRDS <- function(net, output_dir, output_filename) {
  outfile <- file.path(output_dir, paste0(output_filename, ".rds"))
  saveRDS(net, file = outfile)
  cat(paste0("Output saved to file:\n  ", outfile, "\n"))
}


.saveNetworkObjects <- function(net, output_dir, output_filename,
                                pdf_width, pdf_height)
{
  # Save node & cluster data
  node_file <-
    file.path(output_dir, paste0(output_filename, "_NodeMetadata.csv"))
  net$node_data <- apply(net$node_data, MARGIN = 2, FUN = as.character)
  utils::write.csv(net$node_data, file = node_file, row.names = FALSE)
  cat(paste0("Node-level meta-data saved to file:\n  ", node_file, "\n"))

  if ("cluster_data" %in% names(net)) {
    cluster_file <-
      file.path(output_dir, paste0(output_filename, "_ClusterMetadata.csv"))
    net$cluster_data <- apply(net$cluster_data, MARGIN = 2, FUN = as.character)
    utils::write.csv(net$cluster_data, file = cluster_file, row.names = FALSE)
    cat(paste0("Cluster-level meta-data saved to file:\n  ", cluster_file, "\n"))
  }

  # Save plots to a single pdf
  if ("plots" %in% names(net)) {
    saveNetworkPlots(
      net$plots, file.path(output_dir, paste0(output_filename, ".pdf")),
      pdf_width, pdf_height)
  }

  # Save igraph
  igraph_outfile <-
    file.path(output_dir, paste0(output_filename, "_EdgeList.txt"))
  igraph::write_graph(net$igraph, file = igraph_outfile, format = "edgelist")
  cat(paste0("Network igraph saved in edgelist format to file:\n  ", igraph_outfile, "\n"))

  # Save adjacency matrix
  if (inherits(net$adjacency_matrix, "matrix")) {
    matrix_outfile <-
      file.path(output_dir, paste0(output_filename, "_AdjacencyMatrix.csv"))
    utils::write.csv(net$adjacency_matrix, matrix_outfile, row.names = FALSE)
    cat(paste0("Adjacency matrix saved to file:\n  ", matrix_outfile, "\n"))
  } else if (inherits(net$adjacency_matrix, "dgCMatrix")) {
    matrix_outfile <-
      file.path(output_dir, paste0(output_filename, "_AdjacencyMatrix.mtx"))
    Matrix::writeMM(net$adjacency_matrix, matrix_outfile)
    cat(paste0("Adjacency matrix saved to file:\n  ", matrix_outfile, "\n"))
    if ("adj_mat_a" %in% names(net)) {
      matrix_outfile <-
        file.path(output_dir, paste0(output_filename, "_AdjacencyMatrix_ChainA.mtx"))
      Matrix::writeMM(net$adj_mat_a, matrix_outfile)
      cat(paste0("Adjacency matrix for first chain saved to file:\n  ", matrix_outfile, "\n"))
      matrix_outfile <-
        file.path(output_dir, paste0(output_filename, "_AdjacencyMatrix_ChainB.mtx"))
      Matrix::writeMM(net$adj_mat_b, matrix_outfile)
      cat(paste0("Adjacency matrix for second chain saved to file:\n  ", matrix_outfile, "\n"))
    }
  }

}




# Filtering and Subsetting ------------------------------------------------

filterInputData <- function(
    data, seq_col, min_seq_length = NULL, drop_matches = NULL,
    subset_cols = NULL) {

  cat(paste0("Input data contains ", nrow(data), " rows.\n"))

  # coerce sequence column(s) to character
  for (i in 1:length(seq_col)) {
    if (!is.character(data[[seq_col[[i]]]])) {
      data[[seq_col[[i]]]] <- as.character(data[[seq_col[[i]]]])
    }
  }

  # filter by min seq length
  if (!is.null(min_seq_length)) {
    cat(paste0("Removing sequences with length fewer than ", min_seq_length, " characters..."))
    drop_rows <- .RowDropsBySequenceLength(data, seq_col, min_seq_length)
    if (sum(drop_rows) > 0) { data <- data[!drop_rows, , drop = FALSE] }
    cat(paste0(" Done. ", nrow(data), " rows remaining.\n"))
  }

  # filter by seq content
  if (!is.null(drop_matches)) {
    cat(paste0("Removing sequences containing matches to the expression '", drop_matches, "'..."))
    drop_rows <- .RowDropsBySequenceContent(data, seq_col, drop_matches)
    if (sum(drop_rows) > 0) { data <- data[!drop_rows, , drop = FALSE] }
    cat(paste0(" Done. ", nrow(data), " rows remaining.\n"))
  }

  # subset data columns
  if (!is.null(subset_cols)) {
    data <- .subsetColumns(data, c(seq_col, subset_cols))
  }

  return(data)
}


getNeighborhood <- function(
    data, seq_col, target_seq, dist_type = "hamming", max_dist = 1)
{
  if (!target_seq %in% data[[seq_col]]) { return(NULL) }
  dist_fun <- hamDistBounded
  if (dist_type == "levenshtein") { dist_fun <- levDistBounded }
  dists_to_targetseq <- sapply(
    X = data[[seq_col]], FUN = dist_fun, b = target_seq, k = max_dist)
  # get data for sequences within the specified radius
  out <- data[dists_to_targetseq != -1, , drop = FALSE]
  return(out)
}

# convert one or more column references from numeric to character
.convertColRef <- function(colrefs, data) {
  if (is.numeric(colrefs)) { colrefs <- names(data)[colrefs] }
  return(colrefs)
}


# FUNCTION: Filter rep-seq data to remove rows for clonotype sequences with
# length below the specified cutoff
.RowDropsBySequenceLength <- function(data, seq_col, min_length = 3) {
  drop_rows <- rep(FALSE, nrow(data))
  for (i in 1:length(seq_col)) {
    drop_rows <- drop_rows | (nchar(data[[seq_col[[i]]]]) < min_length)
  }
  return(drop_rows)
}

# FUNCTION: Filter rep-seq data to remove rows for clonotype sequences with
# length below the specified cutoff
.RowDropsBySequenceContent <- function(data, seq_col, drop_matches) {
  drop_rows <- rep(FALSE, nrow(data))
  for (i in 1:length(seq_col)) {
    drop_rows <- drop_rows | grepl(drop_matches, data[[seq_col[[i]]]])
  }
  return(drop_rows)
}

.processSubsetCols <- function(subset_cols, other_cols) {
  if (is.null(subset_cols)) { return(NULL) }
  return(c(subset_cols, other_cols))
}

.subsetColumns <- function(data, cols_to_keep) {
  cols_to_keep <- intersect(unique(cols_to_keep), names(data))
  if (is.null(cols_to_keep)) {
    warning("'cols_to_keep' is NULL: returning all columns"); return(data)
  }
  if (length(cols_to_keep) == 0) {
    out <- data; warning("'cols_to_keep' is empty: returning all columns")
  } else { out <- data[ , cols_to_keep, drop = FALSE] }
  return(out)
}


# .filterClonesBySequenceLength <- function(data, seq_col, min_length = 3) {
#   drop_rows <- rep(FALSE, nrow(data))
#   for (i in 1:length(seq_col)) {
#     drop_rows <- drop_rows | nchar(data[ , seq_col[[i]]]) < min_length
#   }
#   if (sum(drop_rows) > 0) {
#     if (ncol(data) == 1) {
#       out <- as.data.frame(data[-drop_rows, ]); colnames(out) <- colnames(data)
#     } else {
#       out <- data[-drop_rows, ]
#     }
#     return(out)
#   } else {
#     return(data)
#   }
# }
#
#
# .filterClonesBySequenceContent <- function(data, seq_col, drop_matches) {
#   cat(paste0("Removing sequences containing matches to the expression '", drop_matches, "'..."))
#   drop_rows <- rep(FALSE, nrow(data))
#   for (i in 1:length(seq_col)) {
#     drop_rows <- drop_rows | grepl(drop_matches, data[ , seq_col[[i]]])
#   }
#   if (sum(drop_rows) > 0) {
#     if (ncol(data) == 1) {
#       out <- as.data.frame(data[-drop_rows, ]); names(out) <- names(data)
#     } else {
#       out <- data[-drop_rows, ]
#     }
#   }
#   cat(paste0(" Done. ", nrow(out), " rows remaining.\n"))
#   return(out)
# }


# return subset of data corresponding to adjacency matrix
.subsetDataForAdjacencyMatrix <- function(data, adjacency_matrix) {
  return(data[as.numeric(dimnames(adjacency_matrix)[[1]]), , drop = FALSE])
}


# # FUNCTION: EXTRACT DATA SUBSET FOR ALL SEQUENCES WITHIN SPECIFIED RADIUS OF
# # TARGET SEQUENCE BY SPECIFIED DISTANCE TYPE
# # If a sample_col is provided, only samples that possess the target sequence
# # will be included
# getSimilarClones <- function(
    #     target_seq, # specified candidate sequence for the neighborhood
#     data, # data frame containing rep seq data, possibly from multiple samples
#     seq_col, # col name/# containing clone sequences
#     sample_col = NULL, # optional col name/# containing sample IDs (only samples possessing target seq will be included)
#     dist_type = "hamming", # options are "hamming" and "levenshtein"
#     max_dist = 1, # Maximum Levenshtein distance allowed for inclusion in neighborhood
#     drop_matches = NULL # regular expression for chars to filter sequences by
# ) {
#   if (!is.data.frame(data)) {
#     data <- as.data.frame(data)
#   }
#   # If sample_id is supplied, subset data keeping only samples with target seq
#   if (is.null(sample_col)) {
#     data_samples_w_targetseq <- data
#   } else {
#     ### SUBSET DATA: SAMPLES WITH TARGET SEQ ###
#     cat("Finding all samples that possess the target sequence...")
#     # Get row ids of merged data corresponding to target seq
#     rows_for_targetseq <- grep(pattern = paste0("^", target_seq, "$"),
#                                x = data[ , seq_col])
#     # Extract rows of merged data corresponding to samples with target seq
#     data_samples_w_targetseq <-
#       data[
#         data[ , sample_col] %in% data[rows_for_targetseq, sample_col], ]
#     cat(" Done.\n")
#   }
#   # Remove sequences that match expression in `drop_matches`
#   if (!is.null(drop_matches)) {
#     drop_matches <- grep(drop_matches, data_samples_w_targetseq[ , seq_col])
#     if (length(drop_matches) > 0) {
#       data_samples_w_targetseq <- data_samples_w_targetseq[-drop_matches, ]
#     }
#   }
#   ### SUBSET DATA: NEIGHBORHOOD OF TARGET SEQUENCE ###
#   # Compute list of bounded distances between target seq and seqs
#   # possessed by samples with target seq (values are -1 where bound is exceeded)
#   # returned vector will be of type integer; names will be the sequences
#   cat("Gathering all cells/clones with receptor sequences similar to the target sequence...")
#   data_targetseq_neighborhood <- .neighborhood(
#     data_samples_w_targetseq, seq_col, target_seq, dist_type, max_dist)
#   cat(paste0(" Done. ", nrow(data_targetseq_neighborhood), " similar cells/clones found.\n"))
#   return(data_targetseq_neighborhood)
# }


# INPUT:
#   rep-seq data w/ specification for count column and one or more grouping cols
#   (clone seq is typically the grouping variable used)
# DO:
#   aggregate the counts by group based on the grouping columns
#   add variable counting the number of reads for each
aggregateIdenticalClones <- function(
    data, # data frame containing columns below
    clone_col,
    count_col, # name or number of column of `data` containing clone counts
    freq_col,
    grouping_cols = NULL # optional integer or character vector specifying additional grouping columns
) {

  # Convert column specifications from numeric to character if not already
  if (is.numeric(clone_col)) { clone_col <- names(data)[clone_col] }
  if (is.numeric(count_col)) { count_col <- names(data)[count_col] }
  if (is.numeric(freq_col)) { freq_col <- names(data)[freq_col] }

  # Define grouping variable(s)
  grouping_variables <- list(data[[clone_col]])
  names(grouping_variables) <- clone_col
  if (!is.null(grouping_cols)) {
    if (is.numeric(grouping_cols)) { grouping_cols <- names(data)[grouping_cols] }
    for (i in 1:length(grouping_cols)) {
      grouping_variables$newvar <- data[[grouping_cols[[i]]]]
      names(grouping_variables)[[length(grouping_variables)]] <-
        grouping_cols[[i]] } }

  # aggregate the reads by group
  cat("Aggregating reads (rows) by unique clone sequence...")
  data_to_aggregate <- list("AggregatedCloneCount" = data[ , c(count_col)],
                            "AggregatedCloneFrequency" = data[ , c(freq_col)])
  agg_counts <- stats::aggregate(data_to_aggregate,
                                 by = grouping_variables, FUN = sum)

  # add variable for num reads (row count)
  groups <- as.data.frame(grouping_variables)
  # names(groups)[[1]] <- "temporary_placeholder_name" # for summarize function
  num_reads <- dplyr::summarize(dplyr::group_by_all(groups),
                                UniqueCloneCount = length({{ clone_col }}))
  names(num_reads)[[1]] <- clone_col # replace placeholder name with orig

  # Merge aggregate counts with num reads
  out <- merge(agg_counts, num_reads, by = c(clone_col, grouping_cols))
  cat(paste0(" Done. ", nrow(out), " unique clone sequences found.\n"))

  return(out)
}

# Adjacency Matrices ------------------------------------------------------


# FUNCTION: COMPUTE ADJACENCY MATRIX FOR LEVENSHTEIN OR HAMMING DISTANCE
# Returns sparse matrix
sparseAdjacencyMatFromSeqs <- function(
    seqs, # List of tcr/clonotype sequences
    dist_type = "hamming", # supports "levenshtein" and "hamming"
    max_dist = 1, # Maximum distance threshold for edge/adjacency between two sequences
    drop_isolated_nodes = TRUE # Drop sequences/nodes with zero degree?
) {
  # attempt to coerce seqs to character vector
  if (length(seqs) == 0) stop("'seqs' has zero length")
  seqs <- as.vector(seqs, mode = "character")
  if (!is.character(seqs)) stop("'seqs' must be cocercible to a character vector")
  if (!is.vector(seqs)) stop("'seqs' must be cocercible to a character vector")

  # Compute adjacency matrix
  if (dist_type %in% c("levenshtein", "Levenshtein, lev, Lev, l, L")) {
    cat(paste0("Computing network edges based on a max ", dist_type, " distance of ", max_dist, "..."))
    out <- .levAdjacencyMatSparse(seqs, max_dist, drop_isolated_nodes)
  } else if (dist_type %in% c("hamming", "Hamming", "ham", "Ham", "h", "H")) {
    cat(paste0("Computing network edges based on a max ", dist_type, " distance of ", max_dist, "..."))
    out <- .hamAdjacencyMatSparse(seqs, max_dist, drop_isolated_nodes)
  } else {
    stop('invalid option for `dist_type`')
  }
  cat(" Done.\n")
  # Number of nodes with positive network degree
  num_nodes <- dim(out)[[1]]
  if (num_nodes == 0) {
    warning("No edges exist using the specified distance cutoff")
  } else {
    if (drop_isolated_nodes) {
      cat(paste("Network contains", num_nodes, "nodes (after removing isolated nodes).\n"))
      # Import record of selected column IDs and use for matrix row names
      clone_ids <- utils::read.table("col_ids.txt")
      dimnames(out)[[1]] <- clone_ids$V1
      dimnames(out)[[2]] <- seqs[clone_ids$V1]
      # cat(paste0("The row names of the adjacency matrix contain the original index values of the corresponding sequences; the column names contain the sequences themselves. They can be accessed using `dimnames()`\n"))
    } else {
      cat(paste("Network contains", num_nodes, "nodes.\n"))
    }
  }

  if (file.exists("col_ids.txt")) { unlink("col_ids.txt") } # cleanup

  return(out)
}


# Adjacency Matrix: Euclidean Distance on Atchley Factor Embedding
# (Only applicable to TCR CDR3 Amino Acid Sequences)
# This function is intended for building the network for a single cluster, where
# the adjacency matrix is typically dense
adjacencyMatAtchleyFromSeqs <- function(
    seqs, # List of TCR CDR3 amino acid sequences corresponding to the seqs
    contig_ids = seq_along(seqs), # used by BriseisEncoder to perform the Atchley-factor embedding of the TCR sequences
    max_dist, # Maximum Euclidean distance threshold for edge/adjacency between two sequences
    return_type = "adjacency_matrix", # can be set to "distance_matrix" to return the distance matrix instead
    outfile_distance_matrix = NULL # savefile for Euclidean distance matrix
) {
  # Embed amino acid seqs in Euclidean 30-space by Atchley factor representation
  embedded_values <- encodeTCRSeqsByAtchleyFactor(seqs, contig_ids)

  # Compute Euclidean distance matrix on embedded sequence values
  cat("Computing Euclidean distances between the embedded values...")
  distance_matrix <- as.matrix(stats::dist(embedded_values[ , -1]))
  cat(" Done.\n")

  if (!is.null(outfile_distance_matrix)) {
    # Save distance matrix to file
    utils::write.csv(distance_matrix, outfile_distance_matrix)
    cat(paste0("Distance matrix saved to file:\n  ", outfile_distance_matrix,
               "\n")) }

  if (return_type == "distance_matrix") {
    return(distance_matrix)
  } else {
    # Convert distance matrix to adjacency matrix using specified bound
    cat(paste0("Generating adjacency matrix based on a maximum distance of ",
               max_dist, "..."))
    adjacency_matrix <-
      matrix(1, nrow = nrow(distance_matrix), ncol = ncol(distance_matrix))
    adjacency_matrix[distance_matrix > max_dist] <- 0
    cat(" Done.\n")
    return(adjacency_matrix)
  }
}





# Network Building --------------------------------------------------------


# input adjacency matrix
# returns igraph network
generateNetworkFromAdjacencyMat <- function(adjacency_matrix) {
  set.seed(9999)
  net <- igraph::graph_from_adjacency_matrix(adjacency_matrix, weighted = TRUE)
  net <- igraph::as.undirected(
    igraph::simplify(net, remove.multiple = T, remove.loops = T))
  return(net)
}


# Input vector of sequences
# return igraph network
.generateNetworkFromSeqs <- function(
    seqs, # character vector of receptor sequences
    dist_type = "hamming", # supports "levenshtein", "hamming", "euclidean_on_atchley"
    dist_cutoff = 1, # max dist threshold for edges
    drop_isolated_nodes = TRUE, # forced to FALSE for dist_type = "euclidean_on_atchley"
    contig_ids = seq_along(seqs), # for dist_type = "euclidean_on_atchley"
    outfile_adjacency_matrix = NULL, # save file for adjacency matrix
    outfile_distance_matrix = NULL, # save file for distance matrix (only for Euclidean on Atchley)
    return_type = "network" # can use "adjacency_matrix" to return the adjacency mat
) {
  ### COMPUTE ADJACENCY MATRIX ###
  if (dist_type %in% c("levenshtein", "hamming")) {
    adjacency_matrix <-
      sparseAdjacencyMatFromSeqs(seqs = seqs,
                                 dist_type = dist_type,
                                 max_dist = dist_cutoff,
                                 drop_isolated_nodes = drop_isolated_nodes)
    if (!is.null(outfile_adjacency_matrix)) {
      Matrix::writeMM(adjacency_matrix, outfile_adjacency_matrix)
    }
    # If no nodes are connected (empty matrix), return NULL
    if (sum(dim(adjacency_matrix)) == 0) { return(NULL) }
  } else if (dist_type == "euclidean_on_atchley") {
    adjacency_matrix <-
      adjacencyMatAtchleyFromSeqs(
        seqs = seqs,
        contig_ids = contig_ids,
        max_dist = dist_cutoff,
        outfile_distance_matrix = outfile_distance_matrix)
  } else { stop("invalid option for argument `dist_type`") }
  if (return_type == "adjacency_matrix") {
    return(adjacency_matrix)
  } else {
    network <- generateNetworkFromAdjacencyMat(adjacency_matrix)
    return(network)
  }
}

# Input data frame
# return list of network objects (single chain)
.generateSingleChainNetwork <- function(
    data, seq_col, dist_type, dist_cutoff, drop_isolated_nodes
) {
  adjacency_matrix <- .generateNetworkFromSeqs(
    data[[seq_col]], dist_type, dist_cutoff, contig_ids = rownames(data),
    return_type = "adjacency_matrix", drop_isolated_nodes = drop_isolated_nodes)
  if (is.null(adjacency_matrix)) { return(NULL) }
  net <- generateNetworkFromAdjacencyMat(adjacency_matrix)
  if (drop_isolated_nodes & dist_type != "euclidean_on_atchley") {
    data <- .subsetDataForAdjacencyMatrix(data, adjacency_matrix)
  }

  return(list("igraph" = net, "adjacency_matrix" = adjacency_matrix,
              "node_data" = as.data.frame(data)))
}

# input data frame
# return list of network objects (dual chain)
.generateDualChainNetwork <- function(
    data, a_col, b_col, dist_type, dist_cutoff, drop_isolated_nodes
) {
  # adjacency matrix for alpha chain
  cat("Computing graph adjacency based on sequences in first chain:\n")
  adj_mat_a <- sparseAdjacencyMatFromSeqs(
    data[[a_col]], dist_type, dist_cutoff, drop_isolated_nodes = FALSE)
  # If no nodes are connected (empty matrix), return NULL
  if (sum(dim(adj_mat_a)) == 0) {
    warning("No edges exist in the network for the first chain using the specified distance type and cutoff")
    return(NULL)
  }
  # adjacency matrix for beta chain
  cat("Computing graph adjacency based on sequences in second chain:\n")
  adj_mat_b <- sparseAdjacencyMatFromSeqs(
    data[[b_col]], dist_type, dist_cutoff, drop_isolated_nodes = FALSE)
  # If no nodes are connected (empty matrix), return NULL
  if (sum(dim(adj_mat_b)) == 0) {
    warning("No edges exist in the network for the second chain using the specified distance type and cutoff")
    return(NULL)
  }

  # Combine adjacency matrices for both chains
  # (only edges present for both chains will become edges in the combined graph)
  cat("Intersecting the adjacencies from both chains...")
  adjacency_matrix <- adj_mat_a + adj_mat_b
  adjacency_matrix[adjacency_matrix == 1] <- 0
  adjacency_matrix[adjacency_matrix == 2] <- 1
  cat(" Done.\n")

  # Generate network from combined adjacency matrix
  cat("Building network based on the combined adjacencies... ")
  net <- generateNetworkFromAdjacencyMat(adjacency_matrix); cat(" Done.\n")

  # Drop isolated nodes from final network if specified
  if (drop_isolated_nodes) { cat("Dropping isolated nodes...")
    nodes_to_keep <- igraph::degree(net) > 0
    # If no edges exist, return NULL with warning
    if (sum(nodes_to_keep) == 0) {
      warning("No edges exist in the combined network for both chains using the specified distance type and cutoff")
      return(NULL)
    }
    adjacency_matrix <- adjacency_matrix[nodes_to_keep, nodes_to_keep]
    data <- data[nodes_to_keep, , drop = FALSE]
    # regenerate network without isolated nodes
    net <- generateNetworkFromAdjacencyMat(adjacency_matrix)
    cat(" Done.\n")
  }
  cat(paste("Network contains", nrow(data), "nodes.\n"))
  return(list("igraph" = net, "adjacency_matrix" = adjacency_matrix,
              "adj_mat_a" = adj_mat_a, "adj_mat_b" = adj_mat_b,
              "node_data" = as.data.frame(data)))
}

# Input data frame
# generate single or dual chain network
# returns list of network objects
generateNetworkObjects <- function(
    data, seq_col, dist_type = "hamming", dist_cutoff = 1,
    drop_isolated_nodes = TRUE
) {
  if (length(seq_col) == 1) {
    return(.generateSingleChainNetwork(
      data, seq_col,
      dist_type, dist_cutoff, drop_isolated_nodes))
  } else if (length(seq_col) == 2) {
    return(.generateDualChainNetwork(
      data, seq_col[[1]], seq_col[[2]],
      dist_type, dist_cutoff, drop_isolated_nodes))
  }
  # if (file.exists("col_ids.txt")) { file.remove("col_ids.txt") } # cleanup
}


# Network Properties ------------------------------------------------------


# input igraph network and data frame;
# return input data frame augmented with node-level properties
addNodeNetworkStats <- function(
    data, # rep-seq data corresponding to the network
    net, # igraph network object
    stats_to_include = node_stat_settings(),
    cluster_fun = cluster_fast_greedy
) {

  if (typeof(stats_to_include) != "list")  {
    if (stats_to_include == "all") {
      stats_to_include <- node_stat_settings(all_stats = TRUE)
    } else if (stats_to_include == "cluster_id_only") {
      stats_to_include <- node_stat_settings(
        degree = FALSE, cluster_id = TRUE, transitivity = FALSE,
        eigen_centrality = FALSE, centrality_by_eigen = FALSE,
        betweenness = FALSE, centrality_by_betweenness = FALSE,
        authority_score = FALSE, coreness = FALSE, page_rank = FALSE)
    } }
  if (stats_to_include$degree | stats_to_include$all_stats) {
    data$degree <- igraph::degree(net) }

  if (stats_to_include$cluster_id | stats_to_include$all_stats) {
    cat("Computing cluster membership within the network...")
    data$cluster_id <- as.factor(as.integer(cluster_fun(net)$membership))
    cat(" Done.\n")
  }

  cat(paste0("Computing node-level network statistics..."))
  if (stats_to_include$transitivity | stats_to_include$all_stats) {
    data$transitivity <- igraph::transitivity(net, type = "local")
  }
  if (stats_to_include$closeness | stats_to_include$all_stats) {
    data$closeness <- igraph::closeness(net, mode = "all", weights = NA)
  }
  if (stats_to_include$centrality_by_closeness | stats_to_include$all_stats) {
    data$centrality_by_closeness <-
      igraph::centr_clo(net, mode = "all", normalized = T)$res
  }
  if (stats_to_include$eigen_centrality | stats_to_include$all_stats) {
    data$eigen_centrality <-
      igraph::eigen_centrality(net, directed = T, weights = NA)$vector
  }
  if (stats_to_include$centrality_by_eigen | stats_to_include$all_stats) {
    data$centrality_by_eigen <-
      igraph::centr_eigen(net, directed = T, normalized = T)$vector
  }
  if (stats_to_include$betweenness | stats_to_include$all_stats) {
    data$betweenness <- igraph::betweenness(net, directed = T, weights = NA) }

  if (stats_to_include$centrality_by_betweenness | stats_to_include$all_stats) {
    data$centrality_by_betweenness <-
      igraph::centr_betw(net, directed = T, normalized = T)$res
  }
  if (stats_to_include$authority_score | stats_to_include$all_stats) {
    data$authority_score <- igraph::authority_score(net, weights = NA)$vector
  }
  if (stats_to_include$coreness | stats_to_include$all_stats) {
    data$coreness <- igraph::coreness(net, mode = "all")
  }
  if (stats_to_include$page_rank | stats_to_include$all_stats) {
    data$page_rank <- igraph::page_rank(net)$vector
  }
  cat(" Done.\n")
  return(data)
}

# return list specifying node-level properties by T/F
# pass to `stats_to_include` argument of `computeNodeNetworkStats()`
node_stat_settings <- function(
    degree = TRUE,
    cluster_id = FALSE,
    transitivity = TRUE,
    closeness = FALSE,
    centrality_by_closeness = FALSE,
    eigen_centrality = TRUE,
    centrality_by_eigen = TRUE,
    betweenness = TRUE,
    centrality_by_betweenness = TRUE,
    authority_score = TRUE,
    coreness = TRUE,
    page_rank = TRUE,
    all_stats = FALSE
) {
  list(degree = degree,
       cluster_id = cluster_id,
       transitivity = transitivity,
       closeness = closeness,
       centrality_by_closeness = centrality_by_closeness,
       eigen_centrality = eigen_centrality,
       centrality_by_eigen = centrality_by_eigen,
       betweenness = betweenness,
       centrality_by_betweenness = centrality_by_betweenness,
       authority_score = authority_score,
       coreness = coreness,
       page_rank = page_rank,
       all_stats = all_stats)
}


# Input igraph network and data frame
# return input data frame augmented with cluster ID
addClusterMembership <- function(
    data, net, fun = cluster_fast_greedy) {
  cat("Computing cluster membership within the network...")
  data$cluster_id <- as.factor(as.integer(fun(net)$membership))
  cat(" Done.\n")
  return(data)
}

# Input adjacency matrix and data frame
# return new data frame containing cluster-level properties
getClusterStats <- function(
    data, # rep-seq data for network, with node-level network stats
    adjacency_matrix, # adjacency matrix for network
    seq_col = NULL, # name or number of column of `data` containing the clone sequences
    count_col = NULL, # name or number of column of `data` containing the clone counts
    cluster_id_col = NULL, # optional name or number of column of `data` containing the cluster IDs
    degree_col = NULL, # optional name or number of column of `data` containing the network degree
    cluster_fun = cluster_fast_greedy
) {

  # Compute Cluster ID and network degree if not provided
  if (is.null(cluster_id_col) | is.null(degree_col)) {
    net <- generateNetworkFromAdjacencyMat(adjacency_matrix)
    if (is.null(cluster_id_col)) {
      cluster_id_col <- "cluster_id"
      data <- addClusterMembership(data, net, cluster_fun)
    }
    if (is.null(degree_col)) {
      degree_col <- "deg"
      data$deg <- igraph::degree(net)
    }
  }

  # 1. case with multiple seq cols (i.e. dual chain)
  if (!is.null(seq_col)) {
    if (length(seq_col) > 1) {
      return(.computeClusterStatsDualChain(
        data, adjacency_matrix, seq_col, count_col, cluster_id_col, degree_col))
    }
  }

  # 2. case with one seq col (or none)
  return(.computeClusterStats(
    data, adjacency_matrix, seq_col, count_col, cluster_id_col, degree_col))

}

# Input adjacency matrix and data frame
# return new data frame containing cluster-level properties
# (single chain only)
.computeClusterStats <- function(data, adjacency_matrix, seq_col, count_col,
                                 cluster_id_col, degree_col) {

  # Compute sequence length
  if (!is.null(seq_col)) { seq_lengths <- nchar(data[[seq_col]]) }

  # Tabulate the number of nodes in each cluster
  out <- as.data.frame(table(data[[cluster_id_col]]))
  colnames(out) <- c("cluster_id", "node_count")
  num_clusters <- nrow(out) # Total number of clusters
  cat(paste0("Computing statistics for the ", num_clusters, " clusters in the network..."))

  ### INITIALIZE VALUES ###

  out$eigen_centrality_eigenvalue <-
    out$eigen_centrality_index <-
    out$closeness_centrality_index <-
    out$degree_centrality_index <-
    out$edge_density <-
    out$assortativity <-
    out$global_transitivity <-
    out$diameter_length <-
    out$seq_w_max_count <-
    out$max_count <-
    out$agg_count <-
    out$seq_w_max_degree <-
    out$max_degree <-
    out$mean_degree <-
    out$mean_seq_length <- NA

  ### COMPUTE STATS FOR EACH CLUSTER ###
  for (i in 1:num_clusters) {

    # current row of cluster data
    cluster_row <- which(out$cluster_id == i)

    # Boolean vec: Rows of node data for current cluster
    node_ids <- data$cluster_id == i

    # Mean sequence length in cluster
    if (!is.null(seq_col)) {
      out$mean_seq_length[[cluster_row]] <-
        round(mean(seq_lengths[node_ids]), 2)
    }

    # Mean degree in cluster
    out$mean_degree[[cluster_row]] <-
      round(mean(data[node_ids, degree_col]), 2)

    # Maximum degree (and corresponding seq) within cluster
    max_deg <- max(data[node_ids, degree_col])
    out$max_degree[[cluster_row]] <- max_deg

    if (!is.null(seq_col)) {
      node_id_max_deg <- which(node_ids & data[[degree_col]] == max_deg)[[1]]
      out$seq_w_max_degree[[cluster_row]] <-
        as.character(data[node_id_max_deg, seq_col])
    }

    if (!is.null(count_col)) {
      # Total aggregate clonotype count in cluster
      out$agg_count[[cluster_row]] <- sum(data[node_ids, count_col])

      # Maximum clonotype count (and corresponding seq) within cluster
      max_count <- max(data[node_ids, count_col])
      out$max_count[[cluster_row]] <- max_count

      if (!is.null(seq_col)) {
        node_id_max_count <- which(node_ids & data[[count_col]] == max_count)[[1]]
        out$seq_w_max_count[[cluster_row]] <-
          as.character(data[node_id_max_count, seq_col])
      }
    }

    # Build cluster network to get network properties for the cluster
    cluster_adjacency_matrix <- as.matrix(adjacency_matrix[node_ids, node_ids])
    cluster <- generateNetworkFromAdjacencyMat(cluster_adjacency_matrix)

    # Diameter (longest geodesic distance)
    out$diameter_length[[cluster_row]] <-
      length(igraph::get_diameter(cluster, directed = T))

    # Assortativity
    out$assortativity[[cluster_row]] <-
      igraph::assortativity_degree(cluster, directed = F)

    # Transitivity
    out$global_transitivity[[cluster_row]] <-
      igraph::transitivity(cluster, type = "global")  # cluster is treated as an undirected network

    # Density: The proportion of present edges from all possible ties.
    out$edge_density[[cluster_row]] <-
      igraph::edge_density(cluster, loops = F)

    # Centralization on degree
    out$degree_centrality_index[[cluster_row]] <-
      igraph::centr_degree(cluster, mode = "in", normalized = T)$centralization

    # Centralization on Closeness (centrality based on distance to others in the graph)
    out$closeness_centrality_index[[cluster_row]] <-
      igraph::centr_clo(cluster, mode = "all", normalized = T)$centralization

    # Centralization on Eigenvector (centrality proportional to the sum of connection centralities)
    #  (values of the first eigenvector of the graph adjacency matrix)
    out$eigen_centrality_index[[cluster_row]] <-
      igraph::centr_eigen(cluster, directed = T, normalized = T)$centralization

    out$eigen_centrality_eigenvalue[[cluster_row]] <-
      igraph::eigen_centrality(cluster, directed = T, weights = NA)$value
  }
  cat(" Done.\n")

  return(out)
}


# Input adjacency matrix and data frame
# return new data frame containing cluster-level properties
# (dual chain only)
.computeClusterStatsDualChain <- function(
    data, adjacency_matrix, seq_col, count_col,
    cluster_id_col, degree_col)
{

  # Compute sequence length
  A_seq_lengths <- nchar(data[[seq_col[[1]]]])
  B_seq_lengths <- nchar(data[[seq_col[[2]]]])

  # Tabulate the number of nodes in each cluster
  out <- as.data.frame(table(data[[cluster_id_col]]))
  colnames(out) <- c("cluster_id", "node_count")
  num_clusters <- nrow(out) # Total number of clusters
  cat(paste0("Computing statistics for the ", num_clusters, " clusters in the network..."))

  ### INITIALIZE VALUES ###

  out$eigen_centrality_eigenvalue <-
    out$eigen_centrality_index <-
    out$closeness_centrality_index <-
    out$degree_centrality_index <-
    out$edge_density <-
    out$assortativity <-
    out$global_transitivity <-
    out$diameter_length <-
    out$B_seq_w_max_count <-
    out$A_seq_w_max_count <-
    out$max_count <-
    out$agg_count <-
    out$B_seq_w_max_degree <-
    out$A_seq_w_max_degree <-
    out$max_degree <-
    out$mean_degree <-
    out$mean_B_seq_length <-
    out$mean_A_seq_length <- NA

  ### COMPUTE STATS FOR EACH CLUSTER ###
  for (i in 1:num_clusters) {

    # current row of cluster data
    cluster_row <- which(out$cluster_id == i)

    # Boolean vec: Rows of node data for current cluster
    node_ids <- data$cluster_id == i

    # Mean sequence length in cluster
    out$mean_A_seq_length[[cluster_row]] <- round(mean(A_seq_lengths[node_ids]), 2)
    out$mean_B_seq_length[[cluster_row]] <- round(mean(B_seq_lengths[node_ids]), 2)

    # Mean degree in cluster
    out$mean_degree[[cluster_row]] <- round(mean(data[node_ids, degree_col]), 2)

    # Maximum degree (and corresponding seq) within cluster
    max_deg <- max(data[node_ids, degree_col])
    out$max_degree[[cluster_row]] <- max_deg
    node_id_max_deg <- which(node_ids & data[[degree_col]] == max_deg)[[1]]
    out$A_seq_w_max_degree[[cluster_row]] <-
      as.character(data[node_id_max_deg, seq_col[[1]]])
    out$B_seq_w_max_degree[[cluster_row]] <-
      as.character(data[node_id_max_deg, seq_col[[2]]])

    if (!is.null(count_col)) {
      # Total aggregate clonotype count in cluster
      out$agg_count[[cluster_row]] <- sum(data[node_ids, count_col])

      # Maximum clonotype count (and corresponding seq) within cluster
      max_count <- max(data[node_ids, count_col])
      out$max_count[[cluster_row]] <- max_count

      node_id_max_count <- which(node_ids & data[[count_col]] == max_count)[[1]]
      out$A_seq_w_max_count[[cluster_row]] <-
        as.character(data[node_id_max_count, seq_col[[1]]])
      out$B_seq_w_max_count[[cluster_row]] <-
        as.character(data[node_id_max_count, seq_col[[2]]])
    }

    # Build cluster network to get network properties for the cluster
    cluster_adjacency_matrix <- as.matrix(adjacency_matrix[node_ids, node_ids])
    cluster <- generateNetworkFromAdjacencyMat(cluster_adjacency_matrix)

    # Diameter (longest geodesic distance)
    out$diameter_length[[cluster_row]] <-
      length(igraph::get_diameter(cluster, directed = T))

    # Assortativity
    out$assortativity[[cluster_row]] <-
      igraph::assortativity_degree(cluster, directed = F)

    # Transitivity
    out$global_transitivity[[cluster_row]] <-
      igraph::transitivity(cluster, type = "global")  # cluster is treated as an undirected network

    # Density: The proportion of present edges from all possible ties.
    out$edge_density[[cluster_row]] <-
      igraph::edge_density(cluster, loops = F)

    # Centralization on degree
    out$degree_centrality_index[[cluster_row]] <-
      igraph::centr_degree(cluster, mode = "in", normalized = T)$centralization

    # Centralization on Closeness (centrality based on distance to others in the graph)
    out$closeness_centrality_index[[cluster_row]] <-
      igraph::centr_clo(cluster, mode = "all", normalized = T)$centralization

    # Centralization on Eigenvector (centrality proportional to the sum of connection centralities)
    #  (values of the first eigenvector of the graph adjacency matrix)
    out$eigen_centrality_index[[cluster_row]] <-
      igraph::centr_eigen(cluster, directed = T, normalized = T)$centralization

    out$eigen_centrality_eigenvalue[[cluster_row]] <-
      igraph::eigen_centrality(cluster, directed = T, weights = NA)$value
  }
  cat(" Done.\n")

  return(out)
}

# wrapper to getClusterStats, for potentially avoiding redoing clustering computation
# after computing node level stats
.getClusterStats2 <- function(data, adjacency_matrix, seq_col, count_col,
                              cluster_fun = cluster_fast_greedy) {
  cluster_id_col <- degree_col <- NULL
  if ("cluster_id" %in% names(data)) { cluster_id_col <- "cluster_id" }
  if ("degree" %in% names(data)) { degree_col <- "degree" }
  return(getClusterStats(data, adjacency_matrix, seq_col,
                         count_col, cluster_id_col, degree_col, cluster_fun))
}


# Visualization -----------------------------------------------------------



plotNetworkGraph <- function(igraph,
                             plot_title = NULL,
                             plot_subtitle = NULL,
                             color_nodes_by = NULL,
                             color_scheme = "default",
                             color_legend = "auto",
                             color_title = "auto",
                             edge_width = 0.1,
                             size_nodes_by = 0.5,
                             node_size_limits = NULL,
                             size_title = "auto",
                             outfile = NULL
) {
  set.seed(9999)
  layout <- igraph::layout_components(igraph)

  graph_plot <-
    ggraph::ggraph(igraph, layout = layout) +
    ggraph::geom_edge_link0(width = edge_width, colour = "grey")

  if (!is.null(color_nodes_by)) {
    # Custom node color scheme
    if (!is.null(size_nodes_by)) {
      # Custom node size scheme
      if (is.numeric(size_nodes_by) & length(size_nodes_by) == 1) {
        # Fixed node sizes
        graph_plot <- graph_plot +
          ggraph::geom_node_point(
            ggplot2::aes(color = color_nodes_by), size = size_nodes_by)
      } else if (length(size_nodes_by) > 1) {
        # size_nodes_by is a numeric vector specifying the size of each node
        graph_plot <- graph_plot +
          ggraph::geom_node_point(
            ggplot2::aes(color = color_nodes_by, size = size_nodes_by))
        if (length(node_size_limits) == 2) {
          # Rescale node sizes if specified
          graph_plot <-
            graph_plot + ggplot2::scale_size(range = node_size_limits)
        }
      }
    } else { # size_nodes_by is null or invalid
      graph_plot <- graph_plot +
        ggraph::geom_node_point(ggplot2::aes(color = color_nodes_by))

    }
  } else { # color_nodes_by is null or invalid
    if (!is.null(size_nodes_by)) {
      if (is.numeric(size_nodes_by) & length(size_nodes_by) == 1) {
        graph_plot <- graph_plot +
          ggraph::geom_node_point(size = size_nodes_by)

      } else if (length(size_nodes_by) > 1) {
        graph_plot <- graph_plot +
          ggraph::geom_node_point(ggplot2::aes(size = size_nodes_by))
        if (!is.null(node_size_limits)) {
          graph_plot <-
            graph_plot + ggplot2::scale_size(range = node_size_limits)
        }
      }
    } else { # size_nodes_by null or invalid
      graph_plot <- graph_plot + ggraph::geom_node_point()
    }
  }

  if (is.null(color_nodes_by)) {
    color_type <- "continuous"
  } else {
    color_type <- ggplot2::scale_type(color_nodes_by)[[1]]
  }

  if (color_legend == "auto") {
    if (is.null(color_nodes_by)) { color_legend <- FALSE
    } else {
      if (color_type == "continuous" || length(unique(color_nodes_by)) <= 20) {
        color_legend <- TRUE
      } else {
        # too many legend values: hide legend; add color varname to subtitle
        color_legend <- FALSE
        color_varname <- ifelse(is.null(color_title),
                                no = color_title,
                                yes = deparse(substitute(color_nodes_by)))
        subtitle_affix <- paste0("Nodes colored by ", color_varname)
        plot_subtitle <- ifelse(is.null(plot_subtitle), yes = subtitle_affix,
                                no = paste0(plot_subtitle, "\n", subtitle_affix))
      }
    }
  }

  graph_plot <- graph_plot +
    ggraph::theme_graph(base_family = "sans") +
    ggplot2::labs(title = plot_title, subtitle = plot_subtitle)

  if (color_legend == TRUE) {
    if (is.null(color_title)) {
      graph_plot <- graph_plot +
        ggplot2::guides(color = ggplot2::guide_legend(title = color_title))
    } else if (color_title != "auto") {
      graph_plot <- graph_plot +
        ggplot2::guides(color = ggplot2::guide_legend(title = color_title))
    }
  } else {
    graph_plot <- graph_plot + ggplot2::guides(color = "none")
  }

  if (is.null(size_title)) {
    graph_plot <- graph_plot +
      ggplot2::guides(size = ggplot2::guide_legend(title = size_title))
  } else if (size_title != "auto") {
    graph_plot <- graph_plot +
      ggplot2::guides(size = ggplot2::guide_legend(title = size_title))
  }


  # Convert node-color variable to factor if discrete
  # if (color_type == "discrete") { color_nodes_by <- as.factor(color_nodes_by) }

  if (color_scheme != "default") {

    if (color_type == "continuous") {
      if (color_scheme %in% c("A", "B", "C", "F", "G",
                              "magma", "inferno", "plasma", "rocket", "mako")) {
        graph_plot <- graph_plot +
          ggraph::scale_color_viridis(option = color_scheme,
                                      begin = 0.2, end = 0.8)
      } else if (color_scheme %in% c("D", "viridis")) {
        graph_plot <- graph_plot +
          ggraph::scale_color_viridis(option = color_scheme,
                                      begin = 0, end = 0.9)
      } else if (color_scheme %in% c("E", "H", "cividis",  "turbo")) {
        graph_plot <- graph_plot +
          ggraph::scale_color_viridis(option = color_scheme)
      } else if (color_scheme %in% c("A-1", "B-1", "C-1", "F-1", "G-1",
                                     "magma-1", "inferno-1", "plasma-1", "rocket-1", "mako-1")) {
        graph_plot <- graph_plot +
          ggraph::scale_color_viridis(option = strsplit(color_scheme, "-1")[[1]],
                                      begin = 0.2, end = 0.8, direction = -1)
      } else if (color_scheme %in% c("D-1", "viridis-1")) {
        graph_plot <- graph_plot +
          ggraph::scale_color_viridis(option = strsplit(color_scheme, "-1")[[1]],
                                      begin = 0, end = 0.9, direction = -1)
      } else if (color_scheme %in% c("E-1", "H-1", "cividis-1",  "turbo-1")) {
        graph_plot <- graph_plot +
          ggraph::scale_color_viridis(option = strsplit(color_scheme, "-1")[[1]], direction = -1)
      } else { warning("Value for 'color_scheme' is not a valid option for continuous variables; using default color scheme instead") }

    } else { # discrete color scheme
      if (color_scheme %in% c("A", "B", "C", "F", "G",
                              "magma", "inferno", "plasma", "rocket", "mako")) {
        graph_plot <- graph_plot +
          ggraph::scale_color_viridis(option = color_scheme, discrete = TRUE,
                                      begin = 0.2, end = 0.8)
      } else if (color_scheme %in% c("D", "viridis")) {
        graph_plot <- graph_plot +
          ggraph::scale_color_viridis(option = color_scheme, discrete = TRUE,
                                      begin = 0, end = 0.9)
      } else if (color_scheme %in% c("E", "H", "cividis",  "turbo")) {
        graph_plot <- graph_plot +
          ggraph::scale_color_viridis(option = color_scheme, discrete = TRUE,)
      } else if (color_scheme %in% c("A-1", "B-1", "C-1", "F-1", "G-1",
                                     "magma-1", "inferno-1", "plasma-1", "rocket-1", "mako-1")) {
        graph_plot <- graph_plot +
          ggraph::scale_color_viridis(option = strsplit(color_scheme, "-1")[[1]], discrete = TRUE,
                                      begin = 0.2, end = 0.8, direction = -1)
      } else if (color_scheme %in% c("D-1", "viridis-1")) {
        graph_plot <- graph_plot +
          ggraph::scale_color_viridis(option = strsplit(color_scheme, "-1")[[1]], discrete = TRUE,
                                      begin = 0, end = 0.9, direction = -1)
      } else if (color_scheme %in% c("E-1", "H-1", "cividis-1",  "turbo-1")) {
        graph_plot <- graph_plot +
          ggraph::scale_color_viridis(option = strsplit(color_scheme, "-1")[[1]], discrete = TRUE,
                                      direction = -1)
      } else if (color_scheme %in% grDevices::hcl.pals()) {
        graph_plot <- graph_plot +
          ggplot2::scale_color_manual(
            values = grDevices::hcl.colors(n = length(unique(color_nodes_by)),
                                           palette = color_scheme))
      } else { warning("Value for 'color_scheme' is not a valid option for discrete variables; using default color scheme instead") } }
  }


  if (!is.null(outfile)) {
    grDevices::pdf(file = outfile, width = 12, height = 8)
    print(graph_plot)
    grDevices::dev.off()
    cat(paste0("Plot of network graph saved to file:\n  ", outfile, "\n"))
  }
  return(graph_plot)
}


generateNetworkGraphPlots <- function(
    igraph, data, print_plots = TRUE,
    plot_title = NULL, plot_subtitle = NULL,
    color_nodes_by = NULL, color_scheme = "default",
    color_legend = "auto", color_title = "auto",
    edge_width = 0.1, size_nodes_by = 0.5,
    node_size_limits = NULL, size_title = "auto")
{
  # harmonize inputs
  new_inputs <- .harmonizePlottingInputs(
    color_nodes_by, color_scheme, color_title, size_nodes_by, size_title)
  color_scheme <- new_inputs$color_scheme
  color_title <- new_inputs$color_title; size_title <- new_inputs$size_title
  if (is.character(size_nodes_by)) { size_nodes_by <- data[[size_nodes_by]] }

  plotlist <- list()
  if (is.null(color_nodes_by)) {
    cat("Generating graph plot...")
    plotlist$uniform_color <-
      plotNetworkGraph(
        igraph, plot_title = plot_title, plot_subtitle = plot_subtitle,
        color_nodes_by = NULL,
        color_title = NULL,
        color_scheme = color_scheme,
        color_legend = FALSE,
        edge_width = edge_width,
        size_nodes_by = size_nodes_by,
        size_title = size_title,
        node_size_limits = node_size_limits)
    if (print_plots) { print(plotlist$uniform_color) }
    cat(" Done.\n")
  } else {
    for (j in 1:length(color_nodes_by)) {
      cat(paste0("Generating graph plot with nodes colored by ",
                 color_nodes_by[[j]], "..."))
      plotlist$newplot <-
        plotNetworkGraph(
          igraph, plot_title = plot_title, plot_subtitle = plot_subtitle,
          color_nodes_by = data[[color_nodes_by[[j]]]],
          color_title = color_title[[j]],
          color_scheme = color_scheme[[j]],
          color_legend = color_legend,
          edge_width = edge_width,
          size_nodes_by = size_nodes_by,
          size_title = size_title,
          node_size_limits = node_size_limits)
      if (print_plots) { print(plotlist$newplot) }
      names(plotlist)[[length(names(plotlist))]] <- color_nodes_by[[j]]
      cat(" Done.\n")
    }
  }
  return(plotlist)
}


# Input graph plot and vector of node labels; return annotated plot
addGraphLabels <- function(plot, node_labels, size = 5, color = "black") {
  plot +
    ggraph::geom_node_text(ggplot2::aes(label = node_labels),
                           size = size, color = color)
}


# Input network list and graph plot;
# return plot with top n clusters labeled
addClusterLabels <- function(plot, net,
                             top_n_clusters = 20, criterion = "node_count",
                             size = 5, color = "black") {
  dat <- net$node_data
  cdat <- net$cluster_data

  # Sort cluster stats by specified criterion
  cdat <- cdat[order(-cdat[[criterion]]) , ]

  # Identify the top n clusters by specified criterion
  clusters <- as.integer(cdat$cluster_id[1:top_n_clusters])

  # Create vector of node labels; only one node per top n cluster is labeled
  node_labels <- as.integer(dat$cluster_id)
  node_labels[duplicated(node_labels) | !node_labels %in% clusters] <- NA

  # Annotate plot with cluster labels and return
  return(addGraphLabels(plot, node_labels, size, color))

}


.generateNetworkGraphPlotsGuarded <- function(
    igraph, data, print_plots,
    plot_title = NULL, plot_subtitle = NULL,
    color_nodes_by = NULL, color_scheme = "default",
    color_legend = "auto", color_title = "auto",
    edge_width = 0.1, size_nodes_by = 0.5,
    node_size_limits = NULL, size_title = "auto")
{
  if (nrow(data) > 1e06) {
    warning("Network contains over 1 million nodes; depending on the number of network edges, this may exceed ggraph limitations. Skipping automatic generation of network graph; you can attempt to generate the graph manually using `plotNetworkGraph()`")
    return(invisible(NULL))
  } else {
    return(generateNetworkGraphPlots(
      igraph, data, print_plots, plot_title, plot_subtitle, color_nodes_by,
      color_scheme, color_legend, color_title, edge_width, size_nodes_by,
      node_size_limits, size_title))
  }
}

.harmonizePlottingInputs <- function(
    color_nodes_by, color_scheme, color_title, size_nodes_by, size_title)
{
  # node color palette and legend title
  if (is.null(color_nodes_by)) { color_title <- NULL
  } else { # color_nodes_by is non NULL
    if (length(color_nodes_by) > 1) { # extend objects to vectors if needed
      if (length(color_scheme) == 1) {
        color_scheme <- rep(color_scheme, length(color_nodes_by)) }
      if (!is.null(color_title)) { if (length(color_title) == 1) {
        color_title <- rep(color_title, length(color_nodes_by)) }
      } else { # color_title is NULL
        color_title <- rep("", length(color_nodes_by)) } } # (hack, since can't have NULL vector entries)
    if (!is.null(color_title)) { # Set default color legend title if applicable
      for (i in 1:length(color_title)) { if (color_title[[i]] == "auto") {
        color_title[[i]] <- color_nodes_by[[i]] } } } }

  # size legend title
  if (!is.null(size_title)) { if (size_title == "auto") { # default size title
    if (is.numeric(size_nodes_by)) { size_title <- NULL } # fixed node sizes
    if (is.character(size_nodes_by)) { size_title <- size_nodes_by } } }
  return(list(color_scheme = color_scheme,
              color_title = color_title,
              size_title = size_title))
}


# check if color_nodes_by is "auto"; if so, look for variables to use
.passColorNodesBy <- function(color_nodes_by, data, count_col) {
  if (length(color_nodes_by) == 1) { if (color_nodes_by == "auto") {
    if ("cluster_id" %in% names(data)) { color_nodes_by <- "cluster_id"
    } else if ("transitivity" %in% names(data)) { color_nodes_by <- "transitivity"
    } else if ("closeness" %in% names(data)) { color_nodes_by <- "closeness"
    } else if ("degree" %in% names(data)) { color_nodes_by <- "degree"
    } else if ("eigen_centrality" %in% names(data)) { color_nodes_by <- "eigen_centrality"
    } else if ("betweenness" %in% names(data)) { color_nodes_by <- "betweenness"
    } else if ("coreness" %in% names(data)) { color_nodes_by <- "coreness"
    } else if ("authority_score" %in% names(data)) { color_nodes_by <- "authority_score"
    } else if ("page_rank" %in% names(data)) { color_nodes_by <- "page_rank"
    } else if (!is.null(count_col)) { color_nodes_by <- count_col
    } else { color_nodes_by <- NULL # default to uniform node colors
    } } }
  return(color_nodes_by)
}

.makePlotTitle <- function(plot_title, type = "standard", network_name = NULL) {
  if (is.null(plot_title)) { return(plot_title) }
  if (type == "standard") {
    if (plot_title == "auto") {
      # if (!is.null(seq_col)) { if (length(seq_col) == 2) { return("Network\nby Dual-Chain Similarity") } }
      if (!is.null(network_name)) { return(network_name)
      } else { return("Network by Receptor Sequence Similarity") }
    }
  } else if (type == "pub_clust_rep") {
    if (plot_title == "auto") { return("Network of Public Clusters by Representative Sequence") }
  }
  return(plot_title)
}

.makePlotSubtitle <- function(plot_subtitle, type = "standard", seq_col = NULL,
                              dist_type = NULL, dist_cutoff = NULL)
{
  if (is.null(plot_subtitle)) { return(plot_subtitle) }
  if (type == "standard") {
    if (plot_subtitle == "auto") {
      if (!is.null(seq_col)) { if (length(seq_col) == 2) { return(paste("Each node denotes a single TCR/BCR cell\nEdges denote a maximum", dist_type, "distance of", dist_cutoff, "between sequences in corresponding chains\n")) } }
      if (dist_type == "euclidean_on_atchley") { return(paste("Each node denotes a single TCR/BCR cell or clone\nSequences encoded numerically using deep learning based on Atchley factor representation\nEdges denote a maximum Euclidean distance of", dist_cutoff, "between encoded values\n")) }
      return(paste("Each node denotes a single TCR/BCR cell or clone\nEdges denote a maximum", dist_type, "distance of", dist_cutoff, "between receptor sequences\n"))
    }
  } else if (type == "pub_clust_rep") {
    if (plot_subtitle == "auto" & seq_col == "seq_w_max_count") { return("Each node corresponds to a cluster\nSimilarity is based on the sequence with the highest count in each cluster") }
    return(paste0("Each node corresponds to a cluster\nSimilarity is based on ", seq_col))
  }
  return(plot_subtitle)
}



# Atchley Factor Embedding ------------------------------------------------

# Embed TCR CDR3 amino acid sequences in Euclidean 30-space based on the Atchley
# factor representations of their elements, using a trained encoding model
encodeTCRSeqsByAtchleyFactor <- function(
    cdr3_AA, # List of TCR CDR3 amino acid sequences
    contig_ids = seq_along(cdr3_AA) # used by BriseisEncoder
) {
  .checkPythonModules()
  if (length(cdr3_AA) != length(contig_ids)) {
    stop("length of `cdr3_AA` and `contig_ids` must match")
  }
  # Write sequences and contig_ids to temporary file
  tempfile_clones <- file.path(getwd(),
                               "Atchley_factor_tcr_only_tmp.csv")
  utils::write.csv(data.frame("contig_id" = contig_ids, "cdr3" = cdr3_AA),
                   file = tempfile_clones, row.names = FALSE)

  # Write files for trained encoder and Atchley factor table to working dir
  file_trained_encoder <-
    system.file(file.path("python", "TrainedEncoder.h5"),
                package = "NAIR")
  file_atchley_table <-
    system.file(file.path("python", "Atchley_factors.csv"),
                package = "NAIR")
  utils::write.csv(data.frame("sysfiles" = c(file_atchley_table,
                                             file_trained_encoder)),
                   file.path(getwd(), "temp_sysfiles.csv"),
                   row.names = FALSE)

  # Run BriseisEncoder python function
  cat("Embedding TCR CDR3 amino acid sequences in 30-dimensional Euclidean space based on Atchley factor representation and a trained encoding model using deep learning routines from the keras & tensorflow python modules. This may produce some warnings...\n")
  reticulate::py_run_file(
    system.file(file.path("python", "BriseisEncoder_modified.py"),
                package = "NAIR"))

  # Import embedded values and remove temp files
  embedded_values <- utils::read.csv("temp_atchley_factors_encoded.csv")
  file.remove("Atchley_factor_tcr_only_tmp.csv", "temp_sysfiles.csv",
              "temp_atchley_factors_encoded.csv")
  cat("Embedding complete.\n")
  warning("The encoder was trained on TCR CDR3 sequences; results not valid for other amino acid sequences")
  return(embedded_values)
}
