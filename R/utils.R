# File Loading ------------------------------------------------------------
loadDataFromFileList <- function(
    file_list,
    input_type,
    data_symbols = NULL,
    header = TRUE, sep = ""
)
{
  .checkargs.InputFiles(
    file_list, input_type, data_symbols, header, sep
  )
  if (input_type == "csv") {
    data <- plyr::ldply(
      file_list, utils::read.csv, header = header, sep = sep
    )
  } else if (input_type %in% c("table", "tsv", "txt")) {
    data <- plyr::ldply(
      file_list, utils::read.table, header = header, sep = sep
    )
  } else if (input_type == "rds") {
    data <- plyr::ldply(file_list, readRDS)
  } else if (input_type %in% c("rda", "Rda", "Rdata", "rdata")) {
    stopifnot(
      "must specify data symbols via argument `data_symbols`" =
        !is.null(data_symbols)
    )
    if (length(file_list) > length(data_symbols)) {
      data_symbols <- rep(data_symbols, length.out = length(file_list))
    }
    temp_symbols <- paste0("data_for_sample", 1:length(file_list))
    for (i in 1:length(file_list)) {
      load(file_list[[i]])
      assign(temp_symbols[[i]], get(data_symbols[[i]]))
      rm(data_symbols[[i]])
    }
    data_list <- mget(temp_symbols)
    rm(list = temp_symbols)
    data <- do.call(rbind, data_list)
  }
  data
}

combineSamples <- function(
    file_list,
    input_type,
    data_symbols = NULL,
    header = TRUE, sep = "",
    seq_col,
    min_seq_length = NULL,
    drop_matches = NULL,
    subset_cols = NULL,
    sample_ids = NULL,
    subject_ids = NULL,
    group_ids = NULL
) {
  .checkargs.CombineSamples(
    file_list, input_type, data_symbols, header, sep,
    seq_col, min_seq_length, drop_matches, subset_cols,
    sample_ids, subject_ids, group_ids
  )
  cat(">>> Loading and compiling data from all samples:\n")

  if (input_type %in% c("rda", "Rda", "Rdata", "rdata")) {
    stopifnot("must specify data symbols via argument `data_symbols`" =
                !is.null(data_symbols)
    )
    if (length(file_list) > length(data_symbols)) {
      data_symbols <- rep(data_symbols, length.out = length(file_list))
    }
  }
  temp_symbols <- paste0("data_for_sample", 1:length(file_list))

  for (i in 1:length(file_list)) {
    cat(paste0("Loading sample ", i, ": "))
    tmp <- .loadDataFromFile(
      file_list[[i]], input_type, data_symbols[[i]], header, sep
    )
    if (i == 1) {
      seq_col <- .convertColRef(seq_col, tmp)
      subset_cols <- .convertColRef(subset_cols, tmp)
    }
    tmp <- filterInputData(tmp, seq_col, min_seq_length, drop_matches,
                           subset_cols
    )
    if (!is.null(sample_ids)) { tmp$SampleID <- sample_ids[[i]] }
    if (!is.null(subject_ids)) { tmp$SubjectID <- subject_ids[[i]] }
    if (!is.null(group_ids)) { tmp$GroupID <- group_ids[[i]] }
    assign(temp_symbols[[i]], tmp)
    rm(tmp)
  }

  data_list <- mget(temp_symbols)
  rm(list = temp_symbols)
  data <- do.call(rbind, data_list)
  cat("All samples loaded.\n")

  data

}

.loadDataFromFile <- function(
    input_file,
    input_type = "rds",
    data_symbol = NULL,
    header = TRUE, sep = ""
)
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
    load(input_file)
    assign(x = "data", value = get(data_symbol))
  } else {
    data <- utils::read.table(input_file, header, sep)
  }
  data
}

# File Saving -------------------------------------------------------------
saveNetwork <- function(
    net,
    output_dir = getwd(),
    output_type = "individual",
    output_filename = "MyRepSeqNetwork",
    pdf_width = 12, pdf_height = 10
)
{
  .checkargs.saveNetwork(
    net, output_dir, output_type, output_filename, pdf_width, pdf_height
  )
  if (is.null(output_dir)) { return(invisible(FALSE)) } # exit if no output dir
  .createOutputDir(output_dir)
  if (output_type == "rds") {
    .saveNetworkRDS(net, output_dir, output_filename)
    if ("plots" %in% names(net)) {
      saveNetworkPlots(
        net$plots, file.path(output_dir, paste0(output_filename, ".pdf")),
        pdf_width, pdf_height
      )
    }
  } else if (output_type == "individual") {
    .saveNetworkObjects(
      net, output_dir, output_filename, pdf_width, pdf_height
    )
  } else {
    .saveNetworkRDA(net, output_dir, output_filename)
    if ("plots" %in% names(net)) {
      saveNetworkPlots(
        net$plots, file.path(output_dir, paste0(output_filename, ".pdf")),
        pdf_width, pdf_height
      )
    }
  }
  invisible(TRUE)
}

saveNetworkPlots <- function(
    plotlist,
    outfile = "MyRepSeqNetwork.pdf",
    pdf_width = 12, pdf_height = 10
) {
  .checkargs.saveNetworkPlots(
    plotlist, outfile, pdf_width, pdf_height
  )
  grDevices::pdf(file = outfile, width = pdf_width, height = pdf_height)
  for (j in 1:length(plotlist)) { print(plotlist[[j]]) }
  grDevices::dev.off()
  cat(paste0("Network graph plots saved to file:\n  ", outfile, "\n"))
  invisible(TRUE)
}

.ensureOutputDir <- function(output_dir) {
  stopifnot("output_dir is required" = !is.null(output_dir))
  if (!dir.exists(output_dir)) { .createOutputDir(output_dir) }
  stopifnot("could not create output_dir" = dir.exists(output_dir))
}

.createOutputDir <- function(dirname) {
  if (!is.null(dirname) && !dir.exists(dirname)) {
    dir.create(dirname, showWarnings = FALSE, recursive = TRUE)
    if (!dir.exists(dirname)) {
      stop(paste0(
        "Unable to create directory ", dirname,
        ". Check to confirm that a valid directory name was provided."
      ))
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
            file = file.path(output_dir, paste0(output_name, ".rds"))
    )
  } else if (output_type == "csv") {
    data <- apply(data, MARGIN = 2, FUN = as.character)
    utils::write.csv(data, row.names = FALSE,
                     file = file.path(output_dir, paste0(output_name, ".csv"))
    )
  } else if (output_type %in% c("tsv", "table")) {
    data <- apply(data, MARGIN = 2, FUN = as.character)
    utils::write.table(data, row.names = FALSE,
                       file = file.path(output_dir, paste0(output_name, ".tsv"))
    )
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
  node_file <-
    file.path(output_dir, paste0(output_filename, "_NodeMetadata.csv"))
  net$node_data <- apply(net$node_data, MARGIN = 2, FUN = as.character)
  utils::write.csv(net$node_data, file = node_file, row.names = FALSE)
  cat(paste0(
    "Node-level meta-data saved to file:\n  ", node_file, "\n"
  ))

  if ("cluster_data" %in% names(net)) {
    cluster_file <-
      file.path(output_dir, paste0(output_filename, "_ClusterMetadata.csv"))
    net$cluster_data <- apply(net$cluster_data, MARGIN = 2, FUN = as.character)
    utils::write.csv(net$cluster_data, file = cluster_file, row.names = FALSE)
    cat(paste0(
      "Cluster-level meta-data saved to file:\n  ", cluster_file, "\n"
    ))
  }

  if ("plots" %in% names(net)) {
    saveNetworkPlots(
      net$plots, file.path(output_dir, paste0(output_filename, ".pdf")),
      pdf_width, pdf_height
    )
  }

  igraph_outfile <-
    file.path(output_dir, paste0(output_filename, "_EdgeList.txt"))
  igraph::write_graph(net$igraph, file = igraph_outfile, format = "edgelist")
  cat(paste0(
    "Network igraph saved in edgelist format to file:\n  ", igraph_outfile, "\n"
  ))

  if (inherits(net$adjacency_matrix, "matrix")) {
    matrix_outfile <-
      file.path(output_dir, paste0(output_filename, "_AdjacencyMatrix.csv"))
    utils::write.csv(net$adjacency_matrix, matrix_outfile, row.names = FALSE)
    cat(paste0(
      "Adjacency matrix saved to file:\n  ", matrix_outfile, "\n"
    ))
  } else if (inherits(net$adjacency_matrix, "sparseMatrix")) {
    matrix_outfile <-
      file.path(output_dir, paste0(output_filename, "_AdjacencyMatrix.mtx"))
    Matrix::writeMM(net$adjacency_matrix, matrix_outfile)
    cat(paste0(
      "Adjacency matrix saved to file:\n  ", matrix_outfile, "\n"
    ))
    if ("adj_mat_a" %in% names(net)) {
      matrix_outfile <-
        file.path(
          output_dir, paste0(output_filename, "_AdjacencyMatrix_ChainA.mtx")
        )
      Matrix::writeMM(net$adj_mat_a, matrix_outfile)
      cat(paste0("Adjacency matrix for first chain saved to file:\n  ",
                 matrix_outfile, "\n"
      ))
      matrix_outfile <-
        file.path(
          output_dir, paste0(output_filename, "_AdjacencyMatrix_ChainB.mtx")
        )
      Matrix::writeMM(net$adj_mat_b, matrix_outfile)
      cat(paste0("Adjacency matrix for second chain saved to file:\n  ",
                 matrix_outfile, "\n"
      ))
    }
  }

}


# Filtering and Subsetting ------------------------------------------------
filterInputData <- function(
    data,
    seq_col,
    min_seq_length = NULL,
    drop_matches = NULL,
    subset_cols = NULL,
    count_col = NULL
) {

  data <- as.data.frame(data)
  .checkargs.filterInputData(
    data, seq_col, min_seq_length, drop_matches, subset_cols, count_col
  )
  cat(paste0("Input data contains ", nrow(data), " rows.\n"))
  for (i in 1:length(seq_col)) {
    if (!is.character(data[[seq_col[[i]]]])) {
      data[[seq_col[[i]]]] <- as.character(data[[seq_col[[i]]]])
    }
    NA_indices <- is.na(data[[seq_col[[i]]]])
    if (sum(NA_indices) > 0) {
      warning(paste("Dropping", sum(NA_indices),
                    "rows containing NA values in sequence column", i
      ))
      data <- data[!NA_indices, ]
    }
  }
  if (!is.null(min_seq_length)) {
    cat(paste0("Removing sequences with length fewer than ",
               min_seq_length, " characters..."
    ))
    drop_rows <- .RowDropsBySequenceLength(data, seq_col, min_seq_length)
    if (sum(drop_rows) > 0) { data <- data[!drop_rows, , drop = FALSE] }
    cat(paste0(" Done. ", nrow(data), " rows remaining.\n"))
  }
  if (!is.null(drop_matches)) {
    cat(paste0(
      "Removing sequences containing matches to the expression '",
      drop_matches, "'..."
    ))
    drop_rows <- .RowDropsBySequenceContent(data, seq_col, drop_matches)
    if (sum(drop_rows) > 0) { data <- data[!drop_rows, , drop = FALSE] }
    cat(paste0(" Done. ", nrow(data), " rows remaining.\n"))
  }
  if (!is.null(subset_cols)) {
    data <- .subsetColumns(data, c(seq_col, count_col, subset_cols))
  }
  if (!is.null(count_col)) {
    if (!is.numeric(data[[count_col]])) {
      data[[count_col]] <- as.numeric(data[[count_col]])
    }
    NaN_indices <- is.na(data[[count_col]])
    if (sum(NaN_indices) > 0) {
      warning(paste("Dropping", sum(NaN_indices),
                    "rows containing NA/NaN values in count column", i
      ))
      data <- data[!NaN_indices, ]
    }
  }
  data
}

getNeighborhood <- function(
    data,
    seq_col,
    target_seq,
    dist_type = "hamming",
    max_dist = 1
) {
  if (!target_seq %in% data[[seq_col]]) { return(NULL) }
  ham_keys <- c("hamming", "Hamming", "ham", "Ham", "h", "H")
  lev_keys <- c("levenshtein", "Levenshtein", "lev", "Lev", "l", "L")
  if (dist_type %in% ham_keys) {
    dist_fun <- hamDistBounded
  } else if (dist_type %in% lev_keys) {
    dist_fun <- levDistBounded
  } else {
    stop("invalid option for dist_type")
  }
  dists_to_targetseq <- sapply(
    X = data[[seq_col]],
    FUN = dist_fun, b = target_seq, k = max_dist
  )
  out <- data[dists_to_targetseq != -1, , drop = FALSE]
  out
}

aggregateIdenticalClones <- function(
    data,
    clone_col,
    count_col,
    freq_col,
    grouping_cols = NULL
) {

  if (is.numeric(clone_col)) { clone_col <- names(data)[clone_col] }
  if (is.numeric(count_col)) { count_col <- names(data)[count_col] }
  if (is.numeric(freq_col)) { freq_col <- names(data)[freq_col] }

  grouping_variables <- list(data[[clone_col]])
  names(grouping_variables) <- clone_col
  if (!is.null(grouping_cols)) {
    if (is.numeric(grouping_cols)) {
      grouping_cols <- names(data)[grouping_cols]
    }
    for (i in 1:length(grouping_cols)) {
      grouping_variables$newvar <- data[[grouping_cols[[i]]]]
      names(grouping_variables)[[length(grouping_variables)]] <-
        grouping_cols[[i]]
    }
  }
  cat("Aggregating reads (rows) by unique clone sequence...")
  data_to_aggregate <- list("AggregatedCloneCount" = data[ , c(count_col)],
                            "AggregatedCloneFrequency" = data[ , c(freq_col)]
  )
  agg_counts <- stats::aggregate(data_to_aggregate,
                                 by = grouping_variables, FUN = sum)
  groups <- as.data.frame(grouping_variables)
  num_reads <- dplyr::summarize(dplyr::group_by_all(groups),
                                UniqueCloneCount = length({{ clone_col }})
  )
  names(num_reads)[[1]] <- clone_col
  out <- merge(agg_counts, num_reads, by = c(clone_col, grouping_cols))
  cat(paste0(" Done. ", nrow(out), " unique clone sequences found.\n"))
  out
}

.convertColRef <- function(colrefs, data) {
  if (is.numeric(colrefs)) { colrefs <- names(data)[colrefs] }
  colrefs
}

.RowDropsBySequenceLength <- function(data, seq_col, min_length = 3) {
  drop_rows <- rep(FALSE, nrow(data))
  for (i in 1:length(seq_col)) {
    drop_rows <- drop_rows | (nchar(data[[seq_col[[i]]]]) < min_length)
  }
  drop_rows
}

.RowDropsBySequenceContent <- function(data, seq_col, drop_matches) {
  drop_rows <- rep(FALSE, nrow(data))
  for (i in 1:length(seq_col)) {
    drop_rows <- drop_rows | grepl(drop_matches, data[[seq_col[[i]]]])
  }
  drop_rows
}

.processSubsetCols <- function(subset_cols, other_cols) {
  if (is.null(subset_cols)) { return(NULL) }
  c(subset_cols, other_cols)
}

.subsetColumns <- function(data, cols_to_keep) {
  cols_to_keep <- intersect(unique(cols_to_keep), names(data))
  if (is.null(cols_to_keep)) {
    warning("'cols_to_keep' is NULL: returning all columns")
    return(data)
  }
  if (length(cols_to_keep) == 0) {
    out <- data
    warning("'cols_to_keep' is empty: returning all columns")
  } else {
    out <- data[ , cols_to_keep, drop = FALSE]
  }
  out
}

# return subset of data corresponding to adjacency matrix
.subsetDataForAdjacencyMatrix <- function(data, adjacency_matrix) {
  return(data[as.numeric(dimnames(adjacency_matrix)[[1]]), , drop = FALSE])
}
