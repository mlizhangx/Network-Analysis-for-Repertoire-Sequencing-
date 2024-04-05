# NAIR: Network Analysis of Immune Repertoire
# Copyright (C) 2023 Li Zhang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# File Loading ------------------------------------------------------------
loadDataFromFileList <- function(
    file_list,
    input_type,
    data_symbols = NULL,
    header, sep, read.args
)
{
  if (missing(header)) {
    header <-  switch(input_type, "txt" = FALSE, "table" = FALSE, TRUE)
  }
  if (missing(sep)) {
    sep <- switch(input_type, "csv" = ",", "csv2" = ";", "tsv" = "\t", "")
  }
  if (missing(read.args)) { read.args <- NULL }
  .checkargs.InputFiles(
    file_list, input_type, data_symbols, header, sep, read.args
  )
  if (input_type %in% c("csv", "csv2", "tsv", "txt", "table")) {
    read.args <- .checkReadArgs(read.args, header, sep)
  }
  if (input_type == "rda" && length(file_list) > length(data_symbols)) {
    data_symbols <- rep(data_symbols, length.out = length(file_list))
  }
  data_list <- sapply(1:length(file_list), .loadFile.i, files = file_list,
                      type = input_type, symbols = data_symbols, head = header,
                      sep = sep, args = read.args, simplify = FALSE
  )
  for (i in 1:length(file_list)) {
    data_list[[i]]$rowid <- paste0("file", i, ".", rownames(data_list[[i]]))
  }
  # data_list <- vector(mode = "list", length = length(file_list))
  # names(data_list) <- paste0("file", 1:length(file_list))
  # for (i in 1:length(file_list)) {
  #   data_list[[i]] <- .loadDataFromFile(
  #     file_list[[i]], input_type, data_symbols[[i]], header, sep, read.args
  #   )
  # }
  data <- do.call(rbind, data_list)
  rownames(data) <- data$rowid
  data$rowid <- NULL
  data
}

combineSamples <- function(
    file_list,
    input_type,
    data_symbols = NULL,
    header, sep, read.args,
    seq_col = NULL,
    min_seq_length = NULL,
    drop_matches = NULL,
    subset_cols = NULL,
    sample_ids = NULL,
    subject_ids = NULL,
    group_ids = NULL,
    verbose = FALSE
) {
  if (missing(header)) {
    header <-  switch(input_type, "txt" = FALSE, "table" = FALSE, TRUE)
  }
  if (missing(sep)) {
    sep <- switch(input_type, "csv" = ",", "csv2" = ";", "tsv" = "\t", "")
  }
  if (missing(read.args)) { read.args <- NULL }
  .checkargs.InputFiles(
    file_list, input_type, data_symbols, header, sep, read.args
  )
  if (input_type %in% c("csv", "csv2", "tsv", "txt", "table")) {
    read.args <- .checkReadArgs(read.args, header, sep)
  }
  .orNull(.MUST.isCharOrIntegerVector, seq_col)
  .orNull(.MUST.hasLength, seq_col, len = c(1, 2))
  if (!is.null(seq_col)) {
    min_seq_length <- .check(min_seq_length, .isNonneg, NULL, ornull = TRUE)
    drop_matches <- .check(drop_matches, .isString, NULL, ornull = TRUE)
  }
  subset_cols <- .check(subset_cols, .isCharOrNumericVector, NULL,
                        ornull = TRUE
  )
  sample_ids <- .checkIDs(sample_ids, length(file_list), ornull = TRUE)
  subject_ids <- .checkIDs(subject_ids, length(file_list), ornull = TRUE,
                           allow_dupes = TRUE
  )
  group_ids <- .checkIDs(group_ids, length(file_list), ornull = TRUE,
                         allow_dupes = TRUE
  )
  if (!is.null(sample_ids)) {
    sample_ids <- as.vector(sample_ids, mode = "character")
  }
  if (!is.null(subject_ids)) {
    subject_ids <- as.vector(subject_ids, mode = "character")
  }
  if (!is.null(group_ids)) {
    group_ids <- as.vector(group_ids, mode = "character")
  }
  if (input_type == "rda" && length(file_list) > length(data_symbols)) {
    data_symbols <- rep(data_symbols, length.out = length(file_list))
  }
  msg <- .makemsg(verbose)
  # temp_symbols <- paste0("sample", 1:length(file_list))
  data_list <- vector(mode = "list", length = length(file_list))
  for (i in 1:length(file_list)) {
    msg("Loading sample ", i, "...", newline = FALSE)
    data_list[[i]] <- .loadDataFromFile(
      file_list[[i]], input_type, data_symbols[[i]], header, sep, read.args
    )
    if (i == 1) {
      seq_col <- .convertColRef(seq_col, data_list[[i]])
      subset_cols <- .convertColRef(subset_cols, data_list[[i]])
    }
    if (!is.null(seq_col)) {
      data_list[[i]] <- filterInputData(data_list[[i]], seq_col,
                                        min_seq_length, drop_matches,
                                        subset_cols, verbose = verbose
      )
    }
    if (!is.null(sample_ids)) {
      data_list[[i]]$SampleID <- sample_ids[[i]]
      data_list[[i]]$rowid <- paste0(sample_ids[[i]], ".",
                                     rownames(data_list[[i]])
      )
    } else {
      data_list[[i]]$rowid <- paste0("file", i, ".", rownames(data_list[[i]]))
    }
    if (!is.null(subject_ids)) { data_list[[i]]$SubjectID <- subject_ids[[i]] }
    if (!is.null(group_ids)) { data_list[[i]]$GroupID <- group_ids[[i]] }
    # assign(temp_symbols[[i]], tmp)
    # rm(tmp)
  }
  # data_list <- mget(temp_symbols)
  # rm(list = temp_symbols)
  data <- do.call(rbind, data_list)
  rownames(data) <- data$rowid
  data$rowid <- NULL
  data
}

.loadDataFromFile <- function(
    input_file,
    input_type,
    data_symbol,
    header, sep, read.args
) {
  if (input_type == "rds") {
    data <- readRDS(input_file)
  } else if (input_type == "rda") {
    load(input_file)
    assign(x = "data", value = get(data_symbol))
  } else {
    if (is.null(read.args)) {
      read.args <- list(header = header, sep = sep)
    }
    read.args$file <- input_file
    read_fun <- switch(input_type,
                       "csv" = utils::read.csv,
                       "csv2" = utils::read.csv2,
                       "tsv" = utils::read.delim,
                       utils::read.table
    )
    data <- do.call(read_fun, read.args)
  }

  data
}

.loadFile.i <- function(i, files, type, symbols, head, sep, args) {
  # Same as .loadDataFromFile, but designed for use with sapply()
  if (type == "rds") {
    data <- readRDS(files[[i]])
  } else if (type == "rda") {
    load(files[[i]])
    assign(x = "data", value = get(symbols[[i]]))
  } else {
    if (is.null(args)) {
      args <- list(header = head, sep = sep)
    }
    args$file <- files[[i]]
    read_fun <- switch(type,
                       "csv" = utils::read.csv,
                       "csv2" = utils::read.csv2,
                       "tsv" = utils::read.delim,
                       utils::read.table
    )
    data <- do.call(read_fun, args)
  }
  data
}

# File Saving -------------------------------------------------------------
saveNetwork <- function(
    net,
    output_dir,
    output_type = "rds",
    output_name = "MyRepSeqNetwork",
    pdf_width = 12,
    pdf_height = 10,
    verbose = FALSE,
    output_filename = deprecated()
) {
  if (lifecycle::is_present(output_filename)) {
    lifecycle::deprecate_soft(
      when = "1.0.1",
      what = "saveNetwork(output_filename)",
      with = "saveNetwork(output_name)"
    )
    output_name <- output_filename
  }
  .MUST.isBaseNetworkOutput(net)
  msg <- .makemsg(verbose)
  if (!.orNull(.isString, output_dir)) {
    warning(sQuote("output_dir"), " is not a character string. ",
            "Output will not be saved"
    )
    return(FALSE)
  }
  .createOutputDir(output_dir)
  if (is.null(output_dir)) {
    msg(sQuote("output_dir"), " is ", dQuote("NULL"),
        ". Output will not be saved."
    )
    return(FALSE)
  } else if (!isTRUE(dir.exists(output_dir))) {
    warning("directory ", dQuote(output_dir),
            " specified for ", sQuote("output_dir"),
            " does not exist and could not be created. Output will not be saved"
    )
    return(FALSE)
  }
  output_type <- .checkOutputType(output_type)
  output_name <- .checkOutputName(output_name, "MyRepSeqNetwork")
  pdf_width <- .check(pdf_width, .isPos, 12)
  pdf_height <- .check(pdf_width, .isPos, 10)
  if (output_type == "individual") {
    .saveNetworkObjects(net, output_dir, output_name, verbose)
  } else if (output_type == "rds") {
    .saveNetworkRDS(net, output_dir, output_name, verbose)
  } else {
    .saveNetworkRDA(net, output_dir, output_name, verbose)
  }
  if ("plots" %in% names(net)) {
    outfile_layout <- NULL
    if (output_type == "individual") {
      outfile_layout <- file.path(output_dir,
                                  paste0(output_name, "_GraphLayout.txt")
      )
    }
    saveNetworkPlots(
      net$plots, file.path(output_dir, paste0(output_name, ".pdf")),
      pdf_width, pdf_height, outfile_layout, verbose
    )
  }
  TRUE
}

saveNetworkPlots <- function(
    plotlist,
    outfile,
    pdf_width = 12,
    pdf_height = 10,
    outfile_layout = NULL,
    verbose = FALSE
) {
  .MUST.isPlotlist(plotlist)
  .MUST.isStringOrConnection(outfile)
  # if (!.isValidFilename(basename(outfile))) {
  #   warning(dQuote(outfile), " may be an unsafe file name")
  # }
  outfile_layout <- .checkOutfileLayout(outfile_layout, plotlist)
  outfile_layout <- .check(outfile_layout, .isStringOrConnection, NULL,
                           ornull = TRUE
  )
  # if (!is.null(outfile_layout) && !.isValidFilename(basename(outfile_layout))) {
  #   warning(dQuote(outfile_layout), " may be an unsafe file name")
  # }
  pdf_width <- .check(pdf_width, .isPos, 12)
  pdf_height <- .check(pdf_width, .isPos, 10)
  msg <- .makemsg(verbose)
  grDevices::pdf(file = outfile, width = pdf_width, height = pdf_height)
  for (j in 1:length(plotlist)) {
    if (inherits(plotlist[[j]], "gg")) {
      print(plotlist[[j]])
    }
  }
  grDevices::dev.off()
  msg("Network graph plots printed to pdf file:\n  ", outfile)
  if (!is.null(outfile_layout) && "graph_layout" %in% names(plotlist)) {
    write(plotlist[[which(names(plotlist) == "graph_layout")]],
          file = outfile_layout,
          ncolumns = 2
    )
    msg("Graph layout for plots saved to file:\n  ", outfile_layout)
  }
  invisible(TRUE)
}

.createOutputDir <- function(dirname) {
  if (!is.null(dirname) && !isTRUE(dir.exists(dirname))) {
    dir.create(dirname, showWarnings = FALSE, recursive = TRUE)
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
    # data <- apply(data, MARGIN = 2, FUN = as.character)
    utils::write.csv(data,
                     file = file.path(output_dir, paste0(output_name, ".csv"))
    )
  } else if (output_type == "csv2") {
    # data <- apply(data, MARGIN = 2, FUN = as.character)
    utils::write.csv2(data,
                      file = file.path(output_dir, paste0(output_name, ".csv"))
    )
  } else if (output_type %in% c("tsv")) {
    # data <- apply(data, MARGIN = 2, FUN = as.character)
    utils::write.table(data, row.names = TRUE, sep = "\t",
                       file = file.path(output_dir, paste0(output_name, ".tsv"))
    )
  } else if (output_type %in% c("txt", "table")) {
    # data <- apply(data, MARGIN = 2, FUN = as.character)
    utils::write.table(data, row.names = TRUE,
                       file = file.path(output_dir, paste0(output_name, ".txt"))
    )
  } else {
    save(data, file = file.path(output_dir, paste0(output_name, ".rda")))
  }
}

.saveNetworkRDA <- function(net, output_dir, output_filename,
                            verbose = FALSE
) {
  msg <- .makemsg(verbose)
  outfile <- file.path(output_dir, paste0(output_filename, ".rda"))
  save(net, file = outfile)
  msg("List 'net' saved to file:\n  ", outfile)
}

.saveNetworkRDS <- function(net, output_dir, output_filename,
                            verbose = FALSE
) {
  msg <- .makemsg(verbose)
  outfile <- file.path(output_dir, paste0(output_filename, ".rds"))
  saveRDS(net, file = outfile)
  msg("List of network objects saved to file:\n  ", outfile)
}

.saveNetworkObjects <- function(net, output_dir, output_filename,
                                verbose = FALSE
) {
  msg <- .makemsg(verbose)
  if ("details" %in% names(net) && isTRUE(inherits(net$details, "list"))) {
    details_outfile <- file.path(
      output_dir, paste0(output_filename, "_Details.txt")
    )
    utils::write.table(net$details, details_outfile, row.names = TRUE,
                       col.names = FALSE, quote = FALSE, sep = "\n\t"
    )
    saveRDS(net$details, file = details_outfile)
    msg("Network details saved to data file:\n  ", details_outfile)
  }
  node_file <- file.path(
    output_dir, paste0(output_filename, "_NodeMetadata.csv")
  )
  utils::write.csv(net$node_data, file = node_file)
  msg("Node metadata saved to file:\n  ", node_file)

  if ("cluster_data" %in% names(net)) {
    cluster_file <- file.path(
      output_dir, paste0(output_filename, "_ClusterMetadata.csv")
    )
    utils::write.csv(net$cluster_data, file = cluster_file, row.names = FALSE)
    msg("Cluster metadata saved to file:\n  ", cluster_file)
  }

  igraph_outfile <- file.path(
    output_dir, paste0(output_filename, "_EdgeList.txt")
  )
  igraph::write_graph(net$igraph, file = igraph_outfile, format = "edgelist")
  msg("Network igraph saved in edgelist format to file:\n  ", igraph_outfile)

  if (inherits(net$adjacency_matrix, "matrix")) {
    matrix_outfile <- file.path(
      output_dir, paste0(output_filename, "_AdjacencyMatrix.csv")
    )
    utils::write.csv(net$adjacency_matrix, matrix_outfile, row.names = FALSE)
    msg("Adjacency matrix saved to file:\n  ", matrix_outfile)
  } else if (inherits(net$adjacency_matrix, "sparseMatrix")) {
    matrix_outfile <- file.path(
      output_dir, paste0(output_filename, "_AdjacencyMatrix.mtx")
    )
    Matrix::writeMM(net$adjacency_matrix, matrix_outfile)
    msg("Adjacency matrix saved to file:\n  ", matrix_outfile)
    if ("adj_mat_a" %in% names(net)) {
      matrix_outfile <- file.path(
        output_dir, paste0(output_filename, "_AdjacencyMatrix_ChainA.mtx")
      )
      Matrix::writeMM(net$adj_mat_a, matrix_outfile)
      msg("Adjacency matrix for first chain saved to file:\n  ", matrix_outfile)
      matrix_outfile <- file.path(
        output_dir, paste0(output_filename, "_AdjacencyMatrix_ChainB.mtx")
      )
      Matrix::writeMM(net$adj_mat_b, matrix_outfile)
      msg(
        "Adjacency matrix for second chain saved to file:\n  ", matrix_outfile
      )
    }
  }

  if ("plots" %in% names(net)) {
    plots_outfile <- file.path(output_dir,
                               paste0(output_filename, "_Plots.rda")
    )
    plots <- net$plots
    save(plots, file = plots_outfile)
    msg("List of network graph plots named 'plots' saved to data file:\n  ",
        plots_outfile
    )
  }
  list_outfile <- file.path(output_dir, paste0(output_filename, ".rda"))
  save(net, file = list_outfile)
  msg("List 'net' saved to file:\n  ", list_outfile)

}


# Filtering and Subsetting ------------------------------------------------
filterInputData <- function(
    data,
    seq_col,
    min_seq_length = NULL,
    drop_matches = NULL,
    subset_cols = NULL,
    count_col = NULL,
    verbose = FALSE
) {
  # if (lifecycle::is_present(count_col)) {
  #   lifecycle::deprecate_soft(
  #     when = "1.0.1",
  #     what = "filterInputData(count_col)",
  #   )
  # }
  data_name <- deparse(substitute(data))
  data <- as.data.frame(data)
  .MUST.isDataFrame(data, data_name)
  .MUST.isSeqColrefs(seq_col, data, deparse(substitute(seq_col)), data_name)
  min_seq_length <- .check(min_seq_length, .isNonneg, NULL, ornull = TRUE)
  drop_matches <- .check(drop_matches, .isString, NULL, ornull = TRUE)
  subset_cols <- .checkDataColrefs(subset_cols, data, NULL)
  count_col <- .checkCountCol(count_col, data, NULL)
  msg <- .makemsg(verbose)
  msg("Input data contains ", nrow(data), " rows.")
  for (i in 1:length(seq_col)) {
    if (!is.character(data[[seq_col[[i]]]])) {
      data[[seq_col[[i]]]] <- as.character(data[[seq_col[[i]]]])
    }
    NA_indices <- is.na(data[[seq_col[[i]]]])
    if (sum(NA_indices) > 0) {
      if (length(seq_col == 2)) {
        msg("dropping ", sum(NA_indices),
            " rows containing NA/NaN values in sequence column ", i
        )
      } else {
        msg("dropping ", sum(NA_indices),
            " rows containing NA/NaN values in sequence column",
        )
      }
      data <- data[!NA_indices, ]
    }
  }
  if (!is.null(min_seq_length)) {
    msg("Removing sequences with length fewer than ",
        min_seq_length, " characters...", newline = FALSE
    )
    drop_rows <- .RowDropsBySequenceLength(data, seq_col, min_seq_length)
    if (sum(drop_rows) > 0) { data <- data[!drop_rows, , drop = FALSE] }
    msg(" Done. ", nrow(data), " rows remaining.")
  }
  if (!is.null(drop_matches)) {
    msg("Removing sequences containing matches to ",
        dQuote(drop_matches), "...", newline = FALSE
    )
    drop_rows <- .RowDropsBySequenceContent(data, seq_col, drop_matches)
    if (sum(drop_rows) > 0) { data <- data[!drop_rows, , drop = FALSE] }
    msg(" Done. ", nrow(data), " rows remaining.")
  }
  if (!is.null(subset_cols)) {
    data <- .subsetColumns(data, c(seq_col, subset_cols))
  }
  if (!is.null(count_col)) {
    # if (!is.numeric(data[[count_col]])) {
    #   data[[count_col]] <- as.numeric(data[[count_col]])
    # }
    NA_indices <- is.na(data[[count_col]])
    if (sum(NA_indices) > 0) {
      warning("dropping ", sum(NA_indices),
              " rows containing NA/NaN values in count column"
      )
      data <- data[!NA_indices, ]
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
  data_name <- deparse(substitute(data))
  data <- as.data.frame(data)
  .MUST.isDataFrame(data, data_name)
  .MUST.isSeqColref(seq_col, data, deparse(substitute(seq_col)), data_name)
  .MUST.isString(target_seq)
  dist_type <- .checkDistType(dist_type, "hamming")
  max_dist <- .check(max_dist, .isNonneg, 1)
  if (!target_seq %in% data[[seq_col]]) {
    return(NULL)
  }
  dist_fun <- ifelse(dist_type == "levenshtein",
                     yes = levDistBounded, no = hamDistBounded
  )
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
    grouping_cols = NULL,
    verbose = FALSE
) {
  data_name <- deparse(substitute(data))
  data <- as.data.frame(data)
  .MUST.isDataFrame(data, data_name)
  .MUST.isSeqColref(clone_col, data, deparse(substitute(clone_col)), data_name)
  .MUST.isCountColref(count_col, data,
                      deparse(substitute(count_col)), data_name
  )
  .MUST.isCountColref(freq_col, data, deparse(substitute(freq_col)), data_name)
  grouping_cols <- .checkDataColrefs(grouping_cols, data, ornull = TRUE)
  msg <- .makemsg(verbose)
  clone_col <- .convertColRef(clone_col, data)
  count_col <- .convertColRef(count_col, data)
  freq_col <- .convertColRef(freq_col, data)
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
  msg("Aggregating reads (rows) by unique clone sequence...", newline = FALSE)
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
  msg(" Done. ", nrow(out), " unique clone sequences found.")
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
  if (is.null(subset_cols)) {
    return(NULL)
  }
  c(subset_cols, other_cols)
}

.subsetColumns <- function(data, cols_to_keep) {
  cols_to_keep <- intersect(unique(cols_to_keep), names(data))
  if (length(cols_to_keep) == 0) {
    out <- data
  } else {
    out <- data[ , cols_to_keep, drop = FALSE]
  }
  out
}

.subsetDataForAdjacencyMatrix <- function(data, adjacency_matrix) {
  # return subset of data corresponding to adjacency matrix
  return(data[as.numeric(dimnames(adjacency_matrix)[[1]]), , drop = FALSE])
}



# Message and Console Output ----------------------------------------------

.makemsg <- function(verbose) {
  verbose <- .checkTF(verbose, FALSE)
  if (verbose) {
    msg <- function(..., newline = TRUE) { message(..., appendLF = newline) }
  } else {
    msg <- function(..., newline = TRUE) { invisible(NULL) }
  }
  msg
}


