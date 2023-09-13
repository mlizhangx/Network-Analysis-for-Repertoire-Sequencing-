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

# Top-Level Functions -----------------------------------------------------

findAssociatedSeqs <- function(
    file_list,
    input_type,
    data_symbols = NULL,
    header, sep, read.args,
    sample_ids = deprecated(),
    subject_ids = NULL,
    group_ids,
    groups = deprecated(),
    seq_col,
    freq_col = NULL,
    min_seq_length = 7,
    drop_matches = "[*|_]",
    min_sample_membership = 5,
    pval_cutoff = 0.05,
    outfile = NULL,
    verbose = FALSE
) {
  .checkDeprecated.findAssociatedSeqs(sample_ids, groups)
  if (missing(header)) {
    header <-  switch(input_type, "txt" = FALSE, "table" = FALSE, TRUE)
  }
  if (missing(sep)) {
    sep <- switch(input_type, "csv" = ",", "csv2" = ";", "tsv" = "\t", "")
  }
  if (missing(read.args)) { read.args <- NULL }
  .checkargs.findAssociatedSeqs(file_list, subject_ids, group_ids, seq_col)
  if (!is.null(outfile) &&
      (!.isString(outfile) || !.isValidFilename(basename(outfile)))
  ) {
    warning("value of argument ", sQuote("outfile"),
            " is invalid (see `?findAssociatedSeqs`).\n",
            "Output will be returned, but will not be saved to file."
    )
    outfile <- NULL
  }
  group_ids <- as.vector(group_ids, mode = "character")
  if (!is.null(subject_ids)) { subject_ids <- as.character(subject_ids) }
  freq_col <- .check(freq_col, .isCharOrNumericScalar, NULL, ornull = TRUE)
  min_seq_length <- .check(min_seq_length, .isNonneg, 7, ornull = TRUE)
  drop_matches <- .check(drop_matches, .isString, "[*|_]", ornull = TRUE)
  min_sample_membership <- .check(min_sample_membership, .isNonneg, 5,
                                  ornull = TRUE
  )
  pval_cutoff <- .check(pval_cutoff, .isPos, 0.05)
  data <- combineSamples(file_list = file_list,
                         input_type = input_type,
                         data_symbols = data_symbols,
                         header = header, sep = sep, read.args = read.args,
                         seq_col = seq_col,
                         min_seq_length = min_seq_length,
                         drop_matches = drop_matches,
                         sample_ids = 1:length(file_list),
                         subject_ids = subject_ids,
                         group_ids = group_ids,
                         subset_cols = c(seq_col, freq_col),
                         verbose = verbose
  )
  .MUST.isSeqColref(seq_col, data, deparse(substitute(seq_col)))
  if (nrow(data) < 2) {
    warning(
      "Returning NULL since fewer than two observations remain after ",
      "loading the data and applying filters for ",
      sQuote("min_sample_membership"), " and ", sQuote("drop_matches")
    )
    return(invisible(NULL))
  }
  msg <- .makemsg(verbose)
  n_samples <- length(group_ids)
  groups <- unique(group_ids)
  ids_g0 <- group_ids == groups[[1]]
  ids_g1 <- group_ids == groups[[2]]
  n_samples_g0 <- sum(ids_g0)
  n_samples_g1 <- sum(ids_g1)
  samples_or_subjects <- "samples"
  if (!is.null(subject_ids)) {
    if (any(duplicated(subject_ids))) { samples_or_subjects <- "subjects" }
  }
  if (samples_or_subjects == "subjects") {
    n_g0 <- length(unique(subject_ids[ids_g0])) # num subjects
    n_g1 <- length(unique(subject_ids[ids_g1])) # num subjects
    n <- length(unique(subject_ids))            # num subjects
    msg("Data contains ", n_samples, " samples and ", n, " subjects, ",
        n_g0, " of which belong to group ", groups[[1]], " and ",
        n_g1, " of which belong to group ", groups[[2]], "."
    )
  } else {
    n_g0 <- n_samples_g0
    n_g1 <- n_samples_g1
    n <- n_samples
    msg("Data contains ", n, " samples, ",
        n_samples_g0, " of which belong to group ", groups[[1]], " and ",
        n_samples_g1, " of which belong to group ", groups[[2]], "."
    )
  }
  out <- .filterSeqsBySampleMembership(data, seq_col, "SampleID",
                                       min_sample_membership, msg
  )
  out <- .filterByFisherPvalue(out, data, seq_col, "GroupID", groups,
                               n_g0, n_g1, pval_cutoff, freq_col,
                               samples_or_subjects, msg
  )
  if (!is.null(outfile)) {
    utils::write.csv(out, outfile, row.names = FALSE)
    msg("Output saved to file:\n  ", outfile)
  }
  invisible(out)
}

findAssociatedSeqs2 <- function(
    data,
    seq_col,
    sample_col,
    subject_col = sample_col,
    group_col,
    groups = deprecated(),
    freq_col = NULL,
    min_seq_length = 7,
    drop_matches = "[*|_]",
    min_sample_membership = 5,
    pval_cutoff = 0.05,
    outfile = NULL,
    verbose = FALSE
) {
  .checkDeprecated.findAssociatedSeqs2(groups)
  data_name <- deparse(substitute(data))
  data <- as.data.frame(data)
  .checkargs.findAssociatedSeqs2(data, seq_col, sample_col, group_col,
                                 data_name, deparse(substitute(seq_col))
  )
  if (!is.null(outfile) &&
      (!.isString(outfile) || !.isValidFilename(basename(outfile)))
  ) {
    warning("value of argument ", sQuote("outfile"),
            " is invalid (see `?findAssociatedSeqs`).\n",
            "Output will be returned, but will not be saved to file."
    )
    outfile <- NULL
  }
  seq_col <- .convertColRef(seq_col, data)
  sample_col <- .convertColRef(sample_col, data)
  subject_col <- .check(subject_col, .isDataColref, default = sample_col,
                        data = data, nse = FALSE, dquote = TRUE
  )
  freq_col <- .check(freq_col, .isDataColref, NULL, ornull = TRUE, data = data)
  freq_col <- .convertColRef(freq_col, data)
  min_seq_length <- .check(min_seq_length, .isNonneg, 7, ornull = TRUE)
  drop_matches <- .check(drop_matches, .isString, "[*|_]", ornull = TRUE)
  min_sample_membership <- .check(min_sample_membership, .isNonneg, 5,
                                  ornull = TRUE
  )
  pval_cutoff <- .check(pval_cutoff, .isPos, 0.05)
  data <- filterInputData(data, seq_col,
                          min_seq_length, drop_matches,
                          subset_cols = NULL,
                          verbose = verbose
  )
  if (nrow(data) < 2) {
    warning(
      "Returning NULL since fewer than two observations remain after ",
      "applying filters for ",
      sQuote("min_sample_membership"), " and ", sQuote("drop_matches")
    )
    return(invisible(NULL))
  }
  msg <- .makemsg(verbose)
  n_subjects <- length(unique(data[ , subject_col]))
  groups <- unique(data[ , group_col])
  rowids_g0 <- data[ , group_col] == groups[[1]]
  n_g0 <- length(unique(data[rowids_g0, subject_col]))
  n_g1 <- n_subjects - n_g0
  samples_or_subjects <- "samples"
  if (subject_col != sample_col) { samples_or_subjects <- "subjects" }
  msg("Data contains ", n_subjects, " ", samples_or_subjects, ", ",
      n_g1, " of which belong to group ", groups[[1]], " and ",
      n_g0, " of which belong to group ", groups[[2]], "."
  )
  out <- .filterSeqsBySampleMembership(data, seq_col, sample_col,
                                       min_sample_membership, msg
  )
  out <- .filterByFisherPvalue(out, data, seq_col, group_col, groups,
                               n_g0, n_g1, pval_cutoff, freq_col,
                               samples_or_subjects, msg
  )
  if (!is.null(outfile)) {
    utils::write.csv(out, outfile, row.names = FALSE)
    msg("Output saved to file:\n  ", outfile)
  }
  invisible(out)
}


findAssociatedClones <- function(
    file_list,
    input_type,
    data_symbols = NULL,
    header, sep, read.args,
    sample_ids = paste0("Sample", 1:length(file_list)),
    subject_ids = NULL,
    group_ids,
    seq_col,
    assoc_seqs,
    nbd_radius = 1,
    dist_type = "hamming",
    min_seq_length = 6,
    drop_matches = NULL,
    subset_cols = NULL,
    output_dir,
    output_type = "rds",
    verbose = FALSE
) {
  if (missing(header)) {
    header <-  switch(input_type, "txt" = FALSE, "table" = FALSE, TRUE)
  }
  if (missing(sep)) {
    sep <- switch(input_type, "csv" = ",", "csv2" = ";", "tsv" = "\t", "")
  }
  if (missing(read.args)) { read.args <- NULL }
  .checkargs.findAssociatedClones(
    file_list, input_type, data_symbols, header, sep, read.args,
    subject_ids, group_ids, seq_col, assoc_seqs, output_dir
  )
  if (input_type %in% c("csv", "csv2", "tsv", "txt", "table")) {
    read.args <- .checkReadArgs(read.args, header, sep)
  }
  .createOutputDir(output_dir)
  .requireOutputDir(output_dir)
  output_type <- .checkOutputType(output_type, "findAssociatedClones")
  dist_type <- .checkDistType(dist_type, "hamming")
  nbd_radius <- .check(nbd_radius, .isNonneg, 1)
  min_seq_length <- .check(min_seq_length, .isNonneg, 6, ornull = TRUE)
  drop_matches <- .check(drop_matches, .isString, NULL, ornull = TRUE)
  subset_cols <- .check(subset_cols, .isCharOrNumericVector, NULL,
                        ornull = TRUE
  )
  sample_ids <- as.character(sample_ids)
  group_ids <- as.character(group_ids)
  if (!is.null(subject_ids)) { subject_ids <- as.character(subject_ids) }
  msg <- .makemsg(verbose)

  msg("Beginning search for associated clones...")
  tmpdirs <- replicate(length(assoc_seqs), tempfile("tmp_assoc_seqs"))
  .createDirectories(tmpdirs)
  for (i in 1:length(file_list)) {
    msg("Processing sample ", i, " of ", length(file_list),
        " (", sample_ids[[i]], ")..."
    )
    .findAssociatedClonesOneSample(input_file = file_list[[i]],
                                   input_type = input_type,
                                   data_symbols = data_symbols,
                                   header = header, sep = sep,
                                   sample_id = sample_ids[[i]],
                                   sample_index = i,
                                   subject_id = subject_ids[[i]],
                                   group_id = group_ids[[i]],
                                   seq_col = seq_col,
                                   assoc_seqs = assoc_seqs,
                                   nbd_radius = nbd_radius,
                                   dist_type = dist_type,
                                   min_seq_length = min_seq_length,
                                   drop_matches = drop_matches,
                                   subset_cols = subset_cols,
                                   output_dirs = tmpdirs,
                                   msg = msg
    )
  }
  msg("Done processing samples. Compiling results...")
  for (i in 1:length(assoc_seqs)) {
    msg("Gathering data from all samples for sequence ", i,
        " (", assoc_seqs[[i]], ")...", newline = FALSE
    )
    .compileNeighborhood(assoc_seq_index = i,
                         input_dir = tmpdirs[[i]],
                         sample_ids = sample_ids,
                         output_dir = output_dir,
                         output_type = output_type,
                         msg = msg
    )
    msg(" Done.")
  }
  unlink(tmpdirs, recursive = TRUE)
  msg("All tasks complete. Output is contained in the following directory:",
      "\n  ", output_dir
  )
  invisible(TRUE)
}



buildAssociatedClusterNetwork <- function(
    file_list,
    input_type = "rds",
    data_symbols = "data",
    header = TRUE, sep,
    read.args = list(row.names = 1),
    seq_col,
    min_seq_length = NULL,
    drop_matches = NULL,
    drop_isolated_nodes = FALSE,
    node_stats = TRUE,
    stats_to_include = chooseNodeStats(cluster_id = TRUE),
    cluster_stats = TRUE,
    color_nodes_by = "GroupID",
    output_name = "AssociatedClusterNetwork",
    verbose = FALSE,
    ...
) {
  drop_isolated_nodes <- .checkTF(drop_isolated_nodes, FALSE)
  node_stats <- .checkTF(node_stats, TRUE)
  cluster_stats <- .checkTF(cluster_stats, TRUE)
  # need to ensure computation of cluster membership
  if (isTRUE(node_stats)) {
    stats_to_include <- .checkStatsToInclude(stats_to_include,
                                             chooseNodeStats(cluster_id = TRUE)
    )
    if (.isLogicalVector(stats_to_include) &&
        .hasElem(stats_to_include, "cluster_id") &&
        !isTRUE(stats_to_include[["cluster_id"]])
    ) {
      stats_to_include[["cluster_id"]] <- TRUE
    }
  } else if (isFALSE(cluster_stats)) {
    node_stats <- TRUE
    stats_to_include <- exclusiveNodeStats(cluster_id = TRUE)
  }
  color_nodes_by <- .check(color_nodes_by, .isCharVector, "GroupID",
                           ornull = TRUE
  )
  output_name <- .checkOutputName(output_name, "AssociatedClusterNetwork")
  msg <- .makemsg(verbose)
  msg("Loading neighborhood data...", newline = FALSE)
  data <- loadDataFromFileList(file_list, input_type, data_symbols,
                               header, sep, read.args
  )
  msg(" Done.")
  msg("Removing duplicates of clones belonging to multiple neighborhoods...",
      newline = FALSE
  )
  tmp_rownames <- rownames(data)
  tmp_rownames <- sapply(
    tmp_rownames,
    # Strip first file*. prefix, which corresponds to assoc nbd file
    function(x) { substr(x, start = 1 + regexpr("\\.", x), stop = nchar(x)) },
    USE.NAMES = FALSE
  )
  # Rownames now identify sample and original row ID
  dupe_idx <- duplicated(tmp_rownames)
  data <- data[!dupe_idx, , drop = FALSE]
  rownames(data) <- tmp_rownames[!dupe_idx]
  msg(" Done.")
  msg("Building global network of associated clusters...")
  buildRepSeqNetwork(data = data, seq_col = seq_col,
                     min_seq_length = min_seq_length,
                     drop_matches = drop_matches,
                     drop_isolated_nodes = drop_isolated_nodes,
                     node_stats = node_stats,
                     stats_to_include = stats_to_include,
                     cluster_stats = cluster_stats,
                     color_nodes_by = color_nodes_by,
                     output_name = output_name,
                     verbose = verbose,
                     ...
  )
}



# Helpers -----------------------------------------------------------------
.findAssociatedClonesOneSample <- function(
    input_file, input_type, data_symbols, header, sep, sample_id,
    sample_index,
    subject_id, group_id, seq_col, assoc_seqs, nbd_radius, dist_type,
    min_seq_length, drop_matches, subset_cols, output_dirs, msg
) {
  data <- .loadDataFromFile(input_file, input_type, data_symbols, header, sep)
  .MUST.isSeqColref(seq_col, data)
  seq_col <- .convertColRef(seq_col, data)
  subset_cols <- .checkDataColrefs(subset_cols, data, NULL)
  subset_cols <- .convertColRef(subset_cols, data)
  data <- filterInputData(data, seq_col,
                          min_seq_length, drop_matches, subset_cols
  )
  data$SampleID <- as.character(sample_id)
  if (!is.null(subject_id)) { data$SubjectID <- as.character(subject_id) }
  data$GroupID <- as.character(group_id)
  msg("Finding clones in a neighborhood of each associated sequence...")
  for (i in 1:length(assoc_seqs)) {
    .getNbdOneSample(seq = assoc_seqs[[i]],
                     data = data, seq_col = seq_col,
                     dist_type = dist_type, nbd_radius = nbd_radius,
                     outfile = file.path(output_dirs[[i]],
                                         paste0(sample_index, ".rds")
                     ),
                     msg = msg,
                     i = i
    )
  }
  msg("Sample complete.")
}


.getNbdOneSample <- function(
    seq, data, seq_col, dist_type, nbd_radius, outfile, msg, i
) {
  nbd <- getNeighborhood(data, seq_col, seq, dist_type, nbd_radius)
  if (is.null(nbd)) {
    msg("Sample does not possess sequence ", i, " (", seq, ").")
  } else {
    msg(nrow(nbd), " clones found in the neighborhood for sequence ", i,
        " (", seq, ")."
    )
    if (nrow(nbd) > 0) {
      nbd$AssocSeq <- seq
      saveRDS(nbd, file = outfile)
    }
  }
}


.compileNeighborhood <- function(
    assoc_seq_index, input_dir, sample_ids, output_dir, output_type, msg
) {
  file_list <- file.path(input_dir, paste0(1:length(sample_ids), ".rds"))
  data <- loadDataFromFileList(file_list[file.exists(file_list)], "rds")
  tmp_rownames <- rownames(data)
  tmp_rownames <- sapply(
    tmp_rownames, # Strip file*. prefix
    function(x) { substr(x, start = 1 + regexpr("\\.", x), stop = nchar(x)) },
    USE.NAMES = FALSE
  )
  rownames(data) <- paste0(data$SampleID, ".", tmp_rownames)
  msg("(", nrow(data), " clones)", newline = FALSE)
  nbd_name <- paste0("assoc_nbd_", assoc_seq_index)
  if (nchar(.sanitizeFilenamePart(data$AssocSeq[[1]])) > 0) {
    nbd_name <- paste0(nbd_name, "_", .sanitizeFilenamePart(data$AssocSeq[[1]]))
  }
  .saveDataGeneric(data, output_dir,
                   output_name = nbd_name,
                   output_type = output_type
  )
}


.filterSeqsBySampleMembership <- function(
    data, seq_col, sample_col, min_sample_membership, msg
) {
  msg("Extracting list of unique sequences... ", newline = FALSE)
  out <- data.frame("ReceptorSeq" = unique(data[[seq_col]]))
  msg("Done. ", nrow(out), " unique sequences present.")
  msg("Computing sample membership (this may take a while)...", newline = FALSE)
  out$shared_by_n_samples <- sapply(
    out$ReceptorSeq,
    function(x) { length(unique(data[data[[seq_col]] == x, sample_col])) }
  )
  msg(" Done.")
  if (is.null(min_sample_membership)) {
    return(out)
  }
  if (min_sample_membership > 0) {
    out <- out[out$shared_by_n_samples >= min_sample_membership, , drop = FALSE]
    msg(nrow(out), " sequences remain after filtering by sample membership.")
    if (nrow(out) == 0) {
      stop("no sequences pass filter for sample membership. ",
           "Try using a lower value of ", sQuote("min_sample_membership")
      )
    }
  }
  out
}

.filterByFisherPvalue <- function(
    unique_seq_data, data, seq_col, group_col, groups,
    n_g0, n_g1, pval_cutoff, freq_col, samples_or_subjects = "subjects", msg
) {
  out <- unique_seq_data
  out$fisher_pvalue <- out$samples_g0 <- out$samples_g1 <- out$label <- NA
  if (samples_or_subjects == "subjects") {
    out$shared_by_n_subjects <- out$subjects_g0 <- out$subjects_g1 <- NA
    cols <- c(
      "ReceptorSeq", "fisher_pvalue",
      "shared_by_n_samples", "samples_g0", "samples_g1",
      "shared_by_n_subjects", "subjects_g0", "subjects_g1", "label"
    )

  } else {
    cols <- c(
      "ReceptorSeq", "fisher_pvalue",
      "shared_by_n_samples", "samples_g0", "samples_g1", "label"
    )
  }
  if (!is.null(freq_col)) {
    out$max_freq <- NA
    cols <- c(cols[1:length(cols)], c("max_freq", "label"))
  }
  out <- out[ , cols]
  msg("Filtering by Fisher's exact test P-value...", newline = FALSE)
  rowids_g0 <- data[[group_col]] == groups[[1]]
  for (i in 1:nrow(out)) {
    rowids_clone <- data[[seq_col]] == out$ReceptorSeq[[i]]
    out$samples_g0[[i]] <- length(
      unique(data[rowids_clone & rowids_g0, "SampleID"])
    )
    out$samples_g1[[i]] <- length(
      unique(data[rowids_clone & !rowids_g0, "SampleID"])
    )
    if (samples_or_subjects == "subjects") {
      out$shared_by_n_subjects[[i]] <- length(
        unique(data[rowids_clone, "SubjectID"])
      )
      out$subjects_g0[[i]] <- length(
        unique(data[rowids_clone & rowids_g0, "SubjectID"])
      )
      out$subjects_g1[[i]] <- length(
        unique(data[rowids_clone & !rowids_g0, "SubjectID"])
      )
    }
    if (samples_or_subjects == "subjects") {
      n_g0_with <- out$subjects_g0[[i]]
      n_g1_with <- out$subjects_g1[[i]]
    } else {
      n_g0_with <- out$samples_g0[[i]]
      n_g1_with <- out$samples_g1[[i]]
    }
    out$fisher_pvalue[[i]] <-
      stats::fisher.test(
        data.frame("g0" = c(n_g1_with, n_g1 - n_g1_with),
                   "g1" = c(n_g0_with, n_g0 - n_g0_with)
        )
      )$p.value
    out$label[[i]] <- paste0(
      "Sequence present in ", out$shared_by_n_samples[[i]], " samples "
    )
    if (samples_or_subjects == "subjects") {
      out$label[[i]] <- paste0(out$label[[i]], "and ",
                               out$shared_by_n_subjects[[i]], " subjects (",
                               n_g0_with, " in group ", groups[[1]], ", ",
                               n_g1_with, " in group ", groups[[2]], ")"
      )
    } else {
      out$label[[i]] <- paste0(out$label[[i]], "(",
                               n_g0_with, " in group ", groups[[1]], ", ",
                               n_g1_with, " in group ", groups[[2]], ")"
      )
    }
    out$label[[i]] <- paste0(out$label[[i]],
                             "\nFisher's exact test P-value: ",
                             signif(out$fisher_pvalue[[i]], digits = 3)
    )
    if (!is.null(freq_col)) {
      out$max_freq[[i]] <- max(data[rowids_clone, freq_col])
      out$label[[i]] <- paste0(
        out$label[[i]], ", Max frequency across all samples: ",
        signif(out$max_freq[[i]], digits = 3)
      )
    }
  }
  out <- out[out$fisher_pvalue < pval_cutoff, , drop = FALSE]
  out <- out[order(out$fisher_pvalue), , drop = FALSE]
  msg(" Done. ", nrow(out), " sequences remain.")
  out
}