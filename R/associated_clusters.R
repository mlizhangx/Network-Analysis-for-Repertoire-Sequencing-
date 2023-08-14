# Top-Level Functions -----------------------------------------------------

findAssociatedSeqs <- function(
    file_list,
    input_type,
    data_symbols = NULL,
    header = TRUE, sep = "",
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
    outfile = "associated_seqs.csv"
) {
  .checkDeprecated.findAssociatedSeqs(sample_ids, groups)
  .checkargs.findAssociatedSeqs(
    file_list, input_type, data_symbols, header, sep,
    subject_ids, group_ids, seq_col, freq_col,
    min_seq_length, drop_matches, min_sample_membership, pval_cutoff, outfile
  )
  group_ids <- as.character(group_ids)
  if (!is.null(subject_ids)) { subject_ids <- as.character(subject_ids) }
  samples_or_subjects <- "samples"
  if (!is.null(subject_ids)) {
    if (any(duplicated(subject_ids))) { samples_or_subjects <- "subjects" }
  }
  groups <- unique(group_ids)
  stopifnot(
    "group_ids must contain exactly two unique values" = length(groups) == 2
  )
  n_samples <- length(group_ids)
  ids_g0 <- group_ids == groups[[1]]
  ids_g1 <- group_ids == groups[[2]]
  n_samples_g0 <- sum(ids_g0)
  n_samples_g1 <- sum(ids_g1)
  if (samples_or_subjects == "subjects") {
    n_g0 <- length(unique(subject_ids[ids_g0])) # num subjects
    n_g1 <- length(unique(subject_ids[ids_g1])) # num subjects
    n <- length(unique(subject_ids))            # num subjects
    cat(paste0("Data contains ", n_samples, " samples and ", n, " subjects, ",
               n_g0, " of which belong to group ", groups[[1]], " and ",
               n_g1, " of which belong to group ", groups[[2]], ".\n"
    ))
  } else {
    n_g0 <- n_samples_g0
    n_g1 <- n_samples_g1
    n <- n_samples
    cat(paste0("Data contains ", n, " samples, ",
               n_samples_g0, " of which belong to group ", groups[[1]], " and ",
               n_samples_g1, " of which belong to group ", groups[[2]], ".\n"
    ))
  }
  data <- combineSamples(
    file_list = file_list, input_type = input_type, data_symbols = data_symbols,
    header = header, sep = sep, seq_col = seq_col,
    min_seq_length = min_seq_length, drop_matches = drop_matches,
    sample_ids = 1:length(file_list), subject_ids = subject_ids,
    group_ids = group_ids, subset_cols = c(seq_col, freq_col)
  )
  if (nrow(data) < 2) {
    warning(paste(
      "insufficient remaining receptor sequences; at least two needed",
      "to proceed. Aborting search for associated sequences"
    ))
    return(invisible(NULL))
  }
  out <- .filterSeqsBySampleMembership(
    data, seq_col, sample_col = "SampleID", min_sample_membership
  )
  out <- .filterByFisherPvalue(
    out, data, seq_col, group_col = "GroupID", groups,
    n_g0, n_g1, pval_cutoff, freq_col, samples_or_subjects
  )
  if (!is.null(outfile)) {
    utils::write.csv(out, outfile, row.names = FALSE)
    cat(paste0("Output saved to file:\n  ", outfile, "\n"))
  }
  cat(
    "All done. Sorting results by Fisher's exact test P-value and returning.\n"
  )
  invisible(out)
}

findAssociatedSeqs2 <- function(
    data,
    seq_col,
    sample_col,
    subject_col = sample_col,
    group_col,
    groups = NULL,
    freq_col = NULL,
    min_seq_length = 7,
    drop_matches = "[*|_]",
    min_sample_membership = 5,
    pval_cutoff = 0.05,
    outfile = "associated_seqs.csv"
) {
  data <- as.data.frame(data)
  .checkargs.findAssociatedSeqs2(
    data, seq_col, sample_col, subject_col, group_col,
    groups, freq_col,
    min_seq_length, drop_matches, min_sample_membership, pval_cutoff, outfile
  )
  seq_col <- .convertColRef(seq_col, data)
  freq_col <- .convertColRef(freq_col, data)
  data <- filterInputData(
    data, seq_col, min_seq_length, drop_matches, subset_cols = NULL
  )
  if (nrow(data) < 2) {
    warning(paste(
      "insufficient remaining receptor sequences; at least two needed",
      "to proceed. Aborting search for associated sequences"
    ))
    return(invisible(NULL))
  }
  samples_or_subjects <- "samples"
  if (subject_col != sample_col) { samples_or_subjects <- "subjects" }
  groups <- unique(data[ , group_col])
  stopifnot(
    "`group_col` of `data` must contain exactly two unique values" =
      length(groups) == 2
  )
  n_subjects <- length(unique(data[ , subject_col]))
  rowids_g0 <- data[ , group_col] == groups[[1]]
  n_g0 <- length(unique(data[rowids_g0, subject_col]))
  n_g1 <- n_subjects - n_g0
  cat(paste0("Data contains ", n_subjects, " ", samples_or_subjects, ", ",
             n_g1, " of which belong to group ", groups[[1]], " and ",
             n_g0, " of which belong to group ", groups[[2]], ".\n"
  ))
  cat("All samples loaded. ")
  out <- .filterSeqsBySampleMembership(
    data, seq_col, sample_col, min_sample_membership
  )
  out <- .filterByFisherPvalue(
    out, data, seq_col, group_col, groups,
    n_g0, n_g1, pval_cutoff, freq_col, samples_or_subjects
  )
  if (!is.null(outfile)) {
    utils::write.csv(out, outfile, row.names = FALSE)
    cat(paste0("Output saved to file:\n  ", outfile, "\n"))
  }
  cat(
    "All done. Sorting results by Fisher's exact test P-value and returning.\n"
  )
  invisible(out)
}


findAssociatedClones <- function(
    file_list,
    input_type,
    data_symbols = NULL,
    header = TRUE,
    sep = "",
    sample_ids = paste0("Sample", 1:length(file_list)),
    subject_ids = NULL,
    group_ids,
    seq_col,
    assoc_seqs,
    nbd_radius = 1,
    dist_type = "hamming",
    min_seq_length = 6,
    drop_matches = "[*|_]",
    subset_cols = NULL,
    output_dir = file.path(getwd(), "associated_neighborhoods"),
    output_type = "csv",
    verbose = FALSE
) {
  .checkargs.findAssociatedClones(
    file_list, input_type, data_symbols, header, sep,
    sample_ids, subject_ids, group_ids, seq_col,
    assoc_seqs, nbd_radius, dist_type, min_seq_length, drop_matches,
    subset_cols, output_dir, output_type, verbose
  )
  .ensureOutputDir(output_dir)
  sample_ids <- as.character(sample_ids)
  group_ids <- as.character(group_ids)
  if (!is.null(subject_ids)) { subject_ids <- as.character(subject_ids) }

  cat(paste0("<<< Beginning search for associated clones >>>\n"))
  tmpdirs <- file.path(tempdir(), assoc_seqs)
  .createDirectories(tmpdirs)
  for (i in 1:length(file_list)) {
    cat(paste0("Processing sample ", i, " of ", length(file_list),
               " (", sample_ids[[i]], "):\n"
    ))
    .findAssociatedClonesOneSample(
      file_list[[i]], input_type, data_symbols, header, sep, sample_ids[[i]],
      subject_ids[[i]], group_ids[[i]], seq_col, assoc_seqs, nbd_radius,
      dist_type, min_seq_length, drop_matches, subset_cols, tmpdirs, verbose
    )
  }
  cat(paste0(">>> Done processing samples. Compiling results:\n"))
  for (i in 1:length(assoc_seqs)) {
    cat(paste0("Gathering data from all samples for sequence ", i,
               " (", assoc_seqs[[i]], ")..."
    ))
    .compileNeighborhood(
      tmpdirs[[i]], sample_ids, output_dir, output_type, verbose
    )
    cat(" Done.\n")
  }
  unlink(tmpdirs, recursive = TRUE)
  cat(paste0(
    ">>> All tasks complete. Output is contained in the following directory:",
    "\n  ", output_dir, "\n"
  ))
  invisible(TRUE)
}


# findAssociatedClones2 <- function(
    #
#   ## Input ##
#   data, sample_col, seq_col,
#
#   ## Search Criteria ##
#   assoc_seqs, nbd_radius = 1, dist_type = "hamming",
#   min_seq_length = 6, drop_matches = "[*|_]",
#
#   ## Output ##
#   subset_cols = NULL,
#   output_dir = file.path(getwd(), "associated_neighborhoods"),
#   output_type = "csv",
#   verbose = FALSE
# ) {
#   .ensureOutputDir(output_dir)
#
#   cat(paste0("<<< Beginning search for associated clones >>>\n"))
#
#   tmpdir <- tempdir()
#   tmpdirs <- file.path(tmpdir, assoc_seqs); .createDirectories(tmpdirs)
#
#   sample_col <- .convertColRef(sample_col, data)
#   sample_list <- unique(data[ , sample_col])
#
#   for (i in 1:length(sample_list)) {
#     cat(paste0("Processing sample ", i, " of ", length(sample_list), " (", sample_list[[i]], "):\n"))
#     rows_current_sample <- data[ , sample_col] == sample_list[[i]]
#     .findAssociatedClonesOneSample2(
#       data[rows_current_sample, ], seq_col, assoc_seqs, nbd_radius,
#       dist_type, min_seq_length, drop_matches, subset_cols, tmpdirs, verbose)
#   }
#
#   cat(paste0(">>> Done processing samples. Compiling results:\n"))
#
#   for (i in 1:length(assoc_seqs)) {
#     cat(paste0("Gathering data from all samples for sequence ", i, " (", assoc_seqs[[i]], ")..."))
#     .compileNeighborhood(tmpdirs[[i]], sample_list, output_dir, output_type, verbose)
#     cat(" Done.\n")
#   }
#
#   unlink(tmpdir)
#   cat(paste0(">>> All tasks complete. Output is contained in the following directory:\n  ", output_dir, "\n"))
# }



buildAssociatedClusterNetwork <- function(
    file_list,
    input_type = "csv",
    data_symbols = "data",
    header = TRUE, sep = ",",
    seq_col,
    min_seq_length = NULL,
    drop_matches = NULL,
    drop_isolated_nodes = FALSE,
    node_stats = TRUE,
    stats_to_include = chooseNodeStats(cluster_id = TRUE),
    cluster_stats = TRUE,
    color_nodes_by = "GroupID",
    output_name = "AssociatedClusterNetwork",
    ...
) {
  .checkargs.buildAssociatedClusterNetwork(
    file_list, input_type, data_symbols, header, sep,
    seq_col, min_seq_length, drop_matches, drop_isolated_nodes,
    node_stats, stats_to_include, cluster_stats, color_nodes_by, output_name
  )
  if (typeof(stats_to_include) %in% c("list", "logical") &&
      !stats_to_include[["cluster_id"]]) {
    stats_to_include[["cluster_id"]] <- TRUE
  }
  data <- loadDataFromFileList(file_list, input_type, data_symbols, header, sep)
  cat("<<< Building network of associated clones >>>\n")
  net <- buildRepSeqNetwork(
    data = data, seq_col = seq_col, min_seq_length = min_seq_length,
    drop_matches = drop_matches, drop_isolated_nodes = drop_isolated_nodes,
    node_stats = node_stats, stats_to_include = stats_to_include,
    cluster_stats = cluster_stats, color_nodes_by = color_nodes_by,
    output_name = output_name, ...
  )
  if (is.null(net)) {
    return(NULL)
  }
  invisible(net)
}



# Helpers -----------------------------------------------------------------
.findAssociatedClonesOneSample <- function(
    input_file, input_type, data_symbols, header, sep, sample_id,
    subject_id, group_id, seq_col, assoc_seqs, nbd_radius, dist_type,
    min_seq_length, drop_matches, subset_cols, output_dirs, verbose
) {
  data <- .loadDataFromFile(input_file, input_type, data_symbols, header, sep)
  seq_col <- .convertColRef(seq_col, data)
  subset_cols <- .convertColRef(subset_cols, data)
  data <- filterInputData(
    data, seq_col, min_seq_length, drop_matches, subset_cols
  )
  data$SampleID <- as.character(sample_id)
  if (!is.null(subject_id)) { data$SubjectID <- as.character(subject_id) }
  data$GroupID <- as.character(group_id)
  cat("Finding clones in a neighborhood of each associated sequence...")
  for (i in 1:length(assoc_seqs)) {
    .getNbdOneSample(
      assoc_seqs[[i]], data, seq_col, dist_type, nbd_radius,
      file.path(output_dirs[[i]], paste0(sample_id, ".rds")), verbose, i
    )
  }
  if (verbose) { cat("\n") }
  cat(" Done.\n")
}

# .findAssociatedClonesOneSample2 <- function(
    #     data, seq_col, assoc_seqs, nbd_radius, dist_type,
#     min_seq_length, drop_matches, subset_cols, output_dirs, verbose)
# {
#   seq_col <- .convertColRef(seq_col, data)
#   subset_cols <- .convertColRef(subset_cols, data)
#   data <- filterInputData(data, seq_col, min_seq_length, drop_matches,
#                           subset_cols)
#
#   cat("Finding clones in a neighborhood of each associated sequence...")
#   for (i in 1:length(assoc_seqs)) {
#     # cat(paste0("Finding clones in a neighborhood of sequence ", i, " (", assoc_seqs[[i]], ")..."))
#     .getNbdOneSample(
#       assoc_seqs[[i]], data, seq_col, dist_type, nbd_radius,
#       file.path(output_dirs[[i]], paste0(sample_id, ".rds")), verbose, i)
#   }
#   if (verbose) { cat("\n") }
#   cat(" Done.\n")
# }


.getNbdOneSample <- function(
    seq, data, seq_col, dist_type, nbd_radius, outfile, verbose, i
) {
  nbd <- getNeighborhood(data, seq_col, seq, dist_type, nbd_radius)
  if (is.null(nbd)) {
    if (verbose) {
      cat(paste0("\nSample does not possess sequence ", i, " (", seq, ")."))
    }
  } else {
    if (verbose) {
      cat(paste0(
        "\n", nrow(nbd), " clones found in the neighborhood for sequence ", i,
        " (", seq, ")."
      ))
    }
    if (nrow(nbd) > 0) {
      nbd$AssocSeq <- seq
      saveRDS(nbd, file = outfile)
    }
  }
}


.compileNeighborhood <- function(
    input_dir, sample_ids, output_dir, output_type, verbose
) {
  file_list <- file.path(input_dir, paste0(sample_ids, ".rds"))
  data <- loadDataFromFileList(file_list[file.exists(file_list)], "rds")
  if (verbose) { cat("(", paste0(nrow(data), " clones)")) }
  .saveDataGeneric(data, output_dir,
                   output_name = data$AssocSeq[[1]],
                   output_type = output_type
  )
}


.filterSeqsBySampleMembership <- function(
    data, seq_col, sample_col, min_sample_membership
) {
  cat(paste0("Extracting list of unique sequences... "))
  out <- data.frame("ReceptorSeq" = unique(data[[seq_col]]))
  cat("Done. ", paste0(nrow(out), " unique sequences present.\n"))
  cat("Computing sample membership (this may take a while)...")
  out$shared_by_n_samples <- sapply(
    out$ReceptorSeq,
    function(x) { length(unique(data[data[[seq_col]] == x, sample_col])) }
  )
  cat(" Done.\n")
  if (is.null(min_sample_membership)) {
    return(out)
  }
  if (min_sample_membership > 0) {
    out <- out[out$shared_by_n_samples >= min_sample_membership, , drop = FALSE]
    cat(paste0(
      nrow(out), " sequences remain after filtering by sample membership.\n"
    ))
    if (nrow(out) == 0) {
      stop(paste("no sequences pass filter for sample membership.",
                 "Try using a lower value of min_sample_membership"
      ))
    }
  }
  out
}

.filterByFisherPvalue <- function(
    unique_seq_data, data, seq_col, group_col, groups,
    n_g0, n_g1, pval_cutoff, freq_col, samples_or_subjects = "subjects"
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
  out <- out[ , cols]
  cat("Filtering by Fisher's exact test P-value...")
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
      out$label[[i]] <- paste0(
        out$label[[i]], ", Max frequency across all samples: ",
        signif(max(data[rowids_clone, freq_col]), digits = 3)
      )
    }
  }
  out <- out[out$fisher_pvalue < pval_cutoff, , drop = FALSE]
  out <- out[order(out$fisher_pvalue), , drop = FALSE]
  cat(paste0(" Done. ", nrow(out), " sequences remain.\n"))
  out
}