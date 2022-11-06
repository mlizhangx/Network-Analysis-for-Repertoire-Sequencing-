
# Input:
#   Combined bulk rep-seq data from multiple samples
# Do:
#   Get potential associated tcr seqs filtered and ranked by Fisher PV:
#   1. Filter clone list:
#       Remove seqs with special characters (*, _, |)
#       Remove duplicate seqs
#       Remove seqs below min length
#       Remove seqs shared by fewer than min # samples
#   2. Compute Fisher PV for filtered clones
#   3. Sort and filter list of potential associated clones by Fisher PV

findAssociatedClones <- function(

  data,
  seq_col,
  freq_col,
  sample_col,
  subject_col = sample_col,
  group_col,
  groups = c("group0", "group1"),  # 0 = reference, 1 = comparison
  min_seq_length = 7,
  min_sample_membership = 5,
  pval_cutoff = 0.05,
  drop_chars = "[*|_]",
  outfile = "potential_associated_clones.csv"

) {


  # clone_seq_col <- amino_col
  # if (clone_seq_type == "nucleotide") { clone_seq_col <- nucleo_col }

  #### FILTER CLONES ####
  # Coerce sequence column to character if needed
  if (!is.character(data[ , seq_col])) {
    data[ , seq_col] <- as.character(data[ , seq_col]) }

  # Determine if data distinguishes subject from sample
  samples_or_subjects <- "samples"
  if (subject_col != sample_col) { samples_or_subjects <- "subjects" }

  # Filter by seq length
  if (!is.null(min_seq_length)) {
    cat(paste0("Removing sequences with length less than ", min_seq_length, "..."))
    data <- filterClonesBySequenceLength(data, seq_col,
                                         min_length = min_seq_length)
    cat(paste0(" Done. ", nrow(data), " rows remaining.\n"))
  }

  # Filter seqs with special chars
  if (!is.null(drop_chars)) {
    cat(paste0("Removing sequences containing matches to the expression '", drop_chars, "'..."))
    drop_matches <- grep(drop_chars, data[ , seq_col])
    if (length(drop_matches) > 0) { data <- data[-drop_matches, ] }
    cat(paste0(" Done. ", nrow(data), " rows remaining.\n")) }

  # Get subject case/control counts
  n_subjects <- length(unique(data[ , subject_col]))
  rowids_case_subjects <- data[ , group_col] == groups[[1]]
  n_case_subjects <-
    length(unique(data[rowids_case_subjects, subject_col]))
  n_control_subjects <- n_subjects - n_case_subjects
  cat(paste0(
    "Data contains ", n_subjects, " ", samples_or_subjects, ", ",
    n_control_subjects, " of which belong to the reference group and ",
    n_case_subjects, " of which belong to the comparison group.\n"))

  # Get unique clone sequences and initialize output
  out <- data.frame("ReceptorSeq" = # list of unique clone seqs
                      unique(data[ , seq_col]))
  # out <- data.frame("AminoAcidSeq" = # list of unique clone seqs
  #                     as.character(unique(data[ , seq_col])))
  # output_clone_seq_col <- "AminoAcidSeq"
  # if (clone_seq_type == "nucleotide") {
  #   names(out)[1] <- output_clone_seq_col <- "NucleotideSeq"
  # }
  out$fisher_pvalue <- rep(NA, nrow(out))
  out$shared_by_n_samples <- rep(NA, nrow(out))
  # out$shared_by_n_case_subjects <- rep(NA, nrow(out))
  # out$max_freq <- rep(NA, nrow(out))
  out$label <- rep(NA, nrow(out))


  #### SAMPLE MEMBERSHIP ####
  cat(paste0(nrow(out), " unique receptor sequences found. Computing sample membership (this could take a while)..."))
  out$shared_by_n_samples <-
    sapply(out$ReceptorSeq,
           # sapply(out[ , output_clone_seq_col],
           function(x) {
             length(unique(data[data[ , seq_col] == x, sample_col])) })
  cat(" Done.\n")

  # Drop sequences shared by fewer samples than the specified minimum
  out <- out[out$shared_by_n_samples >= min_sample_membership, ]
  cat(paste0(nrow(out), " sequences remain after filtering by sample membership.\n"))
  stopifnot("no sequences pass filter for sample membership" = nrow(out) > 0)

  #### FISHER'S EXACT TESTS ####
  cat("Performing Fisher's exact tests...")
  # Iterate over filtered list of clones
  for (i in 1:nrow(out)) {
    # logical vector for rows of merged data corresponding to current clone seq
    # grepl(pattern = paste0("^", clone, "$"), x = data[ , seq_col])
    rowids_clone <- data[ , seq_col] == out$ReceptorSeq[[i]]
    # rowids_clone <- data[ , seq_col] == out[ , output_clone_seq_col][[i]]

    # Number of case subjects with target sequence
    shared_by_n_case_subjects <-
      length(unique(
        data[rowids_clone & rowids_case_subjects, subject_col]))
    # Number of control subjects with target sequence
    shared_by_n_control_subjects <-
      length(unique(
        data[rowids_clone & !rowids_case_subjects, subject_col]))

    # Compute fisher p-value for target sequence
    out$fisher_pvalue[[i]] <-
      stats::fisher.test(data.frame(
        "control" =
          c(shared_by_n_control_subjects,
            n_control_subjects - shared_by_n_control_subjects),
        "case" =
          c(shared_by_n_case_subjects,
            n_case_subjects - shared_by_n_case_subjects))
      )$p.value

    # Compute max clone fraction of target seq among samples containing it
    # out$max_freq[[i]] <- max(data[rowids_clone, freq_col])

    # Description of current target sequence (e.g., for plot subtitle)
    out$label[[i]] <- paste0(
      "Sequence present in ", out$shared_by_n_samples[[i]], " samples ")
    if (samples_or_subjects == "subjects") {
      out$label[[i]] <- paste0(
        out$label[[i]],
        "and ", shared_by_n_case_subjects + shared_by_n_control_subjects, " subjects ")
    }
    out$label[[i]] <- paste0(
      out$label[[i]], "(of which ", shared_by_n_case_subjects, " are in the comparison group)",
      "\nFisher's exact test P-value: ", signif(out$fisher_pvalue[[i]], digits = 3),
      ", Max frequency across all samples: ", signif(max(data[rowids_clone, freq_col]), digits = 3))
  }
  cat(" Done.\n")


  #### FILTER/SORT BY P-VALUE ####

  # Drop sequences with Fisher P-value above specified cutoff
  out <- out[out$fisher_pvalue < pval_cutoff, ]
  cat(paste0(nrow(out), " sequences remain after filtering by Fisher's exact test P-value.\n"))

  # Sort candidate sequence metadata by fisher P-value and return
  cat("Sorting the data for these sequences by Fisher's exact test P-value and returning.\n")
  out <- out[order(out$fisher_pvalue), ]
  if (!is.null(outfile)) {
    utils::write.csv(out, outfile, row.names = FALSE)
    cat(paste0("Output saved to file:\n  ", outfile, "\n"))
  }
  return(out)
}
