
# Input:
#   combined bulk rep-seq data from multiple samples
# Do:
#   Get potential associated tcr seqs filtered and ranked by Fisher PV:
#   1. Filter clone list:
#       Remove seqs with special characters (*, _, |)
#       Remove duplicate seqs
#       Remove seqs below min length
#       Remove seqs shared by fewer than min # samples
#   2. Compute Fisher PV for filtered clones
#   3. Sort and filter list of potential associated clones by Fisher PV

getPotentialAssociatedClones <- function(

  data,
  clone_col,
  freq_col,
  sample_col,
  subject_col = sample_col,
  group_col,
  treatment_groups,
  min_seq_length = 7,
  min_sample_membership = 5,
  pval_cutoff = 0.05,
  drop_chars = "[*|_]"

) {

  ### 1. FILTER CLONE LIST ###

  # Determine if data distinguishes subject from sample
  samples_or_subjects <- "samples"
  if (subject_col != sample_col) { samples_or_subjects <- "subjects" }

  # Drop sequences with specified chars
  if (!is.null(drop_chars)) {
    data <- data[-grep(drop_chars, data[ , clone_col]), ] }

  # Drop sequences below specified length
  data <- data[nchar(data[ , clone_col]) >= min_seq_length, ]
  cat(paste0(
    "Data contains ", nrow(data), " clones after filtering for minimum sequence length and special characters.\n"))

  # Get subject treatment/control counts
  n_subjects <- length(unique(data[ , subject_col]))
  rowids_treatment_subjects <- data[ , group_col] %in% treatment_groups
  n_treatment_subjects <-
    length(unique(data[rowids_treatment_subjects, subject_col]))
  n_control_subjects <- n_subjects - n_treatment_subjects
  cat(paste0(
    "Data contains a total of ", n_subjects, " ", samples_or_subjects,
    ", of which ", n_treatment_subjects, " are labeled as treatment and ",
    n_control_subjects, " are labeled as control.\n"))

  # Get unique clone sequences and initialize output
  out <- data.frame("cloneSeq" = # list of unique clone seqs
                      as.character(unique(data[ , clone_col])))
  out$pv_fisher <- rep(NA, nrow(out))
  out$shared_by_n_samples <- rep(NA, nrow(out))
  # out$shared_by_n_treatment_subjects <- rep(NA, nrow(out))
  # out$max_freq <- rep(NA, nrow(out))
  out$label <- rep(NA, nrow(out))


  ### 2. FILTER CLONES BY SAMPLE MEMBERSHIP ###
  cat(paste0(nrow(out), " unique clononotype sequences identified. Computing sample membership for these sequences...\n"))
  out$shared_by_n_samples <-
    sapply(out$cloneSeq,
           function(x) {
             length(unique(data[data[ , clone_col] == x, sample_col])) })

  # Drop sequences shared by fewer samples than the specified minimum
  out <- out[out$shared_by_n_samples >= min_sample_membership, ]
  cat(paste0(nrow(out), " unique clone sequences meet the requirement of appearing in at least ", min_sample_membership, " samples. Performing Fisher's exact tests for these sequences...\n"))


  ### 3. PERFORM FISHER'S EXACT TESTS FOR FILTERED CLONES ###
  # Iterate over filtered list of clones
  for (i in 1:nrow(out)) {
    # logical vector for rows of merged data corresponding to current clone seq
    # grepl(pattern = paste0("^", clone, "$"), x = data[ , clone_col])
    rowids_clone <- data[ , clone_col] == out$cloneSeq[[i]]

    # Number of treatment subjects with target sequence
    shared_by_n_treatment_subjects <-
      length(unique(
        data[rowids_clone & rowids_treatment_subjects, subject_col]))
    # Number of control subjects with target sequence
    shared_by_n_control_subjects <-
      length(unique(
        data[rowids_clone & !rowids_treatment_subjects, subject_col]))

    # Compute fisher p-value for target sequence
    out$pv_fisher[[i]] <-
      stats::fisher.test(data.frame(
        "control" =
          c(shared_by_n_control_subjects,
            n_control_subjects - shared_by_n_control_subjects),
        "treatment" =
          c(shared_by_n_treatment_subjects,
            n_treatment_subjects - shared_by_n_treatment_subjects))
      )$p.value

    # Compute max clone fraction of target seq among samples containing it
    # out$max_freq[[i]] <- max(data[rowids_clone, freq_col])

    # Description of current target sequence (e.g., for plot subtitle)
    out$label[[i]] <- paste0(
      "Sequence present in ", out$shared_by_n_samples[[i]], " samples ")
    if (samples_or_subjects == "subjects") {
      out$label[[i]] <- paste0(
        out$label[[i]],
        "and ", shared_by_n_treatment_subjects + shared_by_n_control_subjects, " subjects ")
    }
    out$label[[i]] <- paste0(
      out$label[[i]], "(of which ", shared_by_n_treatment_subjects, " are categorized as treatment", samples_or_subjects, ")",
      "\nFisher's exact test P-value: ", signif(out$pv_fisher[[i]], digits = 3),
      ", Max clone frequency across all samples: ", signif(max(data[rowids_clone, freq_col]), digits = 3))
  }

  ### 4. SORT AND FILTER POTENTIAL CLONE LIST BY FISHER P-VALUE ###
  # Drop sequences with Fisher P-value above specified cutoff
  out <- out[out$pv_fisher < pval_cutoff, ]
  cat(paste0("Tests complete. Found ", nrow(out), " sequences with Fisher's exact test P-value below ", pval_cutoff,
             ". Sorting these sequences by P-value and returning them in a data frame containing the P-value, total sample membership and description for each."))
  # Sort candidate sequence metadata by fisher P-value and return
  return(out[order(out$pv_fisher), ])
}
