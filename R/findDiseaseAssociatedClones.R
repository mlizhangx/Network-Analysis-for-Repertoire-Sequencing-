
# Input:
#     merged sample data
# Do:
#     1. Make filtered clone list:
#         Remove seqs with special characters (*, _, |)
#         Remove duplicate seqs
#         Remove seqs below min length
#     2. Compute "shared by n samples" for filtered clone list
#     3. Filter clones using "shared by n samples" >= specified value
#     4. Compute Fisher P-value and metadata for filtered clones
#     5. Filter clones using Fisher P-value >= specified cutoff
#     6. Sort rows of metadata for filtered clones by Fisher P-value and return
findDiseaseAssociatedClones <- function(
  data,
  clone_col, clone_frac_col,
  pat_id_col, sample_id_col, disease_col,
  disease_encoding,
  min_seq_length = 7, # min length for candidate clone seqs
  min_sample_membership = 5, # min num samples to which a candidate must belong
  pval_cutoff = 0.05, # only keep candidates with fisher PV below threshold
  drop_chars = "[*|_]" # ignore seqs with these chars (regular expression)
) {

  # Drop sequences with specified chars from internal copy of merged sample data
  if (!is.null(drop_chars)) {
    data <- data[-grep(drop_chars, data[ , clone_col]), ] }

  # Drop sequences with length below the specified minimum
  data <- data[nchar(data[ , clone_col]) >= min_seq_length, ]

  cat(paste0(
    "Data contains ", nrow(data), " clones for which the designated clonotype sequence in column '", clone_col,
    "' has length at least ", min_seq_length, " and is free of special characters '*', '|' and '_'.\n"))

  # Get counts of number of total/disease/healthy patients
  n_pt_total <- length(unique(data[ , pat_id_col]))
  rowids_disease_pt <- data[ , disease_col] %in% disease_encoding
  rowids_healthy_pt <- !rowids_disease_pt
  n_pt_disease_total <-
    length(unique(data[rowids_disease_pt, pat_id_col]))
  n_pt_healthy_total <- n_pt_total - n_pt_disease_total
  cat(paste0("Data contains a total of ", n_pt_total, " patients: ",
             n_pt_disease_total, " disease patients and ",
             n_pt_healthy_total, " healthy patients.\n"))

  # Get unique clone sequences
  out <- data.frame("clone_seq" = # list of unique clone seqs
                      as.character(unique(data[ , clone_col])))
  cat(paste0(nrow(out), " unique clononotype sequences identified. Computing sample membership for these sequences...\n"))
  out$seq_length <- nchar(out$clone_seq)

  # Compute variable "shared by n samples" for remaining sequences
  out$shared_by_n_samples <-
    sapply(out$clone_seq,
           function(x) {
             length(unique(data[data[ , clone_col] == x, sample_id_col])) })

  # Drop sequences shared by fewer samples than the specified minimum
  out <- out[out$shared_by_n_samples >= min_sample_membership, ]
  cat(paste0(nrow(out), " unique clone sequences meet the requirement of appearing in at least ", min_sample_membership, " samples. Performing Fisher's exact tests for these sequences...\n"))

  ### GET FISHER PVALUE AND OTHER METADATA FOR REMAINING SEQUENCES ###
  # Add variable for "shared by n healthy patients" to candidate sequence data
  out$shared_by_n_pt_healthy <- rep(NA, nrow(out))
  # Add variable for "shared by n disease patients" to candidate sequence data
  out$shared_by_n_pt_disease <- rep(NA, nrow(out))
  # Add variable for Fisher P-value to candidate sequence data
  out$pv_fisher <- rep(NA, nrow(out))
  # Add variable for max clone fraction to candidate sequence data
  out$max_clonefrac <- rep(NA, nrow(out))
  # Add variable for description of sequence/meta
  out$label <- rep(NA, nrow(out))

  for (i in 1:nrow(out)) {
    # logical vector for rows of merged data corresponding to current clone seq
    # grepl(pattern = paste0("^", clone, "$"), x = data[ , clone_col])
    rowids_clone <- data[ , clone_col] == out$clone_seq[[i]]

    # Number of disease patients with candidate sequence
    out$shared_by_n_pt_disease[[i]] <-
      length(unique(
        data[rowids_clone & rowids_disease_pt, pat_id_col]))
    # Number of healthy patients with candidate sequence
    out$shared_by_n_pt_healthy[[i]] <-
      length(unique(
        data[rowids_clone & rowids_healthy_pt, pat_id_col]))
    # Compute fisher p-value for candidate sequence
    out$pv_fisher[[i]] <-
      stats::fisher.test(data.frame(
        "healthy" = c(out$shared_by_n_pt_healthy[[i]],
                      n_pt_healthy_total - out$shared_by_n_pt_healthy[[i]]),
        "disease" = c(out$shared_by_n_pt_disease[[i]],
                      n_pt_disease_total - out$shared_by_n_pt_disease[[i]]))
      )$p.value
    # Compute max clone fraction of candidate seq among samples containing it
    out$max_clonefrac[[i]] <- max(data[rowids_clone, clone_frac_col])

    # Description of current candidate sequence (e.g., for plot subtitle)
    out$label[[i]] <- paste0(
      "Sequence present in ", out$shared_by_n_samples[[i]], " samples and ",
      out$shared_by_n_pt_disease[[i]] + out$shared_by_n_pt_healthy[[i]],
      " patients (", out$shared_by_n_pt_disease[[i]], " of which have positive disease status)",
      "\nFisher P-value = ", signif(out$pv_fisher[[i]], digits = 3),
      ", Max Clone Fraction in Network = ", signif(out$max_clonefrac[[i]], digits = 3))
  }
  # Drop sequences with Fisher P-value above specified cutoff
  out <- out[out$pv_fisher < pval_cutoff, ]
  cat(paste0("Tests complete. Found ", nrow(out), " sequences with Fisher's exact test P-value below ", pval_cutoff,
  ". Returning metadata for these sequences, with rows sorted by P-value.\n"))
  # Sort candidate sequence metadata by fisher P-value and return
  return(out[order(out$pv_fisher), ])
}
