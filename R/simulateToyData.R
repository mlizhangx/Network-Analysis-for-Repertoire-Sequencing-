simulateToyData <- function(
    samples = 2,
    chains = 1,
    sample_size = 100,
    prefix_length = 7,
    prefix_chars = c("G", "A", "T", "C"),
    prefix_probs = rbind(
      "sample1" = c(12, 4, 1, 1),
      "sample2" = c(4, 12, 1, 1)),
    affixes = c("AATTGG", "AATCGG", "AATTCG",
                "AATTGC", "AATTG", "AATTC"),
    affix_probs = rbind(
      "sample1" = c(10, 4, 2, 2, 1, 1),
      "sample2" = c(1, 1, 1, 2, 2.5, 2.5)),
    num_edits = 0,
    edit_pos_probs = function(seq_length) {
      stats::dnorm(seq(-4, 4, length.out = seq_length))
    },
    edit_ops = c("insertion", "deletion", "transmutation"),
    edit_probs = c(5, 1, 4),
    new_chars = prefix_chars,
    new_probs = prefix_probs,
    include_counts = FALSE,
    output_dir = NULL,
    no_return = FALSE,
    seed_value = 42
) {

  set.seed(seed_value)
  seqs <- character(samples * sample_size)
  ranges <- (matrix(seq_len(samples * sample_size), nrow = samples))
  if (num_edits > 0) {
    edit_op <- function(seq, pos, char) {
      op <- sample(edit_ops, size = 1, prob = edit_probs)
      pos2 <- pos + 1; pos3 <- nchar(seq)
      if (op == "deletion") { char <- NULL }
      if (op != "insertion") { pos <- pos - 1}
      seq <- paste0(substr(seq, 1, pos), char, substr(seq, pos2, pos3))
    }
  }

  for (i in 1:samples) {
    prefix <- apply(matrix(sample(prefix_chars,
                                  size = prefix_length * sample_size,
                                  replace = TRUE, prob = prefix_probs[i, ]),
                           ncol = sample_size),
                    MARGIN = 2,
                    FUN = function(x) paste0(x, collapse = ""))
    seqs[ranges[ , i]] <-
      paste0(prefix, sample(affixes, size = sample_size,
                            replace = TRUE, prob = affix_probs[i, ]))
    if (num_edits > 0) {
      for (j in 1:num_edits) {
        char_draws <- sample(new_chars, size = sample_size, replace = TRUE,
                             prob = new_probs[i, ])
        edit_pos <- sapply(seqs[ranges[ , i]],
                           function(x) {
                             sample(seq_len(nchar(x)), size = 1,
                                    prob = edit_pos_probs(nchar(x)))
                           })
        tmp <- c(seqs[ranges[ , i]], as.character(edit_pos), char_draws)
        tmp <- matrix(tmp, ncol = 3)
        seqs[ranges[ , i]] <-
          apply(tmp, MARGIN = 1,
                FUN = function(x) { edit_op(x[[1]], as.integer(x[[2]]), x[[3]]) })

      }
    }
  }

  if (include_counts) {
    counts <- stats::rbinom(samples * sample_size,
                            size = 300, prob = 0.1)
  }

  if (chains == 1) {
    dat <- data.frame(CloneSeq = seqs)
    if (include_counts) {
      dat$CloneCount <- dat$CloneFrequency <- counts
      for (i in 1:samples) {
        dat$CloneFrequency[ranges[ , i]] <-
          counts[ranges[ , i]] / sum(counts[ranges[ , i]])
      }
    }
  } else if (chains == 2) {
    alpha_seqs <- beta_seqs <- seqs
    modify_indices <- sample(
      c(TRUE, FALSE), size = sample_size * samples, replace = TRUE)
    beta_seqs[modify_indices] <-
      paste0(beta_seqs[modify_indices], "G")
    dat <- data.frame(AlphaSeq = alpha_seqs, BetaSeq = beta_seqs)
    if (include_counts) {
      dat$Count <- counts
      dat$UMIs <- pmax(1, stats::rbinom(sample_size * samples,
                                        size = 100, prob = 0.03))
    }
  } else {
    stop("`chains` must be 1 or 2")
  }
  dat$SampleID <- as.character(rep(1:samples, each = sample_size))

  if (!is.null(output_dir)) {
    for (i in 1:samples) {
      saveRDS(dat[ranges[ , i], , drop = FALSE],
              file = file.path(output_dir, paste0(i, ".rds")))
    }
  }

  if (!no_return) { return(dat) } else { return(TRUE) }
}


